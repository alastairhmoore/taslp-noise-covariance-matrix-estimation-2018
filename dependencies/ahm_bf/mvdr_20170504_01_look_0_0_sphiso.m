function[] = mvdr_20170504_01_look_0_0_sphiso(...
    in_wav_file_path ...          % path to received signal
    ,out_wav_file_path ...        % path where enhanced signal(s) should be written
    ,ht_file_path ...             % ignored - path to head tracker signal
    ,listener_characteristics ... % ignored - path to listener characteristics profile
    ,in_params ...                % structure with config information for test bench
    ,oracle_data ...              % ignored - structure containing information and/or paths to information which would not be available in real-world application but may be useful in development
    ,saved_data_dir ...           % path to folder where any supplementary/intermediate results should be written
    ,temp_data_dir ...            % path to folder where temporary files should be written - normally empty but may be useful for debugging
    )
%mvdr_20170504_01_look_0_0_sphiso implements an mvdr
%beamformer where the PSD/covariance matrix is determined numerically under
%an assumption that noise is diffuse and spherically isotropic. The look 
%direction is fixed as frontal with respect to the (possibly rotated) head
%
%Usage:
%  mvdr_20170504_01_look_0_0_sphiso(...
%      in_wav_file_path, out_wav_file_path, ht_file_path, listener_characteristics, ...
%      in_params, oracle_data, saved_data_dir, temp_data_dir)
%
%Outputs:
%  None
%
%Inputs:
%  in_wav_file_path:
%      received signal
%  out_wav_file_path:
%      path where enhanced signal(s) should be written
%  ht_file_path:
%      head tracker signal
%  listener_characteristics:
%      [not yet implemented]
%  in_params:
%      struct
%  oracle_data:
%      sruct containing information which would not be available in a
%      real-world application but may be useful in development
%  saved_data_dir:
%      folder for storing intermediate data for exchange between modules or
%      supplementary results which should not be included in metrics.
%      May be empty if test bench dictates that this data should not be retained
%      after the experiment, in which case intermediate data should be
%      saved to temporary storage (see tempname).
%  temp_data_dir:
%      folder for storing data which is not required once processing is
%      complete. Normally empty, in which case use standard temporary directory
%      (see tempname) but useful for debugging.
%
%Alastair Moore, May 2017

%% validate the input
% not really required in this simple function but will use it as a template

% in_wav_file_path
if ~exist(in_wav_file_path,'file'), error('Couldn''t find %s',in_wav_file_path), end

% final_data_dir
final_data_dir = fileparts(out_wav_file_path);
check_output_dir_exists(final_data_dir)

% ht_file_path
%if ~exist(ht_file_path,'file'), error('Couldn''t find %s',ht_file_path), end

% listener_characteristics (haven't specified this yet)
% if not present, assume default listener?

% oracle_data (struct)
% check the required field(s) are present
required_oracle_data = {...
    'ema_fcn'...
    };
for ireq = 1:length(required_oracle_data)
    if ~isfield(oracle_data,required_oracle_data{ireq})
        error('oracle_data is missing %s field',required_oracle_data{ireq})
    end
end

% saved_data_dir (string)
if nargin<7 || isempty(saved_data_dir)
    error('saved_data_dir should not be empty - want to retain the data')
end
check_output_dir_exists(saved_data_dir);

% temp_data_dir (string)
% if nargin<8 || isempty(temp_data_dir)
%     temp_data_dir = tempname;
% end
% check_output_dir_exists(temp_data_dir);

%% defualt parameters
params.fs = [];
params.c = soundspeed();
params.frame_duration = 0.01;
params.frame_overlap_frac = 0.5;

%% update with input parameters
if ~isempty(in_params)
    params = override_valid_fields(params,in_params);
end


%% read in data and do pre-processing
[x, fs_in] = audioread(in_wav_file_path); %[nSamples,nChans]
if isempty(params.fs)
    params.fs = fs_in;
else
    x = resample(x,params.fs,fs_in);
end
nReqOutSamples = size(x,1);
nChans = size(x,2);

% setup
fs = params.fs;
mvdr_fcn = @fcn_20170222_01_mvdr_weights;


% get a single impulse response so we can decide the filter length
ema = oracle_data.ema_fcn();
ema.prepareData(fs);
H_origin_t = ema.H0.';
t0m1 = ema.t0-1;
ir = ema.getImpulseResponseForSrc(0,degtorad(90)); % fixed look direction (assumed known)

len_filt = size(ir,1);

% frame lengths come directly from params
sig_frame_len = round((1/params.frame_overlap_frac) * params.frame_duration*fs);
sig_frame_inc = round(params.frame_duration * fs);

% stft properties for filter design of length len_filt
% the size of the fft is fixed by len_filt
% frame increment is fixed by the frame rate of the signal in the filtering
% stage
win_anal{1} = hamming(len_filt,'periodic');
win_anal{1} = win_anal{1}./(sum(win_anal{1}(1:sig_frame_inc:len_filt)));
win_anal{2} = ones(len_filt,1); %actually we never do the inverse fft

% stft properties for cola fast convolution
n_fft = sig_frame_len+len_filt-1;
win_cola{1} = hamming(sig_frame_len,'periodic');
win_cola{1} = win_cola{1}./(sum(win_cola{1}(1:sig_frame_inc:sig_frame_len)));
win_cola{2} = ones(n_fft,1);

% extend signal to allow for full convolution
x= [x; zeros(len_filt-1,nChans)];

% pad to ensure that all samples are processed (and incomplete blocks can
% be ignored
len_pad_pre = sig_frame_len-sig_frame_inc;
len_pad_post = sig_frame_len-sig_frame_inc;
x = [zeros(len_pad_pre,nChans);x;zeros(len_pad_post,nChans)];

nOutSamplesMax = size(x,1);
minOutFrames = ceil(nOutSamplesMax/sig_frame_inc)+ceil((sig_frame_len-sig_frame_inc)/sig_frame_inc);

% intialise the analysis/design stft
[~,x_tail_anal,pm_anal] = stft_v2('init',nChans,win_anal,sig_frame_inc,len_filt,fs);

% initialise the filtering stft
[~,x_tail_cola,pm_cola] = stft_v2('init',nChans,win_cola,sig_frame_inc,n_fft,fs);

%allocate mememory for output structure
y_fixed_mono = zeros(nOutSamplesMax,1);
y_fixed_mono_tail = [];


%% psd matrix is independent of the input signal or the array orientation
% but we do need to know the fft size

% define an equally spaced grid in azimuth and inclination
% set separation to ~10 degrees - based on voicebox's sphrharm.m 
n_inc = round(180/10);
n_az = round(360/10);
az_v = (0:n_az-1) * (2*pi/n_az); %vector of azimuths
inc_v=(1:2:2*n_inc)*pi*0.5/n_inc;%vector of inclinations
vq=(n_inc-abs(n_inc+1-2*(1:n_inc))).^(-1).*exp(-1i*(inc_v+0.5*pi));
weight_inc=(-2*sin(inc_v).*real(fft(vq).*exp(-1i*(0:n_inc-1)*pi/n_inc))/n_inc);
%----
az_frac = 2*pi/n_az;

% get a full grid
[inc_g,az_g] = ndgrid(inc_v(:),az_v(:));
quad_g = repmat(1/(4*pi) * az_frac * weight_inc(:),1,length(az_v));
if sum(quad_g(:))-1 > 1e-14,error('quadrature weights don''t sum to one'),end

% for spherically isotropic noise the power is equal in all directions
pow_g = ones(size(inc_g));

% get the psd matrix
R = fcn_20170421_01_diffuse_pow_to_psd_matrix(oracle_data.ema_fcn(),...
    fs,len_filt,az_g,inc_g,pow_g,quad_g); % [nSenors,nSensors,nFreq]

% define the fixed beamformers using the steering vector obtained earlier
H_fixed_t = rfft(ir,[],1).'; % [nChans,nFreq]

% relative transfer function wrt origin
rH_fixed_mono_t = bsxfun(@rdivide,H_fixed_t,H_origin_t);

% get filter weights in freq domain
W_fixed_mono_t = conj(mvdr_fcn(R,rH_fixed_mono_t));

% filter in time domain
w_fixed_mono_t = irfft(W_fixed_mono_t,[],2);

% needs to be causal
w_fixed_mono_t = circshift(w_fixed_mono_t,t0m1,2);

% filter in freq domain ready for cola fast convolution
W_fixed_mono = rfft(w_fixed_mono_t,n_fft,2).';





%% everything else needs to be done 'online'
fprintf('Entering outer loop - processing in chunks\n');
block_size = 1000;
nBlocks = ceil(length(x)/block_size);
idc_in = 0;
idc_out = 0;
t_fr = [];
nFreq_cola = length(pm_cola.f);
for iblock = 1:nBlocks
    idc_in_start = idc_in(end)+1;
    idc_in_end = min(idc_in(end)+block_size,length(x));
    idc_in = idc_in_start:idc_in_end;
    
    fprintf('Fwd STFTs  ');
    [X_cola,x_tail_cola,pm_cola] = stft_v2('fwd',x(idc_in,:),pm_cola,x_tail_cola);
        
    % work out how many frames we can process
    nFrames = size(X_cola,3);
    
    fprintf('Apply  ');    
    Y_fixed_mono = sum(bsxfun(@times,X_cola,W_fixed_mono),2);
    
    fprintf('Inv STFTs  ');
    [tmp_y_fixed_mono,y_fixed_mono_tail] = stft_v2('inv',Y_fixed_mono,pm_cola,y_fixed_mono_tail);
    
    % indices to fill
    idc_out_start = idc_out(end)+1;
    idc_out_end = idc_out_start+size(tmp_y_fixed_mono,1)-1;
    idc_out = idc_out_start:idc_out_end;
    
    % fill them
    y_fixed_mono(idc_out,:) = tmp_y_fixed_mono;
    
    fprintf('\n');
end
%% append the last frames' tails
idc_out_start = idc_out(end)+1;
idc_out_end = idc_out_start+size(y_fixed_mono_tail,1)-1;
idc_out = idc_out_start:idc_out_end;

y_fixed_mono(idc_out,:) = y_fixed_mono_tail;

%% remove the padding
% from beginning
y_fixed_mono(1:len_pad_pre,:) = [];

% filters in this case are non-causal
% to ensure the output is the same length and time aligned with the input,
% need to remove the pre- and post- ringing
%t0m1 = t0-1; %is the number of samples to remove for the original
y_fixed_mono(1:t0m1,:) = [];

% can now just crop to the required length
y_fixed_mono(nReqOutSamples+1:end,:) = [];

%% write out
audiowrite(out_wav_file_path,y_fixed_mono,fs);
