function[] = fsb_20170419_03_binaural_look_0_0(...
    in_wav_file_path ...          % path to received signal
    ,out_wav_file_path ...        % path where enhanced signal(s) should be written
    ,ht_file_path ...             % ignored - path to head tracker signal
    ,listener_characteristics ... % ignored - path to listener characteristics profile
    ,in_params ...                % structure with config information for test bench
    ,oracle_data ...              % ignored - structure containing information and/or paths to information which would not be available in real-world application but may be useful in development
    ,saved_data_dir ...           % path to folder where any supplementary/intermediate results should be written
    ,temp_data_dir ...            % path to folder where temporary files should be written - normally empty but may be useful for debugging
    )
%tb_20170321_ahm_01_bf_fixed_and_tracking gives the results of two basic beamformers
%one which has a fixed look direction straight ahead the other uses the
%ground truth relative doa of the desired source
%
%Usage:
%  fsb_20170419_03_binaural_look_0_0(...
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
%Alastair Moore, April 2017

%% validate the input
% not really required in this simple function but will use it as a template

% in_wav_file_path
if ~exist(in_wav_file_path,'file'), error('Couldn''t find %s',in_wav_file_path), end

% final_data_dir
final_data_dir = fileparts(out_wav_file_path);
check_output_dir_exists(final_data_dir)

% ht_file_path
if ~exist(ht_file_path,'file'), error('Couldn''t find %s',ht_file_path), end

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

% stft properties for cola fast convolution
n_fft = sig_frame_len+len_filt-1;
win{1} = hamming(sig_frame_len,'periodic');
win{1} = win{1}./(sum(win{1}(1:sig_frame_inc:sig_frame_len)));
win{2} = ones(n_fft,1);

% define the fixed beamformers using the steering vector obtained earlier
H_fixed_t = rfft(ir,[],1).'; % [nChans,nFreq]

% relative transfer function wrt reference channels
rH_fixed_left_t = bsxfun(@rdivide,H_fixed_t,H_fixed_t(ema.refChanLeft,:));
rH_fixed_right_t = bsxfun(@rdivide,H_fixed_t,H_fixed_t(ema.refChanRight,:));

% get filter weights in freq domain
W_fixed_left_t = conj(fcn_20170419_01_fsb_weights(rH_fixed_left_t));
W_fixed_right_t = conj(fcn_20170419_01_fsb_weights(rH_fixed_right_t));

% filter in time domain
w_fixed_left_t = irfft(W_fixed_left_t,[],2);
w_fixed_right_t = irfft(W_fixed_right_t,[],2);

% needs to be causal
w_fixed_left_t = circshift(w_fixed_left_t,t0m1,2);
w_fixed_right_t = circshift(w_fixed_right_t,t0m1,2);

% filter in freq domain ready for cola fast convolution
W_fixed_left =rfft(w_fixed_left_t,n_fft,2).';
W_fixed_right = rfft(w_fixed_right_t,n_fft,2).';

% extend signal to allow for full convolution
x= [x; zeros(len_filt-1,nChans)];

% pad to ensure that all samples are processed (and incomplete blocks can
% be ignored
len_pad_pre = sig_frame_len-sig_frame_inc;
len_pad_post = sig_frame_len-sig_frame_inc;
x = [zeros(len_pad_pre,nChans);x;zeros(len_pad_post,nChans)];

nOutSamplesMax = size(x,1);

% initialise the filtering stft
[~,x_tail,pm] = stft_v2('init',nChans,win,sig_frame_inc,n_fft,fs);

%allocate mememory for output structure
y_fixed_bin = zeros(nOutSamplesMax,2);
y_fixed_left_tail = [];
y_fixed_right_tail = [];

%% everything else needs to be done 'online'
fprintf('Entering outer loop - processing in chunks\n');
block_size = 1000;
nBlocks = ceil(length(x)/block_size);

idc_in = 0;
idc_out = 0;
for iblock = 1:nBlocks
    idc_in_start = idc_in(end)+1;
    idc_in_end = min(idc_in(end)+block_size,length(x));
    idc_in = idc_in_start:idc_in_end;
    
    fprintf('Fwd STFTs  ');
    [X,x_tail,pm] = stft_v2('fwd',x(idc_in,:),pm,x_tail);
    
    % initialise inner loop over each frame
    [nFreq,~,nFrames] = size(X);
    
    fprintf('Apply  ');    
    Y_fixed_left = sum(bsxfun(@times,X,W_fixed_left),2);
    Y_fixed_right = sum(bsxfun(@times,X,W_fixed_right),2);
    
    
    fprintf('Inv STFTs  ');
    [tmp_y_fixed_left,y_fixed_left_tail] = stft_v2('inv',Y_fixed_left,pm,y_fixed_left_tail);
    [tmp_y_fixed_right,y_fixed_right_tail] = stft_v2('inv',Y_fixed_right,pm,y_fixed_right_tail);
    
    % indices to fill
    idc_out_start = idc_out(end)+1;
    idc_out_end = idc_out_start+size(tmp_y_fixed_left,1)-1;
    idc_out = idc_out_start:idc_out_end;
    
    % fill them
    y_fixed_bin(idc_out,1) = tmp_y_fixed_left;
    y_fixed_bin(idc_out,2) = tmp_y_fixed_right;

    fprintf('\n');
end
%% append the last frames' tails
idc_out_start = idc_out(end)+1;
idc_out_end = idc_out_start+size(y_fixed_left_tail,1)-1;
idc_out = idc_out_start:idc_out_end;

y_fixed_bin(idc_out,1) = y_fixed_left_tail;
y_fixed_bin(idc_out,2) = y_fixed_right_tail;

%% remove the padding
% from beginning
y_fixed_bin(1:len_pad_pre,:) = [];

% filters in this case are non-causal
% to ensure the output is the same length and time aligned with the input,
% need to remove the pre- and post- ringing
%t0m1 = t0-1; %is the number of samples to remove for the original
y_fixed_bin(1:t0m1,:) = [];

%can now just crop to the required length
y_fixed_bin(nReqOutSamples+1:end,:) = [];

%% write out
audiowrite(out_wav_file_path,y_fixed_bin,fs);
