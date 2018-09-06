function [y,w_mat] = beamforming_hr(sig,air,pp)
% BEAMFORMING denoises signal using frequency-domain beamformers.
%
% INPUTS:
% sig:  time-domain multichannel signal, with one column as one channel.
% air:	acoutic environment
%       air.fs: sampling rate
%       air.h: room impulse responses(RIRs). The RIRs are only used in[3]
%       to compute the relative transfer function of early RIRs.
%       OPTIONAL air.noi: if air contains this field, the algorithm use
%       oracle noise information for mvdr/mwf.
%       Other fields of air may change depending on whether the scenario is
%       simulated or generated using the Oldenburg database.
% pp:   parameters of the algorithm. Default values are shown in brackets.
%       pp.win_len(256): Window size of short-time-Fourier-transform
%       (STFT). 
%       pp.overlap(0.5): Overlap between frames in range [0,1]
%       pp.f_debug(0):   Whether operates in the debug mode 
%       pp.f_rtf(0):     The steering vectors are computed using RIRs if f_rtf == 1.
%                        Otherwise, use the geometry information contained in "air" instead. 
%       pp.beamformer('mwf'): choose the beamformer type
%       pp.noi_est('s'): If 's', use souden's method [4] for noise
%       covariance matrix estimation. If 'h', use [3] instead.
%       pp.channel(1): the reference channel
% OUTPUTS:
% y:    time-domain single-channel beamformer output
% w_mat:n_mic*n_fft_2*n_frame matrix consisting of frequency-domain filter
%       coefficients in each frame.
%
% REFERENCES: 
% [1] Benesty, Jacob, Jingdong Chen, and Yiteng Huang.
% Microphone Array Signal Processing. Springer Science & Business Media,
% 2008.
% [2] Pan, Chao, Jingdong Chen, and J. Benesty. Performance Study
% of the MVDR Beamformer as a Function of the Source Incidence Angle.
% IEEE/ACM Transactions on Audio, Speech, and Language Processing 22, no. 1
% (2014): 67-79.
% [3] Hendriks, R.C., and T. Gerkmann. Noise Correlation Matrix Estimation for
% Multi-Microphone Speech Enhancement. IEEE Transactions on Audio, Speech, and Language 
% Processing 20, no. 1 (2012): 223-233. doi:10.1109/TASL.2011.2159711.
% [4] Souden, M., Jingdong Chen, J. Benesty, and S. Affes. an Integrated Solution for Online Multichannel Noise Tracking and Reduction.
% IEEE Transactions on Audio, Speech, and Language Processing 19, no. 7 (September 2011).
% 
% Wei Xue, Imperial College, Feb 14, 2017


%--- Check and unpack parameters ---
if nargin<3
    pp = [];
end

if isstruct(pp) || isempty(pp)
    default_fields = {'win_len','overlap','f_debug','f_rtf'};
    default_values = [256, 0.5, 0, 0];
    n_fields = length(default_fields);
	for i_field = 1:n_fields
		if ~isfield(pp,default_fields{i_field})
			eval([default_fields{i_field} '=' num2str(default_values(i_field)) ';']);
		else
			eval([default_fields{i_field} '=pp.' default_fields{i_field} ';']);
		end
	end

    if ~isfield(pp,'channel')
        pp.channel = 1;
    end
	
    if ~isfield(pp,'beamformer')
        pp.beamformer = 'mwf';
    end

    if ~isfield(pp,'noi_est')
        pp.noi_est = 's';
    end 

else
	error('please specify the parameters.');	
end

if ~isfield(air,'rollInterpolant') || ~isfield(air,'pitchInterpolant') || ~isfield(air,'yawInterpolant') || ~isfield(air,'ema')
    error('please specify the trajectory of head rotation.')
end

if strcmp(pp.beamformer,'mwf')    
    ppbf.channel = pp.channel;
    ppbf.w = 3;
end

sig_f_multi = rfft_multi(sig,win_len,overlap);
[n_mic,n_fft_2,n_frame] = size(sig_f_multi);


%--------------------perform mvdr/wmf------------------------
%Initialize ryy and rvv, set parameters
if strcmp(pp.noi_est,'h') || strcmp(pp.beamformer,'mvdr')
    if isempty(air)
        error('missing the acoustic information.');
    else
        if ~f_rtf
            steering_mat = pos2stv(air,n_fft_2);
        else
            steering_mat = 0;
        end
    end
end

sig_f = rfft_multi(sig,win_len,overlap);
ryy_mat = stft2rxx(sig_f,0.92);

if isfield(air,'noi') %non-blind, the noise signal is a field of pp
	noi_f = rfft_multi(air.noi,win_len,overlap);
	rvv_mat = stft2rxx(noi_f,0.92);
	noi_flg = 1;
else %blind
    rvv = mean(ryy_mat(:,:,:,1:12),4);
    switch pp.noi_est %set parameters for noise estimation
        case 'h'
            pyy = abs(sig_f).^2;        
            noi_psd_mat = zeros(n_mic,n_fft_2,n_frame);  
            for i_mic = 1:n_mic
                noi_tmp = estnoiseg(permute(pyy(i_mic,:,:),[3 2 1]),0.016);
                noi_psd_mat(i_mic,:,:) = permute(noi_tmp,[3 2 1]);
            end       
            rpy = zeros(n_mic,n_mic,n_fft_2);
            qq = struct('steering_mat',steering_mat,'rpy',rpy,'eq',11);
        otherwise
            spp = zeros(n_fft_2,1);
            qq = struct('spp',spp,'phi_0',-2e3);
    end
	noi_flg = 0;
end

n_ref = length(pp.channel);
h_bf = zeros(n_mic,n_ref);
clean_sig_f = zeros(n_ref,n_fft_2,n_frame);
w_mat = zeros(n_mic,n_ref,n_fft_2,n_frame);
h_bf_mat = zeros(n_mic,n_ref,n_fft_2);



if f_debug
    h_bf_frames = zeros(n_mic,n_fft_2,n_frame);
end

%framewise processing
for i_frame = 1:n_frame
    
    %_______________rir of current frame__________
    [air.h,t0] = inst_h(i_frame,air,win_len,overlap);
    len_filt = length(air.h(:,1));
    ir_origin = zeros(len_filt,1);
    ir_origin(t0) = 1;
    air.h(:,end+1) = ir_origin;
    steering_mat = rir2stv(air.h,n_fft_2);



	sig_f = sig_f_multi(:,:,i_frame);
    ryy = ryy_mat(:,:,:,i_frame);
	%_______________Noise estimation_______________
    if noi_flg > 0.5
        rvv = rvv_mat(:,:,:,i_frame); %use oracle noise       
    else     
        switch pp.noi_est  %noise estimation by different methods, default:souden's method
            case 'h'
                qq.noi_mat = noi_psd_mat(:,:,i_frame);
                [rvv,qq] = estnoicov_h(sig_f,rvv,qq);
            otherwise
                qq.ryy = ryy;
                [rvv,qq] = estnoicov_s(sig_f,rvv,qq);
        end    
    end

    %_______________beamforming_______________
    for i_freq = 1:n_fft_2
           
        switch pp.beamformer %computing beamformer coefficients, default: mwf
            case 'mvdr'
                for i_ref = 1:n_ref
                    h_bf(:,i_ref) = h_mvdr(rvv(:,:,i_freq),steering_mat(1:n_mic,i_freq)./steering_mat(i_ref,i_freq)); %mvdr
                end
            otherwise
                h_bf = h_mwf(ryy(:,:,i_freq),rvv(:,:,i_freq),ppbf); %wmf
        end

        clean_sig_f(:,i_freq,i_frame) = permute(h_bf'*sig_f(:,i_freq),[3 1 2]); %spatial filtering
        h_bf_mat(:,:,i_freq) = h_bf;
    end
    
    w_mat(:,:,:,i_frame) = h_bf_mat;
    h_bf_frames(:,:,i_frame) = permute(h_bf_mat(:,1,:),[1 3 2]);
end
y = irfft_multi(clean_sig_f,win_len,overlap);
y = y(1:length(sig(:,1)),:);


function [ir,t0] = inst_h(i_frame,air,win_len,overlap)
% INST_H returns the ir of one frame.
%

tt = (i_frame-1)*win_len*(1-overlap)/air.fs;
roll_fr = air.rollInterpolant(tt);
pitch_fr = air.pitchInterpolant(tt);
yaw_fr = air.yawInterpolant(tt);
ax_vec_rot = rotate_cart_axes_ypr(yaw_fr,pitch_fr,roll_fr);
[src_az_rot,src_inc_rot] = dir_wrt_rotated_axes(ax_vec_rot,air.src_az,air.src_inc);
[ir,t0] = air.ema.getIRForSrcDoa(src_az_rot,src_inc_rot);