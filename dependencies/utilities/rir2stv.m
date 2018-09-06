function steering_mat = rir2stv(rir,n_fft_2)
% RIR2STV computes the steering vectors for according to the room impulse responses (RIRs).
%
% INPUTS:
% rir:  	multichannel RIRs where each column represents one channel.
% n_fft_2:	(FFT length)/2+1
% 
% OUTPUT:
% steering_mat: n_mic*n_fft_2 steering matrix
%
% REFERENCES: 
% [1] Benesty, Jacob, Jingdong Chen, and Yiteng Huang.
% Microphone Array Signal Processing. Springer Science & Business Media,
% 2008.
% 
% Wei Xue, Imperial College, Feb 14, 2017

n_mic = size(rir,2);
n_fft = 2*(n_fft_2-1);

rir_f = fft(rir,n_fft).';
steering_mat = rir_f./repmat(rir_f(1,:)+eps,n_mic,1);
steering_mat = steering_mat(:,1:n_fft_2);