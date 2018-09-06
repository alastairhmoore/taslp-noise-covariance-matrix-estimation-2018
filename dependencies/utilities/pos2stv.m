function steering_mat = pos2stv(air,n_fft_2)
% POS2STV computes the steering vectors for one position.
%
% INPUTS:
% air:  	acoustic environment
%			air.rcv: 3*n_mic microphone positions
%			air.src: 3*1 source position 
%			air.c:	 speed of sound in the air
% 			air.fs:		sampling rate
% n_fft_2:	FFT length
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

n_mic = size(air.rcv,2);

distance_vec = sqrt(sum((air.rcv-repmat(air.src,1,n_mic)).^2));
distance_diff = (distance_vec-distance_vec(1))/air.c*air.fs;

steering_mat = exp(2i*pi*distance_diff(:)*(0:n_fft_2-1)/(2*(n_fft_2-1)));