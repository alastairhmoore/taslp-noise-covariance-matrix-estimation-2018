function X2 = stft2psd(X,alpha_s)
% STFT2PSD returns the psd using exponential smoothing.
%
% INPUTS:
% X: 3-D n_mic*n_fft*n_frame STFT matrix.
% alpha_s: exponential smoothing factor.
%
% OUTPUT:
% X2:  3-D n_mic*n_fft*n_frame PSD matrix
%
% Wei Xue, Imperial College, Feb 14, 2017

[~,~,n_frame] = size(X);
X2 = abs(X).^2;
for i_frame = 2:n_frame	
	X2(:,:,i_frame) = alpha_s*X2(:,:,i_frame-1) + (1-alpha_s)*abs(X(:,:,i_frame)).^2;
end