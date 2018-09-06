function X2 = stft2cpsd(X,Y,alpha_s,n_average)
% STFT2CPSD returns the cross psd using exponential smoothing.
%
% INPUTS:
% X: 3-D n_mic*n_fft*n_frame STFT matrix.
% Y: 3-D n_mic*n_fft*n_frame STFT matrix.
% alpha_s: exponential smoothing factor.
%
% OUTPUT:
% X2:  3-D n_mic*n_fft*n_frame PSD matrix
%
% Wei Xue, Imperial College, Mar 08, 2017

k=0;
if ndims(X)<3
	X = permute(X,[3 1 2]);
	Y = permute(Y,[3 1 2]);
    k=1;
end

[~,~,n_frame] = size(X);
X2 = X.*conj(Y);

if nargin < 4
	for i_frame = 2:n_frame	
		X2(:,:,i_frame) = alpha_s*X2(:,:,i_frame-1) + (1-alpha_s)*X(:,:,i_frame).*conj(Y(:,:,i_frame));	
	end
else
	for i_frame = 2:n_frame
		if i_frame >= n_average
			X2(:,:,i_frame) = mean(X(:,:,i_frame-n_average+1:i_frame).*conj(Y(:,:,i_frame-n_average+1:i_frame)),3);
		else
			X2(:,:,i_frame) = mean(X(:,:,1:i_frame).*conj(Y(:,:,1:i_frame)),3);
		end
	end
end

if k==1
	X2 = permute(X2,[2 3 1]);
end