function x = irfft_multi(x_f,w_len,ov)
% IRFFT_MULTI gives the inverse DFT of multichannel data.
%
% INPUTS:
% x_f:  	n_mic*n_fft*n_frame DFT
% w_len: 	window length 
% ov:		overlap
%
% OUTPUT:
% x:   	 	n_sig*n_mic time-domain multichannel signal
%
% REFERENCES: 
% [1] D. M. Brookes, VOICEBOX: A speech processing toolbox for MATLAB, http://www.ee.ic.ac.uk/hp/staff/dmb/
% voicebox/voicebox.html, 1997-2017.

% Wei Xue, Imperial College, Feb 14, 2017

if ndims(x_f)<3
	x_f = permute(x_f,[3 1 2]);
end

w=sqrt(hamming(w_len+1))'; 
w(end)=[]; % for now always use sqrt hamming window
inc_len = fix(w_len*(1-ov));
w = w/sqrt(sum(w(1:inc_len:w_len).^2));

[n_mic,n_fft,n_frame] = size(x_f);
n_sig = n_frame*inc_len-inc_len+w_len;
x = zeros(n_sig,n_mic);

for i_mic = 1:n_mic
	x_f_mic = permute(x_f(i_mic,:,:),[3 2 1]);
	x_mic = overlapadd(irfft(x_f_mic,w_len,2),w,inc_len);
	x(:,i_mic) = x_mic(:);
end