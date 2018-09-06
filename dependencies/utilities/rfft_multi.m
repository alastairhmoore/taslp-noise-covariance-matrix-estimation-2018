function x_f = rfft_multi(x,w_len,ov)
% RFFT_MULTI gives the STFT of multichannel real data, based on the rfft function in [1].
%
% INPUTS:
% x:  	 	time-domain multichannel signal, with one column as one channel.
% w_len: 	window length 
% ov:		overlap
%
% OUTPUT:
% x_f:   	n_mic*n_fft_2*n_frame STFT
%
% REFERENCES: 
% [1] D. M. Brookes, VOICEBOX: A speech processing toolbox for MATLAB, http://www.ee.ic.ac.uk/hp/staff/dmb/
% voicebox/voicebox.html, 1997-2017.

% Wei Xue, Imperial College, Feb 14, 2017

w=sqrt(hamming(w_len+1))'; 
w(end)=[]; % for now always use sqrt hamming window
inc_len = fix(w_len*(1-ov));
w=w/sqrt(sum(w(1:inc_len:w_len).^2)); 
[n_sig,n_mic] = size(x);
n_frame = max(fix((n_sig-w_len+inc_len)/inc_len),0)+1;
n_fft = w_len;
y_f = zeros(n_frame,n_fft/2+1,n_mic);

for i_mic = 1:n_mic
	y=enframe(x(:,i_mic),w,inc_len,'z');
	y_f(:,:,i_mic)=rfft(y,n_fft,2);
end

if n_mic>1
	x_f = permute(y_f,[3 2 1]);
else
	x_f = y_f.';
end