function[out,out_tail] = istft(X,pm,in_tail)
%ISTFT takes 3D short time fourier transform matrix at n_freq=fix(1+n_fft/2)
%frequencies and returns 2D multichannel signal matrix
%
% Inputs
%       X  STFT matrix [n_freq, n_channels, n_frames]
%      pm  structure returned from stft
%
% Outputs
%       x  output signal [n_samples n_channels]
%
% 

if nargin < 3
    block_based = 0;
else
    block_based = 1;
end

% dim = size(X);
% n_freqs = dim(1);
% if length(dim)==2
%     %need to check whether we have only 1 channel or 1 frame
%     if length(pm.fr_st)~=1;
%         %definitel should have more than one frame so must have only one
%         %channel
%         n_chans = 1;
%         n_frames = dim(2);
%         X = permute(X,[1 3 2]);
%     else
%         n_chans = dim(2);
%         n_frames = 1;
%     end       
% else
%     n_chans = dim(2);
%     n_frames = dim(3);
% end
[n_freqs,n_chans,n_frames] = size(X);

% do ifft
%x = irfft(X);
x = ifft(X, pm.n_fft, 'symmetric'); %this is ~5 times faster that irfft

% remove fft padding
x = x((pm.fft_pre_pad+1):(end-pm.fft_post_pad),:,:);

% apply window
x = x .* repmat(pm.w(:),[1, n_chans, n_frames]);

%% overlap add frames
% using sparse and clever indexing
%ii = repmat(pm.fr_idc,1,n_chans);
%ij = repmat(1:n_chans,numel(pm.fr_idc),1);
%x = permute(x, [1 3 2]);
%out = full(sparse(ii(:),ij(:),x(:),pm.len_x_pad,n_chans));

% using old school for loop is ~10 times faster
n_w = length(pm.w);
fr_st = 1 + [0:n_frames-1]*pm.inc;          %frame start indices
fr_idc = bsxfun(@plus,fr_st,[0:n_w-1]');    %frame indices fully expanded
n_samples = fr_idc(end,end);
out = zeros(n_samples,n_chans);
for n = 1:n_frames
    out(fr_idc(:,n),:) = out(fr_idc(:,n),:) + x(:,:,n);
end

%deal with tail / padding
if block_based
    %add the previous tail to the leading samples
    if ~isempty(in_tail)
        out(1:size(in_tail,1),:) = out(1:size(in_tail,1),:) + in_tail;
    end
    %retain the unfilled samples for the next block - starting from ninc
    %samples after last frame start
    out_tail = out(fr_st(end)+pm.inc:end,:);
    out(fr_st(end)+pm.inc:end,:) = [];
else
    out(end-pm.post_pad_len+[1:pm.post_pad_len],:) = [];
    out([1:pm.pre_pad_len],:) = [];
end
