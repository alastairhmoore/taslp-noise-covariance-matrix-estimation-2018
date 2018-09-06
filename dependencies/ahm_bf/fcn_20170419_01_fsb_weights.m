function[w] = fcn_20170419_01_fsb_weights(d)
%d is steering vector for each TF bin [nChans,nFreq,nFrames]
%debug_level controls how to handle singular or ill-conditioned matrices

% validate sizes
[nChans,nFreq,nFrames] = size(d);

%preallocate
w = zeros(nChans,nFreq,nFrames);

%for clarity and simplicity, loop over time and frequency indices
for ifreq = 1:nFreq
    for iframe = 1:nFrames
        w(:,ifreq,iframe) = d(:,ifreq,iframe)./(d(:,ifreq,iframe)'*d(:,ifreq,iframe));
    end
end
