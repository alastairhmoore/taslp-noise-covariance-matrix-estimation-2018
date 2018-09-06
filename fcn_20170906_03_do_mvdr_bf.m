function[z,w] = fcn_20170906_03_do_mvdr_bf(x,d,Rx)

debug_level = 0; % suppress all warnings
stop_on_singular = 0; %on returning catch singular values


[nFrames, nChans] = size(x);
[nCoefs,nExtras] = size(d);
if nCoefs~=nChans || nExtras~=1
    error('badly formed input')
end
[nFramesRx,nCovTerms] = size(Rx);
if nCovTerms~=nChans^2 || nFramesRx~=nFrames
    error('Rx dimensions don''t match')
end

% we know that d is constant
% check to see if Rx is actually constant
diff_Rx = diff(Rx,[],1);
if all(diff_Rx(:) == 0)
    % can just calculate the filter once
    %pause()
    Rx = reshape(Rx(1,:),nChans,nChans);
    [w,issingular] = fcn_20170222_01_mvdr_weights(Rx,d,debug_level);
    if stop_on_singular && any(issingular(:))
        keyboard
    end
    z = x * conj(w);
else
    % need to expand d to match Rx
    
    Rx = permute(reshape(Rx,nFrames,nChans,nChans),[2 3 4 1]);
    d = repmat(d,1,1,nFrames);
    [w,issingular] = fcn_20170222_01_mvdr_weights(Rx,d,debug_level);
    if stop_on_singular && any(issingular(:))
        %figure;plot(issingular,'x')
        keyboard
    end
    z = zeros(nFrames,1);
    for iframe = 1:nFrames
        z(iframe) = x(iframe,:) * conj(w(:,1,iframe));
    end
end
        
        