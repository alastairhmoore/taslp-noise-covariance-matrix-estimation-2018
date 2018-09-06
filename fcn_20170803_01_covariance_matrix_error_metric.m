function[J,scale_factor_saved] = fcn_20170803_01_covariance_matrix_error_metric(R_est,R_target,doScaleInvariant)
%   R_est: [nFrames, nChans^2, nSmooth]
%R_target: [nFrames, nChans^2, 1]

[nFrames, nChansSq, nSmooth] = size(R_est);

if doScaleInvariant
% horrible for loop just to check the idea
scale_factor_saved = zeros(nFrames,1,nSmooth);
for iframe = 1:nFrames
    for ismooth = 1:nSmooth
        this_vec = permute(R_est(iframe,:,ismooth),[2,1,3]); %[nChansSq,1]
        scale_factor_saved(iframe,1,ismooth) = real(R_target(iframe,:) * this_vec / (this_vec.' * this_vec));
    end
end

% apply the scale factor
R_est = bsxfun(@times,R_est,scale_factor_saved);
end

J = zeros(nFrames,nSmooth);
for ismooth = 1:nSmooth
    J(:,ismooth) = sum(abs(R_est(:,:,ismooth) - R_target).^2,2)./sum(abs(R_target).^2,2);
end

if nargout>1 && ~doScaleInvariant
    scale_factor_saved = ones(nFrames,1,nSmooth);
end