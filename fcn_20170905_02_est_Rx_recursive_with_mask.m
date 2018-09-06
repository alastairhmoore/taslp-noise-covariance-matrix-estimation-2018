function[Rx_est,model_params] = fcn_20170905_02_est_Rx_recursive_with_mask(xxH,smooth,Rx_init,mask)
% Estimate of Rx obtained using recursive smoothing
% inputs
%     xxH: Per frame microphone covariance matrix [nFrames,nCovTerms]
%  smooth: Smoothing factor
% Rx_init: Optional intialisation of estimate [1, nCovTerms]
%    mask: Binary speech abscence indicator - only update estimate when 1 [nFrames,1]
%
% N.B. If R_init is empty and mask(1)==0 there will be an error

sz_Rx = size(xxH);
if numel(sz_Rx)~=2, error('Rx_in should be 2D matrix'),end
nFrames = sz_Rx(1);
nCovTerms = sz_Rx(2);
Rx_est = zeros(nFrames,nCovTerms);

if nargin < 4 || isempty(mask)
    mask = ones(nFrames,1);
else
    mask = mask(:);
    if size(mask)~=[nFrames, 1]
        error('mask does not match dimensions of xxH')
    end
end

if nargin<3 || isempty(Rx_init)
    if mask(1)==0,error('Cannot initialise Rx_est for first frame'),end
    Rx_est(1,:) = xxH(1,:);
else
    if mask(1)
        Rx_est(1,:) = smooth * Rx_init + (1-smooth) * xxH(:,1);
    else
        Rx_est(1,:) = Rx_init;
    end
end

for iframe = 2:nFrames
    if mask(iframe)
        Rx_est(iframe,:) = smooth * Rx_est(iframe-1,:) ...
                           + (1-smooth) * xxH(iframe,:);
    else
        Rx_est(iframe,:) = Rx_est(iframe-1,:);
    end
end

model_params.estimator_name = 'recursive smoothing'; 
model_params.smooth_factor = smooth;