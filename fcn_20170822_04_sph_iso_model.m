function[Rx_est,model_params] = fcn_20170822_04_sph_iso_model(xxH,smooth,Rx_init,tildeHnm,Acell,per_frame_pose)
% Estimate of Rx obtained using spherically isotropic model of diffuse field
%
% !!! N.B. No scaling is currently implemented so absolute error will be high !!!
%
% inputs
%            xxH: Per frame microphone covariance matrix [nFrames,nCovTerms]
%         smooth: forgetting factor
%        Rx_init: Initialisation vector for Rx_est
%       tildeHnm: SH model of complex conjugate of array manifold
%          Acell: SH covariance of diffuse sound field with power distribution
%                 described by sum of SHs
% per_frame_pose: struct containing column vectors 'roll','pitch' and 'yaw'
%                 each of which is [nFrames 1]

sz_Rx = size(xxH);
if numel(sz_Rx)~=2, error('Rx_in should be 2D matrix'),end
nFrames = sz_Rx(1);
nCovTerms = sz_Rx(2);
Rx_est = zeros(nFrames,nCovTerms);

% Use Acell to initialise SH dimensions
nSHDiffuse = numel(Acell);
sz_SHCov = size(Acell{1});
nSHPWD = unique(sz_SHCov);
if ~isequal(sz_SHCov,[nSHPWD nSHPWD]),error('Elements of Acell should be square'),end

sphOrdDiffuse = sqrt(nSHDiffuse)-1;
sphOrdPWD = sqrt(nSHPWD)-1;
if rem(sphOrdDiffuse,1)~=0,error('Number of diffuse SHs is incorrect'),end
if rem(sphOrdPWD,1)~=0,error('Number of pwd SHs is incorrect'),end

% Validate size of tildeHnm
sz_Hnm = size(tildeHnm);
nSHPWD_Hnm = sz_Hnm(1);
nSensors = sz_Hnm(2);
if prod(sz_Hnm)~=nSHPWD_Hnm*nSensors,error('tildeHnm should be a 2D matrix'),end
if nSHPWD_Hnm~=nSHPWD,error('Number of SHs in tildeHnm is different to Acell'),end
if nCovTerms~=nSensors^2,error('Number of sensors in tildeHnm is different to Rx_in'),end

% Mapping of diffuse SH components to sensor covariance
B = zeros(nSensors,nSensors,nSHDiffuse);
for ish = 1:nSHDiffuse
    B(:,:,ish) = tildeHnm' * Acell{ish} * tildeHnm;
end
Bt = reshape(B,nCovTerms,nSHDiffuse).'; % [nSHDiffuse nCovTerms]

% initialise to spherically isotropic
fnm_est = zeros(nSHDiffuse,1);
fnm_est(1) = 1;

% if initialisation is given then use it to scale
Rx_sph_iso = fnm_est' * Bt;

Rx_est = repmat(Rx_sph_iso,nFrames,1);

model_params.estimator_name = 'spherically isotropic model';
