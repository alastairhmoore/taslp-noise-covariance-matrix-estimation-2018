function[Rv_est,model_params] = fcn_20180327_01_ewls_ncm_est(vvH,lambda,rv_init,tildeHnm,Gcell,per_frame_pose,mask,in_params)
% Estimate of Rx obtained using SH model of noise field power distribution
% inputs
%            vvH: Per frame microphone covariance matrix [nFrames,nCovTerms]
%         lambda: exponential weighting factor
%        rv_init: Initialisation vector for rv_est
%       tildeHnm: SH model of complex conjugate of array manifold
%          Gcell: SH covariance of diffuse sound field with power distribution
%                 described by sum of SHs
% per_frame_pose: struct containing column vectors 'roll','pitch' and 'yaw'
%                 each of which is [nFrames 1]
%           mask: Binary speech abscence indicator - only update model estimate when 1 [nFrames,1]

show_C_updates = 0;

params.sensor_noise_model = 2;
allowed_vals.sensor_noise_model = [0 ...  % assume no noise
                                   ,1 ... % assume diagonal noise
                                   ,2 ... % independent parameter for each sensor
                                   ];   
                               
if nargin>7 && ~isempty(in_params)
    params = override_valid_fields(params,in_params,allowed_vals);
end

if ~isempty(rv_init)
    warning('rv_init is ignored')
end

sz_vvH = size(vvH);
if numel(sz_vvH)~=2, error('vvH should be 2D matrix'),end
nFrames = sz_vvH(1);
nCovTerms = sz_vvH(2);
vvH = vvH.'; %[nCovTerms,nFrames]
Rv_est = zeros(nCovTerms,nFrames);
autoCorrEst = zeros(nCovTerms,nCovTerms);
crossCorrEst = zeros(nCovTerms,1);

if nargin < 7 || isempty(mask)
    mask = ones(nFrames,1);
else
    mask = mask(:);
    if size(mask)~=[nFrames, 1]
        error('mask does not match dimensions of xxH')
    end
end

% Use Gcell to initialise SH dimensions
nSHDiffuse = numel(Gcell);
sz_SHCov = size(Gcell{1});
nSHPWD = unique(sz_SHCov);
if ~isequal(sz_SHCov,[nSHPWD nSHPWD]),error('Elements of Gcell should be square'),end

sphHarmOrdDiffuse = sqrt(nSHDiffuse)-1;
sphHarmOrdPWD = sqrt(nSHPWD)-1;
if rem(sphHarmOrdDiffuse,1)~=0,error('Number of diffuse SHs is incorrect'),end
if rem(sphHarmOrdPWD,1)~=0,error('Number of pwd SHs is incorrect'),end

% Validate size of tildeHnm
sz_Hnm = size(tildeHnm);
nSHPWD_Hnm = sz_Hnm(1);
nSensors = sz_Hnm(2);
if prod(sz_Hnm)~=nSHPWD_Hnm*nSensors,error('tildeHnm should be a 2D matrix'),end
if nSHPWD_Hnm~=nSHPWD,error('Number of SHs in tildeHnm is different to Gcell'),end
if nCovTerms~=nSensors^2,error('Number of sensors in tildeHnm is different to Rx_in'),end

% Mapping of diffuse SH components to sensor covariance
Gmat = zeros(nSHPWD^2,nSHDiffuse);
for ish = 1:nSHDiffuse
    Gmat(:,ish) = Gcell{ish}(:);
end
B = kron(tildeHnm.',tildeHnm') * Gmat; %[nCovTerms nSHDiffuse]


% Use ElobesMicArray object to handle conversion of rotations from roll,
% pitch, yaw to Euler angles required by Wigner D
ema = RigidSphereHearingAidArray4; %N.B. Choice is arbitrary since we don't get impulse response from it


% mapping of spatially white paramteters to microphone covariance terms
% is direct
switch params.sensor_noise_model
    case 0
        M = [];
        nWhiteParams = 0;
    case 1
        M = reshape(eye(nSensors),nCovTerms,1);
        nWhiteParams = 1;
    case 2
        M = zeros(nSensors,nSensors,nSensors);
        for isensor = 1:nSensors
            M(isensor,isensor,isensor) = 1;
        end
        M = reshape(M,nCovTerms,nSensors);
        nWhiteParams = nSensors;
end

% initialise paramters to spherically isotropic with no spatially white noise
nParams = nSHDiffuse+nWhiteParams;
%fnm_est = zeros(nParams,1);
%fnm_est(1) = 1;

% % if initialisation is given then just use it to scale the isotropic
% % (rotation invariant) estimate
% rv_sph_iso = [B,M] * conj(fnm_est);
% if ~isempty(rv_init)
%     scale_factor = real(rv_init.' * rv_sph_iso.' / (rv_sph_iso * rv_sph_iso.'));
%     fnm_est = scale_factor * fnm_est;
% end

fnm_est_hist = zeros(nFrames,nParams);
        
% - loop over frames
reverseStr='';
for iframe = 1:nFrames
    if rem(iframe,100)==1
        msg = sprintf('Processing %d/%d', iframe, nFrames);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    if iframe==1 ...
            || per_frame_pose.roll(iframe) ~= per_frame_pose.roll(iframe-1) ...
            || per_frame_pose.pitch(iframe) ~= per_frame_pose.pitch(iframe-1) ...
            || per_frame_pose.yaw(iframe) ~= per_frame_pose.yaw(iframe-1)
        
        % ** need to update our model **
        if show_C_updates
            msg = sprintf('\nNew segment - update C\n');
            fprintf(msg);
            reverseStr = sprintf('');
        end
        
        
        % apply rotation to obtain euler angles
        ema.setPoseRollPitchYaw(per_frame_pose.roll(iframe),...
            per_frame_pose.pitch(iframe),...
            per_frame_pose.yaw(iframe))
        [alpha,beta,gamma] = fcn_20180326_02_body_rotation_qr_to_wignerD_euler(ema.poseQuaternion);
        D = sparse(fcn_20170616_02_myWignerDMatrix_Jacobi(sphHarmOrdDiffuse,...
            alpha,beta,gamma));
        C = [B*D.',M];
        
        %CCH = C*C';
        %CHC = C' * C;
        CTconjC = C.' * conj(C);
        if iframe==1
            autoCorrEst = zeros(nParams,nParams);
            crossCorrEst = zeros(nParams,1);
            warning('off','MATLAB:nearlySingularMatrix')
        elseif iframe==2
            warning('on','MATLAB:nearlySingularMatrix')
        end
    end
    
    if mask(iframe)
        %only update model parameters if mask is non-zero
        d = vvH(:,iframe); %[nSensors^2 x 1] vector

        autoCorrEst = lambda * autoCorrEst + CTconjC;
        crossCorrEst = lambda * crossCorrEst + C.' * conj(d); %C*conj(d);
        paramEst = autoCorrEst\crossCorrEst;
    
    end
    % the microphone PSD matrix can always be updated to reflect change in C
    Rv_est(:,iframe) =  C * conj(paramEst);
    fnm_est_hist(iframe,:) = paramEst.';
end
fprintf('\nDone!\n');

Rv_est = Rv_est.'; %convert back to [nFrames, nCovTerms]

model_params.estimator_name = 'EWLS SH+W coeffs';
model_params.smooth_factor = lambda;
model_params.fnm_est = paramEst;
model_params.fnm_est_hist = fnm_est_hist;
