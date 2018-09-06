function[x,xxH,posePerFrame,Rx_model_sh] = fcn_20180326_03_generate_acoustic_noise(fnm,posePerFrame,nFramesPerSegment,H,doaQuad)
% use the array manifold H for each planewave direction in doa_quad to
% determine the spatial covariance due to a unit plane wave
% then use pose to rotate fnm and evaluate at the planewave directions

nSHDiffuse = size(fnm,1);
if numel(fnm)~=nSHDiffuse, error('fnm should be a column vector'), end
sphHarmOrdDiffuse = sqrt(nSHDiffuse)-1;


% precompute the basis functions for inverse SFT
Y = sphBasis(doaQuad.az(:),doaQuad.inc(:),sphHarmOrdDiffuse);

% microphone covariance per doa
[nDOA,nSensors] = size(H);
nCovTerms = nSensors^2;
Rx_per_doa = reshape(...
    bsxfun(@times,H,conj(permute(H,[1 3 2]))),... % [nDOA,nSensors,nSensors]
    nDOA,nCovTerms);
quad_weighted_Rx_per_doa = bsxfun(@times,doaQuad.weights(:),Rx_per_doa);


% Use ElobesMicArray object to handle conversion of rotations from roll,
% pitch, yaw to Euler angles required by Wigner D
ema = RigidSphereHearingAidArray4; %N.B. Choice is arbitrary since we don't get impulse response from it

%% generate the random fluctuations to nominal poses

nFrames = size(posePerFrame.yaw,1);
nSegments = ceil(nFrames/nFramesPerSegment);

%% generate the random fluctuations to noise level from each direction
unit_variance_noise_per_doa = (1/sqrt(2)) * (randn(nFrames,nDOA) + 1i * randn(nFrames,nDOA)); %[nDOA,nFrames]
% weight according to sampling scheme
quad_weighted_noise_per_doa = unit_variance_noise_per_doa .* sqrt(doaQuad.weights.');

%assign empty variables
Rx_model_sh = zeros(nFrames,nCovTerms);
x = zeros(nFrames,nSensors);   %diffuse noise component

for iframe=1:nFrames
    if iframe==1 || ...
            posePerFrame.roll(iframe)~=posePerFrame.roll(iframe-1) || ...
            posePerFrame.pitch(iframe)~=posePerFrame.pitch(iframe-1) || ...
            posePerFrame.yaw(iframe)~=posePerFrame.yaw(iframe-1)
        
        % apply rotation to obtain euler angles
        ema.setPoseRollPitchYaw(posePerFrame.roll(iframe),...
            posePerFrame.pitch(iframe),...
            posePerFrame.yaw(iframe));
        [alpha,beta,gamma] = fcn_20180326_02_body_rotation_qr_to_wignerD_euler(ema.poseQuaternion);
        %for mic rotation apply the inverse rotation to sound field
        Dinv = sparse(fcn_20170616_02_myWignerDMatrix_Jacobi(sphHarmOrdDiffuse,...
            -gamma,-beta,-alpha));
        fnm_rot = Dinv * fnm; 
        
        pow_per_doa = Y * fnm_rot; %[nDoa,1]
        if max(imag(pow_per_doa))>3*eps
            error('non real doa power');
        end
        %pow_per_doa = doaQuad.weights .* real(pow_per_doa);
        pow_per_doa = real(pow_per_doa);
        doa_gain_factor = sqrt(pow_per_doa);
                
        % total Rx is weighted sum over DOAs       
        this_Rx = pow_per_doa.' * quad_weighted_Rx_per_doa; %[1,nCovTerms]
    end
    
    Rx_model_sh(iframe,:) = this_Rx;
    
    noise_sig_per_doa = quad_weighted_noise_per_doa(iframe,:) .* doa_gain_factor.';
    
    %[1 nSensors] = [1 nDOA]*[nDOA nSensors]
    x(iframe,:) = noise_sig_per_doa * H;
end

xxH = reshape(bsxfun(@times,x,conj(permute(x,[1 3 2]))),nFrames,nCovTerms);