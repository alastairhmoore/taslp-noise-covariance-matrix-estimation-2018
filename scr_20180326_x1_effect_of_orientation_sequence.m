exp_id = '20180326_x1';
sphHarmOrdPWDModel = 15;
sphHarmOrdDiffuseGt = 2;
sphHarmOrdDiffuseMax = sphHarmOrdDiffuseGt;
pose_measurement_noise_var_rad = deg2rad(1);
f_target = 2200;
ema = RigidSphere_20171120_01_horiz_up_20deg_circle4();

methods_to_run = [ 1 ... % gt
                   2 ... % white
                   3 ... % sph iso
                   4];   % EWLS ord x0.99999
               
nTests = 1;
nPoses = 50;
nFramesPerPose = 100;
idcPosesToEval = []; % only evaluate beamformer on last segment
maskOffset = 0; % number of iterations before the 'target' that source is assumed to be active (irrelevant here because we do not mask anywhere)

scr_20180326_02_precompute_variables

for iposeseq = 1:3
    switch iposeseq
        case 1
            yaw = linspace(0,6*pi,nPoses).'; %equally spaced over whole circle
            yaw_full_set = repmat(yaw,1,nTests);
            pitch_full_set = zeros(size(yaw_full_set));
            roll_full_set = zeros(size(yaw_full_set));
        case 2
            yaw = linspace(0,6*pi,nPoses).'; %equally spaced over whole circle
            yaw_full_set = repmat(yaw,1,nTests);
            rng(2355,'twister');
            pitch_full_set = deg2rad(20)*randn(size(yaw_full_set));
            roll_full_set = deg2rad(10)*randn(size(yaw_full_set));
        case 3
            rng(52116,'twister');
            yaw_full_set = -pi+2*pi*rand(nPoses,nTests); % everything is random
            pitch_full_set = -pi+2*pi*rand(nPoses,nTests);
            roll_full_set =  -pi+2*pi*rand(nPoses,nTests);
    end
    
    scr_id = sprintf('%s_gt_o%d_%d_Hz_poseseq_%d',...
        exp_id,sphHarmOrdDiffuseGt,...
        f_target,iposeseq);
    
    scr_20180326_01_main_loop
end