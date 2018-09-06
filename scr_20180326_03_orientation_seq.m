
nPoses = 5;
nFramesPerPose = 250; 
idcPosesToEval = 5; % only evaluate beamformer on last segment
maskOffset = 50; % number of iterations before the 'target' that source is assumed to be active
% define the array orientations
% - yaw is deterministic
yaw = linspace(0,pi/2,nPoses-1).'; %equally spaced over 1/4 circle
yaw(end+1) = 0; % last step is same as original pose
yaw_full_set = repmat(yaw,1,nTests);
% - pitch and roll add randomnsess
rng(927361,'twister');
pitch_full_set = deg2rad(20)*randn(nPoses,nTests);
pitch_full_set(end,:) = 0;
roll_full_set = deg2rad(10)*randn(nPoses,nTests);
roll_full_set(end,:) = 0;