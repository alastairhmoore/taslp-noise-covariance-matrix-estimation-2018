%% check for empty workspace
vars = who;
if ~isempty(vars)
    error('This script requires an empty workspace. Please clear your workspace before running this script.')
end

%% intialise everything
clear
clear classes
restoredefaultpath

%% add paths
addpath(genpath('dependencies'))

data_out_dir = 'saved_results';
check_output_dir_exists(data_out_dir);

%% universal constants
GLOBAL_SEED = 998212; %ensure tests are repeatable
fs = 10000;
sphHarmOrdPWDMax = 16; % use this to determine quad grid for the ground truth sound field
zero_pad_len = 0; % option to extend array manifold in time domain
noise_field_steering_fcn = @fcn_20160304_01_dodecahedron_vertices_sph;    
nFields = 20;
diffuse_to_sensor_noise_ratio_db = 20;
sensor_noise_intersensor_variance = 0.1; %variance in variances as fraction of unity
doScaleInvariant = 1;

%% one script per experiment
%%
scr_20180326_x1_effect_of_orientation_sequence
%%
scr_20180326_x2_convergence_vs_other_methods
%%
scr_20180326_x3_effect_of_PWD_order
%%
scr_20180326_x4_effect_of_est_ord_and_array



