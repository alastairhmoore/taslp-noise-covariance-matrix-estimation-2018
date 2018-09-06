exp_id = '20180326_x2';
sphHarmOrdPWDModel = 15;
sphHarmOrdDiffuseGt = 2;
sphHarmOrdDiffuseMax = sphHarmOrdDiffuseGt;
pose_measurement_noise_var_rad = deg2rad(1);
f_target = 2200;
ema = RigidSphere_20171120_01_horiz_up_20deg_circle4();

methods_to_run = [ 1  2  3 ... % gt, white, sph iso
                   5  6  7  18 8 ... % RS 0.9 0.95 0.99 0.995 0.999
                   15 16 17 9  4 10 ... % EWLS ord x 0.9 0.99 0.999 0.9999 0.99999 0.999999
                 ];


nTests = 1;
scr_20180326_03_orientation_seq

scr_20180326_02_precompute_variables

scr_id = sprintf('%s_gt_o%d_%d_Hz',...
        exp_id,sphHarmOrdDiffuseGt,...
        f_target);
    
scr_20180326_01_main_loop