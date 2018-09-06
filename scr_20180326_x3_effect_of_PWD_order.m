exp_id = '20180326_x3';
sphHarmOrdPWDModel = [];                        % variable under test
sphHarmOrdDiffuseGt = 2;
sphHarmOrdDiffuseMax = sphHarmOrdDiffuseGt;
pose_measurement_noise_var_rad = deg2rad(1);
f_target = [];                                  % variable under test
ema = RigidSphere_20171120_01_horiz_up_20deg_circle4();

methods_to_run = [ 1 ... % gt
                   4 ... % EWLS ord x 0.99999
                 ];  

f_target_seq = [100 220 470 1000 2200 4700] ;
pwd_order_seq{1} = [1 2 3 4 14];
pwd_order_seq{2} = [1 2 3 4 14];
pwd_order_seq{3} = [1 2 3 4 14];
pwd_order_seq{4} = [1 2 3 4 5 6 7 8 14];
pwd_order_seq{5} = [1 2 3 4 5 6 7 8 10 12 14];
pwd_order_seq{6} = [1 2 3 4 5 6 7 8 10 12 14];

nTests = 1;
scr_20180326_03_orientation_seq

for ifreqseq = length(f_target_seq):-1:1
    f_target = f_target_seq(ifreqseq);
    for iordseq = 1:length(pwd_order_seq{ifreqseq})
        sphHarmOrdPWDModel = pwd_order_seq{ifreqseq}(iordseq);
        scr_20180326_02_precompute_variables

        scr_id = sprintf('%s_gt_o%d_%d_Hz_pwd_ord_%d',...
            exp_id,sphHarmOrdDiffuseGt,...
            f_target,sphHarmOrdPWDModel);
        scr_20180326_01_main_loop
    end
end