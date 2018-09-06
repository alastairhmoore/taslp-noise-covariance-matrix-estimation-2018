exp_id = '20180326_x4';
sphHarmOrdPWDModel = 15;
sphHarmOrdDiffuseGt = [];                       % variable under test
sphHarmOrdDiffuseMax = 4;
pose_measurement_noise_var_rad = deg2rad(1);
f_target = [];                                  % variable under test
ema = [];                                       % variable under test

methods_to_run = [ 1  2  3 ... % gt, white, sph iso
                   6  7 18 8 ... % RS 0.95 0.99 0.995 0.999
                  11 12 13 14 ... % EWLS ord 1 2 3 4 0.99999
                 ];

f_target_seq = [2200];
sph_ord_gt_seq = [1 2 3 4];
mic_arrays = {@RigidSphere_20171120_01_horiz_up_20deg_circle4,...
    @RigidSphere_20171102_02_horiz_up_20deg_circle16};
mic_lab = {'circ4','circ16'};


nTests = 1;
scr_20180326_03_orientation_seq

for ifreqseq = length(f_target_seq):-1:1
    f_target = f_target_seq(ifreqseq);

    for imic = 1:length(mic_arrays)    
        ema = mic_arrays{imic}();
        
        for iordseq = 1:length(sph_ord_gt_seq)
            sphHarmOrdDiffuseGt = sph_ord_gt_seq(iordseq);
            
            scr_20180326_02_precompute_variables
            
            scr_id = sprintf('%s_gt_o%d_%d_Hz_%s_array',...
                exp_id,sphHarmOrdDiffuseGt,...
                f_target,mic_lab{imic});
            
            scr_20180326_01_main_loop
        end
    end
end