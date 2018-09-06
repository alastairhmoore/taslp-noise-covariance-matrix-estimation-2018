% function handles for estimation methods
sph_iso_fcn = @fcn_20170822_04_sph_iso_model;
recursive_smooth_fcn = @fcn_20170905_02_est_Rx_recursive_with_mask;
ewls_fcn = @fcn_20180327_01_ewls_ncm_est;

% get random seeds
rng(GLOBAL_SEED,'twister');
acoustic_seeds = randi(2^32-1,nFields,nTests);
sensor_seeds = randi(2^32-1,nFields,nTests);
pose_err_seeds = randi(2^32-1,nFields,nTests);

% get SH coefficients of ground truth noise power distributions
fnm = fcn_20171110_01_load_steered_ord_N_cardioids(sphHarmOrdDiffuseGt,noise_field_steering_fcn);
if nFields>size(fnm,2)
    error('More fields required than supplied');
end
fnm = fnm(:,1:nFields);


% predefine output fields
saved_gt.xxH = cell(nFields,nTests);
saved_gt.Rx_diffuse =  cell(nFields,nTests);
saved_gt.Rx_total =  cell(nFields,nTests);
saved_gt.bf_noise_ref = cell(nFields,nTests);
saved_gt.fnm = cell(nFields,nTests);
saved_gt.sensor_var = cell(nFields,nTests);
saved_gt.per_frame_pose = cell(nFields,nTests);
saved_gt.speech_absence_mask = cell(nFields,nTests);

% dereived parameters
nMethods = length(methods_to_run);
nIterationsToMask = nFramesPerPose + maskOffset; %all of the last pose + the lead in


%% loop over the fields
for ifield = nFields:-1:1
    fprintf('\n\nProcessing field %d/%d\n\n', nFields-ifield+1, nFields);
        
    for itest = nTests:-1:1
        rot_var = [];
        nFramesPerSegment = 10;
        
        pose.yaw = reshape(repmat(yaw_full_set(:,itest).',nFramesPerPose,1),[],1);
        pose.pitch = reshape(repmat(pitch_full_set(:,itest).',nFramesPerPose,1),[],1);
        pose.roll = reshape(repmat(roll_full_set(:,itest).',nFramesPerPose,1),[],1);
        
        %% make the data
        this_fnm = fnm(:,ifield);
        rng(acoustic_seeds(ifield,itest),'twister')
        [x,xxH,per_frame_pose,Rx_model_sh] = fcn_20180326_03_generate_acoustic_noise(this_fnm,pose,nFramesPerSegment,H,doa_quad);
        
        % realise the sensor noise with nominally unity variance
        nFrames = size(x,1);
        rng(sensor_seeds(ifield,itest),'twister')
        [v,vvH,this_sensor_noise_var] = fcn_20171122_04_generate_sensor_noise(nFrames,nSensors,sensor_noise_intersensor_variance);
        
        % get the average power in diffuse and white signals
        mean_x_pow_db = 10*log10(mean(mean(abs(x).^2,1),2));
        mean_y_pow_db = 10*log10(mean(mean(abs(v).^2,1),2));
        y_scale_db = (mean_x_pow_db-mean_y_pow_db) - diffuse_to_sensor_noise_ratio_db;
        signal_gain = 10^(y_scale_db/20);
        pow_gain  = 10^(y_scale_db/10);
        
        % add it in
        x = x + signal_gain * v;
        xxH = xxH + pow_gain * vvH;
        
        % add ground truth spatially white variance as flattened diagonal
        % to each frame
        this_sensor_noise_var = this_sensor_noise_var * pow_gain;
        Rx_gt = bsxfun(@plus,Rx_model_sh,reshape(diag(this_sensor_noise_var),1,[]));
        
        
        % make a mask where estimates are not updated - this is where we
        % evaluate the noise reduction performance
        mask = ones(nFrames,1);
        for ii=1:length(idcPosesToEval)
            idc = idcPosesToEval(ii);
            mask(nFramesPerPose*(idc-1) + [1:nIterationsToMask]-maskOffset) = 0;
        end
        
        
        % create noise to add to pose
        rng(pose_err_seeds(ifield,itest),'twister');
        noisy_per_frame_pose.yaw = per_frame_pose.yaw + pose_measurement_noise_var_rad * randn(size(per_frame_pose.yaw));
        noisy_per_frame_pose.pitch = per_frame_pose.pitch + pose_measurement_noise_var_rad * randn(size(per_frame_pose.pitch));
        noisy_per_frame_pose.roll = per_frame_pose.roll + pose_measurement_noise_var_rad * randn(size(per_frame_pose.roll));
        
        % generate target signal - relative to head
        % - pick the allowed DOAs
        %         [doa_vec(:,1),doa_vec(:,2),doa_vec(:,3)] = mysph2cart(doa_quad.az,doa_quad.inc,ones(size(doa_quad.az)));
        %        [O_vec(:,1),O_vec(:,2),O_vec(:,3)] = mysph2cart(0,pi/2,1);
        %        angle_from_O = distcos(doa_vec,O_vec)*pi;
        %        i_near = find(angle_from_O<deg2rad(10));
        %x_sig = zeros(size(x));
        %x_sig(mask==0,:) = fcn_20170906_02_generate_unit_variance_target_source_random_doa(...
        %    sig_seeds(ifield),sum(mask==0),H(i_near,:));
        
        % steering vector is always dead ahead with respect to the
        % (rotated) array
        %        [min_val,i_O] = min(angle_from_O);
        %        d = H(i_O,:).';
        
        % target signal is the beamformed clean signal using
        %        w_ref = d/(d' * d);
        %bf_sig_ref = x_sig * conj(w_ref); %this orientation is avoids lots of permutations
        
        % only need to measure the power in the beamformed noise
        %        bf_noise_ref = x * conj(w_ref);
        
        % store the ground truth to be saved
        saved_xxH{ifield,itest} = xxH;
        saved_gt.xxH{ifield,itest} = xxH;
        saved_gt.Rx_diffuse{ifield,itest} = Rx_model_sh;
        saved_gt.Rx_total{ifield,itest} = Rx_gt;
        %        saved_gt.bf_noise_ref{ifield,itest} = bf_noise_ref;
        saved_gt.fnm{ifield,itest} = this_fnm;
        saved_gt.sensor_var{ifield,itest} = this_sensor_noise_var;
        saved_gt.per_frame_pose{ifield,itest} = per_frame_pose;
        saved_gt.speech_absence_mask = mask;
        
        
        %% loop over methods
        for ii = nMethods:-1:1
            imethod = methods_to_run(ii);
            
            %% do the estimation - some methods below are the same but with different parameter values
            switch imethod
                case 1
                    est(ifield,itest,ii).legend_lab = sprintf('gt');
                    est(ifield,itest,ii).Rx_est = Rx_gt;
                case 2
                    est(ifield,itest,ii).legend_lab = sprintf('white');
                    Rx_one_frame = reshape(eye(nSensors),1,nSensors^2);
                    est(ifield,itest,ii).Rx_est = repmat(Rx_one_frame,nFrames,1);
                case 3
                    est(ifield,itest,ii).legend_lab = 'sph iso';
                    fprintf('Doing %s...\n',est(ifield,itest,ii).legend_lab);
                    [Rx_est,model_params] = sph_iso_fcn(xxH,[],[],tildeHnm,Gcell(1),[]);
                    est(ifield,itest,ii).Rx_est = Rx_est;
                    est(ifield,itest,ii).model_params = model_params;
                case 4
                    sphHarmOrdDiffuseEst = sphHarmOrdDiffuseGt;
                    nSHDiffuseEst = (sphHarmOrdDiffuseEst+1)^2;
                    smooth = 0.99999;
                    pm.sensor_noise_model = 1;
                    est(ifield,itest,ii).legend_lab = sprintf('ewls w%d - %2.6f',...
                        pm.sensor_noise_model,...
                        smooth);
                    fprintf('Doing %s...\n',est(ifield,itest,ii).legend_lab);
                    [Rx_est,model_params] = ewls_fcn(xxH,smooth,[],tildeHnm,Gcell(1:nSHDiffuseEst),noisy_per_frame_pose,mask,pm);
                    est(ifield,itest,ii).Rx_est = Rx_est;
                    est(ifield,itest,ii).model_params = model_params;
                case 5
                    smooth = 0.9;
                    est(ifield,itest,ii).legend_lab = sprintf('recursive smooth - %2.5f',smooth);
                    fprintf('Doing %s...\n',est(ifield,itest,ii).legend_lab);
                    [Rx_est,model_params] = recursive_smooth_fcn(xxH,smooth,[],mask);
                    est(ifield,itest,ii).Rx_est = Rx_est;
                    est(ifield,itest,ii).model_params = model_params;
                case 6
                    smooth = 0.95;
                    est(ifield,itest,ii).legend_lab = sprintf('recursive smooth - %2.5f',smooth);
                    fprintf('Doing %s...\n',est(ifield,itest,ii).legend_lab);
                    [Rx_est,model_params] = recursive_smooth_fcn(xxH,smooth,[],mask);
                    est(ifield,itest,ii).Rx_est = Rx_est;
                    est(ifield,itest,ii).model_params = model_params;
                case 7
                    smooth = 0.99;
                    est(ifield,itest,ii).legend_lab = sprintf('recursive smooth - %2.5f',smooth);
                    fprintf('Doing %s...\n',est(ifield,itest,ii).legend_lab);
                    [Rx_est,model_params] = recursive_smooth_fcn(xxH,smooth,[],mask);
                    est(ifield,itest,ii).Rx_est = Rx_est;
                    est(ifield,itest,ii).model_params = model_params;
                case 8
                    smooth = 0.999;
                    est(ifield,itest,ii).legend_lab = sprintf('recursive smooth - %2.5f',smooth);
                    fprintf('Doing %s...\n',est(ifield,itest,ii).legend_lab);
                    [Rx_est,model_params] = recursive_smooth_fcn(xxH,smooth,[],mask);
                    est(ifield,itest,ii).Rx_est = Rx_est;
                    est(ifield,itest,ii).model_params = model_params;
                case 9
                    sphHarmOrdDiffuseEst = sphHarmOrdDiffuseGt;
                    nSHDiffuseEst = (sphHarmOrdDiffuseEst+1)^2;
                    smooth = 0.9999;
                    pm.sensor_noise_model = 1;
                    est(ifield,itest,ii).legend_lab = sprintf('ewls w%d - %2.6f',...
                        pm.sensor_noise_model,...
                        smooth);
                    fprintf('Doing %s...\n',est(ifield,itest,ii).legend_lab);
                    [Rx_est,model_params] = ewls_fcn(xxH,smooth,[],tildeHnm,Gcell(1:nSHDiffuseEst),noisy_per_frame_pose,mask,pm);
                    est(ifield,itest,ii).Rx_est = Rx_est;
                    est(ifield,itest,ii).model_params = model_params;
                case 10
                    sphHarmOrdDiffuseEst = sphHarmOrdDiffuseGt;
                    nSHDiffuseEst = (sphHarmOrdDiffuseEst+1)^2;
                    smooth = 0.999999;
                    pm.sensor_noise_model = 1;
                    est(ifield,itest,ii).legend_lab = sprintf('ewls w%d - %2.6f',...
                        pm.sensor_noise_model,...
                        smooth);
                    fprintf('Doing %s...\n',est(ifield,itest,ii).legend_lab);
                    [Rx_est,model_params] = ewls_fcn(xxH,smooth,[],tildeHnm,Gcell(1:nSHDiffuseEst),noisy_per_frame_pose,mask,pm);
                    est(ifield,itest,ii).Rx_est = Rx_est;
                    est(ifield,itest,ii).model_params = model_params;
                case 11
                    sphHarmOrdDiffuseEst = 1;
                    nSHDiffuseEst = (sphHarmOrdDiffuseEst+1)^2;
                    smooth = 0.99999;
                    pm.sensor_noise_model = 1;
                    est(ifield,itest,ii).legend_lab = sprintf('ewls w%d - %2.6f',...
                        pm.sensor_noise_model,...
                        smooth);
                    fprintf('Doing %s...\n',est(ifield,itest,ii).legend_lab);
                    [Rx_est,model_params] = ewls_fcn(xxH,smooth,[],tildeHnm,Gcell(1:nSHDiffuseEst),noisy_per_frame_pose,mask,pm);
                    est(ifield,itest,ii).Rx_est = Rx_est;
                    est(ifield,itest,ii).model_params = model_params;
                case 12
                    sphHarmOrdDiffuseEst = 2;
                    nSHDiffuseEst = (sphHarmOrdDiffuseEst+1)^2;
                    smooth = 0.99999;
                    pm.sensor_noise_model = 1;
                    est(ifield,itest,ii).legend_lab = sprintf('ewls w%d - %2.6f',...
                        pm.sensor_noise_model,...
                        smooth);
                    fprintf('Doing %s...\n',est(ifield,itest,ii).legend_lab);
                    [Rx_est,model_params] = ewls_fcn(xxH,smooth,[],tildeHnm,Gcell(1:nSHDiffuseEst),noisy_per_frame_pose,mask,pm);
                    est(ifield,itest,ii).Rx_est = Rx_est;
                    est(ifield,itest,ii).model_params = model_params;
                case 13
                    sphHarmOrdDiffuseEst = 3;
                    nSHDiffuseEst = (sphHarmOrdDiffuseEst+1)^2;
                    smooth = 0.99999;
                    pm.sensor_noise_model = 1;
                    est(ifield,itest,ii).legend_lab = sprintf('ewls w%d - %2.6f',...
                        pm.sensor_noise_model,...
                        smooth);
                    fprintf('Doing %s...\n',est(ifield,itest,ii).legend_lab);
                    [Rx_est,model_params] = ewls_fcn(xxH,smooth,[],tildeHnm,Gcell(1:nSHDiffuseEst),noisy_per_frame_pose,mask,pm);
                    est(ifield,itest,ii).Rx_est = Rx_est;
                    est(ifield,itest,ii).model_params = model_params;
                case 14
                    sphHarmOrdDiffuseEst = 4;
                    nSHDiffuseEst = (sphHarmOrdDiffuseEst+1)^2;
                    smooth = 0.99999;
                    pm.sensor_noise_model = 1;
                    est(ifield,itest,ii).legend_lab = sprintf('ewls w%d - %2.6f',...
                        pm.sensor_noise_model,...
                        smooth);
                    fprintf('Doing %s...\n',est(ifield,itest,ii).legend_lab);
                    [Rx_est,model_params] = ewls_fcn(xxH,smooth,[],tildeHnm,Gcell(1:nSHDiffuseEst),noisy_per_frame_pose,mask,pm);
                    est(ifield,itest,ii).Rx_est = Rx_est;
                    est(ifield,itest,ii).model_params = model_params;
                case 15
                    sphHarmOrdDiffuseEst = sphHarmOrdDiffuseGt;
                    nSHDiffuseEst = (sphHarmOrdDiffuseEst+1)^2;
                    smooth = 0.9;
                    pm.sensor_noise_model = 1;
                    est(ifield,itest,ii).legend_lab = sprintf('ewls w%d - %2.6f',...
                        pm.sensor_noise_model,...
                        smooth);
                    fprintf('Doing %s...\n',est(ifield,itest,ii).legend_lab);
                    [Rx_est,model_params] = ewls_fcn(xxH,smooth,[],tildeHnm,Gcell(1:nSHDiffuseEst),noisy_per_frame_pose,mask,pm);
                    est(ifield,itest,ii).Rx_est = Rx_est;
                    est(ifield,itest,ii).model_params = model_params;
                case 16
                    sphHarmOrdDiffuseEst = sphHarmOrdDiffuseGt;
                    nSHDiffuseEst = (sphHarmOrdDiffuseEst+1)^2;
                    smooth = 0.99;
                    pm.sensor_noise_model = 1;
                    est(ifield,itest,ii).legend_lab = sprintf('ewls w%d - %2.6f',...
                        pm.sensor_noise_model,...
                        smooth);
                    fprintf('Doing %s...\n',est(ifield,itest,ii).legend_lab);
                    [Rx_est,model_params] = ewls_fcn(xxH,smooth,[],tildeHnm,Gcell(1:nSHDiffuseEst),noisy_per_frame_pose,mask,pm);
                    est(ifield,itest,ii).Rx_est = Rx_est;
                    est(ifield,itest,ii).model_params = model_params;
                case 17
                    sphHarmOrdDiffuseEst = sphHarmOrdDiffuseGt;
                    nSHDiffuseEst = (sphHarmOrdDiffuseEst+1)^2;
                    smooth = 0.999;
                    pm.sensor_noise_model = 1;
                    est(ifield,itest,ii).legend_lab = sprintf('ewls w%d - %2.6f',...
                        pm.sensor_noise_model,...
                        smooth);
                    fprintf('Doing %s...\n',est(ifield,itest,ii).legend_lab);
                    [Rx_est,model_params] = ewls_fcn(xxH,smooth,[],tildeHnm,Gcell(1:nSHDiffuseEst),noisy_per_frame_pose,mask,pm);
                    est(ifield,itest,ii).Rx_est = Rx_est;
                    est(ifield,itest,ii).model_params = model_params;
                case 18
                    smooth = 0.995;
                    est(ifield,itest,ii).legend_lab = sprintf('recursive smooth - %2.5f',smooth);
                    fprintf('Doing %s...\n',est(ifield,itest,ii).legend_lab);
                    [Rx_est,model_params] = recursive_smooth_fcn(xxH,smooth,[],mask);
                    est(ifield,itest,ii).Rx_est = Rx_est;
                    est(ifield,itest,ii).model_params = model_params;                otherwise
                    error('Not specified!');
            end
            
            %% get error in noise estimation
            est(ifield,itest,ii).est_cov_mean_sq_err = fcn_20170803_01_covariance_matrix_error_metric(est(ifield,itest,ii).Rx_est,Rx_gt,doScaleInvariant);
            
            % get beamformer and resulting bf noise
            [est(ifield,itest,ii).bf_noise, est(ifield,itest,ii).w] = ...
                fcn_20170906_03_do_mvdr_bf(x,d,est(ifield,itest,ii).Rx_est);
            
        end
    end
end

save(sprintf('%s/dat_%s.mat',data_out_dir,scr_id),...
    '-v7.3','est','saved_gt')

%% compute summary data
mse_mat = zeros(nFields,nTests,nMethods);
bf_mat = zeros(nFields,nTests,nMethods);
nr_mat = zeros(nFields,nTests,nMethods);
legstr = {};
for ifield = 1:nFields
    for itest = 1:nTests
        sum_x_pow = 0;
        for isensor = 1:nSensors
            icov = (isensor-1)*nSensors + isensor;
            sum_x_pow = sum_x_pow + mean(saved_gt.xxH{ifield,itest}(mask==0,icov),1,'omitnan');
        end
        mean_x_pow = sum_x_pow/nSensors;
        
        for imethod = 1:nMethods
            
            mse_mat(ifield,itest,imethod) = 10*log10(mean(est(ifield,itest,imethod).est_cov_mean_sq_err(mask==0),1,'omitnan'));
            bf_mat(ifield,itest,imethod) = 10*log10(mean(abs(est(ifield,itest,imethod).bf_noise(mask==0)).^2,1,'omitnan'));
            nr_mat(ifield,itest,imethod) = 10*log10(mean(abs(est(ifield,itest,imethod).bf_noise(mask==0)).^2,1,'omitnan')/mean_x_pow);
            
        end
        if ifield==1 && itest==1
            legstr{imethod} = est(ifield,itest,imethod).legend_lab;
        end
    end
end

check_output_dir_exists(data_out_dir);
save(sprintf('%s/dat_%s_summary.mat',data_out_dir,scr_id),...
    'mse_mat','bf_mat','nr_mat','sphHarmOrdDiffuseGt','sphHarmOrdDiffuseEst',...
    'legstr');
