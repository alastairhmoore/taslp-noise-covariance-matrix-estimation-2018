function[ fnm] = fcn_20171110_01_load_steered_ord_N_cardioids(sphHarmOrd,steering_directions_fcn)

% caridoid weights for order sphHarmOrd (only deg 0 so applies to real or
% complex SHs)
addpath('/Users/amoore1/third_party/polarch/Spherical-Array-Processing');
b_n = beamWeightsCardioid2Spherical(sphHarmOrd);
fnm_prototype = zeros((sphHarmOrd+1)^2,1);
ord_list = 0:sphHarmOrd;
idc_deg_0 = ord_list.^2 + ord_list + 1;
fnm_prototype(idc_deg_0) = b_n;

% directions in which to steer them
[az,inc,~] = steering_directions_fcn();
nFields = size(az,1);

% do SH rotation
% inc is the angle formed with z axis - beta rotation
% ax is the angle around z axis - gamma rotation
fnm = zeros(size(fnm_prototype,1),nFields);
for ifield = 1:nFields
    D = fcn_20170616_02_myWignerDMatrix_Jacobi(sphHarmOrd,...
        az(ifield),...
        inc(ifield),...
        0);
    fnm(:,ifield) = D * fnm_prototype;
end