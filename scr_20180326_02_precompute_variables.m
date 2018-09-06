% precompute variables that can be

% get directions for quadrature grid - directions from which incident plane
% waves arrive
[az,inc,weights] = fcn_20171102_01_quad_grid(sphHarmOrdPWDMax);
doa_quad.az = az(:);
doa_quad.inc = inc(:);
doa_quad.weights = weights(:);


% create an array
ema.prepareData(fs);

%% get steering vector
az_look = deg2rad(0);
inc_look = deg2rad(90);
ir0 = ema.getImpulseResponseForSrc(az_look,inc_look);
len_filt = size(ir0,1);

nfft = len_filt + zero_pad_len;
if nfft<size(ir0,1)
    error('the FFT is shorter than the data length')
end
hrtf = rfft(ir0,nfft,1); %[nFreq, nSensors, nDir=1]
fscale = (0:size(hrtf,1)).' * ema.fs/nfft;

% select the frequency of interest
[~,ifreq] = find_nearest(f_target,fscale);
H = permute(hrtf(ifreq,:,:),[2 1 3]); %[nSensors, nFreq=1, nDir=1]

% assign array manifold to steering vector 
d = squeeze(H); %[nSensors, 1];

%% get array manifold for all directions on quadrature grid
hrir = ema.getImpulseResponseForSrc(doa_quad.az,doa_quad.inc);
hrtf = rfft(hrir,nfft,1); %[nFreq, nSensors, nDir]
H = permute(hrtf(ifreq,:,:),[3 2 1]);

%get HRTF in SH domain
fprintf('Calculating HRTF in SH domain...');
[tildeHnm, fscale_check] = fcn_20170721_01_get_hrtf_as_SH(ema,sphHarmOrdPWDModel,1,zero_pad_len); %[nFreq,nSensors,1,nSH]
if ~isequal(fscale,fscale_check)
    error('Need same frequency bins in space domain HRTF and SH domain HRTF')
end
[nFreq,nSensors,~,~] = size(tildeHnm);
tildeHnm = squeeze(tildeHnm(ifreq,:,1,:)).'; %[nSHPWD, nSensors]
fprintf('Done!\n')

%get Gaunt coefficients
fprintf('Calculating Gaunt coefficients...');
% [nSHPWD nSHPWD nSHDiffuse]
Gmat = fcn_20170626_01_triple_sph_harm_integral(sphHarmOrdDiffuseMax,sphHarmOrdPWDModel,sphHarmOrdPWDModel);
%form 2 - keep SH components separate and sparse
nSHDiffuseMax = (sphHarmOrdDiffuseMax+1)^2;
Gcell = cell(nSHDiffuseMax,1);
for ish = 1:nSHDiffuseMax
    Gcell{ish} = sparse(squeeze(Gmat(ish,:,:)));
end
fprintf('Done!\n')


