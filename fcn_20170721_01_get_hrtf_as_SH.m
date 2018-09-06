function[Hnm,fscale] = fcn_20170721_01_get_hrtf_as_SH(ema,sphHarmOrder,doConj,numZerosPad)

if nargin < 4 || isempty(numZerosPad)
    numZerosPad = 0;
end
if nargin < 3 || isempty(doConj)
    doConj = 0;
end

% precompute the SFT of the HRTF
nIncReg = sphHarmOrder+2;
nAzReg = 2*sphHarmOrder + 2;
[inc,az,weights]=sphrharm('cg',nIncReg,nAzReg);
[inc,az] = ndgrid(inc,az);
weights = weights*(2*pi/nAzReg);
weights = repmat(weights(:),1,nAzReg);
inc = inc(:);
az = az(:);
weights =weights(:);


hrir = ema.getImpulseResponseForSrc(az,inc);
hrir = permute(weights,[2 3 1]) .* hrir; %apply weighting in time domain for fewest multiplications
len_filt = size(hrir,1);
nfft = len_filt + numZerosPad;
hrtf = rfft(hrir,nfft,1); %[nFreq, nSensors, nDir]
if doConj
    hrtf = conj(hrtf); 
end

conjYlm = conj(sphBasis(az,inc,sphHarmOrder));

%sphrharm_conjYlm: [nDir,nSH]
% sum over directions, with nSH in 4th dimension
% conjHconjlm: [nFreq, nSensors, 1, nSH]
Hnm = sum(bsxfun(@times,hrtf,...
    permute(conjYlm, [3 4 1 2])),...
    3);
fscale = (0:size(hrtf,1)).' * ema.fs/nfft;
