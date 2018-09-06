function[R] = fcn_20170421_01_diffuse_pow_to_psd_matrix(ema,fs,nfft,az,inc,pow,quad)
% given the power arriving from each specified direction (and an optional
% quadrature scalar for each direction) simulate the microphone signals 
% which would be experienced in a diffuse field and from these compute an
% mvdr beamformer


% ema: handle of instantiated ElobesMicArray object
az = az(:);
inc = inc(:);
pow = pow(:);
if nargin>6 && ~isempty(quad)
    quad = quad(:);
else
    quad = ones(size(az));
end

duration = 0.3;
nSigSamples = round(duration*fs);
nDoa = length(az);


ema.prepareData(fs);
ir = ema.getImpulseResponseForSrc(az(:),inc(:)); %[len_filt,nSensors,nDirections]
per_doa_weight = permute(sqrt(pow).*quad,[2 1]); %[1 nDoa]
% have to loop becuase fftfilt cant do 3dimensional array
per_doa_sig = bsxfun(@times,randn(nSigSamples,nDoa),per_doa_weight); %[nSigSamples,nDoa]
mic_sigs = zeros(nSigSamples,ema.nSensors);
ir = permute(ir,[1 3 2]);
for ichan = 1:ema.nSensors
    mic_sigs(:,ichan) = sum(fftfilt(ir(:,:,ichan),per_doa_sig),2);
end
    
% stft
nwin = nfft;
ninc = round(0.5*nwin);
win{1} = sqrt(hamming(nwin,'periodic'));     % wola window for overlap factor of 2
win{1} = win{1}/sqrt(sum(win{1}(1:ninc:nwin).^2));      % normalize window
win{2} = win{1};

padding = zeros(nwin-ninc,ema.nSensors);
mic_sigs = [padding;mic_sigs;padding];
X = stft_v2('fwd',mic_sigs,win,ninc,nfft,fs);

% make PSD matrix
X = permute(X,[2 4 1 3]); % [nSensors, 1, nFreq, nFrames];
R = mean(bsxfun(@times,X,conj(permute(X,[2 1 3 4]))),4); % [nSenors,nSensors,nFreq]
