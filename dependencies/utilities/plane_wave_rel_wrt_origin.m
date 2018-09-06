function[nc_ir, t0] = plane_wave_rel_wrt_origin(mic_xyz,doas_azinc,fs,c,nSamples)

halfMinSincSamples = 30; %how many samples before and after the middle of the sinc that are needed

if nargin<4 || isempty(c)
    c = soundspeed;
end

% get unit vectors pointing *towards* doas
% these are normal to plane of the propagating wave
nDOAs = size(doas_azinc,1);
[doa_x,doa_y,doa_z]  = mysph2cart(doas_azinc(:,1),doas_azinc(:,2),ones(nDOAs,1));

% Distance to plane: D = \hat{n} \cdot \hat{x}_0 + p
% (http://mathworld.wolfram.com/Point-PlaneDistance.html)
% p is 0 for plane intersection origin
% x_0 is our point(s) of interest
% cdot indicates dot product, or simply inner product

%check orientation of mic_xyz data
[nMics,nAx,nDummy] = size(mic_xyz);
if nAx~=3 || nDummy~=1
    error('mic_xyz must be an Nx3 matrix')
end

%check nSamples is odd, if specified
if nargin > 4 && ~isempty(nSamples)
    if rem(nSamples,2)~=1
        error('nSamples must be an odd number')
    end
end

% grid of distances
d = sum(bsxfun(@times,permute(mic_xyz,[2 1 3]),... %[3,nMics,1]
                      permute([doa_x,doa_y,doa_z],[2 3 1])),... %[3,1,nDOAs]
        1); %[1,nMics, nDOAs] 

% determine required impulse response length from dimensions of array
tFlightMin = max(abs(d(:))) / c;
%halfNSamplesMin = tFlightMin * fs + halfMinSincSamples;
%nReqSamples = 2*halfNSamplesMin + 1; %ensure it is odd to have delta at origin
halfNSamplesMin = max(ceil(tFlightMin * fs),halfMinSincSamples);
%make sure there delta is at least 1/4 into impulse response
nReqSamples = 4*halfNSamplesMin + 1; %ensure it is odd to have delta at origin 

if nargin<5 || isempty(nSamples)
    nSamples = nReqSamples;
else
    if nSamples < nReqSamples
        error('Requested nSamples is too low')
    end
end 
halfNSamples = (nSamples-1)/2;
    
% calculate in the frequency domain
N_fft = nSamples; %assignment for clarity

f = (0:((N_fft+2)/2 - 1)).' * fs/N_fft; %[nf 1] where nf is number of bins up to nyquist

mag = ones(length(f),nMics,nDOAs);
phase = 2*pi/c * bsxfun(@times, d, f);

t0 = halfNSamples+1; %this is the sample index of true zero time
nc_ir = circshift(ifft(mag .* exp(1i*phase),N_fft,'symmetric'),[t0-1 0 0]); %avoid wraparound - now non-causal
