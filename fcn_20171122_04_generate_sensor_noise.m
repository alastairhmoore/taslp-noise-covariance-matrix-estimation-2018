function[y,yyH,sensor_noise_var] = fcn_20171122_04_generate_sensor_noise(nFrames,nSensors,intersensor_variance_frac)

nCovTerms = nSensors^2;

%get the actual values for variance in each microphone
sensor_noise_var = 1 + (intersensor_variance_frac * randn(1,nSensors));

% convert these to the gain which must be applied to zero mean complex gaussian
complex_sensor_noise_gain_factor = (1/sqrt(2)) * sqrt(sensor_noise_var);

% sensor noise - should not stop xxH from being hermitian
y = complex_sensor_noise_gain_factor .* (randn(nFrames,nSensors) + 1i * randn(nFrames,nSensors));

% get the outer product as well just for convenience
yyH = reshape(bsxfun(@times,y,conj(permute(y,[1 3 2]))),nFrames,nCovTerms);