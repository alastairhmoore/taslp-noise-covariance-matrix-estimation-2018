function[] = save_ht_rpy_as_wav(fn,roll,pitch,yaw,fs,rad_or_deg)

% want yaw from -180 to 180 to represent full scale
% here we use degrees, but elsewhere we might use radian, in whcih case ist
% just from -pi to pi
if nargin<6 || isempty(rad_or_deg)
    rad_or_deg = 'rad';
end
switch rad_or_deg
    case 'rad'
        denom = pi;
    case 'deg'
        denom = 180;
    otherwise
        error('rad_or_deg must be ''rad'' or ''deg''')
end

sz_r = size(roll);
sz_p = size(pitch);
sz_y = size(yaw);

if ~isequal(sz_r,sz_p,sz_y) || any([sz_r(2:end),sz_p(2:end),sz_y(2:end)]>1)
    error('roll pitch and yaw must be column vectors of the same length')
end

if ~strcmp(fn(end-4+(1:4)),'.wav')
    fn = [fn '.wav'];
end

audiowrite(fn,[roll,pitch,yaw]./denom,fs,'BitsPerSample',16);