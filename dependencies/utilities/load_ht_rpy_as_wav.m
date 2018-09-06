function[varargout] = load_ht_rpy_as_wav(fn,rad_or_deg)
if ~ismember(nargout,[0 2 4])
    error('invalid calling syntax - incorrect number of outputs')
end
varargout = cell(nargout,1);

% want yaw from -180 to 180 to represent full scale
% here we use degrees, but elsewhere we might use radian, in whcih case ist
% just from -pi to pi
if nargin<2 || isempty(rad_or_deg)
    rad_or_deg = 'rad';
end
switch rad_or_deg
    case 'rad'
        nom = pi;
    case 'deg'
        nom = 180;
    otherwise
        error('rad_or_deg must be ''rad'' or ''deg''')
end

[in,fs] = audioread(fn);
if size(in,2)~=3
    error('Input file does not contain 3 channels - cannot be roll,pitch and yaw')
end

%  scale to rad or deg
in = in*nom;

switch nargout
    case 0
        figure;set(gcf,'name',fn);
        tscale = (0:size(in)-1).'./fs;
        sp(1)=subplot(311);
        plot(tscale,in(:,1));xlabel('Time [s]');ylabel(sprintf('Roll [%s]',rad_or_deg)); 
        sp(2)=subplot(312);
        plot(tscale,in(:,2));xlabel('Time [s]');ylabel(sprintf('Pitch [%s]',rad_or_deg));
        sp(3)=subplot(313);
        plot(tscale,in(:,3));xlabel('Time [s]');ylabel(sprintf('Yaw [%s]',rad_or_deg));
        setappdata(gcf,'lp',linkprop(sp,'xlim'));
    case 2
        varargout{1} = in;
        varargout{2} = fs;
    case 4
        for i=1:3
            varargout{i} = in(:,i);
        end
        varargout{4} = fs;
end        
