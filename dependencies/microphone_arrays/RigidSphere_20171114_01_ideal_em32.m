classdef RigidSphere_20171114_01_ideal_em32 < RigidSphereArray
    properties (SetAccess=protected)
        % properties from ElobesMicArray already implemented in SphericalHarmonicSoundFieldArray
        
    end
    methods
        function[obj] = RigidSphere_20171114_01_ideal_em32()
            
            % use superclass to create the object with input parameter
            obj = obj@RigidSphereArray(0.042);
            
            % override the default properties of the superclass
            obj.sensorCartesianPositionsDefault = predefinedSensorPositions();
            obj.refChan = 0; % reference is the origin
            
        end
    end
    
end

function[sensor_pos] = predefinedSensorPositions()
% evaluate to determine the postitions of the elements relative to
% the origin
nMics = 32;
radius = 0.042;     % radius on which microphones lie [metres]


sensor_angles_deg = [69 0;...
    90 32;...
    111 0;...
    90 328;...
    32 0;...
    55 45;...
    90 69;...
    125 45;...
    148 0;...
    125 315;...
    90 291;...
    55 315;...
    21 91;...
    58 90;...
    121 90;...
    159 89;...
    69 180;...
    90 212;...
    111 180;...
    90 148;...
    32 180;...
    55 225;...
    90 249;...
    125 225;...
    148 180;...
    125 135;...
    90 111;...
    55 135;...
    21 269;...
    58 270;...
    122 270;...
    159 271];
inc = sensor_angles_deg(:,1)*pi/180;
az = sensor_angles_deg(:,2)*pi/180;

sensor_pos = radius * [cos(az).*sin(inc), sin(az).*sin(inc), cos(inc)]; % [x,y,z] offsets of sensors


%% Don't currently implement a sampling scheme that uses quadrature weights
%  put these here for future reference

% Quadrature weight for microphone at centre of pentagonal face
quad_pent = 5*pi/42;
% Quadrature weight for microphone at centre of hexagonal face
quad_hex = 9*pi/70;

quad = zeros(nMics, 1);
% All microphones are at centre of hexagonal face, except...
quad(1:nMics) = quad_hex;
% 12 microphones which are at centre of pentagonal face
quad([2 4 5 9 14 15 18 20 21 25 30 31]) = quad_pent;


end