classdef FreeField_20180202_01_ULA_configurable < FreeFieldArray
        properties (SetAccess=protected)
            % properties from ElobesMicArray
            sensorCartesianPositionsDefault
            refChan
    
            % properties from BinauralArray
            refChanLeft
            refChanRight
            channelsLeft
            channelsRight
        end
    methods
        function[obj] = FreeField_20180202_01_ULA_configurable(nSensors,sensorSpacing)
       
            % Use superclass to create the object
            obj = obj@FreeFieldArray();
            
            % Populate the parameters
            obj.sensorCartesianPositionsDefault = predefinedSensorPositions(nSensors,sensorSpacing);
            obj.refChan = 1;  % reference is the origin
        end
    end
    
end

function[sensor_pos] = predefinedSensorPositions(nSensors,sensorSpacing)
% evaluate to determine the postitions of the elements relative to
% the origin
% to number from left to right from listener perspective align with y-axis

%halfN = (nSensors-1)/2;
%if rem(halfN,1)~=0, error('nSensors must be odd'),end
x_pos = sensorSpacing * (0:-1:-(nSensors-1)).';
y_pos = zeros(nSensors,1);
z_pos = zeros(nSensors,1);

sensor_pos = [x_pos,y_pos,z_pos]; % [x,y,z] offsets of sensors
end