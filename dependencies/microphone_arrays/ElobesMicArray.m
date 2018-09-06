%ElobesMicArray creates an abstract class which defines the essential methods
%of a microphone arrat
%Idea is to provide a default constructor that allows us to wrap different
%types of data.
classdef ElobesMicArray < handle
    properties (SetAccess=protected)
        %supportsSpherical               = 0; % [0]/1 - concept of obtaining SH representation of received signals
        %supportsBinaural                = 0; % [0]/1 - concept of left/right channels)
        supportsRotation                = 0; % [0]/1 - can sensible data be obtained if array is rotated from its default position
        fs                              = []; % sample rate in Hertz
        availableInterpolationMethods   = {}; % defines the available interpolation methods
        poseAxes                        = []; % [x axis vec, y axis vec, z axis vec]
        poseQuaternion                  = []; % quartenion represening the rotation required to get from the default orientation to the current one
        sensorCartesianPositions        = []; % M x 3 vector of cartesian coordinates
        arrayGeometry                   = struct('vertices',[],'faces',[]); % data for plotting the shape of the physical array
        t0                              = []; % index in samples of time zero (i.e. arrival time of direct path signal) at the origin
        nfft                            = []; % size of fft to use for internal time/frequency transformation
        H0                              = []; % [fix(nfft/2)+1 1] vector specifying the freqeuncy domain response at the origin - sets bandwidth and bulk delay
        f                               = []; % [fix(nfft/2)+1 1] vector specifying the bin frequencies of the data in H0.    
       
        % -- fields for analytical arrays --
        c                   = soundspeed; % speed of sound [Hz]
        
    end
    properties (Dependent)
        nSensors                               % depends on sensorCartesianPositionsDefaults
        hasDefinedSensorPositions              % depends on sensorCartesianPositionsDefaults
    end
    
    % properties which users can set directly
    % TODO: Add validation to these so they can't be set to incorrect
    % values
    properties
        interpolationMethod             = ''; % defines which type of interpolation to be used in getImpulseResponseForSrc method
        doFastPlot                      = 1;  % 0/[1] - turns off fancy graphics, lighting, etc in the default plot method to improve speed
    end
 
    properties (Abstract, SetAccess=protected)
        % subclass must define these 
        sensorCartesianPositionsDefault % M x 3 vector of cartesian coordinates before any rotations
        refChan                         % channel to treat as the reference        
    end
    
    methods (Abstract)
        prepareData(obj,req_fs,varargin)
        % runs any preprocessing required to ready the object beofore getting the data
        % minimally the sample rate req_fs is needed
        % use additional inputs in key/value pairs to pass more parameters
        % which are specific to a particular subclass
        
        [minRadius, maxRadius] = getValidSrcRadiusRange(obj)
        % provides bounds on the acceptable range of source distances that
        % can be supported by a particular array/dataset
        
        [ir,rel_src_az,rel_src_inc] = getFrequencyResponseForSrc(obj,src_az,src_inc,varargin)
        % returns the multichannel frequency response due to the defined
        % source(s) positions relative to the current orientation. The relative
        % directions are also returned [nSamples,nChans,nSrc]
        % use additional inputs in key/value pairs to pass more parameters
        % which are specific to a particular subclass
        % e.g. defining the source distance, directivity, orientation
        % all implementations should work without any additional inputs
       
    end
    
    
    methods
        function[obj] = ElobesMicArray()
            % not sure there is anything to do here
            % would be nice to use default pose to intialise the value of
            % sensorCartesianPositions but we need the value of
            % sensorCartesianPositionsDefault to be set first
            obj.sensorCartesianPositions = obj.sensorCartesianPositionsDefault;
            obj.poseAxes = eye(3);
            obj.poseQuaternion = roteu2qr('xyz',[0,0,0]);
        end         
        function[] = setPoseRollPitchYawDegrees(obj,roll_deg,pitch_deg,yaw_deg)
            % Convenient helper function
            % Relies on subclass's implementation of setPoseRollPitchYaw
            setPoseRollPitchYaw(obj,degtorad(roll_deg),degtorad(pitch_deg),degtorad(yaw_deg));
        end
        function[] = setPoseRollPitchYaw(obj,roll,pitch,yaw)
            if ~obj.supportsRotation
                % Not all arrays support rotation. Nevertheless incuding
                % this functionality in the base class is more convenient
                % than re-implementing it
                error('This microphone array does not support rotation')
            end
            if ~isscalar(roll) || ~isscalar(pitch) || ~isscalar(yaw)
                error('roll,pitch and yaw must be scalars')
            end
            % rotates the array by the specified Euler rotations
            % based on code from rotate_cart_axes_ypr.m
            % requires voicebox
            obj.poseQuaternion = roteu2qr('xyz',[roll,pitch,yaw]);
            % rotqrvec expects positions as column vectors
            obj.poseAxes = rotqrvec(obj.poseQuaternion,eye(3)).';
            obj.sensorCartesianPositions = rotqrvec(...
                obj.poseQuaternion,...
                obj.sensorCartesianPositionsDefault.').';
            
            % N.B. Don't rotate obj.arrayGeometry - potentially huge so
            % only do it if we need to plot it
            
            % No need to store these for now
            %obj.roll = roll;
            %obj.pitch = pitch;
            %obj.yaw = yaw;
        end
        function[rel_src_az,rel_src_inc] = planeWaveDoaWrtArray(obj,src_az,src_inc)
            % planeWaveDoaWrtArray finds the input source(s)' directions with respect 
            % to the (rotated) array
            % plane wave means that the DOA is independent of source
            % distance and of any array translation
            % based on code from dir_wrt_rotated_axes.m
            % requires voicebox
            if obj.poseAxes==eye(3)
                % array is not rotated so nothing to do
                rel_src_az = src_az;
                rel_src_inc = src_inc;
            else
                % DOAs as unit vectors pointing towards DOA
                [x,y,z] = mysph2cart(src_az,src_inc,ones(size(src_az)));

                % inclination is the angle between each vector and the (rotated) z axis
                rel_src_inc = pi*distcos(obj.poseAxes(3,:),[x,y,z]);
                rel_src_inc = rel_src_inc(:);

                % project src vectors into rotated xy plane
                x_proj = (obj.poseAxes(1,:) * [x.';y.';z.']).';
                y_proj = (obj.poseAxes(2,:) * [x.';y.';z.']).';
                rel_src_az = atan2(y_proj,x_proj);
            end
        end
        function[ir,rel_src_az,rel_src_inc] = getImpulseResponseForSrc(obj,src_az,src_inc,varargin)
            [H,rel_src_az,rel_src_inc] = getFrequencyResponseForSrc(obj,src_az,src_inc,varargin{:});
            ir = ifft(H,obj.nfft,'symmetric');
        end
        
        
        function[sensorAzimuth,sensorInclination,sensorRadius] = getSensorPositionsSphericalCoordinates(obj)
            [sensorAzimuth,sensorInclination,sensorRadius] = mycart2sph(obj.sensorCartesianPositions);
        end
        
        %% Plotting functions
        function[] = plot(obj,axh)
            if nargin<2 || isempty(axh)
                axh = gca;
            end
            % old logic
%             if obj.doFastPlot
%                 %only plot the sensor positions
%                 obj.plotSensorPositions(axh);
%             else
%                 %plotting multiple elements so need to control hold state
%                 prevHoldState = axh.NextPlot;
%                 % first element uses existing hold state
%                 obj.plotRotatedAxes(axh);
%                 axh.NextPlot = 'add';
%                 obj.plotGeometry(axh);
%                 obj.plotSensorPositions(axh);
%                 % restore hold state
%                 axh.NextPlot = prevHoldState;
%             end
            % new logic
            
                %plotting multiple elements so need to control hold state
                prevHoldState = axh.NextPlot;
                % first element uses existing hold state
                obj.plotRotatedAxes(axh);
                axh.NextPlot = 'add';
                obj.plotBoresight(axh);
            if ~obj.doFastPlot
                obj.plotGeometry(axh);
            end
                obj.plotSensorPositions(axh);
                % restore hold state
                axh.NextPlot = prevHoldState;
            xlabel(axh,'x [m]')
            ylabel(axh,'y [m]')
            zlabel(axh,'z [m]')
        end
        
        function[] = plotRotatedAxes(obj,axh)
            if nargin<2 || isempty(axh)
                axh = gca;
            end
            scale_factor = 1.5*max(max(0.1,obj.maxSensorRadius));
            colors = {'r','g','b'};
            prevHoldState = axh.NextPlot;
            for iax = 1:3
                pos_vec = obj.poseAxes(iax,:);
                new_vec = scale_factor.*[pos_vec;-pos_vec];
                plot3(axh,new_vec(:,1),new_vec(:,2),new_vec(:,3),'color',colors{iax});
                axh.NextPlot = 'add';
            end
            axh.NextPlot = prevHoldState;
            axh.DataAspectRatio = [1 1 1];
        end

        function[] = plotBoresight(obj,axh)
            if nargin<2 || isempty(axh)
                axh = gca;
            end
            scale_factor = 1.5*max(max(0.1,obj.maxSensorRadius));
            prevHoldState = axh.NextPlot;
            pos_vec = obj.poseAxes(1,:);
            new_vec = scale_factor.*[0 0 0;pos_vec];
            plot3(axh,new_vec(:,1),new_vec(:,2),new_vec(:,3),...
                'color',[0 0 0],...
                'linewidth',1.5);
            axh.NextPlot = 'add';
            axh.NextPlot = prevHoldState;
            axh.DataAspectRatio = [1 1 1];
        end
        
        function[] = plotSensorPositions(obj,axh)
            if nargin<2 || isempty(axh)
                axh = gca;
            end
            if obj.hasDefinedSensorPositions
                scatter3(axh,obj.sensorCartesianPositions(:,1)...
                    , obj.sensorCartesianPositions(:,2)...
                    , obj.sensorCartesianPositions(:,3)...
                    , 20,[1 0 0],'o','filled')
                axh.DataAspectRatio = [1 1 1];
            end
        end
        
        function[] = plotGeometry(obj,axh)
            if nargin<2 || isempty(axh)
                axh = gca;
            end
            if ~isempty(obj.arrayGeometry.vertices)
                patch(axh,'vertices',rotqrvec(obj.poseQuaternion,obj.arrayGeometry.vertices.').',...
                    'faces',obj.arrayGeometry.faces,...
                    'facecolor',[0.7 0.8 0.7],...
                    'edgeColor','black',...
                    'facelighting','phong');
                %axis vis3d equal;% off;
                axh.DataAspectRatio = [1 1 1];
                %warning('Need to make this take an axis handle')
                %view ([144 4]);
                light(axh,'Position',[0 0 1],'Style','infinite');
            end
        end
        
        
        %% getters for dependent variables
        function[nSensors] = get.nSensors(obj)
            if obj.hasDefinedSensorPositions
                nSensors = size(obj.sensorCartesianPositionsDefault,1);
            else
                nSensors = 0;
            end
        end
        function[result] = get.hasDefinedSensorPositions(obj)
            if isempty(obj.sensorCartesianPositionsDefault)
               result = false;
            else
                result = true;
            end
        end
    end
    methods (Static)
        function[x,y,z] = mysph2cart(az,inc,r)
            % mysph2cart converts from azimuth, inclination and radius to cartesian x,y,z
            z = r .* cos(inc);
            rcosinc = r .* sin(inc);
            x = rcosinc .* cos(az);
            y = rcosinc .* sin(az);
        end
    end
end


