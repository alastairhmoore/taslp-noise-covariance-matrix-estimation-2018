classdef RigidSphereArray < SphericalHarmonicSoundFieldArray
    % RigidSphereArray provides an abstract class for implementing arrays
    % centred on a rigid spherical baffle
    % No restrictions are placed on the locations of the sensors so it does
    % not make sense to return the spherical harmonic decomposition of the
    % microphone signals here. Instead a subclass should also subclass
    % SphericalArray, which provides that funcionality
    % Slightly awkward/confusing thing is that the superclass
    % SphericalHarmonicSoundFieldArray has defined sensor postions of [0 0
    % 0]. In this class they are not defined but any subclasses must have
    % them defined and not [0 0 0].
    
    
%     properties (Abstract)
%         
%         
%     end
    properties (SetAccess=protected)
        % define the radius of the scattering sphere (often called 'a')
        sphereRadius;
        % sensorModeStrength - Mx1 vector populated during prepareData based on
        % sphereRadius and sensorCartesianPositionsDefault
        sensorModeStrength = [];
        % sphHarmOrderMax - The spherical harmonic truncation order used in
        % the inverse spherical Fourier transform
        %sphHarmOrderMax = [];
        % sensorInverseSftMatrix - The precomputed element-wise product of
        % the SH mode strength and SH function at each sensor
        sphHarmFunc = [];
        sensorInverseSftMatrix = [];
        
    end
    methods
        function[obj] = RigidSphereArray(radius)
            % RigidSphereArray encapsulates the analytical response of microphones
            % on the surface of, or surrounding, a rigid sphere
            % It is an abstract class which realisations of specific microphone
            % arrangements should inherit from
            %
            if nargin == 0
                radius = 0.1;
            end
            % call superclass's constructor
            obj=obj@SphericalHarmonicSoundFieldArray();
            
            % set properties
            obj.sphereRadius = radius;
            
            [x,y,z] = sphere;
            x = x * obj.sphereRadius;
            y = y * obj.sphereRadius;
            z = z * obj.sphereRadius;
            [obj.arrayGeometry.faces,obj.arrayGeometry.vertices,~] = ...
                surf2patch(x,y,z);
        end
            
            
        function[] = prepareData(obj,req_fs,varargin)
            % prepareData does all the precalculations to initialise variables
            % and speed up subsequent calls to getImpulseResponseForSrc()
            
            % parse the inputs to this function
            p = inputParser;
            p.addRequired('req_fs',@isnumeric);
            p.addParameter('sphHarmOrderMax',[],@isnumeric);
            p.addParameter('freqMax',[],@isnumeric)
            p.parse(req_fs,varargin{:});
            
            % going to use superclass to do most of the work. Need to
            % convert some parameters first
            
            %copy parsed parameters for ease of reading
            req_fs = p.Results.req_fs;
            req_sphHarmOrderMax = p.Results.sphHarmOrderMax;
            freqMax = p.Results.freqMax;
            
            
            % validate sphHarmOrderMax
            % the inverse SFT to be used is actually an infinite sum
            % we need to specify/validate that truncation order is adequate
            % the error introduced by a particular truncation order is frequency
            % dependent and also source-directivity dependent
            % in some circumstances we may choose to limit the order deliberately to avoid
            % spatial aliasing for a particular sensor layout
            if isempty(freqMax)
                freqMax = req_fs/2;
            end
            
            % use freqMax and array geometry to determine the required
            % spherical harmonic order
            % k = 2*pi*f/c
            % calculate truncation error at highest frequency and with sensor
            % furthest from the origin based on Jin2014
            SPHERICAL_ORDER_UPPER_LIMIT = 500; %ridiculuously high, so error is relative to this
            MAX_ALLOWED_TRUNCATION_ERROR = 10^(-80/20);
            b = modeStrength('rigid' ....
                ,obj.sphereRadius ...
                ,obj.maxSensorRadius ...
                ,freqMax,SPHERICAL_ORDER_UPPER_LIMIT,obj.c);
            i_first_invalid = find(isnan(b),1,'first'); % values get so small they become nans
            b(i_first_invalid:end) = [];
            tot = cumsum(b .* (2*(0:length(b)-1)+1));
            truncationError = abs(1-tot./tot(end));
            if 0
            figure
            plot(0:length(b)-1,20*log10(truncationError));
            xlabel('Truncation order')
            ylabel('Truncation error [dB]')
            end
            
            % find order which meets predefined threshold ('1+..' because 0 order is index 1)
            sphHarmOrderMaxDefault = -1+find(truncationError<MAX_ALLOWED_TRUNCATION_ERROR,1,'first');
            
            if isempty(req_sphHarmOrderMax)
                obj.maxSphericalHarmonicOrder = sphHarmOrderMaxDefault;
            else               
                obj.maxSphericalHarmonicOrder = req_sphHarmOrderMax;
                if req_sphHarmOrderMax < sphHarmOrderMaxDefault
                    % issue a warning but do it anyway
                    warning('With sphHarmOrderMax set to %d the truncation error at %2.2f Hz is %2.2f dB',...
                        req_sphHarmOrderMax, freqMax, 20*log10(truncationError(req_sphHarmOrderMax+1)));
                end
                
            end
            
            % call superclass's method to populate required properties, especially
            % f (the cetre frequencies of the FFT bins)
            obj.prepareData@SphericalHarmonicSoundFieldArray(req_fs);
            
            % The inverse SFT requires the sensor positions in spherical
            % coordinates
            % radii are used to calculate the mode strength
            % directions are used to evaluate the spherical harmonics
            % N.B. It is equivalent to rotate the underlying sound field as to
            % rotate the sensors. From a computational performance point of
            % we would prefer to evaluate fewer spherical harmonics on each
            % call to getFrequencyResponseForSrc
            % - the soundfield requires nSH evaluations
            % - the sensors require nSensors * nSH evaluations
            % so better to evaluate the SHs for fixed sensor positions
            % in the default array orientation
            
            [sensorAzimuth,sensorInclination,sensorRadius] = mycart2sph(obj.sensorCartesianPositionsDefault);
            nSH = (obj.maxSphericalHarmonicOrder+1)^2;
            obj.sensorModeStrength = zeros(length(obj.f),nSH,obj.nSensors);
            % could possibly speed this up by checking for unique values of
            % radius, but probably not worth the complexity since this
            % function only gets called once
            %fix mode strength at DC
            fft_freq_vec = [0.01; obj.f(2:end)];
            for iSensor = 1:obj.nSensors
                [b,i_b] = modeStrength('rigid' ....
                    ,obj.sphereRadius ...
                    ,sensorRadius(iSensor) ...
                    ,fft_freq_vec,obj.maxSphericalHarmonicOrder,obj.c);
                
                
                obj.sensorModeStrength(:,:,iSensor) = b(:,i_b); % use indexing variable i_b to expand each order by the number of SHs in that order (i.e. for each degree)
            end
            
            % evaluateSphericalHarmonics is inherited from SphericalHarmonicSoundFieldArray
            % true indicates useComplexSH
            % results is [nSensors,nSH]
            obj.sphHarmFunc = obj.evaluateSphericalHarmonics(sensorAzimuth,sensorInclination,true);
            
            % precompute the matrix of weights required to evaluate the
            % inverse SFT
            % result: [nFreq,nSH, 1, nSensors]
            obj.sensorInverseSftMatrix = bsxfun(@times, ...
                permute(obj.sensorModeStrength,[1 2 4 3]), ...  % [nFreq,nSH, 1, nSensors]
                permute(obj.sphHarmFunc,[3 2 4 1]));            % [1    ,nSH, 1, nSensors]
            
        end
        function[H,rel_src_az,rel_src_inc] = getFrequencyResponseForSrc(obj,src_az,src_inc,varargin)
            % getFrequencyResponseForSrc returns the frequency response as a 3D array [nFreq,nSensors,nDOAs]
            %
            % uses superclass to first finds DOAs with respect to the 
            % (possibly rotated) array SphericalHarmonicSoundFieldArray
            % then the modeStrength and sensor angles are used to compute
            % the inverse SFT
            
            [Anm,rel_src_az,rel_src_inc] = getFrequencyResponseForSrc@SphericalHarmonicSoundFieldArray(...
                obj,src_az,src_inc,'useComplexSh',true); %[nFreq,nSH]
            
            % Mulitply the soundfield by the mode-strength adjusted SH function
            % for each sensor direction and sum over all SH
            % permute matrix to get result for each sensor in second
            % dimension
            H = permute(sum(bsxfun(@times,Anm,obj.sensorInverseSftMatrix),2),...
                   [1,4,3,2]); %profile: 52.06
            %H = permute(sum( Anm .* obj.sensorInverseSftMatrix,2),...
            %        [1,4,3,2]); %profile: 54.12
        end
    end
end
    