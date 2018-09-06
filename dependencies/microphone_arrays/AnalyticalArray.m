classdef AnalyticalArray < ElobesMicArray
    % AnalyticalArray provides a base class for arrays whose manifold can
    % be calculated anyalytically
    % No special work is required to handle different types of rotation
    properties (Dependent)
        maxSensorRadius; % distance of the microphone furthest from origin
    end
    methods
        function[obj] = AnalyticalArray()
            % call superclass's constructor
            obj = obj@ElobesMicArray();
            
            % override the default properties of the superclass
            obj.supportsRotation=1;
            obj.availableInterpolationMethods = {'analytical'};
            obj.interpolationMethod = 'analytical';
            
        end
        function[] = prepareData(obj,req_fs,varargin)
            % analytical array so we get to choose the sample rate
            obj.fs = req_fs;
            
            % define the response at the origin (if it hasn't been
            % specified) which will form the basis of the other
            % microphone responses
            
            % f0 of sinc determines the bandwidth
            % len_filt determines the resolution/roll off
            f0=0.95*obj.fs/2;
            len_filt = 99;
            half_len = (len_filt-1)/2;
            % the pulse at the origin
            h = hamming(len_filt) .* (sinc(2*f0*(-half_len:half_len).'/obj.fs));
            
            % need to zero pad both sides to allow for relative propagation
            % delay to furtherst sensor
            pad_len = ceil(obj.maxSensorRadius / obj.c * obj.fs);
            h = [zeros(pad_len,1); h; zeros(pad_len,1)];
            obj.nfft = size(h,1);
            obj.t0 = pad_len+half_len+1; % sample offset of main peak
            obj.H0 = rfft(h);
            obj.f = (0:((obj.nfft+2)/2 - 1)).' * obj.fs/obj.nfft; %[nf 1] where nf is number of bins up to nyquist
            
        end
        function[r] = get.maxSensorRadius(obj)
            % derive some stuff
            if obj.hasDefinedSensorPositions
                sensorRadii = cartnorm(obj.sensorCartesianPositionsDefault);
                r = max(sensorRadii);
            else
                r = 0;
            end
        end
    end
end