classdef FreeFieldArray < AnalyticalArray
    % FreeFieldArray provides the analytical formulation for an array in
    % free space
    methods
        function[H,rel_src_az,rel_src_inc] = getFrequencyResponseForSrc(obj,src_az,src_inc,varargin)
            % getFrequencyResponseForSrc returns the frequency response as a 3D array [nFreq,nSensors,nDOAs]
            %
            % first finds DOAs with respect to the (possibly rotated) array
            %
            % based on code from free_field_plane_wave_rel_wrt_origin
            
            if nargin>3
                warning('Inputs other than src_az and src_inc are not yet supported - assuming plane wave incidence')
            end
            szAz = size(src_az);
            szInc = size(src_inc);
            nDOAs = szAz(1);
            if ~isequal(szAz,szInc) || ~isequal(prod(szAz),nDOAs)
                error('src_az and src_inc must be column vectors of the same size');
            end
            
            % determine the source directions relative to the rotated array
            [rel_src_az,rel_src_inc] = obj.planeWaveDoaWrtArray(src_az,src_inc);
            
            % get unit vectors pointing *towards* doas
            % these are normal to plane of the propagating wave
            [doa_x,doa_y,doa_z]  = mysph2cart(rel_src_az,rel_src_inc,ones(nDOAs,1));
            
            % Distance to plane: D = \hat{n} \cdot \hat{x}_0 + p
            % (http://mathworld.wolfram.com/Point-PlaneDistance.html)
            % p is 0 for plane intersection origin
            % x_0 is our point(s) of interest
            % cdot indicates dot product, or simply inner product
            
            % grid of distances
            d = sum(bsxfun(@times,...
                permute(obj.sensorCartesianPositionsDefault,[2 1 3]),... %[3,nMics,1]
                permute([doa_x,doa_y,doa_z],[2 3 1])),... %[3,1,nDOAs]
                1); %[1,nSensors, nDOAs]
            
            % calculate magnitude and phase relative to response at origin
            % and combine with the response at the origin
            % mag = ones(length(obj.f),obj.nSensors,nDOAs); % relative response
            % mag = bsxfun(@times,mag,abs(obj.H0));         % combined
            mag = repmat(abs(obj.H0),1,obj.nSensors,nDOAs);      % since mag at origin is just ones, save ourselves the multiplication
            
            phase = 2*pi/obj.c * bsxfun(@times, d, obj.f);
            phase = bsxfun(@plus,angle(obj.H0),phase);
            
            H = mag .* exp(1i*phase);           
        end              
        function[rMin, rMax] = getValidSrcRadiusRange(obj)
            % TODO: Add support for point sources with a defined distance
            % For now the only valid range is infinity (i.e. a plane wave)
            % so we don't actually use obj
            rMin = Inf;
            rMax = Inf;
        end        
    end
end