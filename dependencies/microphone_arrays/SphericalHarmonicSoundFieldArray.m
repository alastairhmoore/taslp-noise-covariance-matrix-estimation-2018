classdef SphericalHarmonicSoundFieldArray < AnalyticalArray
    % SphericalHarmonicSoundFieldArray provides the analytical formulation for
    % the sound field in free space
    % By default real spherical harmonics are used such that the time
    % domain response is real valued
    properties
        maxSphericalHarmonicOrder = [];             % must specify this before doing prepareData
    end
    properties (SetAccess = protected)
        % properties from ElobesMicArray
        sensorCartesianPositionsDefault;
        refChan;
    end
    properties (Hidden)
        % Store some precomputed values
        nSH = [];
        k_sphOrd = [];
        fact = [];
    end
    methods
        function[obj] = SphericalHarmonicSoundFieldArray()
            
            % call superclass's constructor
            obj=obj@AnalyticalArray();
            
            % populate the parameters
            obj.sensorCartesianPositionsDefault = []; % Make it empty since there are no actual sensors
            obj.refChan = 1;                          % Corresponds to 0 order SH
        end
        function[H,rel_src_az,rel_src_inc] = getFrequencyResponseForSrc(obj,src_az,src_inc,varargin)
            % getFrequencyResponseForSrc returns the freqeuncy response of each real spherical harmonic [nFreq,nSensors,nDOAs]
            %
            % first finds DOAs with respect to the (possibly rotated) array
            % TODO: Find out whether the inputParser has significant impact
            % on performance
            p = inputParser;
            p.addRequired('src_az',@isnumeric);
            p.addRequired('src_inc',@isnumeric);
            useComplexShDefault = false;
            p.addParameter('useComplexSh',useComplexShDefault,@islogical);
            p.parse(src_az,src_inc,varargin{:});
            
            %copy parsed parameters for ease of reading
            src_az = p.Results.src_az;
            src_inc = p.Results.src_inc;
            useComplexSh = p.Results.useComplexSh;
            
%             if useComplexSh
%                 warning('Complex spherical harmonics cannot be used to directly return an impulse response')
%                 warning('We could/should add a check to identify the calling function')
%             end
            
            szAz = size(src_az);
            szInc = size(src_inc);
            nDOAs = szAz(1);
            if ~isequal(szAz,szInc) || ~isequal(prod(szAz),nDOAs)
                error('src_az and src_inc must be column vectors of the same size');
            end
            
            % determine the source directions relative to the rotated array
            [rel_src_az,rel_src_inc] = obj.planeWaveDoaWrtArray(src_az,src_inc);
            
            % evaluate spherical harmomics in the source directions
            Y = obj.evaluateSphericalHarmonics(rel_src_az,rel_src_inc,useComplexSh);
            if useComplexSh
                Y = conj(Y);
            end
            
            % by defintion everything is real so we can just multiply by the
            % response at the origin
            % H0: [nFreq,1]
            %  Y: [nDOAs,nSH]
            H = bsxfun(@times,obj.H0,permute(Y,[3,2,1])); %profile: 6.69
            %H = obj.H0 .* permute(Y,[3,2,1]); %profile: 9.45
        end
        
        function[Y] = evaluateSphericalHarmonics(obj,az,inc,useComplexSh)
            % enforce column vectors
            az = az(:);
            inc = inc(:); 
            
            % code based on sphBasisReal.m and sphBasis.m)
            Y = zeros(length(az),obj.nSH);
            Y(:,1) = sqrt(obj.k_sphOrd(1)); % zeroth component is constant (omnidirectional)
            for sphOrd = 1:obj.maxSphericalHarmonicOrder
                
                P = legendre(sphOrd,cos(inc)).';      %[nDirections, sphOrd+1];
                % - only defined for positive SH orders
                % - use index sphDeg+1 to get sphOrd component
                
                % 0 degree
                col_i = sphOrd^2 + sphOrd + 1;        %use 0 degree index into output matrix as reference
                Y(:,col_i) = sqrt(obj.k_sphOrd(sphOrd+1)) * P(:,1);
                
                % non-0 degrees
                if useComplexSh
                    % Uses a normalisation factor which
                    % - indexes into the main precomputed normalisations k_sphOrd
                    % - indexes into precomputed factorial functions
                    for sphDeg = 1:sphOrd
                        NormFactor = sqrt(obj.k_sphOrd(sphOrd+1) * ...
                            obj.fact(sphOrd-sphDeg+1)/obj.fact(sphOrd+sphDeg+1)) * P(:,sphDeg+1);
                        Y(:,col_i+sphDeg) = NormFactor .* exp(1i*sphDeg*az); %positive sphDeg;
                        Y(:,col_i-sphDeg) = (-1)^-sphDeg .* conj(Y(:,col_i+sphDeg)); %negative sphDeg
                    end
                else
                    % Uses a normalisation factor which
                    % - undoes the Cordon-Shortley phase included in legendre
                    % - indexes into the main precomputed normalisations k_sphOrd
                    % - indexes into precomputed factorial functions
                    % NormFactor = (-1).^sphDeg .* sqrt(2*k_sphOrd(sphOrd+1)*(fact(1+sphOrd-sphDeg)./fact(1+sphOrd+sphDeg)))
                    for sphDeg = 1:sphOrd
                        NormFactor = (-1)^sphDeg * sqrt(2*obj.k_sphOrd(sphOrd+1) * ...
                            obj.fact(sphOrd-sphDeg+1)/obj.fact(sphOrd+sphDeg+1)) * P(:,sphDeg+1);
                        Y(:,col_i-sphDeg) = NormFactor .* sin(sphDeg*az); %negative sphDeg
                        Y(:,col_i+sphDeg) = NormFactor .* cos(sphDeg*az); %positive sphDeg
                    end
                end
            end
        end
            
        
        function[rMin, rMax] = getValidSrcRadiusRange(obj)
            % TODO: Add support for point sources with a defined distance
            % For now the only valid range is infinity (i.e. a plane wave)
            % so we don't actually use obj
            rMin = Inf;
            rMax = Inf;
        end
        function[] = prepareData(obj,req_fs,varargin)
            %             p = inputParser;
            %             addRequired(p,'req_fs',@isnumeric);
            %             addRequired(p,'order',@isnumeric);
            %             parse(p,req_fs,varargin{:});
            %             obj.maxSphericalHarmonicOrder = p.order;
            %obj.maxSphericalHarmonicOrder = order;
            % is this really the best way to set this parameter?
            if isempty(obj.maxSphericalHarmonicOrder)
                error('Must set value of maxSphericalHarmonicOrder before calling prepareData')
            end
            
            % populate constants that are used in calculation of spherical
            % harmonics
            obj.nSH = (obj.maxSphericalHarmonicOrder+1)^2;                    %[1 1]; - total number of spherical harmonics
            obj.k_sphOrd = (2*(0:obj.maxSphericalHarmonicOrder).' +1)/(4*pi); %[sphOrdMax+1 1]; - use index sphOrd+1 to get sphOrd component
            obj.fact = factorial((0:2*obj.maxSphericalHarmonicOrder).');      %[2*sphOrdMax 1]; - use index i to get factorial(i)
            
            % call superclass's version of method
            % to populate the omni response at the origin and set fs
            obj.prepareData@AnalyticalArray(req_fs);
        end
        
    end
    %% static method - uses
end