classdef (HandleCompatible) BinauralArray
    % BinauralArray defines abstract properties (which must be
    % defined in a subclass) and constant properties (which all objects of
    % this types share) of arrays which are specifically arranged into left
    % and right parts
    properties (Abstract, SetAccess=protected)
        refChanLeft
        refChanRight
        channelsLeft
        channelsRight
    end
    properties (Constant)
        supportsBinaural = 1
    end  
end
        
       