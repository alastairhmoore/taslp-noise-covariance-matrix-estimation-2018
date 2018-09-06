function [x,y,z] = mysph2cart(az,inc,r)

%   Spherical Coordinates to 3D Cartesian
%
%   [x,y,z] = mysph2cart(az,inc,r)
%
%   Inputs:
%       az      Azimuth angle in radians, counterclockwise from xy plane 
%               from the positive x axis (otherwise referred to as phi)
%       inc     Inclination angle in radians, from positive z axis 
%               (otherwise referred to as theta) 
%       r       Radius
%
%   Outputs:
%       x       x-coordinate
%       y       y-coordinate
%       z       z-coordinate
%
%   Notes:
%       The MATLAB function cart2sph reverses phi and theta.
%
%**************************************************************************
% Author:           E. A. P. Habets, M. R. P. Thomas and D. P. Jarrett
% Date:             27 July 2010
% Revised:          A. H. Moore, 6 February 2018
% Version: $Id: mysph2cart.m 351 2011-07-21 15:10:09Z dpj05 $
%**************************************************************************
if nargin<3
    r = ones(size(az));
end

z = r .* cos(inc);
rcosinc = r .* sin(inc);
x = rcosinc .* cos(az);
y = rcosinc .* sin(az);

if nargout==1
    [nDir nEx] = size(az);
    if nEx>1, error('Can''t use single output argument for multidimensional data'),end
    x = [x,y,z];
end