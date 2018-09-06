function[az_rot,inc_rot] = dir_wrt_rotated_axes(ax_vec_rot,az,inc)
%dir_wrt_rotated_axes
%
%Inputs:
% ax_vec_rot: cartesian axes [x,y,z]=eye(3) after undergoing some rotation 
%         az: vector of azimuths wrt original cartesian axes
%        inc: vector of inclinations wrt original cartesian axes
%
%Outputs
%     az_rot: vector of azimuths wrt rotated axes
%    inc_rot: vector of inclinations wrt rotated axes

if ~isequal(size(az),size(inc))
    error('az and inc must be the same size');
end

[x,y,z] = mysph2cart(az,inc,ones(size(az)));

% I think there is probably an easy way to do this - probably using the
% inverse of the quaternion. To get a sense of it though lets find the
% angles explicitly.



% inclination is the angle between the vector and the (rotated) z axis
inc_rot = pi*distcos(ax_vec_rot(:,3).',[x,y,z]);
inc_rot = inc_rot(:);

% project src vectors into rotated xy plane
x_proj = (ax_vec_rot(:,1).' * [x.';y.';z.']).';
y_proj = (ax_vec_rot(:,2).' * [x.';y.';z.']).';

az_rot = atan2(y_proj,x_proj); 