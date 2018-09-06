function[alpha, beta, gamma] = body_rotation_qr_to_wignerD_euler_v2(qr)
% Given the quaternion which rotates a body with respect to fixed world coordinates
% obtain the euler angles needed to obtain WignerD matrix which will
% similarly rotate a vector expressed as spherical harmonic coefficients
%
% http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/rotation.html
%
[rot_angles] = rotqr2eu('zyz',qr);


gamma = rot_angles(1);
beta = rot_angles(2);
alpha = rot_angles(3);
% alpha = neg_angles(1);
% beta = neg_angles(2);
% gamma = neg_angles(3);

