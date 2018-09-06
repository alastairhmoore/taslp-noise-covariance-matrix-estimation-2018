function[ax_vec_rot] = rotate_cart_axes_ypr(yaw, pitch, roll)

ax_vec_rot = eye(3); %[x,y,z] - before rotation

%apply rotations in turn around the head centred axes
% 1 and 3 give the same result as do 2 and 4.
switch 4
    case 1
        qr_roll = rotax2qr(ax_vec_rot(:,1),roll);
        ax_vec_rot = rotqrvec(qr_roll,ax_vec_rot);
        qr_pitch = rotax2qr(ax_vec_rot(:,2),pitch);
        ax_vec_rot = rotqrvec(qr_pitch,ax_vec_rot);
        qr_yaw = rotax2qr(ax_vec_rot(:,3),yaw);
        ax_vec_rot = rotqrvec(qr_yaw,ax_vec_rot);
    case 2
        qr_yaw = rotax2qr(ax_vec_rot(:,3),yaw);
        ax_vec_rot = rotqrvec(qr_yaw,ax_vec_rot);
        qr_pitch = rotax2qr(ax_vec_rot(:,2),pitch);
        ax_vec_rot = rotqrvec(qr_pitch,ax_vec_rot); 
        qr_roll = rotax2qr(ax_vec_rot(:,1),roll);
        ax_vec_rot = rotqrvec(qr_roll,ax_vec_rot);
    case 3
        qr = roteu2qr('zyx',[yaw,pitch,roll]);
        ax_vec_rot = rotqrvec(qr,ax_vec_rot);
    case 4
        qr = roteu2qr('xyz',[roll,pitch,yaw]);
        ax_vec_rot = rotqrvec(qr,ax_vec_rot);
end