%% function of az and inc
myfun = @(az,inc) 0.4 + 0.2*(cos(az)).^2.*sin(inc) + 0.1*sin(az).*sin(inc) + 0.2*cos(inc);
plotfun = @(ax,xyz,color) plot3(xyz(:,1),xyz(:,2),xyz(:,3),'o','color',color,'markerfacecolor',color);

% Rafaely 2015 version
Rz = @(theta)[cos(theta), -sin(theta),           0;
    sin(theta),  cos(theta),           0;
    0,           0,           1 ];
Ry = @(theta)[cos(theta),           0,  sin(theta);
    0,           1,           0;
    -sin(theta),          0,  cos(theta) ];

N_s = 3;

%% array rotation
rot.rpy.yaw = deg2rad(30);
rot.rpy.pitch = deg2rad(10);
rot.rpy.roll = deg2rad(0);
ema = RigidSphere_20171120_01_horiz_up_20deg_circle4;
ema.setPoseRollPitchYaw(rot.rpy.roll,rot.rpy.pitch,rot.rpy.yaw);

%% world view
[world.grid.az,world.grid.inc,world.grid.weight] = fcn_20171102_01_quad_grid(5,'az_sep_deg');

world.grid.xyz = mysph2cart(world.grid.az,world.grid.inc);

world.s.f = myfun(world.grid.az,world.grid.inc);
if any(world.s.f<0),error('Negative functions make 3d plots messy'),end
world.s.xyz = bsxfun(@times,world.s.f,world.grid.xyz);

figure;plotfun(gca,world.s.xyz,'r');fcn_20180328_01_fix3daxes(gca);

world.y = sphBasis(world.grid.az,world.grid.inc,N_s); %[nDOA, nSH] 
world.s.fnm = world.y' * (world.grid.weight .* world.s.f);

%% array view
[array.grid.az,array.grid.inc,array.grid.weight] = fcn_20171102_01_quad_grid(5,'az_sep_deg');
array.grid.xyz = mysph2cart(array.grid.az,array.grid.inc);

[rot.eu.alpha,rot.eu.beta,rot.eu.gamma] = fcn_20180326_02_body_rotation_qr_to_wignerD_euler(ema.poseQuaternion);

% rotate array grid
array.grid_wrt_world.xyz = (Rz(rot.eu.alpha)*Ry(rot.eu.beta)*Rz(rot.eu.gamma)* array.grid.xyz.').';
% convert to spherical coordinates
[array.grid_wrt_world.az,array.grid_wrt_world.inc,~] = mycart2sph(array.grid_wrt_world.xyz);
% evaluate function
array.s.f = myfun(array.grid_wrt_world.az,array.grid_wrt_world.inc);
array.s.xyz = bsxfun(@times,array.s.f,array.grid.xyz);

figure;plotfun(gca,array.s.xyz,'k');fcn_20180328_01_fix3daxes(gca);

% rotate SH representation of world - inverse rotation
Dinv = sparse(fcn_20170616_02_myWignerDMatrix_Jacobi(N_s,...
            -rot.eu.gamma,-rot.eu.beta,-rot.eu.alpha));
array.s.fnm = Dinv * world.s.fnm;
array.y = sphBasis(array.grid.az,array.grid.inc,N_s); %[nDOA, nSH]
array.s.f_from_fnm = array.y * array.s.fnm;
array.s.xyz_from_fnm = bsxfun(@times,array.s.f_from_fnm,array.grid.xyz);

array.s.fnm_from_f = array.y' * (array.grid.weight .* array.s.f);

figure;
plotfun(gca,array.s.xyz,'k');
hold all
plotfun(gca,array.s.xyz_from_fnm,'b');
fcn_20180328_01_fix3daxes(gca);


peak_doa.world.az = 0;
peak_doa.world.inc = pi/2;
peak_doa.world.xyz = mysph2cart(peak_doa.world.az,peak_doa.world.inc);

grid.world.angle_with_peak_doa = acos(grid.world.xyz * peak_doa.world.xyz.');

%N_cardioid = 3;
%s.world.mag = 0.5^N_cardioid * (1+cos(grid.world.angle_with_peak_doa)).^N_cardioid;





grid.array.xyz = rotqrvec(ema.poseQuaternion,grid.world.xyz.').';

N_cardioid = 3;
s.world.mag = 0.5^N_cardioid * (1+cos(grid.world.angle_with_peak_doa)).^N_cardioid;
s.world.xyz = bsxfun(@times,s.world.mag,grid.world.xyz);




figure;
%ax=gca;
for iax = 1:2
ax(iax) = subplot(1,2,iax);
switch iax
    case 1
plot3(grid.world.xyz(:,1),grid.world.xyz(:,2),grid.world.xyz(:,3),'ok','markerfacecolor',[0 0 0]);
hold all
    case 2
plot3(grid.array.xyz(:,1),grid.array.xyz(:,2),grid.array.xyz(:,3),'ok','markerfacecolor',[1 0 0]);
    case 3
        
end
xlabel('x');
ylabel('y');
zlabel('z');
set(gca,'DataAspectRatio',[1 1 1])
axis([-1 1 -1 1 -1 1])
end
setappdata(gcf,'lp',linkprop(ax,{'CameraPosition','CameraUpVector','CameraViewAngle','CameraTarget'}));

plot3(s.world.xyz(:,1),s.world.xyz(:,2),s.world.xyz(:,3),'ok','markerfacecolor',[0 0 0]);

[rot.eu.alpha,rot.eu.beta,rot.eu.gamma] = fcn_20180326_02_body_rotation_qr_to_wignerD_euler(ema.poseQuaternion);
