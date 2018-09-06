function[] = fcn_20180328_01_fix3daxes(ax)
axes(ax)
xlabel('x');
ylabel('y');
zlabel('z');
set(ax,'DataAspectRatio',[1 1 1])
axis([-1 1 -1 1 -1 1])
