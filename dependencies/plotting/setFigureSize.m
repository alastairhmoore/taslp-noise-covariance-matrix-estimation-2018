function[] = setFigureSize(figno,dimensions,units)

px = 0.01; %fudge factors to avoid clipping, set to zero to restore old behaviour
py = 0.04;

%resize figure window
set(figno,'Units',units);
pos = get(figno,'position');
set(figno,'Position',[pos(1:2) dimensions]);

%check that they were set correctly (large or small values may not be allowed)
new_pos = get(figno,'Position');
if abs(new_pos(3)-dimensions(1)) > 0.1 || abs(new_pos(4)- dimensions(2)) > 0.1
    disp('Warning: On screen figure may not be actual size')
end
set(figno,'PaperUnits',units);
set(figno,'PaperPosition',[px*dimensions(1), py*dimensions(2), dimensions]);
set(figno,'PaperSize',[(1+2*px) * dimensions(1), (1+2*px) * dimensions(2)]);