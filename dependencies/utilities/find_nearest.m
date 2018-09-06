function[acheived,i] = find_nearest(desired,available)

[~,i] = min(abs(desired-available));
acheived = available(i);