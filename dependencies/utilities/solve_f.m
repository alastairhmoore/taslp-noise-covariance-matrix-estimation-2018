function y = solve_f(func,target,tor)
% Intro: SOLVE_F computes the solution of a function.
%
%Input:
% func:     Handle of function
% target:   Domain value
% tor:      Threashold
%Output:
% y:     	Solution
%

% Copyright (C) Wei Xue, Imperial College London
% Date: 2016-11-10

x = 1;
bool_label = 0;

while 1
	if abs(func(x)-target)< tor
		break;
	end

	if bool_label == 0
        if func(x) < target
            x_l = x;
            x_r = x;
            while 1
                x_r = x_r*2;
                if func(x_r)>target
                    break;
                end
            end
        else
            x_l = x;
            x_r = x;
            while 1
                x_l = x_l/2;
                if func(x_l)<target
                    break;
                end
            end
        end
        bool_label = 1;
    end
    
    x = (x_l + x_r)/2;
    if func(x) < target
        x_l = x;
    else
        x_r = x;
    end
end
y = x;