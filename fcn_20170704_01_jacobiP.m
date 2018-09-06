function[P] = fcn_20170704_01_jacobiP(n,a,b,x)

% use toolbox to evaluate the recurrence relations
% to avoid duplication of effort, only evaluate the largest n for each
% unique combination of a,b,x;

%TODO: more flexible variable expansion
n = n(:);
a = a(:);
b = b(:);
x = x(:);
szn = size(n);
sza = size(a);
szb = size(b);
szx = size(x);
if ~isequal(szn,sza) || ~isequal(szn,szb) 
    error('Inputs must all be the same size')
end
if ~isequal(szn,szx)
    if isequal(szx,[1 1])
        x = repmat(x,szn(1));
    else
        error('x input must all be the same size or scalar')
    end
end


% matlab jacobiP - tic time: 4.30 s

%version 1 - horrible loop - tic time: 0.35 s
nvals = szn(1);
P = zeros(nvals,1);
for ii = 1:nvals
    tmp = j_polynomial ( 1, n(ii), a(ii), b(ii), x(ii) );
    P(ii) = tmp(end);
end

%version 2 - short loop - tic time: 1.04 s
% profiler suggests that j_polynomial is called the same number of times in
% both cases and I suspect the assignment operation requires that all the
% outputs are written, whereas in version 1 there is perhaps some internal
% optimisation
% cat_a_b_x = [a b x];
% [cat_a_b_x_can, ~,i_can] = unique(cat_a_b_x,'rows');
% maxn = max(n);
% nvals = szn(1);
% ncanvals = size(cat_a_b_x_can,1);
% Ptmp = zeros(ncanvals,maxn+1);
% for ii = 1:ncanvals
%     Ptmp(ii,:) = j_polynomial ( 1, maxn, cat_a_b_x_can(ii,1),...
%         cat_a_b_x_can(ii,2), cat_a_b_x_can(ii,3) );
% end
% P2 = Ptmp(sub2ind([ncanvals maxn+1],i_can,n+1));

