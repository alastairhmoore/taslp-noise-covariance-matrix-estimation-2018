function[A] = fcn_20170619_01_triple_sph_harm_integral(sphOrdMax,sphOrdMax1,sphOrdMax2)
% sphOrd3 refers to the conjugate term

% generate indices for each SH
% use standard 00,1(-1),10,11,...,NN vecorisation

% as rough guide, highest factorial is perhaps 3 * the highest spherical
% harmonic order
sOmax = max([sphOrdMax,sphOrdMax1,sphOrdMax2]);
%fact = factorial((0:3*sOmax).');
local_fact(3*sOmax); %call this to ensure persistent value is populated with a large number


[sOv,sDv] = local_index_gen(sphOrdMax);
[sO1v,sD1v] = local_index_gen(sphOrdMax1);
[sO2v,sD2v] = local_index_gen(sphOrdMax2);

[sO, sO1, sO2] = ndgrid(sOv,sO1v,sO2v);
[sD, sD1, sD2] = ndgrid(sDv,sD1v,sD2v);


% Clebsch-Gordan coefficients have leading \delta_{m2,m+m1} so the matrix
% of terms is very sparse

% find the non-zero terms
idc = find(sD2==sD+sD1);

% pick them out and flatten grid into vectors for ease
n = sO(idc);
n1 = sO1(idc);
n2 = sO2(idc);
m = sD(idc);
m1 = sD1(idc);
m2 = sD2(idc);

% summation over v for all values of v where argument to factorials are >=0
%
% A = k * C(n,n1,n2,0,0,0) * C(n,n1,n2,m,m1,m2

% Caluclate all valid value of C first
% C = sqrt(t0) * sum_v(t(v))
C = zeros(size(n));

for ii = 1:length(n)
    C(ii) = local_C(n(ii),n1(ii),n2(ii),m(ii),m1(ii),m2(ii));
end

% find the entries corresponding to m=m1=m2 = 0
idc_000 = find(all(bsxfun(@eq,[m m1 m2],[0 0 0]),2));

C000 = zeros(size(C));

% populate C000 for corresponding n n1 n2
for ii = 1:length(idc_000)
    C000( n==n(idc_000(ii)) & n1==n1(idc_000(ii)) & n2==n2(idc_000(ii))) = C(idc_000(ii));
end

A = zeros(size(sO)); %preallocate ouptut grid

% final calculation
A(idc) = sqrt( (2*n+1) .* (2*n1+1) ./ ((4 * pi) .* (2*n2 + 1)) ) .* C .* C000;
%A(idc) = C;
%A(idc) = C000;


function[C] = local_C(n,n1,n2,m,m1,m2)
% find range on v (nu in Shabtai2014 but lots of nx already)
lower_bounds = [0;...
                m - n;...
                n1 + m2 - n];
upper_bounds = [n1 + n2 + m;...
                n2 - n + n1;...
                n2 + m2];

v = min(lower_bounds):max(upper_bounds); %full possible range

% t0 is independent of v, [len(idc) 1]
fd1 = local_fact( n + n1 + n2 + 1);
fd2 = local_fact(n - m);
fd3 = local_fact(n + m);
fd4 = local_fact(n1 - m1);
fd5 = local_fact(n1 + m1);
t0_denom = fd1 .* fd2 .* fd3 .* fd4 .* fd5;

if t0_denom==0
    C = 0;
    return
end

fn1 = local_fact( n2 + n - n1);
fn2 = local_fact(n2 - n + n1);
fn3 = local_fact( n + n1 - n2);
fn4 = local_fact(n2 + m2);
fn5 = local_fact(n2 - m2);
t0_num = (2*n2 + 1) .* fn1 .* fn2 .* fn3 .* fn4 .* fn5; % numerator

if t0_num==0
    C = 0;
    return
end
%t0 = zeros(size(fn1));
%t0(t0_denom~=0) = t0_num(t0_denom~=0)./t0_denom(t0_denom~=0);
%
t0 = t0_num/t0_denom;


% expand terms in v, [len(idc) len(v)]
fn1 = local_fact(n1 + n2 + m - v);
fn2 = local_fact(n - m + v);
fd1 = local_fact(v);
fd2 = local_fact(n2 - n + n1 - v);
fd3 = local_fact(n2 + m2 - v);
fd4 = local_fact(v + n - n1 - m2);
tv_num = (-1).^(v + n1 + m1) .* fn1 .* fn2;
tv_den =  fd1 .* fd2 .* fd3 .* fd4;
%remove zeros from denominator to avoid nans
den_zeros = find(tv_den==0);
tv_num(den_zeros) = [];
tv_den(den_zeros) = [];

C = sqrt(t0) .* sum(tv_num./tv_den,2);





%
%
% A = zeros(size(sO));
% for ii = 1:length(idc)
%     % pick out the elements
%     n=sO(idc(ii));
%     n1=sO1(idc(ii));
%     n2=sO2(idc(ii));
%     m=sD(idc(ii));
%     m1=sD1(idc(ii));
%     m2=sD2(idc(ii));
%
%     A(idc(ii)) = sqrt( (2*n+1) * (2*n1+1) / (4 * pi * (2*n2 + 1)) ) * ...
%         nested_C(n,n1,n2,0,0,0) * nested_C(n,n1,n2,m,m1,m2);
% end
%
%
%
%     function[C] = nested_C(n,n1,n2,m,m1,m2)
%         if m2~=(m+m1)
%             C = 0;
%             warning('local_C called when delta condition not met')
%         else
%             % C = sqrt(t0) * sum_v(t(v))
%
%             % find range on v (nu in Shabtai2014 but lots of nx already)
%             lower_bounds = [0;...
%                 m - n;...
%                 n1 + m2 - n];
%             upper_bounds = [n1 + n2 + m;...
%                 n2 - n + n1;...
%                 n2 + m2];
%             v = (max(lower_bounds):min(upper_bounds))';
%             if isempty(v)
%                 C = nan;
%             else
%
%             t0 = (2*n2 + 1) * local_fact( n2 + n - n1) * local_fact(n2 - n + n1) * ...
%                 local_fact( n + n1 - n2) * local_fact(n2 + m2) * local_fact(n2 - m2) / ...
%                 ( local_fact( n + n1 + n2 + 1) * local_fact(n - m) * ...
%                 local_fact(n + m) * local_fact(n1 - m1) * local_fact(n1 + m1) );
%
%
%             tv_num = -1.^(v + n1 + m1) .* local_fact(n1 + n2 + m - v) .* local_fact(n - m + v);
%             tv_den = local_fact(v) .* local_fact(n2 - n + n1 - v) .* local_fact(n2 + m2 - v) ...
%                 .* local_fact(v + n - n1 - m2);
%             C = sqrt(t0) * sum(tv_num./tv_den,1);
%             end
%         end
%     end
%
% %     function[C] = nested_C(n,n1,n2,m,m1,m2)
% %         if m2~=(m+m1)
% %             C = 0;
% %             warning('local_C called when delta condition not met')
% %         else
% %            % C = sqrt(t0) * sum_v(t(v))
% %
% %             % find range on v (nu in Shabtai2014 but lots of nx already)
% %             lower_bound = max([0;...
% %                 n1 - n2 - m;...
% %                 n + m1 - n2]);
% %             upper_bound = min([n + n1 - n2;... %n + n1 - m2;... in Shabtai
% %                 n-m;...
% %                 n1 + m1]);
% %             v = (lower_bound:upper_bound)';
% %             if isempty(v)
% %                 error('no valid values for v')
% %             else
% %                C = 1;
% %             end
% % %
% % %             t0 = (2*n2 + 1) * fact( n2 + n - n1 + 1) * fact(n2 - n1 + n + 1) * ...
% % %                 fact( n + n1 - n2 + 1) * fact(n2 + m2 + 1) * fact(n2 - m2 + 1) / ...
% % %                 ( fact( n + n1 + n2 + 1) * fact(n - m + 1) * ...
% % %                 fact(n + m + 1) * fact(n1 - m1 + 1) * fact(n1 + m1 + 1) );
% % %
% % %
% % %
% % %             tv_num = -1.^(v + n1 + m1) .* fact(n1 + n2 + m - v + 1) .* fact(n - m + v + 1);
% % %             tv_den = fact(v+1) .* fact(n2 - n + n1 - v + 1) .* fact(n2 + m2 - v + 1) ...
% % %                 .* fact(v + n - n1 - m2 + 1);
% % %             C = sqrt(t0) * sum(tv_num./tv_den,1);
% %         end
% %     end
% end

function[sphOrd,sphDeg] = local_index_gen(sphOrdMax)
% could do this faster but does the job
sphOrd = [];
sphDeg = [];

for ord = 0:sphOrdMax
    tmp_deg = (-ord:ord).';
    sphOrd = [sphOrd; ord*ones(size(tmp_deg))];
    sphDeg = [sphDeg; tmp_deg];
end


function[f] = local_fact(x)
% returns the factorial of x
% prod(1:x)  for x>0
% 1          for x==0
% 0          for x<0

persistent fv
f = zeros(size(x));
max_x = max(x(:));
if max_x>length(fv)
    fv = factorial((1:max_x).');
end
f(x==0) = 1;
f(x>0) = fv(x(x>0));
