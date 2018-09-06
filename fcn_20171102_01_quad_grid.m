function[az,inc,weights,sphOrdMax] = fcn_20171102_01_quad_grid(spec,modestr)
%
if nargin<2 || isempty(modestr)
    modestr = 'sh_ord';
end
switch modestr
    case 'sh_ord'
        sphOrdMax = spec;
        if rem(sphOrdMax,1)~=0
            error('SH order must be a whole number')
        end
        n_az = 2*sphOrdMax+2; % +1 is the minimum condition but safer with +2
    case 'az_sep_deg'
        az_sep_deg = spec;
        n_az = round(360/az_sep_deg); %allows sampling up order (naz-1)/2 (15,20,24,30)
        sphOrdMax=floor((n_az-1)/2);

        
end

n_inc = 2*ceil((sphOrdMax)./2)+1; %to achieve same order in inclination need N+1, but force it to be odd 
[inc_v,az_v,quad_v] = sphrharm('g',n_inc,n_az);

% get a full grid
[inc_g,az_g] = ndgrid(inc_v(:),az_v(:));
az_frac = 2*pi/n_az;
% quad_g = repmat(1/(4*pi) * az_frac * quad_v(:),1,length(az_v));
% if sum(quad_g(:))-1 > 1e-14,error('quadrature weights don''t sum to one'),end
quad_g = repmat(az_frac * quad_v(:),1,length(az_v));
if abs(sum(quad_g(:))-4*pi) > 1e-12
    warning('quadrature weights don''t sum to 4pi')
    fprintf('difference is %2.3e\n',abs(sum(quad_g(:))-4*pi))
    pause()
end


az = az_g(:);
inc = inc_g(:);
weights = quad_g(:);