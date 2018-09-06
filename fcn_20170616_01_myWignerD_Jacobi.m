function[D,d] = fcn_20170616_01_myWignerD_Jacobi(n,alpha,beta,gamma)

%n: spherical harmonic order
%m: spherical harmonic degree of transformed variable (i.e. output)
%mp: spherical harmonic degree of original variable (i.e. input)

fact = factorial((0:3*n).');      %[2*n+1 1]; - use index i+1 to get factorial(i)

m = -n:n;
mp = -n:n;

[mp_grid,m_grid] = ndgrid(mp,m);

mu_grid = abs(mp_grid-m_grid);
nu_grid = abs(mp_grid+m_grid);
s_grid = n-(mu_grid+nu_grid)/2;

cat_grid = [mu_grid(:) nu_grid(:) s_grid(:)];
[cat_grid_can, ~,i_can] = unique(cat_grid,'rows');

cosb = cos(beta);
cosb2 = cos(beta/2);
sinb2 = sin(beta/2);

% populate Jacobi polynomial values
%P = jacobiP(cat_grid_can(:,3),cat_grid_can(:,1),cat_grid_can(:,2),cosb);
P = fcn_20170704_01_jacobiP(cat_grid_can(:,3),cat_grid_can(:,1),cat_grid_can(:,2),cosb);

P = reshape(P(i_can),size(mp_grid));
zeta_mpm = ones(size(mp_grid));
idc = find(m_grid<mp_grid);
zeta_mpm(idc) = (-1).^(m_grid(idc)-mp_grid(idc));
d = zeta_mpm .* sqrt( ...
        (fact(s_grid + 1) .* fact(s_grid + mu_grid + nu_grid + 1)) ./ ...
        (fact(s_grid + mu_grid + 1) .* fact(s_grid + nu_grid + 1)) ...
    ) .* ...
    sinb2.^mu_grid .* cosb2.^nu_grid .* P;
D = exp(-1i * mp_grid * alpha) .* d .* exp(-1i * m_grid * gamma);

        