function[Dmat] = fcn_20170616_02_myWignerDMatrix_Jacobi(sphOrd,alpha,beta,gamma)

Dmat = zeros((sphOrd+1)^2);
idc=0;
for l = 0:sphOrd
    D_l = fcn_20170616_01_myWignerD_Jacobi(l,alpha,beta,gamma);
    idc = idc(end)+(1:size(D_l,1));
    Dmat(idc,idc) = D_l;
end