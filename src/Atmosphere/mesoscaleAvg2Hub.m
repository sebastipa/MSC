function [ubk, vbk] = mesoscaleAvg2Hub(u, v, abl)

% reconstruction based on average match and angle perturbation (Stipa 20203)
kappa = 0.4;
nx    = size(u,1);
ny    = size(u,2);

ubk   = zeros(nx,ny);
vbk   = zeros(nx,ny);

for i=1:nx
    for j=1:ny

        uAvgMag     = sqrt((abl.U1 + u(i,j))^2+(abl.V1 + v(i,j))^2);
        uStarLocal  = (kappa*(abl.H1 - abl.z0)*uAvgMag)/ ( abl.H1*(log(abl.H1/abl.z0) - 1) + abl.z0 + 4.7/(2*abl.L)*(abl.H1^2 - abl.z0^2));
        phiPrime    = atan2d((abl.V1 + v(i,j)),(abl.U1 + u(i,j))) - atan2d(abl.V1,abl.U1);
        phiZero     = interp1(abl.zVec, abl.uAngle, abl.Href);

        ubk_ij      = uStarLocal/kappa * (log(abl.Href/abl.z0) + 4.7*abl.Href/abl.L)*cosd(phiZero + phiPrime);
        vbk_ij      = uStarLocal/kappa * (log(abl.Href/abl.z0) + 4.7*abl.Href/abl.L)*sind(phiZero + phiPrime);
        
        ubk(i,j)    = ubk_ij;
        vbk(i,j)    = vbk_ij;

    end
end