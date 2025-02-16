function meso = updateDeepArrayEffect(meso, abl)
%%
tstart     = tic;
fprintf('Extended top-down Model: ');

nx         = meso.nx;
ny         = meso.ny;
x          = meso.x;
y          = meso.y;

limit      = x(meso.farmBuf.imin) - 500;

[ubak, vbak] = mesoscaleAvg2Hub(meso.u1, meso.v1, abl);
cft        = meso.ctf;

z0         = abl.z0;
H          = abl.H;
Href       = abl.Href;
avgD       = abl.avgD;

L          = abl.L;

ubk_da     = zeros(nx,ny);
vbk_da     = zeros(nx,ny);

uStarFree  = sqrt(meso.ubk.^2 + meso.vbk.^2) * 0.4 / (log(Href/z0) + 4.7*Href/abl.L);

parfor i=1:nx
    for j=1:ny

        xp   = x(i);
        yp   = y(j);
        
        [udeep_p, vdeep_p] = topDownModelPoint(xp, yp, x, y, ubak, vbak, limit, cft, uStarFree, z0, H, L, Href, avgD);

        ubk_da(i,j) = udeep_p;
        vbk_da(i,j) = vdeep_p;
    end
end

meso.ubk_da = ubk_da;
meso.vbk_da = vbk_da;

fprintf('Elapsed Time %.2f s\n', toc(tstart));

end