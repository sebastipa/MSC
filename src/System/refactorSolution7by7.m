function meso = refactorSolution7by7(s, meso)

nx       = meso.nx;
ny       = meso.ny;

u1_hat   = zeros(nx,ny);
u2_hat   = zeros(nx,ny);
v1_hat   = zeros(nx,ny);
v2_hat   = zeros(nx,ny);
p_hat    = zeros(nx,ny);
eta1_hat    = zeros(nx,ny);
eta2_hat    = zeros(nx,ny);

id       = 0;

for i=1:nx
    for j=1:ny
        u1_hat(i,j)      = s(id+1);
        u2_hat(i,j)      = s(id+2);
        v1_hat(i,j)      = s(id+3);
        v2_hat(i,j)      = s(id+4);
        p_hat (i,j)      = s(id+5);
        eta1_hat (i,j)   = s(id+6);
        eta2_hat (i,j)   = s(id+7);
        id               = id + 7;
    end
end

meso.u1     = real(ifft2(ifftshift(u1_hat)));
meso.u2     = real(ifft2(ifftshift(u2_hat)));
meso.v1     = real(ifft2(ifftshift(v1_hat)));
meso.v2     = real(ifft2(ifftshift(v2_hat)));
meso.p      = real(ifft2(ifftshift(p_hat)));
meso.eta    = real(ifft2(ifftshift(eta1_hat))) + ...
              real(ifft2(ifftshift(eta2_hat)));

end