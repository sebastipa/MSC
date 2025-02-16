function meso = refactorSolution5by5(s, meso, abl)

nx       = meso.nx;
ny       = meso.ny;

u1_hat   = zeros(nx,ny);
u2_hat   = zeros(nx,ny);
v1_hat   = zeros(nx,ny);
v2_hat   = zeros(nx,ny);
p_hat    = zeros(nx,ny);
eta_hat  = zeros(nx,ny);

id       = 0;

for i=1:nx
    for j=1:ny
        if(i ~= nx/2 + 1 || j ~= ny/2+1)
            u1_hat(i,j)   = s(id+1);
            u2_hat(i,j)   = s(id+2);
            v1_hat(i,j)   = s(id+3);
            v2_hat(i,j)   = s(id+4);
            p_hat (i,j)   = s(id+5);
            id            = id + 5;
    
            W             = -(meso.k(i)*abl.U3 + meso.l(j)*abl.V3);
            K             = abl.gPrime + complexStratificationEvaluate(abl.N, W, meso.k(i), meso.l(j));
            eta_hat(i,j)  = 1/(abl.rho*K) * p_hat (i,j);
        else
            u1_hat(i,j)   = 0.0;
            u2_hat(i,j)   = 0.0;
            v1_hat(i,j)   = 0.0;
            v2_hat(i,j)   = 0.0;
            p_hat (i,j)   = 0.0;
            eta_hat(i,j)  = 0.0;
        end
    end
end

% for i=nx/2+1:nx
%     for j=1:ny
%         ii = mod(-(i-1),nx) + 1;
%         jj = mod(-(j-1),ny) + 1;
%         u1_hat(i,j)   = u1_hat(ii,jj)';
%         u2_hat(i,j)   = u2_hat(ii,jj)';
%         v1_hat(i,j)   = v1_hat(ii,jj)';
%         v2_hat(i,j)   = v2_hat(ii,jj)';
%         p_hat (i,j)   = p_hat (ii,jj)';
%         eta_hat(i,j)  = eta_hat(ii,jj)';
%     end
% end

meso.u1     = real(ifft2(ifftshift(u1_hat )));
meso.u2     = real(ifft2(ifftshift(u2_hat )));
meso.v1     = real(ifft2(ifftshift(v1_hat )));
meso.v2     = real(ifft2(ifftshift(v2_hat )));
meso.p      = real(ifft2(ifftshift(p_hat  )));
meso.eta    = real(ifft2(ifftshift(eta_hat)));

end