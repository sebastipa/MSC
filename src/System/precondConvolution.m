function [L,U, Amat] = precondConvolution(meso, abl)

% Performs incomplete LU factorization of the linear system matrix
nx = meso.nx;
ny = meso.ny;

% layer velocities
U1 = abl.U1;
V1 = abl.V1;
U2 = abl.U2;
V2 = abl.V2;
U3 = abl.U3;
V3 = abl.V3;
fc = abl.fc;

gPrime   = abl.gPrime;       
N        = abl.N;
H1       = abl.H1;
H2       = abl.H2;
rho      = 1.225;

idmat    = 0;
mtid_i   = zeros(25*nx*ny,1);
mtid_j   = zeros(25*nx*ny,1);
values   = zeros(25*nx*ny,1);
id       = 0;

for q=1:nx
    for j=1:ny

        kq =  meso.k(q);
        lj =  meso.l(j);

        [sigma1, sigma2] = advectionCoeffEvaluate(U1, V1, U2, V2, kq, lj);
        [nut1, nut2]     = nutEvaluate(abl.nut1, abl.nut2, kq, lj);
        W                = -1.0*(kq*U3 + lj*V3);
        K                = gPrime + complexStratificationEvaluate(N, W, kq, lj);
        k1               = 1i*K*H1/sigma1;
        k2               = 1i*K*H2/sigma2;

        A  = [sigma1+nut1,  0,             -fc,           0,                  1i*kq/rho;...
              fc,           0,             sigma1+nut1,   0,                  1i*lj/rho;...
              0,            sigma2+nut2,   0,             -fc,                1i*kq/rho;...
              0,            fc,            0,             sigma2+nut2,        1i*lj/rho;...
              k1*kq,        k2*kq,         k1*lj,         k2*lj,              1/rho   ];
        
        for ii=1:5
            for jj=1:5
                id = id + 1;
                mtid_i(id) = idmat+ii;
                mtid_j(id) = idmat+jj;
                values(id) = A(ii,jj);
            end
        end

        idmat = idmat + 5;
    end
end

Amat   = sparse(mtid_i,mtid_j,values);

[L, U] = ilu(Amat,struct('type','ilutp','udiag',1,'droptol',1e-10));

end