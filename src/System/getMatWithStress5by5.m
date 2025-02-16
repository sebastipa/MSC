function Amat = getMatWithStress5by5(meso, abl)

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
DU = U2 - U1;
DV = V2 - V1;
fc = abl.fc;

gPrime   = abl.gPrime;       
N        = abl.N;
H1       = abl.H1;
H2       = abl.H2;
rho      = abl.rho;

c_coeff = abl.tauMag1 / (U1^2 + V1^2) / sqrt(U1^2 + V1^2);
d_coeff = abl.tauMag2 / (DU^2 + DV^2) / sqrt(DU^2 + DV^2);

C11      = ( 2*U1^2 + V1^2 ) * c_coeff;
C12      = ( U1 * V1       ) * c_coeff;
C21      = ( V1 * U1       ) * c_coeff;
C22      = ( 2*V1^2 + U1^2 ) * c_coeff;

D11      = ( 2*DU^2 + DV^2 ) * d_coeff;
D12      = ( DU * DV       ) * d_coeff;
D21      = ( DV * DU       ) * d_coeff;
D22      = ( 2*DV^2 + DU^2 ) * d_coeff;

idmat    = 0;
mtid_i   = zeros(25*(nx*ny-1),1);
mtid_j   = zeros(25*(nx*ny-1),1);
values   = zeros(25*(nx*ny-1),1);
id       = 0;

for q=1:nx
    for j=1:ny
        if(q ~= nx/2+1 || j ~= ny/2+1)
            kq =  meso.k(q);
            lj =  meso.l(j);
    
            [reSigma1, reSigma2] = advectionCoeffEvaluate(U1, V1, U2, V2, kq, lj);
            sigma1           = 1i*reSigma1;
            sigma2           = 1i*reSigma2;
            [nut1, nut2]     = nutEvaluate(abl.nut1, abl.nut2, kq, lj);
            W                = -(kq*U3 + lj*V3);
            K                = gPrime + complexStratificationEvaluate(N, W, kq, lj);
            k1               = K*H1/reSigma1;
            k2               = K*H2/reSigma2;
    
            A  = [sigma1+nut1+(D11 + C11)/H1,  -D11/H1,             -fc+(D12+C12)/H1,           -D12/H1 ,                  1i*kq/rho;...
                  fc+(D21+C21)/H1,             -D21/H1,             sigma1+nut1+(D22+C22)/H1,   -D22/H1,                   1i*lj/rho;...
                  -D11/H2,                     sigma2+nut2+D11/H2,  -D12/H2 ,                   -fc+D12/H2,                1i*kq/rho;...
                  -D21/H2,                     fc+D21/H2,           -D22/H2,                    sigma2+nut2+D22/H2,        1i*lj/rho;...
                  k1*kq,                       k2*kq,               k1*lj,                      k2*lj,                     1/rho   ];
            
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
end

Amat   = sparse(mtid_i,mtid_j,values);

%[L, U] = ilu(Amat,struct('type','ilutp','udiag',1,'droptol',1e-10));

end