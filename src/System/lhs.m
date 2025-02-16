function Ax = lhs(s, meso, abl)

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
rho      = 1.225;

% fields for convolutional terms 
u1_hat   = zeros(nx,ny);
v1_hat   = zeros(nx,ny);
u2_hat   = zeros(nx,ny);
v2_hat   = zeros(nx,ny);

Ax       = zeros(5*nx*ny,1);
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

        A  = [sigma1+nut1,  0,            -fc        ,   0,           1i*kq/rho;...
              fc,           0,            sigma1+nut1,   0,           1i*lj/rho;...
              0,            sigma2+nut2,  0,             -fc,         1i*kq/rho;...
              0,            fc,           0,             sigma2+nut2, 1i*lj/rho;...
              k1*kq,        k2*kq,        k1*lj,         k2*lj        1/rho   ];
        
        Ax(id+1:id+5) = A*s(id+1:id+5);

        % set convolutional term for fft 
%         if(abs(kq)<nx/3 && abs(lj)<ny/3)
%             u1_hat(q,j)   = s(id+1);
%             u2_hat(q,j)   = s(id+2);
%             v1_hat(q,j)   = s(id+3);
%             v2_hat(q,j)   = s(id+4);
%         else
%             u1_hat(q,j)   = 0.0;
%             u2_hat(q,j)   = 0.0;
%             v1_hat(q,j)   = 0.0;
%             v2_hat(q,j)   = 0.0;
%         end

        u1_hat(q,j)   = s(id+1);
        u2_hat(q,j)   = s(id+2);
        v1_hat(q,j)   = s(id+3);
        v2_hat(q,j)   = s(id+4);

        id = id + 5;
    end
end

% form the convolutional term in physical space 
u1     = real(ifft2(ifftshift(u1_hat)));
v1     = real(ifft2(ifftshift(v1_hat)));
u2     = real(ifft2(ifftshift(u2_hat)));
v2     = real(ifft2(ifftshift(v2_hat)));

t1     = zeros(nx,ny);
t2     = zeros(nx,ny);
t3     = zeros(nx,ny);
t4     = zeros(nx,ny);

for q=1:nx
    for j=1:ny
        c_coeff = abl.tau1Mag / (U1^2 + V1^2) / sqrt(U1^2 + V1^2);
        d_coeff = abl.tau2Mag / (DU^2 + DV^2) / sqrt(DU^2 + DV^2);

        C11      = ( 2*U1^2 + V1^2 ) * c_coeff;
        C12      = ( U1 * V1       ) * c_coeff;
        C21      = ( V1 * U1       ) * c_coeff;
        C22      = ( 2*V1^2 + U1^2 ) * c_coeff;
        
        D11      = ( 2*DU^2 + DV^2 ) * d_coeff;
        D12      = ( DU * DV       ) * d_coeff;
        D21      = ( DV * DU       ) * d_coeff;
        D22      = ( 2*DV^2 + DU^2 ) * d_coeff;

        t1(q,j)  = (D11+C11)/H1*u1(q,j) - D11/H1*u2(q,j) + (D12+C12)/H1*v1(q,j) - D12/H1*v2(q,j);
        t2(q,j)  = (D21+C21)/H1*u1(q,j) - D21/H1*u2(q,j) + (D22+C22)/H1*v1(q,j) - D22/H1*v2(q,j);
        t3(q,j)  = - D11/H2*u1(q,j)     + D11/H2*u2(q,j) - D12/H2*v1(q,j)       + D12/H2*v2(q,j);
        t4(q,j)  = - D21/H2*u1(q,j)     + D21/H2*u2(q,j) - D22/H2*v1(q,j)       + D22/H2*v2(q,j);
        
    end
end

% transform back and add to dot product 
t1     = fftshift(fft2(t1));
t2     = fftshift(fft2(t2));
t3     = fftshift(fft2(t3));
t4     = fftshift(fft2(t4));

id     = 0;

for q=1:nx
    for j=1:ny
        Ax(id+1) = Ax(id+1) + t1(q,j);
        Ax(id+2) = Ax(id+2) + t2(q,j);
        Ax(id+3) = Ax(id+3) + t3(q,j);
        Ax(id+4) = Ax(id+4) + t4(q,j);
        id = id + 5;
    end
end
