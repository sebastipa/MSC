function b = rhsGetStressPerturb5by5(s, meso, abl)

% adds the shear stress convolutional term 

% domain data 
nx = meso.nx;
ny = meso.ny;

% layer velocities
U1 = abl.U1;
V1 = abl.V1;
U2 = abl.U2;
V2 = abl.V2;
DU = U2 - U1;
DV = V2 - V1;
H1 = abl.H1;
H2 = abl.H2;

% shear stress jacobian
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

% fields for convolutional terms 
u1_hat   = zeros(nx,ny);
v1_hat   = zeros(nx,ny);
u2_hat   = zeros(nx,ny);
v2_hat   = zeros(nx,ny);

id       = 0;
for q=1:nx
    for j=1:ny

        u1_hat(q,j)   = s(id+1);
        u2_hat(q,j)   = s(id+2);
        v1_hat(q,j)   = s(id+3);
        v2_hat(q,j)   = s(id+4);

        id = id + 5;
    end
end

% transform in physical space 
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

        % only add if perturbation differs from zero 
        if( abs(meso.taux1(q,j)) > 1e-10 )

            t1(q,j)  = (D11+C11)/H1*u1(q,j) - D11/H1*u2(q,j) + (D12+C12)/H1*v1(q,j) - D12/H1*v2(q,j);
            t2(q,j)  = (D21+C21)/H1*u1(q,j) - D21/H1*u2(q,j) + (D22+C22)/H1*v1(q,j) - D22/H1*v2(q,j);
            t3(q,j)  = - D11/H2*u1(q,j)     + D11/H2*u2(q,j) - D12/H2*v1(q,j)       + D12/H2*v2(q,j);
            t4(q,j)  = - D21/H2*u1(q,j)     + D21/H2*u2(q,j) - D22/H2*v1(q,j)       + D22/H2*v2(q,j);
    
            epsx1     =  (meso.taux2(q,j) - meso.taux1(q,j)) / H1;
            epsy1     =  (meso.tauy2(q,j) - meso.tauy1(q,j)) / H1;
            epsx2     = - meso.taux2(q,j) / H2;
            epsy2     = - meso.tauy2(q,j) / H2;
    
            t1(q,j)  = epsx1 - t1(q,j);
            t2(q,j)  = epsx2 - t2(q,j);
            t3(q,j)  = epsy1 - t3(q,j);
            t4(q,j)  = epsy2 - t4(q,j);
        end
    end
end

% transform back 
t1     = fftshift(fft2(t1));
t2     = fftshift(fft2(t2));
t3     = fftshift(fft2(t3));
t4     = fftshift(fft2(t4));

b      = zeros(5*nx*ny,1);
id     = 0;
for q=1:nx
    for j=1:ny
        b(id+1) = t1(q,j);
        b(id+2) = t2(q,j);
        b(id+3) = t3(q,j);
        b(id+4) = t4(q,j);
        id = id + 5;
    end
end

end