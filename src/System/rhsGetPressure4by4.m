function b = rhsGetPressure4by4(meso, abl)

% discretization
nx = meso.nx;
ny = meso.ny;

p_hat  = fftshift(fft2(meso.p));

b       = zeros(4*nx*ny,1);
id      = 0;
for q=1:nx
    for j=1:ny

        kq =  meso.k(q);
        lj =  meso.l(j);

        b(id+1) = - 1i*kq/abl.rho * p_hat(q,j);
        b(id+2) = - 1i*lj/abl.rho * p_hat(q,j);
        b(id+3) = - 1i*kq/abl.rho * p_hat(q,j);
        b(id+4) = - 1i*lj/abl.rho * p_hat(q,j);
        id = id + 4;
    end
end

end