function b = rhsGetTurbines5by5(meso, abl)

% discretization
nx = meso.nx;
ny = meso.ny;

H1       = abl.H1;

fx_hat  = fftshift(fft2(meso.fx));
fy_hat  = fftshift(fft2(meso.fy));

b       = zeros(5*(nx*ny-1),1);
id      = 0;
for q=1:nx
    for j=1:ny
        if(q ~= nx/2 + 1 || j ~= ny/2+1)
            b(id+1) = fx_hat(q,j) / H1;
            b(id+2) = fy_hat(q,j) / H1;
            id = id + 5;
        end
    end
end

end