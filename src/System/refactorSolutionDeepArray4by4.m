function meso = refactorSolutionDeepArray4by4(s, meso)

nx       = meso.nx;
ny       = meso.ny;

ubk_hat   = zeros(nx,ny);
vbk_hat   = zeros(nx,ny);

id       = 0;

for i=1:nx
    for j=1:ny
        ubk_hat(i,j)      = s(id+1);
        vbk_hat(i,j)      = s(id+3);
        id                = id + 4;
    end
end

% for i=nx/2+1:nx
%     for j=1:ny
%         ii = mod(-(i-1),nx) + 1;
%         jj = mod(-(j-1),ny) + 1;
%         ubk_hat(i,j)      = ubk_hat(ii,jj)';
%         vbk_hat(i,j)      = vbk_hat(ii,jj)';
%     end
% end

meso.ur1_da    = real(ifft2(ifftshift(ubk_hat)));
meso.vr1_da    = real(ifft2(ifftshift(vbk_hat)));

end