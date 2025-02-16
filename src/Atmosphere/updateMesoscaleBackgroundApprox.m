function meso = updateMesoscaleBackgroundApprox(meso, sol, abl)

fprintf('3LM-Reconstruction: ');
tstart = tic;

% Performs incomplete LU factorization of the linear system matrix
nx = meso.nx;
ny = meso.ny;

% layer velocities
U1 = abl.U1;
V1 = abl.V1;
U2 = abl.U2;
V2 = abl.V2;
fc = 0;%abl.fc; % it is not necessary here to include fc

H1       = abl.H1;
c_coeff  = abl.tauMag1 / (U1^2 + V1^2) / sqrt(U1^2 + V1^2);
C11      = ( 2*U1^2 + V1^2 ) * c_coeff;
C12      = ( U1 * V1       ) * c_coeff;
C21      = ( V1 * U1       ) * c_coeff;
C22      = ( 2*V1^2 + U1^2 ) * c_coeff;

rho      = abl.rho;

p_hat    = fftshift(fft2(meso.p));

for q=1:nx
    for j=1:ny

        kq =  meso.k(q);
        lj =  meso.l(j);

        if(q ~= nx/2 + 1 || j ~= ny/2+1)

            [nut1, ~]     = nutEvaluate(abl.nut1, abl.nut2, kq, lj);
            [reSigma1, ~] = advectionCoeffEvaluate(U1, V1, U2, V2, kq, lj);
            sigma1        = 1i*reSigma1;
    
            bx            = - 1/rho * 1i * kq *  p_hat(q,j);
            by            = - 1/rho * 1i * lj *  p_hat(q,j);
    
            A             = sigma1 + nut1 + C11/H1;
            B             = sigma1 + nut1 + C22/H1;
            C             = C12/H1 - fc;
            D             = C21/H1 + fc;
            
            meso.ur1(q,j) = (B*bx - C*by) / (A*B - C*D);
            meso.vr1(q,j) = (A*by - D*bx) / (A*B - C*D);
        else
            meso.ur1(q,j) = 0.0;
            meso.vr1(q,j) = 0.0;
        end

    end
end

% hi-pass filter 
% for q=1:nx
%     stddev = 1/20000;
%     f = 1.0 - exp(-0.5*meso.k(q)^2/stddev^2);
%     meso.ur1(q,:) = meso.ur1(q,:)*f;
%     meso.vr1(q,:) = meso.vr1(q,:)*f;
% end

meso.ur1    = real(ifft2(ifftshift(meso.ur1)));
meso.vr1    = real(ifft2(ifftshift(meso.vr1)));

% update background hub velocity with pressure effect
[ubk, vbk] = mesoscaleAvg2Hub(meso.ur1, meso.vr1, abl);
meso.ubk   = ubk; 
meso.vbk   = vbk; 

fprintf('Elapsed Time %.2f s\n',toc(tstart));

end