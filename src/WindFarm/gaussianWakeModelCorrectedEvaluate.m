function wake = gaussianWakeModelCorrectedEvaluate(x, y, z, D, Ct, TI)

if(x>0)


    % coefficients from Bastankhah and Porte Agel
    [as, bs, cs] = findGaussianCoeffsPorteAgel();

    % beta coefficient
    beta        = 0.5 * (1 + sqrt(1 - Ct)) / (sqrt(1 - Ct));

    % merging parameters: merging is done at 2D and has a buffer of 4D
    mergeCenter = 2*D;
    delta       = 4*D;

    % near wake weight (supergaussian)
    nNW          = 2.0*exp(-0.68*x/D) + 2;
    wNW          = 0.5*(1-tanh(7*(x-mergeCenter)/delta));

    % far wake weight (gaussian)
    nFW          = 2;
    wFW          = 0.5*(1+tanh(7*(x-mergeCenter)/delta));

    % merge exponent smoothly
    n            = wNW*nNW + wFW*nFW;
    
    % wake expansion
    sigmaByD     = (as*TI + bs)*x/D + cs*sqrt(beta);

    % deficit function 
    alpha      =  2^(2/n-1) - sqrt( 2^(4/n-2) - (n*Ct)/(16*real(mygamma(2/n))*sigmaByD^(4/n)) ); 

    % wake function
    r          = sqrt(z^2 + y^2);
    wake       = real(alpha * exp( - 0.5 * r^n / sigmaByD^2 / D^2 * D^(2-n)));
else
    wake = 0;
end

end