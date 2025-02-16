function wake = superGaussianWakeModelEvaluate(x, y, z, D, Ct, TI, coeffs)

if(heaviside(x)>0)

    if(nargin < 7)
        % coefficients from the second model tuning from Cathelain
        %[as, bs, cs, bf, cf, af] = findSuperGaussianCoeffsCathelain2(Ct, TI, D);

        % coefficients from the first model tuning from Cathelain
        [as, bs, cs, bf, cf, af] = findSuperGaussianCoeffsCathelain1();

        % coefficients from Bastankhah and Porte Agel
        %[as, bs, cs] = findGaussianCoeffsPorteAgel();
    else
        % coefficients from the second model tuning from Cathelain
        [as, bs, cs, bf, cf, af] = findSuperGaussianCoeffsCathelain2(Ct, TI, D);

        % coefficients from the first model tuning from Cathelain
        %[as, bs, cs, bf, cf, af] = findSuperGaussianCoeffsCathelain1();

        % coefficients from Bastankhah and Porte Agel
        %[as, bs, cs] = findGaussianCoeffsPorteAgel();
        
        % provided from tuning optimization algorithm
        as = coeffs(1);
        bs = coeffs(2);
    end

    n          = af*exp(bf*x/D) + cf;
    beta       = 0.5 * (1 + sqrt(1 - Ct)) / (sqrt(1 - Ct));
    
    % wake expansion proposed by Cathelain
    sigmaByD   = (as*TI + bs)*x/D + cs*sqrt(beta);

    alpha      =  2^(2/n-1) - sqrt( 2^(4/n-2) - (n*Ct)/(16*real(mygamma(2/n))*sigmaByD^(4/n)) ); 

    r          = sqrt(z^2 + y^2);
    wake       = alpha * exp( - 0.5 * r^n / sigmaByD^2 / D^2 * D^(2-n));
else
    wake = 0;
end

end