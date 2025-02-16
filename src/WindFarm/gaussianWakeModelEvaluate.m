function wake = gaussianWakeModelEvaluate(x, y, z, D, Ct, TI, coeffs)

if(heaviside(x)>0)

    if(nargin < 7)
        % coefficients from Bastankhah and Porte Agel
        [as, bs, cs] = findGaussianCoeffsPorteAgel();
    else
        % coefficients from Bastankhah and Porte Agel
        [as, bs, cs] = findGaussianCoeffsPorteAgel();

        % provided from tuning optimization algorithm
        as = coeffs(1);
        bs = coeffs(2);
    end

    % fixes the problem closer than 2D (PorteAgel)
    Ct = Ct*(1 + erf(x/D))/2;

    % wake width at x = 0
    beta      = 0.5 * (1 + sqrt(1 - Ct)) / (sqrt(1 - Ct));

    % wake width at current x
    sigmaByD  = (as*TI + bs) * x / D + cs * sqrt(beta);

    % gaussian coefficient
    alpha     = 1 - sqrt(max(1 - Ct / (8 * sigmaByD^2),0));
    wake      = alpha * exp( - 0.5 / sigmaByD^2 * (z^2 + y^2) / D^2);
else
    wake = 0;
end

end