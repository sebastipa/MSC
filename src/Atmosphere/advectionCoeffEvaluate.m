function [reSigma1, reSigma2] = advectionCoeffEvaluate(U1, V1, U2, V2, k, l)

    % advection coefficients
    A1 = k*U1;
    B1 = l*V1;

    A2 = k*U2;
    B2 = l*V2;

    reSigma1 = (A1 + B1);
    reSigma2 = (A2 + B2);

end