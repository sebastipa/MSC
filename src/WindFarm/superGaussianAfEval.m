function val = superGaussianAfEval(af_p, cf, cs, Ct, D)
    n_0        = af_p + cf;
    beta       = 0.5 * (1 + sqrt(1 - Ct)) / (sqrt(1 - Ct));
    sigma_0    = cs*sqrt(beta)*D;
    sigmaByD_0 = sigma_0/D;
    val        =  2^(2/n_0-1) - sqrt( 2^(4/n_0-2) - (n_0*Ct)/(16*real(mygamma(2/n_0))*sigmaByD_0^(4/n_0)) ); 
end