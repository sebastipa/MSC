function [as, bs, cs, bf, cf, af] = findSuperGaussianCoeffsCathelain2(Ct, TI, D)

as = 0.18;
bs = 0.0119;
cs = 0.0564*Ct + 0.13;
bf = 1.59*exp(-23.31*TI)-2.15;
cf = 2.98;

a  = 0.5 * (1-sqrt(1 - Ct));

% bisection method to find af
maxiter = 50; 
tol     = 1e-5;
iter    = 0;

af_a      = 2;
af_b      = 7;
af        = 0.5*(af_a+af_b);

val      = superGaussianAfEval(af, cf, cs, Ct, D)   - a;
val_a    = superGaussianAfEval(af_a, cf, cs, Ct, D) - a;

while(abs(val) > tol && iter < maxiter)
    if sign(val) == sign(val_a) 
        af_a = af;
    else
        af_b = af;
    end

    af        =  0.5*(af_a+af_b);

    val       = superGaussianAfEval(af, cf, cs, Ct, D)   - a;
    val_a     = superGaussianAfEval(af_a, cf, cs, Ct, D) - a;
  
    iter = iter + 1;

end

if(iter == maxiter)
    fprintf('--> warning: maxiter reached in supergaussian wake model when computing af\n');
end

end