af = -100:100;

Ct = 0.77;
TI = 0.05;
cs = 0.0564*Ct + 0.13;
cf = 2.98;
D  = 126;

a  = 0.5 * (1-sqrt(1 - Ct));

for i=1:length(af)
    val(i) = superGaussianAfEval(af(i), cf, cs, Ct, D) - a;
end

plot(af, val);