function y = mygamma(x)

z = x;

c = [   1.000000000000000174663;
        5716.400188274341379136;
        -14815.30426768413909044;
        14291.49277657478554025;
        -6348.160217641458813289;
        1301.608286058321874105;
        -108.1767053514369634679;
        2.605696505611755827729;
       -0.7423452510201416151527e-2;
        0.5384136432509564062961e-7;
       -0.4023533141268236372067e-8];
   
if(real(x) < 0.0)
    z = -z;
end

g = 9;
t = z + g;
s = 0;
for k = g+2:-1:2
    s = s + c(k) / t;
    t = t - 1.0;
end

s  = s + c(1);
ss = (z + g - 0.5);
s  = log(s * sqrt(pi+pi)) + (z - 0.5) * log(ss) - ss;
f = exp(s);

% recursive formula    
if(real(x) < 0.0)
    f = - pi / (x * f * sin(pi*x));
end

% negative integers 
if(round(x) == x && imag(x) == 0 && real(x)  <= 0)
    f = Inf;
end

% exact results for integer arguments 
if(round(x) == x && imag(x) == 0 && real(x) > 0)
    pp   = 1;
    for zint=1:(x-1)
        pp = pp*zint;
    end
    f = pp;
end

y = f;

end