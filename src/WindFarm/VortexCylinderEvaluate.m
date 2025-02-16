function [ur, ua] = VortexCylinderEvaluate(D, x, y, z, U, Ct)

% Evaluate induced velocities from the tangential vorticity of the vortex
% cylinder as defined by Branlard Gaunaa (2014). That is the only vorticity
% component to have an in-plane induction. Other components are discarded.

if(x>=-1 && x <= 0)
    x = -1;
end

r = sqrt(y^2 + z^2);

% this is an image, should exclude that 
if(x>-1 && z >= 0.55*D)
    ur = 0.0;
    ua = 0.0;
% this is inside this turbine cylinder
elseif(x>-1 && r <= 0.55*D)
    ur = 0.0;
    ua = 0.0;
else
    % compute radius
    R     = D / 2;

    % compute tangential vorticity 
    gamma = U * ( sqrt(1 - Ct) - 1);

    % adjust parameters for singularities  
    if abs(x) < 1e-5
        x = x - 1e-5;
    end
    if(abs(r) < 1e-5)
        r = r + signNoZero(r) *  1e-5;
    end
    if(abs(R - r) < 1e-5)
        r = R + signNoZero(R - r) *  1e-5;
    end

    % compute elliptic integrals 
    kSq   = 4 * r * R / ((R + r)^2 + x^2);
    kOrig = 4 * r * R / (R + r)^2;

    % This is my own implementation (superfast)
    [K, E, Pi] = elliptic123(kSq,kOrig);

    % Matlab elliptic integrals are too slow
    %K     = ellipticK(kSq);
    %E     = ellipticE(kSq);
    %Pi    = ellipticPi(kOrig, kSq);

    % compute radial velocity 
    ur = - gamma / (2* pi) * sqrt(R/r) * ...
         ( (2 - kSq) / sqrt(kSq) * K - 2 / sqrt(kSq) * E);

    % compute axial velocity 
    ua =   gamma / 2 * ...
           (...
              (R - r + abs(R - r)) / (2 * abs(R - r)) + ...
              x * sqrt(kSq) / (2 * pi * sqrt(r * R)) * ...
              (K + (R - r) / (R + r) * Pi) ...
           );
end
   
    function s = signNoZero(x)
        if(x > 0)
            s = 1.0;
        elseif (x < 0)
            s = -1.0;
        else
            s = 1.0;
        end
    end
end