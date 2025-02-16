function bi = rhsGetInflowWithShear(qpoints, hubWind, abl)

% this function computes the average wind speed among the turbine
% quadrature points, considering a logarithmic inflow profile

% 16 quadrature points and uniform weights always used 
Nqp    = 16;      
w      = 1/Nqp; 
bi     = 0;

zRef   = abl.Href;
z0     = abl.z0;
L      = abl.L;

for qp=1:Nqp

    % get quad. point z component (use abs so that is also holds for mirrored turbines)
    z = abs(qpoints(3,qp));

    % get local uStar
    uStar = 0.4 * hubWind / (log(zRef/z0) + 4.7*zRef/L);

    % get velocity interpolating log law
    u = uStar/0.4*(log(z/z0) + 4.7*zRef/L);

    % add to right hand side 
    bi = bi + w*u;
end

end