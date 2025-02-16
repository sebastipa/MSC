function meso = updateBodyForce(meso, farm, abl)

Tx_tot     = 0;
Tx_int_tot = 0;

tstart     = tic;

meso.fx = zeros(meso.nx, meso.ny);
meso.fy = zeros(meso.nx, meso.ny);

% loop over wind turbines and cumulate body force at this mesh
% divide by density to be consistent 
for t=1:farm.Nti

    thrust = [farm.turbines(t).Tx; farm.turbines(t).Ty; 0] / abl.rho;
    
    % add to total thrust
    Tx_tot = Tx_tot + thrust(1);

    imin   = farm.turbines(t).project.imin;
    imax   = farm.turbines(t).project.imax;
    jmin   = farm.turbines(t).project.jmin;
    jmax   = farm.turbines(t).project.jmax;
    rotor  = farm.turbines(t).rotor;

    for i=imin:imax
        for j=jmin:jmax

            % coordinates of 2D the mesh point with z = Href
            point_p = [meso.x(i); meso.y(j); farm.turbines(t).zHub];

            % compute distance 
            dist    = norm(point_p - rotor,2);
                
            % projection function (gaussian)
            g       = 1 / (pi * farm.Lbf^2) * exp( - dist^2 / farm.Lbf^2);

            meso.fx(i,j) = meso.fx(i,j) - g*thrust(1);
            meso.fy(i,j) = meso.fy(i,j) - g*thrust(2);

            % error checking 
            Tx_int_tot = Tx_int_tot + g*thrust(1)*meso.dx*meso.dy;
        end
    end
    
end

meso.f = sqrt(meso.fx.^2 + meso.fy.^2);

fprintf('Force Projection: error on thrust %.4f%%, Elapsed Time %.2f s\n',abs(Tx_tot-Tx_int_tot)/Tx_tot*100, toc(tstart));

end