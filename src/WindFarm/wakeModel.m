function [farm] = wakeModel(farm, sol, abl)

fprintf('Wake Model');
if(sol.localInduction)
    fprintf(' + vortex cylinder')
end
tstart     = tic;

% update wake model runs counter
farm.wakeModelRuns = farm.wakeModelRuns + 1;

% find turbulence intensity at turbine location
% 1: NayafarPorteAgel, 2: from shear stress model, 3: ambient TI
modelTI = 1; % hardcoded to 1 because it is the only relyable one 
farm    = findTurbineTI(farm, abl, modelTI);

% build the linear system to find the turbine inflow
A      = zeros(farm.Ntp,farm.Ntp);
b      = zeros(farm.Ntp,1);

Nqp    = 16;      % number of quadrature points 
w      = 1/Nqp;   % uniform quadrature points weighting 

% loop over wind turbines 
for i=1:farm.Ntp  

    % sheared inflow 
    if(sol.wakeModelLogLaw)
        b(i) = rhsGetInflowWithShear(farm.turbines(i).qpoints, farm.turbines(i).backgroundWindMagP, abl);
    else
        b(i) = farm.turbines(i).backgroundWindMagP;
    end
    
    % compute the result of wake turbine_j on turbine_i
    for j=1:farm.Ntp
              
        % set wind vectors
        x_hat = farm.turbines(j).windDirHat;
        z_hat = [0; 0; 1];
        y_hat = cross(z_hat, x_hat);

        if(dot((farm.turbines(i).rotor - farm.turbines(j).rotor),x_hat)<0)
            continue;
        end
        if(i==j)
            % diagonal terms are 1
            A(i,j) = 1.0;
        else
            % initialize to zero 
            A(i,j) = 0;

            % add contributions from all quadrature points 
            center = farm.turbines(j).rotor;

            for qp=1:Nqp

                % gather quad. point of turbine i
                point_qp   = farm.turbines(i).qpoints(:,qp);

                % gather distance from rot. center of turbine j
                dist_j2qpi = point_qp - center;

                % project in wind coordinates 
                x = dot(dist_j2qpi, x_hat);
                y = dot(dist_j2qpi, y_hat);
                z = dot(dist_j2qpi, z_hat);

                if(farm.tuningWakeModel)
                    coeffs = farm.coeffs;
                    wake   = gaussianWakeModelEvaluate(x, y, z, farm.turbines(j).D, farm.turbines(j).Ct, farm.turbines(j).TI, coeffs);
                else
                    wake = gaussianWakeModelCorrectedEvaluate(x, y, z, farm.turbines(j).D, farm.turbines(j).Ct, farm.turbines(j).TI);
                end
                
                heaviside = x>0;
                a_qp = w * wake * heaviside; 

                A(i,j) = A(i,j) + a_qp;
            end
        end
    end
end

% solve the linear system 
s = A \ b;

% s is the velocity at each turbine position, under the effect of the
% wakes and the varying background flow contained in the right-hand side. 

% now we have to store this velocity, and optionally add local induction
% effects. The algorithm is as follows, for each turbine i
% 1. compute induction from all other turbines except i
% 2. Update the thrust coefficient of turbine i
% 3. Add induction from turbine i 

% Point 2 excludes turbine i for the computation of CT as power curves for
% isolated rotors are usually computed using the "freestream" velocity,
% which does not contain self induction. Then we add the effect of
% induction from turbine i because that is physically what the wind tubrine
% is experiencing.

% store the solution (inflow velocity at each turbine)
for i=1:farm.Ntp

    % initialize induction at turbine i to zero
    i_induction_self = 0;
    i_induction_all  = 0;

    % compute induction from all turbines 
    if(sol.localInduction)

        % induction from all j turbines 
        for j=1:farm.Ntp

			% turbine j center
			center_j = farm.turbines(j).rotor;

            % induction from turbine j
            induction_j = 0;

            for q=1:Nqp

                % gather quad. point of turbine i
                point_qp   = farm.turbines(i).qpoints(:,qp);

			    % compute point to rotor center distance 
			    dist = point_qp - center_j;
			    
			    % set wind vectors (aligned with turbine j which is inducting)
			    x_hat = farm.turbines(j).windDirHat;
			    z_hat = [0; 0; 1];
			    y_hat = cross(z_hat, x_hat);
    
			    % project in wind coordinates 
			    x = dot(dist, x_hat);
			    y = dot(dist, y_hat);
                z = dot(dist, z_hat);
			    
			    % get background velocity at turbine location
                inflow = real(s(j));
			    
			    % only get axial velocity
			    [~, ua] = VortexCylinderEvaluate(farm.turbines(j).D, x, y, z, inflow, farm.turbines(j).Ct);
			    induction_j = induction_j + w*real(ua);
            end

            % add the induction of turbine j to the total induction (keep self induction separate)
            if(j==i || j==(i+farm.Nti))
                i_induction_self = i_induction_self + induction_j;
            else
                i_induction_all  = i_induction_all  + induction_j;
            end

        end
    end

    % update turbine inflow using max-deficit between wake model with
    % variable background and background wind with deep array effect
    if(sol.deepArray)

        % contains all effects but self induction 
        farm.turbines(i).inflow = min(s(i) + i_induction_all, farm.turbines(i).backgroundWindMagPDA);

        % interpolate new CT and CP
        if(farm.useTurbineCurves)

            % make sure interpolation is successful
            if(farm.turbines(i).inflow>farm.curves.wind(end) || farm.turbines(i).inflow<farm.curves.wind(1))

                fprintf(' --> Warning: inflow out of bounds for turbine %d: inflow = %.4f m/s...\n', i, farm.turbines(i).inflow);
                farm.turbines(i).inflow = max(min(farm.turbines(i).inflow, farm.curves.wind(end)), farm.curves.wind(1));
                fprintf(' -->          bounding: inflow = %.4f m/s\n',farm.turbines(i).inflow);
            end

            %fprintf('turbine %d: inflow = %.4f\n', i, farm.turbines(i).inflow);
            farm.turbines(i).Ct      = interp1(farm.curves.wind, farm.curves.ct, farm.turbines(i).inflow);
            farm.turbines(i).Cp      = interp1(farm.curves.wind, farm.curves.cp, farm.turbines(i).inflow);
        end

        farm.turbines(i).uDiskMT = farm.turbines(i).inflow * (1 - (1-sqrt(1-farm.turbines(i).Ct))/2);
        farm.turbines(i).uDiskVC = farm.turbines(i).inflow + i_induction_self;
    else
        % contains all effects but self induction 
        farm.turbines(i).inflow = s(i) + i_induction_all;

        % interpolate new CT and CP
        if(farm.useTurbineCurves)

            % make sure interpolation is successful
            if(farm.turbines(i).inflow>farm.curves.wind(end) || farm.turbines(i).inflow<farm.curves.wind(1))

                fprintf(' --> Warning: inflow out of bounds for turbine %d: inflow = %.4f m/s...\n', i, farm.turbines(i).inflow);
                farm.turbines(i).inflow = max(min(farm.turbines(i).inflow, farm.curves.wind(end)), farm.curves.wind(1));
                fprintf(' -->          bounding: inflow = %.4f m/s\n',farm.turbines(i).inflow);
            end

            farm.turbines(i).Ct      = interp1(farm.curves.wind, farm.curves.ct, farm.turbines(i).inflow);
            farm.turbines(i).Cp      = interp1(farm.curves.wind, farm.curves.cp, farm.turbines(i).inflow);
        end

        farm.turbines(i).uDiskMT = farm.turbines(i).inflow * (1 - (1-sqrt(1-farm.turbines(i).Ct))/2);
        farm.turbines(i).uDiskVC = farm.turbines(i).inflow + i_induction_self;
            
    end

    
    % compute CtPrime if not read from LES 
    if(sol.readTurbineCtPrime)
        % already have CtPrime
        farm.turbines(i).CpPrime = farm.turbines(i).Cp / (1 - (1-sqrt(1-farm.turbines(i).Ct))/2)^3;
    else
        farm.turbines(i).CtPrime = farm.turbines(i).Ct / (1 - (1-sqrt(1-farm.turbines(i).Ct))/2)^2;
        farm.turbines(i).CpPrime = farm.turbines(i).Cp / (1 - (1-sqrt(1-farm.turbines(i).Ct))/2)^3;
    end

    %thrust                  = 0.5 * abl.rho * farm.turbines(i).CtPrime * farm.turbines(i).Ar * farm.turbines(i).uDiskMT^2 * farm.turbines(i).windDirHat;
    %power                   = 0.5 * abl.rho * farm.turbines(i).CpPrime * farm.turbines(i).Ar * farm.turbines(i).uDiskMT^3                              ;
    thrust                   = 0.5 * abl.rho * farm.turbines(i).Ct * farm.turbines(i).Ar * farm.turbines(i).inflow^2 * farm.turbines(i).windDirHat;
    power                    = 0.5 * abl.rho * farm.turbines(i).Cp * farm.turbines(i).Ar * farm.turbines(i).inflow^3                              ;
    farm.turbines(i).T       = norm(thrust,2);
    farm.turbines(i).Tx      = thrust(1);
    farm.turbines(i).Ty      = thrust(2);
    farm.turbines(i).P       = power;

    % save this iteration values 
    farm.turbines(i).iterT(farm.wakeModelRuns) = norm(thrust,2);
    farm.turbines(i).iterP(farm.wakeModelRuns) = power;
end

fprintf(': Elapsed Time %.2f s\n',toc(tstart));

end