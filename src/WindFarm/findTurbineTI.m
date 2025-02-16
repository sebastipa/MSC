function farm = findTurbineTI(farm, abl, modelTI)

% Nayafar and Porte-Agel
if(modelTI == 1)
   
    I0     = abl.TI;
    
    % Local TI (Niayifar and Portè-Agel (2016)) model
    for i=1:farm.Ntp
        IPlus = 0.0;
        for j=1:farm.Ntp
            if(i~=j)
    
                % set wind vectors
                x_hat = farm.turbines(j).windDirHat;
                z_hat = [0; 0; 1];
                y_hat = cross(z_hat, x_hat);
    
                % TI provided from wake of turbine j
                center_j = farm.turbines(j).rotor;
                
                % get rotor center of turbine i
                center_i = farm.turbines(i).rotor;
    
                % gather distance from rot. center of turbine j
                dist_j2i = center_i - center_j;
    
                % project in wind coordinates 
                x = dot(dist_j2i, x_hat);
                y = dot(dist_j2i, y_hat);
                
                % compute TI shed by turbine j at point i
                induction = 0.5  * (1.0 - sqrt(1.0 - farm.turbines(j).Ct));
                
				% original coefficient 
				% IPlus_j   = 0.73 * induction^0.8325 * I0^0.0325 * (x/farm.turbines(j).D)^(-0.32);
                
                % changed coefficient after tuning with LES
                IPlus_j   = 0.8798 * induction^0.8325 * I0^0.0325 * (x/farm.turbines(j).D)^(-0.32);
                
                Ai        = farm.turbines(i).Ar;
                beta      = 0.5  * (1 + sqrt(1 - farm.turbines(j).Ct)) / (sqrt(1 - farm.turbines(j).Ct));
                epsilon   = 0.2  * sqrt(beta);
                ILocal_j  = sqrt(IPlus_j^2 + I0^2);
                kStar_j   = 0.3837 * ILocal_j + 0.003678;
                sigmaWake = (kStar_j * x) + epsilon * farm.turbines(j).D;
                yAbs      = abs(y) + 0.5 * farm.turbines(i).D;
                % fully waked
                if((yAbs - sigmaWake) < 0.0)
                    Awaked = Ai;
                % non-waked
                elseif((yAbs - sigmaWake) > farm.turbines(i).D)
                    Awaked = 0.0;
                % partial waking
                else
                    RPrime = 0.5 * (yAbs - sigmaWake);
                    Awaked = pi * RPrime^2;
                end
                IPlus  = max(IPlus, heaviside(x) * IPlus_j * Awaked / Ai);
            end
        end
        
        % compute local TI and kStar (wake expansion coeff at i-location)
        ILocal                  = sqrt(IPlus^2 + I0^2);
        farm.turbines(i).TI     = ILocal;
    end
% ambient TI
else
    for i=1:farm.Ntp
        farm.turbines(i).TI    = abl.TI;
    end
end

end