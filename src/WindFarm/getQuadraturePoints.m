function qpoints = getQuadraturePoints(wind, rotorCenter, D)

    % wind is a vector (not normalized) directed as the wind. We use CROSS
    % TYPE quadrature points, which puts more point horizontally, so that
    % wake interaction is better captured. Moreover, the CROSS TYPE also
    % has points near the rotor center, where the velocity is low, while
    % the STAR TYPE misses that part of the rotor. It would be interesting
    % to try with more points (something like ADM) and see what happes 

    % STAR QUADRATURE POINTS 
    % quadrature point radiuses (odd or even qp)
    %     nq      = 16;
    %     Rodd    = sqrt((3 + sqrt(3))/6) * D / 2;              
    %     Reven   = sqrt((3 - sqrt(3))/6) * D / 2;  
    %     qpoints = zeros(3,nq);
    % 
    %     x_hat   = wind/norm(wind,2);
    %     z_hat   = [0; 0; 1];
    %     y_hat   = cross(z_hat, x_hat);
    % 
    %     % create quadrature points 
    %     for qp=1:nq
    % 
    %         if(mod(qp,2) == 1)
    %             % odd quadrature point
    %             R = Rodd;
    %         elseif(mod(qp,2) == 0)
    %             % even quadrature point
    %             R = Reven;
    %         end
    % 
    %         % azimuth angle of this quadrature point  
    %         azmth_qp = 2*pi*(qp-1)/nq;
    % 
    %         % z of this quadrature point 
    %         z_qp     = R * sin(azmth_qp) * z_hat;
    % 
    %         % y of this quadrature point 
    %         y_qp     = R * cos(azmth_qp) * y_hat;
    % 
    %         % set this quadrature point coordinates 
    %         qpoints(:, qp) = rotorCenter + y_qp + z_qp;
    %     end


    % CROSS QUADRATURE POINTS 
    % quadrature point radiuses (odd or even qp)
    nq      = 16;
    Deff    = 2*D/3;
    Reff    = Deff/2;
    qpoints = zeros(3,nq);

    x_hat   = wind/norm(wind,2);
    z_hat   = [0; 0; 1];
    y_hat   = cross(z_hat, x_hat);

    % create quadrature points 
    for qp=1:nq

        % z of this quadrature point 
        if(qp>8)
            z_qp     = (-Reff  + (qp-9)*Deff/7) * z_hat;
        else
            z_qp     = 0 * z_hat;
        end

        % y of this quadrature point 
        if(qp>8)
            y_qp     = 0 * y_hat;
        else
            y_qp     = (-Reff  + (qp-1)*Deff/7) * y_hat;
        end

        % set this quadrature point coordinates 
        qpoints(:, qp) = rotorCenter + y_qp + z_qp;
    end
end