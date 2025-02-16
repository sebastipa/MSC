function psi = complexStratificationEvaluate(N, Omega, k, l)

    epsilon = eps;

    % make sure that sign(Omega) != 0
    if(abs(Omega) <= epsilon)
        if(Omega >= 0)
            Omega = epsilon;
        else 
            Omega = -epsilon;
        end
    end

    Nsq      = N^2;
    Osq      = Omega^2;
    wvNmbrSq = (k^2 + l^2);

    if(Osq<Nsq)
        % vertically propagating: sign(m) = - sign(Omega) for radiation condition
        m = - sign(Omega) * sqrt(max(wvNmbrSq, eps)*(Nsq/Osq - 1));

        % complex stratification coeff. 
        psi   = (1i * (Nsq - Osq)) / m;

    elseif(Osq>Nsq)
        % evanescent: positive complex root is choosen for sol. boundness
        m = sqrt(max(wvNmbrSq, eps)*(Nsq/Osq - 1)); 

        % complex stratification coeff. 
        psi   = (1i * (Nsq - Osq)) / m;

    else   
        psi = 0.0;
    end

end