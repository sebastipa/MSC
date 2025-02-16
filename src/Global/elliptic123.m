function [K, E, Pi] = elliptic123(m,n)

    if((abs(m-1.0)) < 1e-10)
        K = 5;
        E = 1.0;
        Pi = 5.0;
    else
        k = sqrt(1-m);

        M = agm(k);
        K = real(pi / (2.0 * M));

        dx  = 0.001;
        k_p = sqrt(1-m+dx);
        k_m = sqrt(1-m-dx);
        dh  = 2*dx;

        M_p = agm(k_p);
        M_m = agm(k_m);
        K_p = real(pi / (2.0 * M_p));
        K_m = real(pi / (2.0 * M_m));
        dKdm = -(K_p - K_m) / dh;

        E = real((1.0 - m)*(2.0*m*dKdm + K));
        
        %fprintf('rf/rj = %f, %f\n',rf,rj);
        
        if(abs(n) > 1e-5)
            if(abs(m) < 1e-5)
                m = m + 0.01;
            end
            rf = RF(0,1.0-m,1);
            rj = RJ(0,1.0-m,1,1.0-n);
            Pi =  real(rf+ 1.0/3.0*n*rj);
        else
            Pi = K;
        end

    end

    % arithmetic geometric mean M(1,k)
    function M = agm(k) 
        a_old = 1.0;
        b_old = k;
        
        while (abs(a_old-b_old) > 1e-5)
            a = (a_old+b_old)/2.0;
            b = sqrt(a_old*b_old);
            
            a_old = a;
            b_old = b;
        end
        
        M = (a_old + b_old) / 2.0;
    end

    % Carlson Rf integral reduction algorithm
    function v = RF(x,y,z)
        maxiter = 100;
        x_old = x;
        y_old = y;
        z_old = z;
        r     = 1e-10;
        
        A0    = (x+y+z)/3.0;
        A     = A0;
        Q     = (3.0*r)^(-1/6.0) * max(max(abs(A0-z),abs(A0-y)),abs(A0-x));
        
        it    = 0;
        res   = 4.0 ^ (-it)*Q;
        
        while(res > abs(A) && it < maxiter)
            lambda = sqrt(x_old)*sqrt(y_old) + sqrt(x_old)*sqrt(z_old) + sqrt(y_old)*sqrt(z_old);
            x_new  = 0.25 * (x_old + lambda);
            y_new  = 0.25 * (y_old + lambda);
            z_new  = 0.25 * (z_old + lambda);
            A      = 0.25 * (A + lambda);
            it     = it+1;
            res    = 4.0 ^ (-it)*Q;
            x_old  = x_new;
            y_old  = y_new;
            z_old  = z_new;
        end
        
        %fprintf('RF ended at iter %d\n',it);
        
        X = (A0 - x) / (4.0^it * A);
        Y = (A0 - y) / (4.0^it * A);
        Z = - X - Y;
        
        E2 = X*Y - Z^2;
        E3 = X*Y*Z;
        
        v  = (1.0 / sqrt(A)) * (1.0 - 1.0/10.0*E2 + 1.0/14.0*E3 + 1.0/24.0*E2^2 - 3.0/44.0*E2*E3);

    end

    % Carlson Rc integral reduction algorithm
    function v = RC(x,y)
        maxiter = 100;
        x_old = x;
        y_old = y;
        r     = 1e-10;
        
        A0    = (x+2.0*y)/3.0;
        A     = A0;
        Q     = (3.0*r)^(-1.0/8.0) * abs(A0-x);
        
        it    = 0;
        res   = 4.0 ^ (-it)*Q;
        
        while(res > abs(A) && it < maxiter)
            lambda = 2.0*sqrt(x_old)*sqrt(y_old) + y_old;
            x_new  = 0.25 * (x_old + lambda);
            y_new  = 0.25 * (y_old + lambda);
            A      = 0.25 * (A + lambda);
            it     = it+1;
            res    = 4.0 ^ (-it)*Q;
            x_old  = x_new;
            y_old  = y_new;
        end
        
        %fprintf('RC ended at iter %d\n',it);
        
        s = (y - A0) / (4^it * A);
        
        v  = (1.0 / sqrt(A)) * (1.0 + 3.0/10.0*s^2 + 1.0/7.0*s^3 + 3.0/8.0*s^4 + 9.0/22.0/s^5 + 159.0/208.0*s^6 + 9.0/8.0*s^7);

    end

    % Carlson Rj integral reduction algorithm
    function v = RJ(x,y,z,p)
        maxiter = 100;
        x_old = x;
        y_old = y;
        z_old = z;
        p_old = p;
        r     = 1e-10;
        
        A0    = (x+y+z+2.0*p)/5.0;
        delta = (p-x)*(p-y)*(p-z);
        A     = A0;
        Q     = (r/4.0)^(-1.0/6.0) * max(max(max(abs(A0-p),abs(A0-y)),abs(A0-z)),abs(A0-x));
        
        it    = 0;
        res   = 4.0 ^ (-it)*Q;
        
        RCsum = 0.0;
        
        while(res > abs(A) && it < maxiter)
            lambda = sqrt(x_old)*sqrt(y_old) + sqrt(x_old)*sqrt(z_old) + sqrt(y_old)*sqrt(z_old);
            x_new  = 0.25 * (x_old + lambda);
            y_new  = 0.25 * (y_old + lambda);
            z_new  = 0.25 * (z_old + lambda);
            p_new  = 0.25 * (p_old + lambda);
            A      = 0.25 * (A + lambda);
            d      = (sqrt(p_old)+sqrt(x_old))*(sqrt(p_old)+sqrt(y_old))*(sqrt(p_old)+sqrt(z_old));
            e      = (4.0^(-3.0*it) * delta) / d^2.0;
            
            % using RF degenerate form works
            RCsum  = RCsum + (4.0^(-it))/d * RF(1,1+e,1+e);
            
            % using RC does not work: produces very high numbers
            % RCsum  = RCsum + (4.0^(-it))/d * RC(1,1+e);
            
            it     = it+1;
            res    = 4.0 ^ (-it)*Q;
            x_old  = x_new;
            y_old  = y_new;
            z_old  = z_new;
            p_old  = p_new;
        end
        
        %fprintf('RJ ended at iter %d\n',it);
        
        X = (A0 - x) / (4.0^it * A);
        Y = (A0 - y) / (4.0^it * A);
        Z = (A0 - z) / (4.0^it * A);
        P = (- X - Y - Z) / 2.0;
        
        E2 = X*Y + X*Z + Y*Z - 3.0*P^2.0;
        E3 = X*Y*Z + 2.0 * E2 * P + 4.0*P^3;
        E4 = (2.0 * X*Y*Z + E2*P + 3.0*P^3) * P;
        E5 = X*Y*Z*P^2.0;
        
        v  = (4.0^(-it) * A^(-3.0/2.0)) * (1.0 - 3.0/14.0*E2 + 1.0/6.0*E3 + 9.0/88.0*E2^2 - 3.0/22.0*E4 - 9.0/52.0*E2*E3 + 3.0/26.0*E5) +...
             6.0 * RCsum;

    end
    
end