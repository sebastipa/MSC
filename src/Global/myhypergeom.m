function y = myhypergeom(a,b,c,z)

    % a,b,c can be complex 
    % z has to be real or complex and >= 1
    
    if(round(c)==c && imag(c) == 0 && real(c) <= 0)
        y = Inf;
        return;
    end

    epsilon = 1e-15;

    term        = 1;
    result      = 1;
    termination = 0;

    for k = 0:30000
        term = term * ((a + k) * (b + k) * z / ((c + k) * (k + 1.0)));
        if(abs(term) < abs(result * epsilon)) 
            termination = termination + 1;
            if(termination >= 3) 
                y = result;
                return;
            end
        else
            termination = 0;
        end
        if(term == 0)
            break;
        end
        result =  result + term;
        if(isnan(result))
            y = NaN;
            return;
        end
    end
    
    y = result;

end