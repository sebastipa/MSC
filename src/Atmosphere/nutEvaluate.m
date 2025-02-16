function [nut1, nut2] = nutEvaluate(nut1, nut2, k, l)

    nut1 = nut1*(k^2 + l^2);
    nut2 = nut2*(k^2 + l^2);  

end