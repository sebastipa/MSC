function [as, bs, cs] = findGaussianCoeffsPorteAgel()


% Stipa tuning on LES no shear stress model (only as and bs)
as = 0.3678;
bs = 0.0028;

% Stipa tuning on LES with shear stress model (only as and bs)
as = 0.9708;
bs = -0.0520;

% Stipa tuning on LES with shear stress model only for first 3LM solution (only as)
as = 0.3731;
bs = 0.0028;

% Niayafar Porte Agel 
as = 0.3837;
bs = 0.003678;
cs = 0.2;

end