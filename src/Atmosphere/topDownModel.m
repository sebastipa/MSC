function [ ...
u,         ...   % velocity mag profile
uStar_hi,  ...   % exit friction vel in high stress layer 
uStar_lo,  ...   % exit friction vel in low stress layer 
deltaFarm, ...   % IBL
deltaWake, ...   % WBL
z0_hi      ...   % wind farm roughness
         ] ...
         = topDownModel(uStar, z0, H, L, cft, zh, D, farmDist, wakeDist)


% implementation of the Meneveau top-down model + wake treatment 
% 1. extra uStar and z0 are required inside the turbine layer, otherwise it 
%    doesn't exist a set of uStar_lo, uStar_hi and z0_hi such that the
%    profile is continuous at H and zh and shows consistent z0_hi decrease
% 2. this is justified by the fact that the WBL grows, hence the upper
%    layer slowly disappears and the wind farm layer recovers while it
%    thinkens. The behavior if uStar_md is such that it is low inside the
%    wind farm (high velocity deficit), while it increases during the WB
%    growth. At this point uStar_md mimics the effect of uStar_hi (as the 
%    upper layer disappeared) hence it starts decreasing. 
% 3. For the shear stress, max(uStar_hi, uStar_md) should be taken as the
%    meaningful uStar in the wind farm wake.  

kappa   = 0.4;

% init nueff
nut     = 0;
beta    = 0;

% init roughnesses
z0_hi   = z0;

% init heights
rbyz      = D/(2*zh);
zTip      = zh+D/2;
deltaFarm = zTip;
deltaWake = zTip;

% init u stars
uStar_hi      = uStar;
uStar_lo      = uStar; 

% compute decay length scale if applicable 
if(wakeDist>0)
    uinf  = uStar/kappa*(log(zh/z0) + 4.7*zh/L);
    phiM  = 1.0 + 4.7*zh/L;
    gamma = (phiM *D^2 * uinf)/(kappa*uStar*zh);
end

% inside the wind farm
if(farmDist > 0)

    % compute effective viscosity 
    nut        = 28*sqrt(0.5*cft);

    % exponential decay in the wake 
    if(wakeDist > 0)
        nut   = nut*exp(-(wakeDist)/gamma);
    end

    beta       = nut/(1+nut);

    % compute stability-dependent corrections to the wind profile 
    S1         = 4.7*(zh+nut*(zh+D/2))/(L*(1+nut));
    S2         = 4.7*(zh+nut*(zh-D/2))/(L*(1+nut));

    % psi        = log(zh/z0*(1-rbyz)^beta);
    psi        = log(zh/z0*(1-rbyz)^beta) + S2;
    
    % z0_hi      = zh*(1+rbyz)^beta * exp(-( cft/2/kappa^2 + psi^(-2) )^(-0.5));
    z0_hi      = zh*(1+rbyz)^beta * exp(-( cft/2/kappa^2 + psi^(-2) )^(-0.5) - S1);

    % compute IBL height 
    deltaFarm  = min(zTip+0.32*z0_hi*(farmDist/z0_hi)^0.8, H);
    deltaWake  = zTip;
    
    % compute u stars
    % uStar_hi   = uStar    * log(deltaFarm/z0) / log(deltaFarm/z0_hi); 
    % uStar_lo   = uStar_hi * ( log(zh/z0_hi*(1+rbyz)^beta) ) / ( log(zh/z0*(1-rbyz)^beta) ); % ensures continuous profile at hub height
    uStar_hi   = uStar    * (log(deltaFarm/z0) + 4.7*deltaFarm/L) / (log(deltaFarm/z0_hi) + 4.7*deltaFarm/L); 
    uStar_lo   = uStar_hi * ( log(zh/z0_hi*(1+rbyz)^beta) + S1) / ( log(zh/z0*(1-rbyz)^beta) + S2); % ensures continuous profile at hub height 
end

% velocity inside the wind farm 
% u = uStar_lo/kappa*log( (zh/z0) * (1 - rbyz)^beta);
u = uStar_lo/kappa*( log((zh/z0) * (1 - rbyz)^beta) + 4.7*(zh+nut*(zh-D/2))/(L*(1+nut)) );

% decay the wind profile (u is velocity at wind farm exit)
if(wakeDist > 0)
    uinf  = uStar/kappa*(log(zh/z0) + 4.7*zh/L);
    phiM  = 1.0 + 4.7*zh/L;
    gamma = (phiM *D^2 * uinf)/(kappa*uStar*zh);
    coeff = 0.5;
    zeta  = max((coeff-1)/(gamma) * (wakeDist) + 1, coeff);

    R0    = u/uinf;
    R     = 1 + (R0-1)*exp(-(wakeDist/gamma)^zeta);
    u     = uinf*R;
end

end