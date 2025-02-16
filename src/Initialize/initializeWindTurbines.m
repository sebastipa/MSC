function turbines = initializeWindTurbines(sol, farm, meso, tData)

% turbines are initialized assuming the wind is aligned to x direction
windDir = [1; 0; 0];

% loop over wind turbine bases (non periodic)
for t=1:farm.Nt
    turbine.base    = [farm.x(t); farm.y(t); 0.0];
    turbine.rotor   = turbine.base + [0.0; 0.0; tData.zHub(t)];
    turbine.qpoints = getQuadraturePoints(windDir, turbine.rotor, tData.D(t));
    
    % save Ct, D, zHub for this turbine 
    turbine.Ct      = tData.Ct(t);
    turbine.Cp      = tData.Cp(t);
    turbine.zHub    = tData.zHub(t);
    turbine.D       = tData.D(t);
    turbine.Ar      = tData.D(t)^2*pi/4.0;
    turbine.CtPrime = 0; % disk thrust coefficient (just output)
    turbine.uDiskMT = 0; % disk velocity from momentum theory (just output)
    turbine.uDiskVC = 0; % disk velocity from vortex cylinder (just output)

    if(turbine.Ct>=1)
        error([application ':io'],'Ct higher ore equal to 1 detected in the wind farm');
    end

    % interpolation from mesoscale domain 
    [ids, weights, iClose, jClose] = findBilinaerInterpParam(turbine.rotor(1), turbine.rotor(2), meso.xs, meso.ys, meso.dx, meso.dy, meso.nx, meso.ny);
    turbine.idsMeso                = ids;
    turbine.weightsMeso            = weights;
    turbine.backgroundWindMagP     = 0.0;    % background wind with pressure effect
    turbine.backgroundWindMagPDA   = 0.0;    % background wind with pressure and deep array effect
    
    % meso/micro scale wind at turbine location (just initialize)
    turbine.windDirHat             = windDir/norm(windDir,2);

    % body force projection (conf. interv > 3-sigma = 99.87% of force recovered)
    iLbf                           = ceil(farm.Lbf/meso.dx);
    jLbf                           = ceil(farm.Lbf/meso.dy);
    nsigma                         = 3;
    project.imin                   = max(iClose-nsigma*iLbf, 1);
    project.imax                   = min(iClose+nsigma*iLbf, meso.nx);
    project.jmin                   = max(jClose-nsigma*jLbf, 1);
    project.jmax                   = min(jClose+nsigma*jLbf, meso.ny);
    turbine.project                = project;

    % rotor-average wind from wake model (just initialize)
    turbine.inflow                 = 0;

    % wake model - specific parameters (just initialize) 
    turbine.kStar                  = 0;
    turbine.TI                     = 0;

    % thrust and power 
    turbine.Tx                     = 0;
    turbine.Ty                     = 0;
    turbine.T                      = 0;
    turbine.P                      = 0;

    % thrust and power by iteration (add slot for initial condition)
    turbine.iterP                  = zeros(sol.iter + 1,1);
    turbine.iterT                  = zeros(sol.iter + 1,1);

    % add this turbine to the list 
    turbines(t) = turbine;
end

% periodize the array, so that interpolation weights are correct
if(sol.periodicSpanwise || sol.mirrorWindFarm)

    if(sol.periodicSpanwise)
        turbines = [turbines, turbines, turbines];
    end

    if(sol.mirrorWindFarm)
        turbines = [turbines, turbines];
    end

    % correct position-dependent quantities 
    for t=1:farm.Ntp
        turbines(t).base    = [farm.xp(t); farm.yp(t); 0.0];
        turbines(t).rotor   = turbines(t).base + [0.0; 0.0; tData.zHub(t)];
        turbines(t).qpoints = getQuadraturePoints(windDir, turbines(t).rotor, tData.D(t));
    end
end
