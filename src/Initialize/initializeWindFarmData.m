function farm = initializeWindFarmData(sol, meso)

fprintf('\nInitializing wind farm data...\n');

fid = fopen('./input/settings.msm', 'r'); 

if (fid<=0)
    application = mfilename;
    error([application ':io'],'Unable to open settings.msm file');
end 

while ~feof(fid)
    tline        = fgetl(fid);
    keyword      = sscanf(tline,'%s',[1 1]);
    tline        = tline(length(keyword)+1:end);
    
    intflgValue  = sscanf(tline,'%d',[1 Inf]);
    scalarValue  = sscanf(tline,'%f',[1 Inf]);
    vectorValue  = sscanf(tline,'%f %f',[1 Inf]);
    word         = sscanf(tline,'%s',[1 1]);

    switch keyword
        case 'xsFarm',            farm.xs           = scalarValue * 1000;
        case 'ysFarm',            farm.ys           = scalarValue * 1000;
        case 'farmInputType',     farm.inputType    = intflgValue;
    end
end 

% read based on inputType flag
frewind(fid);
if(farm.inputType == 3)
    while ~feof(fid)
        tline        = fgetl(fid);
        keyword      = sscanf(tline,'%s',[1 1]);
        tline        = tline(length(keyword)+1:end);
        intflgValue  = sscanf(tline,'%d',[1 Inf]);
        word         = sscanf(tline,'%s',[1 1]);
        scalarValue  = sscanf(tline,'%f',[1 Inf]);
        switch keyword
            case 'excelName',        xlsName          = word;
            case 'useTurbineCurves', useTurbineCurves = intflgValue;    
            case 'Cp',               Cp               = scalarValue;
        end
    end

    % now see if must read turbine curves
    farm.useTurbineCurves = useTurbineCurves;
    if(useTurbineCurves == 1)
        load('input/turbineCurves.mat');
        farm.curves           = turbineCurves;
    end

    % create wind farm layout 
    [~,name] = xlsfinfo(strcat('./input/',xlsName));
    
    t              = 0;
    clusterStartId = 0;

    fprintf('\nCreating generalized cluster\n');
    
    % loop over clusters 
    for n = 2:length(name)
        
        data     = xlsread(strcat('./input/',xlsName), name{n});
        
        % loop over wind turbine clusters
        for ti = 1:size(data,1)
            t     = t + 1;
            
            % get coordinates from the excel
            xt(t) = data(ti,2);
            yt(t) = data(ti,3);
            
            % get inidividual Ct, D, zHub from the excel
            Cp_tmp(t)   = Cp;
            Ct_tmp(t)   = data(ti,8);
            D_tmp(t)    = data(ti,6);
            zHub_tmp(t) = data(ti,7);
        end 
        fprintf(' > %s (%d turbines)\n', name{n}, t-clusterStartId);

        clusterStartId = t;
    end 

    xt = xt - min(xt) + farm.xs;
    yt = yt - min(yt) + farm.ys;

    % rotate coordinates w.r.t min-coords turbine (temporary)
    rotation = 60;
    for i=1:length(xt)
        vx = xt(i) - farm.xs;
        vy = yt(i) - farm.ys;
    
        vxrot =  vx*cosd(rotation) + vy*sind(rotation);
        vyrot = -vx*sind(rotation) + vy*cosd(rotation);
    
        xt(i) = vxrot + farm.xs;
        yt(i) = vyrot + farm.ys;
    end

elseif(farm.inputType == 1 || farm.inputType == 2)
    while ~feof(fid)
        tline        = fgetl(fid);
        keyword      = sscanf(tline,'%s',[1 1]);
        tline        = tline(length(keyword)+1:end);
        intflgValue  = sscanf(tline,'%d',[1 Inf]);
        scalarValue  = sscanf(tline,'%f',[1 Inf]);
        switch keyword
            case 'D',                D                = scalarValue;
            case 'zHub',             zHub             = scalarValue;
            case 'Ct',               Ct               = scalarValue;
            case 'Cp',               Cp               = scalarValue;
            case 'useTurbineCurves', useTurbineCurves = intflgValue;
            case 'Sx',               Sx               = scalarValue;
            case 'Sy',               Sy               = scalarValue;
            case 'Ntx',              Ntx              = intflgValue;
            case 'Nty',              Nty              = intflgValue;
        end
    end

    % now see if must read turbine curves
    farm.useTurbineCurves = useTurbineCurves;
    if(useTurbineCurves == 1)
        load('input/turbineCurves.mat');
        farm.curves           = turbineCurves;
    end
    
    % dimensionalize turbine spacing 
    Sx      = Sx * D;
    Sy      = Sy * D;

    % aligned layout
    if(farm.inputType == 1)
        fprintf('Creating aligned layout:\n');
    
        t  = 0;
        
        dtx = Sx;
        dty = Sy;
    
        for i=1:Ntx
            for j=1:Nty
                x = (i-1)*dtx + farm.xs;
                y = (j-1)*dty + farm.ys;
    
                t = t + 1;
                xt(t) = x;
                yt(t) = y;
               
            end
        end
  
        Ct_tmp   = ones(length(xt),1) * Ct;
        Cp_tmp   = ones(length(xt),1) * Cp;
        D_tmp    = ones(length(xt),1) * D;
        zHub_tmp = ones(length(xt),1) * zHub;

    % staggered layout 
    else
        fprintf('\nCreating staggered layout:\n');
        
        t  = 0;
        
        dtx = Sx;
        dty = Sy;
    
        for i=1:Ntx
            for j=1:Nty
                x = (i-1)*dtx + farm.xs;
                y = (j-1)*dty + farm.ys + mod(i,2)*dty/2;
                
                t = t + 1;
                xt(t) = x;
                yt(t) = y;
               
            end
        end
        
        Ct_tmp   = ones(length(xt),1) * Ct;
        Cp_tmp   = ones(length(xt),1) * Cp;
        D_tmp    = ones(length(xt),1) * D;
        zHub_tmp = ones(length(xt),1) * zHub;
    end
    
else
    error('unknown farmInputType');
end

fclose(fid);

% save non-periodic turbine positions
farm.x  = xt;
farm.y  = yt;
farm.Nt = length(xt);

% correct wind farm mesh if periodic
if(sol.periodicSpanwise)
    
    % correct coordinates 
    xt = [xt, xt, xt];
    yt = [yt, yt - meso.Ly, yt + meso.Ly];
    
    % correct Ct, D, zHub
    Ct_tmp   = [Ct_tmp  , Ct_tmp  , Ct_tmp  ];
    Cp_tmp   = [Cp_tmp  , Cp_tmp  , Cp_tmp  ];
    D_tmp    = [D_tmp   , D_tmp   , D_tmp   ];
    zHub_tmp = [zHub_tmp, zHub_tmp, zHub_tmp];
end

% correct wind farm if mirroring is active 
if(sol.mirrorWindFarm)
    % correct coordinates 
    xt = [xt, xt];
    yt = [yt, yt];
    
    % correct Ct, D, zHub
    Ct_tmp   = [Ct_tmp  ,  Ct_tmp];
    Cp_tmp   = [Cp_tmp  ,  Cp_tmp];
    D_tmp    = [D_tmp   ,  D_tmp ];
    zHub_tmp = [zHub_tmp, -zHub_tmp];
end

% save periodic turbine positions
farm.xp  = xt;
farm.yp  = yt;
farm.Ntp = length(xt);

% save actual turbines number used for body force (actual + periodic, not mirrored)
farm.Nti = farm.Ntp;
if(sol.mirrorWindFarm) 
    farm.Nti = farm.Nti / 2; 
end

% projection width 
farm.Lbf = (max(meso.dx, meso.dy));

% save wake model run counter 
farm.wakeModelRuns = 0;

% save the temporary variable containing Ct, D, zHub for each turbine
tData_tmp.Ct       = Ct_tmp;
tData_tmp.Cp       = Cp_tmp;
tData_tmp.D        = D_tmp;
tData_tmp.zHub     = zHub_tmp;

% init wind turbines (including the periodized turbines)
farm.turbines = initializeWindTurbines(sol, farm, meso, tData_tmp);

% set this flag always to zero (it is an handle to modify wake model coeffs in a tuning loop)
farm.tuningWakeModel = 0;

% print turbine array info 
fprintf(' - periodic : %*d \n'    , -15, sol.periodicSpanwise);
fprintf(' - Nt       : %*d \n'    , -15, farm.Nt);
fprintf(' - xmin     : %*.5f \n'  , -15, min(farm.x));
fprintf(' - xmax     : %*.5f \n'  , -15, max(farm.x));
fprintf(' - ymin     : %*.5f \n'  , -15, min(farm.y));
fprintf(' - ymax     : %*.5f \n\n', -15, max(farm.y));

% override turbine interpolation weights for original one-point coupling
if(sol.singlePointCoupling)

    % compute min x position and average y of the wind farm 
    wind    = [1; 0; 0];
    x_hat   = wind/norm(wind,2);
    z_hat   = [0; 0; 1];
    y_hat   = cross(z_hat, x_hat);

    % baricentric wind farm point
    p_mean  = [mean(farm.x); mean(farm.y); 0];

    % quantities of interest 
    xmin    = 1e20;
    yavg    = 0;
    Davg    = 0;

    for t=1:farm.Nt
        dist   = farm.turbines(t).base - p_mean;
        dist_x = dot(dist, x_hat);
        dist_y = dot(dist, y_hat);

        yavg = yavg + dist_y;

        if(dist_x < xmin)
            xmin = dist_x;
        end

        Davg = Davg + farm.turbines(t).D;
    end

    yavg = yavg / farm.Nt + p_mean(2);
    xmin = xmin + p_mean(1);
    Davg = Davg / farm.Nt;

    % point located at the start of the wind farm 
    p_start  = [xmin; yavg; 0.0];

    % point located 10D upstream for velocity sampling
    farm.p_sample = p_start - 10*Davg*x_hat;

    % interpolation from mesoscale domain 
    [ids, weights, ~, ~] = findBilinaerInterpParam(farm.p_sample(1), farm.p_sample(2), meso.xs, meso.ys, meso.dx, meso.dy, meso.nx, meso.ny);
    
    % all turbine have the same interpolation point for freestream velocity
    for t=1:farm.Ntp
        farm.turbines(t).idsMeso                = ids;
        farm.turbines(t).weightsMeso            = weights;
    end

end


end