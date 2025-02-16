function sol = initializeSolutionData()

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
        case 'iterations',          sol.iter                = intflgValue;
        case 'deepArray',           sol.deepArray           = intflgValue;
        case 'localInduction',      sol.localInduction      = intflgValue;
        case 'readTurbineThrust',   sol.readTurbineThrust   = intflgValue;
        case 'readTurbineCtPrime',  sol.readTurbineCtPrime  = intflgValue;
        case 'periodicSpanwise',    sol.periodicSpanwise    = intflgValue;
        case 'mirrorWindFarm',      sol.mirrorWindFarm      = intflgValue;
        case 'ablModel',            sol.ablModel            = intflgValue;
        case 'displayMicroScale',   sol.displayMicroScale   = intflgValue;
        case 'wakeModelLogLaw',     sol.wakeModelLogLaw     = intflgValue;
        case 'singlePointCoupling', sol.singlePointCoupling = intflgValue;
        case 'excludeBackground',   sol.excludeBackground   = intflgValue;
    end
end

% set iteration counter
sol.i       = 1;
sol.solNorm = zeros(sol.iter,1);

% read background input mode if applicable 
if(sol.ablModel == 1)
elseif(sol.ablModel == 2)
    
    frewind(fid);
    while ~feof(fid)
        tline        = fgetl(fid);
        keyword      = sscanf(tline,'%s',[1 1]);
        tline        = tline(length(keyword)+1:end);
        word         = sscanf(tline,'%s',[1 1]);
        switch keyword
            case 'bkStateName', sol.bkStateName = word;
        end
    end
else
    application = mfilename;
    error([application ':io'],'unknown ablModel in settings.msm');
end

fclose(fid);

% do some checks
% ======================================================================= %

% if turbine thrust is provided there is no need to iterate as thrust is fixed 
if(sol.readTurbineThrust && sol.iter > 1)
    fprintf('--> Warning: reducing iterations to 1 as turbine thrust is fixed\n');
    sol.iter = 1;
end

% if turbine thrust is provided there is no need to read also CtPrime
if(sol.readTurbineThrust && sol.readTurbineCtPrime)
    fprintf('--> Warning: setting readTurbineCtPrime to 0 as turbine thrust is fixed\n');
    sol.readTurbineCtPrime = 0;
end

% if original coupling flag is set to 1 set to zero the following flags: 
% localInduction, deepArray, as they are all part of the
% new model formulation 
if(sol.singlePointCoupling && sol.localInduction == 1)
    fprintf('--> Warning: setting localInduction to 0 as original coupling is active\n');
    sol.localInduction = 0;
end
if(sol.singlePointCoupling && sol.deepArray == 1)
    fprintf('--> Warning: setting deepArray to 0 as original coupling is active\n');
    sol.deepArray = 0;
end

% if exclude background flag is set to 1 set to zero all the following flags
% deepArray, singlePointCoupling as the model only solves
% the micro scale. Moreover, iterations must be zero as there is no
% coupling. 
if(sol.excludeBackground); sol.iter = 0; end
if(sol.excludeBackground && sol.deepArray == 1)
    fprintf('--> Warning: setting deepArray to 0 as MSM is running in wake model mode\n');
    sol.deepArray = 0;
end
if(sol.excludeBackground && sol.singlePointCoupling == 1)
    fprintf('--> Warning: setting singlePointCoupling to 0 as MSM is running in wake model mode\n');
    sol.singlePointCoupling = 0;
end

% disable warnings
warning('off','all');

% ======================================================================= %
