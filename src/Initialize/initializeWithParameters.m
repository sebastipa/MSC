function msm = initializeWithParameters(deltat_i, lapse_i)

% msm is the top level struct of the MSM and contains 
% msm.sol  : solution parameters 
% msm.meso : meso scale domain data 
% msm.micro: micro scale domain data 
% msm.farm : wind farm data 
% msm.abl  : atmospheric boundary layer data

% Init Solution Settings 
sol = initializeSolutionData();

% Init ABL Settings
abl = initializeAblDataWithParameters(sol, deltat_i, lapse_i);

% Create solution state and output directory 
sol.state = strcat(...
    'output_H',num2str(abl.H),...
    '_DeltaT', num2str(abl.dTH),...
    '_Gamma',num2str(abl.dTdz));

% Init Mesoscale Domain
meso = initializeMesoData(sol, abl);

% Init Wind Farm Data
farm = initializeWindFarmData(sol, meso);

% Init Microscale Data 
micro = initializeMicroData(sol, meso, farm);

% Compute wind farm buffer indices 
meso  = setWindFarmBufferOnDomain(meso, farm);

% do some computations for initialization 
if(sol.readTurbineThrust)
    fprintf('\nReading wind farm thrust from input...\n');
    % update wind dir, qpoints, pwr and thr based LES results. Mesoscale
    % inflow calculated based on the mesoscale reconstructed bk hub wind
    startTime = 1e5;
    farm      = readTurbineThrustTOSCA(farm, meso, startTime);
else 

    % update CtPrime if necessary 
    if(sol.readTurbineCtPrime)
        fprintf('\nReading wind farm Ct'' from input...\n');
        startTime = 1e5;
        farm = readTurbineCtPrimeTOSCA(farm, startTime);
    end

    % update wind dir, qpoints and mesoscale inflow using the mesoscale 
    % reconstructed bk hub wind
    farm      = updateTurbinesWindAndDirection(farm, meso, sol, abl);

    % update thr and pwr using wake model
    farm      = wakeModel(farm, sol, abl);
end

% compute average turbine diameter for shear stress model
abl.avgD = averageTurbineDiameter(farm);

% initialize body force
meso = updateBodyForce(meso, farm, abl);

% compute local turbine spacing field 
meso = updateTurbineSpacingField(meso, farm);

% assign output 
msm.sol   = sol;
msm.abl   = abl;
msm.meso  = meso;
msm.micro = micro;
msm.farm  = farm;

end