function msm = solve(msm)

% top level function for multiscale model (MSM) Solution
% ----------------------------------------------------------------------- %
% 1: compute body force on mesoscale domain
% 2: compute localized shear stress using shear stress model 
% 3: run large scale blockage model (3LM)
% 4: reconstruct background hub-wind (3LMR)
% 5: run wake model using reconstructed background wind (Bastankhah + vortex cylinder)
% 6: repeat from 1 for the required number of iterations (usually 2/3)
% ----------------------------------------------------------------------- %

niter   = msm.sol.iter;
solnorm = 0;

while(msm.sol.i<=niter)
  
    tstart = tic;
    fprintf('\nMSM iteration %d\n\n',msm.sol.i);

    % body force 
    msm.meso = updateBodyForce(msm.meso, msm.farm, msm.abl);
    
    % three layer model (Allaerts Meyers)
    msm.meso = updateMesoscalePerturbations(msm.meso,  msm.sol, msm.abl);
    
    % three layer model reconstruction (Stipa)
    msm.meso = updateMesoscaleBackgroundApprox(msm.meso,  msm.sol, msm.abl);

    % deep array model (Stipa)
    if(msm.sol.deepArray)
        msm.meso = updateTurbineCftField(msm.meso, msm.farm);
        msm.meso = updateDeepArrayEffect(msm.meso, msm.abl);
    end
    
    % update background wind speed at turbine position (only if thrust with wake model)
    if(~msm.sol.readTurbineThrust)
        msm.farm = updateTurbinesWindAndDirection(msm.farm, msm.meso, msm.sol, msm.abl);
    
        % wake model 
        msm.farm = wakeModel(msm.farm, msm.sol, msm.abl);
    end

    % compute 3LM residual
    relRes                     = abs(norm(msm.meso.p,2)-solnorm)/norm(msm.meso.p,2);
    solnorm                    = norm(msm.meso.p,2);
    msm.sol.solNorm(msm.sol.i) = solnorm;
    msm.sol.i                  = msm.sol.i + 1; 
    fprintf('\nMSM: residual on pressure perturbation = %e, Elapsed Time = %.2f s\n\n', relRes, toc(tstart));
end

% compute microscale wind velocity if applicable
if(msm.sol.displayMicroScale)
    msm.micro = updateMicroscaleWind(msm.micro, msm.meso, msm.farm, msm.sol, msm.abl);
end

end