function meso = updateMesoscaleBackground(meso, sol, abl)

% initialize 4x4 bloxk-diagonal matrix and rhs
Amat = getRedMatWithStress4by4(meso, abl);
b    = rhsGetPressure4by4     (meso, abl);

fprintf('3LM-Reconstruction: ');
tstart = tic;

% compute initial guess 
s       = Amat \ b;

% save initial guess (only contains pressure effect)
meso    = refactorSolution4by4(s, meso);

% update background hub velocity with pressure effect
[ubk, vbk] = mesoscaleAvg2Hub(meso.ur1, meso.vr1, abl);
meso.ubk   = ubk; 
meso.vbk   = vbk; 

fprintf('Elapsed Time %.2f s\n',toc(tstart));

% if(sol.deepArray && sol.localShearStress)
%     tolerance   = 1e-4;
%     relRes      = 1e20;
%     iter        = 1;
%     maxiter     = 25;
%     solnorm     = 0;
%     
%     % update convolutional shear stress perturbation 
%     while (relRes>tolerance && iter < maxiter)
%         tstart = tic;
%         fprintf(' > iteration %d: ', iter);
%         
%         
%         bb = b + rhsGetStressPerturb4by4(s, meso, abl);
%         s = Amat \ bb;
% 
%     
%         relRes = abs(norm(s,2)-solnorm)/norm(s,2);
%         fprintf('relRes = %e, Elapsed Time = %.2f s\n', relRes, toc(tstart));
%        
%         solnorm = norm(s,2);
%         iter    = iter + 1;
%     end
% 
%     % update meso scale fields based on solution vector s (contains pressure and deep array effects)
%     meso = refactorSolutionDeepArray4by4(s, meso);
% 
%     % update background hub velocity with pressure and deep array effect
%     [ubk_da, vbk_da] = mesoscaleAvg2Hub(meso.ur1_da, meso.vr1_da, abl);
%     meso.ubk_da      = ubk_da; 
%     meso.vbk_da      = vbk_da; 
% 
% else
% 
% end