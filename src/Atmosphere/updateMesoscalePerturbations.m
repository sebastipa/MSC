function meso = updateMesoscalePerturbations(meso, sol, abl)

fprintf('3LM Solution:\n');
tstart = tic;
fprintf(' > assembled matrix in ');

% initialize 5x5 bloxk-diagonal matrix and rhs
Amat = getMatWithStress5by5(meso, abl);
b    = rhsGetTurbines5by5  (meso, abl);

fprintf('%.2f s\n',toc(tstart));

tstart = tic;
fprintf(' > calculated solution in ');

% get some info on the system
% spparms('spumoni',2)

% compute initial guess 
s       = Amat \ b;

fprintf('%.2f s\n',toc(tstart));

% update meso scale fields based on solution vector s
meso = refactorSolution5by5(s, meso, abl);

end