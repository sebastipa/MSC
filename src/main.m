% Multi-Scale Model (MSM) 
% 
% ----------------------------------------------------------------------- %
%
% This is the UBC version of the multi-scale model, developed by Sebastiano
% Stipa. The model is a mix of different models (meso and micro scale models)
% and can model large-scale and small-scale blockage effects, as well as
% deep array effects and wake interactions (intra-farm and farm-farm).
%
% Input files should be contained in the ./input/ directory. 
% Outputs are thrown in the msm struct which is saved in the output directory
%
% Contacts: Sebastiano Stipa (sebastiano.stipa@vki.ac.be)
%
% ----------------------------------------------------------------------- % 

clear  all; close all; clc;
set    (0,'defaulttextinterpreter','latex');
addpath('.\Initialize');
addpath('.\Global');
addpath('.\System');
addpath('.\Atmosphere');
addpath('.\WindFarm');

%% Initialize models 
msm = initialize();

%% Solve models
msm = solve(msm);

%% Save Data
save(strcat(msm.sol.state,'/',msm.sol.state,'.mat'), 'msm','-v7.3');

%% Plot Data 
msmPlotData(msm);
