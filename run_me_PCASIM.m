%%
addpath(genpath('.\functions\'))

%% generate KSD and PSTHs from spike times
createKSDs; % trial-averaged KSD (for PCA) and individual-trial KSD for model
createPSTHs; % PSTHs for display and PCA

%% perform diff PCA
performDiffPCA; % perform diff PCA with different inputs 

%% Model

% simulate response latencies with different thresholds and ND times,
% store the results in a table "SimRes", and save to "SimResTable_Sal.mat"
performSimulation; 

% load behavioral data "salicylate_results.mat", compare results in SimRes,
% and update SimRes with goodness-of-fit data, save to ""SimResTableEval_Sal.mat"
evaluateSimulation;

%% plot figures
plotSimResults;