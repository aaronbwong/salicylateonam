%% This script generates panels of figures from source data
% Please make sure the folder "Data" from the source data is extracted 
% directly at the root directory of the github code, so that it contains:
% github_dir\Data\
% github_dir\Data\Neuro\
% github_dir\Data\Behavior\
%  etc.
% and from the code repository: 
% github_dir\run_me_AllFigures.m
% github_dir\src\paper
%  etc.

addpath('functions\')
addpath(genpath('src\'))

%% Main figures

% Figure 1
Fig1_Spontaneous_ephys
% to check

%% Figure 2
Fig2_FRA_minThr

%% Figure 3
Fig3_FigS2_AM_example_units_ABW

%% Figure 4
Fig4_Population_salicylate_result

%% Figure 5
Fig5_Iceberg_effect

%% Figure 6
Fig6_BehaviouralResults_Raw_Calc

%% Figure 7
Fig7_DRC_STRF_SalicylateCompare

%% Figure 8
Fig8_DRC_signalPower_CGFAverage % some non main figure code to be relocated

%% Figure 9
Fig9_DRC_STRF_PP_Linearity

%% Figure 10
Fig10_plotSimResults

%% Supplementary Figures
% Figure S1
FigS1_SalConcVSTime_SpontRateVsCF

% Figure S2
% see Fig3_FigS2_AM_example_units_ABW

% Figure S3
% XXX data analyzed in Igor

% Figure S4
FigS4_DRC_STRF_Illustration