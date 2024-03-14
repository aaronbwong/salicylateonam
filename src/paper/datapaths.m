%% This script specifies the path to the data
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

% look for root directory that contains this file `datapaths.m`
rootpath = fileparts(which('datapaths.m'));
foldername = strsplit(rootpath,filesep);
rootpath = fullfile(foldername{1:end-2});
clear foldername

% specify all other paths relative to the root
FRAPath = fullfile(rootpath,'Data\Neuro\FRA_analysis');
AMPath = fullfile(rootpath,'Data\Neuro\AMn_analysis');
PopPath = fullfile(rootpath,'Data\Neuro\Population');
SumPath = fullfile(rootpath,'Data\Neuro');
SpkPath = fullfile(rootpath,'Data\Neuro\Spks');
BehPath = fullfile(rootpath,'Data\Behavior');
FigPath = fullfile(rootpath,'Figures\Links');
FullFigPath = fullfile(rootpath,'Figures');
DRCPath = fullfile(rootpath,'Data\Neuro\DRC');
PCAPath = fullfile(rootpath,'Data\PCA_SIM');
SimPath = fullfile(rootpath,'Data\PCA_SIM');
ConcPath = fullfile(rootpath,'Data\salicylate_concentration');
ImgPath = fullfile(rootpath,'Data\Imaging');