DataPath = '.\Data\Neuro\';
OutPath = '.\Data\Neuro\';

%% Generate and collect KSDs
load([DataPath,'AMDSpks_Sal.mat'],'SpkT_all','cids_all','cpos_all','Mouse_all','SortedStim');

% == KSD parameters ===
minT = -0.2; 
maxT = 2.2;
tStep = 0.5e-3;
bws = [1e-3,2e-3,5e-3,10e-3,25e-3,100e-3];
% bws = [10e-3];

tic
    fprintf('Generating average KSDs\r');
Res = struct;

uAM = unique(SortedStim(:,[3,1,2,20]),'rows'); % col 1: Mf; col 2: Md; col 3: Intensity
nUStm = size(uAM,1);
nClu = size(SpkT_all,2);

for b = 1:length(bws)
    bandwidth = bws(b);
    KSD_all = [];
    fprintf('== bandwidth %d ms ==\r',bandwidth*1000);
    
    fprintf('(%d clus ): ',nClu);
    for cc = 1:nClu
        KSD_clu = [];
        fprintf('%d ',cc);
        for st = 1:nUStm
            sel = (SortedStim.Mf == uAM.Mf(st) & SortedStim.Md == uAM.Md(st) & SortedStim.Intensity == uAM.Intensity(st)...
                    & strcmp(SortedStim.condition, uAM.condition(st)));
            [KSD,tt] = avgKSD(SpkT_all(sel,cc),minT,maxT,bandwidth,tStep);
            KSD_clu = [KSD_clu; KSD]; %#ok<AGROW>
        end
        KSD_all = cat(1,KSD_all,reshape(KSD_clu,[1,size(KSD_clu)]));
    end
    fprintf('\r');
KSD_all = permute(KSD_all,[2,3,1]);
suffix = replace(num2str(bandwidth*1000),'.','x');
varName = ['KSD_all_' suffix];
Res.(varName) = KSD_all;
varName = ['tt_' suffix];
Res.(varName) = tt;
end
toc
Res.uAM = uAM;
Res.bws = bws;Res.tStep = tStep;Res.minT = minT;Res.maxT = maxT;
Res.cids_all = cids_all; Res.cpos_all = cpos_all; Res.Mouse_all = Mouse_all;

% == save data ==
save([OutPath,'AMKSDs_Sal.mat'],'-struct','Res');disp('Data saved.');
beep

%% single trial KSD
% Output should be a nTrials x nBin x nClu matrix
load([DataPath, 'AMDSpks_Sal.mat'],'SpkT_all','cids_all','cpos_all','Session_all','SortedStim');

% == KSD parameters ===
minT = -0.2; 
maxT = 2.2;
tStep = 0.5e-3;
bws = [1e-3,2e-3,5e-3,10e-3,25e-3,100e-3];

tic
    fprintf('Generating single trial KSDs\r');
    nTrials = size(SortedStim,1);
    nClu = size(SpkT_all,2);
Res = struct;

for b = 1:length(bws)
    bandwidth = bws(b);
    KSD_all = [];
    fprintf('== bandwidth %d ms ==\r',bandwidth*1000);
    fprintf('(%d clus ): ',nClu);
    for cc = 1:nClu
        KSD_clu = [];
        fprintf('%d ',cc);
        for st = 1:nTrials
            [KSD,tt] = avgKSD(SpkT_all(st,cc),minT,maxT,bandwidth,tStep);
            KSD_clu = [KSD_clu; KSD];
        end
        KSD_all = cat(1,KSD_all,reshape(KSD_clu,[1,size(KSD_clu)]));
    end
    fprintf('\r');
    
    KSD_all = permute(KSD_all,[2,3,1]);
    suffix = replace(num2str(bandwidth*1000),'.','x');
    varName = ['KSD_trls_' suffix];
    Res.(varName) = KSD_all;
    varName = ['tt_' suffix];
    Res.(varName) = tt;
end
toc

Res.SortedStim = SortedStim;
Res.bws = bws;Res.tStep = tStep;Res.minT = minT;Res.maxT = maxT;
Res.cids_all = cids_all; Res.cpos_all = cpos_all; Res.Mouse_all = Mouse_all;
% == save data ==
save([OutPath,'AMKSDs_individual_Sal.mat'],'-struct','Res','-v7.3');disp('Data saved.');

beep


%% local functions
function [KSD,xx] = avgKSD(tempSpk,minT,maxT,bw,tStep)
 
if nargin < 4;     tStep = 0.5e-3; end %s
    NStim = length(tempSpk);
    tempSpk = vertcat(tempSpk{:});
    NSpk = length(tempSpk);
    xx = minT:tStep:maxT;
    if NSpk == 0
        KSD = zeros(size(xx)); return;
    end
    [KSD,~] = ksdensity(tempSpk,xx,'Bandwidth',bw);
    KSD = KSD * NSpk / NStim;
end
