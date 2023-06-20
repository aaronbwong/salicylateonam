DataPath = '.\Data\Neuro\';
OutPath = '.\Data\Neuro\';

%% PSTH generation
% Output should be a nTrials x nBin x nClu matrix
% ==== PSTH parameters (TRANSITION) =====
minT = 0.5; maxT = 1.5;
bins = 1./(2.^(1:10));

load([DataPath,'AMDSpks_Sal.mat'],'SpkT_all','cids_all','cpos_all','Mouse_all','SortedStim');

tic
    fprintf('Generating single trial PSTHs (TRANSITION)\r');

% ===== initialize variables =====
    Res = struct;

% ===== loop through bin sizes ====
for b = 1:length(bins)
    binWidth = bins(b);%1e-3;
    binEdges = minT:binWidth:maxT;
    nBins = length(binEdges)-1;
    fprintf('== bin %d ms ==\r',binWidth*1000);

    PSTH_all = [];
% for m = 1:nMice
    nTrials = size(SortedStim,1);
    nClu = size(SpkT_all,2);
    fprintf('(%d clus ): ',nClu);
    for cc = 1:nClu
        PSTH_clu = NaN(nTrials,nBins);
        fprintf('%d ',cc);
        for st = 1:nTrials
            [PSTH_clu(st,:),~] = histcounts(SpkT_all{st,cc},binEdges);
        end
%         fprintf('\r');
        PSTH_all = cat(3,PSTH_all,PSTH_clu);
    end
    fprintf('\r');

suffix = replace(num2str(binWidth*1000),'.','x');
varName = ['PSTH_all_' suffix];
Res.(varName) = PSTH_all;
varName = ['binEdges_' suffix];
Res.(varName) = binEdges;
end
toc

Res.SortedStim = SortedStim;
Res.cids_all = cids_all; Res.cpos_all = cpos_all; Res.Mouse_all = Mouse_all;

Res.minT = minT; Res.maxT = maxT; Res.bins = bins;
save([OutPath,'AMPSTHs_Sal.mat'],'-struct','Res','-v7.3');
beep


% %%  ==== PSTH parameters FULL LENGTH ===== 
minT = -0.5; maxT = 2.5;

load([OutPath,'AMDSpks_Sal.mat'],'SpkT_all','cids_all','cpos_all','Mouse_all','SortedStim');

tic
fprintf('Generating single trial PSTHs (FULL LENGTH)\r');

% ===== initialize variables =====
    Res = struct;

% ===== loop through bin sizes ====
for b = 1:length(bins)
    binWidth = bins(b);%1e-3;
    binEdges = minT:binWidth:maxT;
    nBins = length(binEdges)-1;
    fprintf('== bin %d ms ==\r',binWidth*1000);

    PSTH_all = [];
% for m = 1:nMice
    nTrials = size(SortedStim,1);
    nClu = size(SpkT_all,2);
    fprintf('(%d clus ): ',nClu);
    for cc = 1:nClu
        PSTH_clu = NaN(nTrials,nBins);
        fprintf('%d ',cc);
        for st = 1:nTrials
            [PSTH_clu(st,:),~] = histcounts(SpkT_all{st,cc},binEdges);
        end
%         fprintf('\r');
        PSTH_all = cat(3,PSTH_all,PSTH_clu);
    end
    fprintf('\r');

suffix = replace(num2str(binWidth*1000),'.','x');
varName = ['PSTH_all_' suffix];
Res.(varName) = PSTH_all;
varName = ['binEdges_' suffix];
Res.(varName) = binEdges;
end
toc

Res.SortedStim = SortedStim;
Res.cids_all = cids_all; Res.cpos_all = cpos_all; Res.Mouse_all = Mouse_all;

Res.minT = minT; Res.maxT = maxT; Res.bins = bins;
save([OutPath,'AMPSTHs_full_Sal.mat'],'-struct','Res','-v7.3');
beep