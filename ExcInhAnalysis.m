% AM_results_salicylate.mat
% AMFR
% 6 x 5 x 2 x 2 x 55
% ModDepth x ModFreq x Int x Set(condition) x clu

% FRA results
% location: W:\AW_repository_transfer\Maurits\FRA_analysis
% size(FRA.FRASR) = 16 x 22 x 2 x X
%                   int x freq x set x cluster

% NeuroDataPath = 'Data\Neuro\';
% load([NeuroDataPath,'AM_results_salicylate.mat'],'units');
% 
% mice = unique(units(:,1));

addpath('.\functions');

salUnits = readtable('.\Data\salicylate_units.xlsx');

FRAPath = 'W:\AW_repository_transfer\Maurits\FRA_analysis';
Mice = unique(salUnits.mouse);
salUnits.ZeroRate_Base = NaN(size(salUnits,1),1);
salUnits.CF_Base = NaN(size(salUnits,1),1);
salUnits.minThr_Base = NaN(size(salUnits,1),1);
salUnits.ZeroRate_Sal = NaN(size(salUnits,1),1);
salUnits.CF_Sal = NaN(size(salUnits,1),1);
salUnits.minThr_Sal = NaN(size(salUnits,1),1);


salUnits.FRAFR_Base = cell(size(salUnits,1),1);
salUnits.FRAFR_Sal = cell(size(salUnits,1),1);


setNum_Base = [1,2]; % 1 or 2
setNum_Sal = 9;
for m = 1:length(Mice)
% m = 2;
    Mouse = num2str(Mice(m),'M%02.0f');
    multiWaitbar('Mouse', m/length(Mice));
    load([FRAPath,'\', Mouse, '\', Mouse '_FRA_both_data.mat'],'FRA','cids','cpos');
    unitIdx = salUnits.mouse == Mice(m);
    setIdx_Base = ismember(FRA.FRASetNum,setNum_Base);
    setIdx_Sal = FRA.FRASetNum == setNum_Sal;
    %if there is no salicylate FRA
    if ~ismember(FRA.FRASetNum,9)
       continue
    end

    UFreq = FRA.UFreq;
    UInt = FRA.UInt;

    cIdx = ismember(cids,salUnits.cids(unitIdx));

    % baseline
    ZeroRate_Base = squeeze(FRA.FRASR(1,1,setIdx_Base,cIdx));
    salUnits.ZeroRate_Base(unitIdx) = ZeroRate_Base;
    salUnits.CF_Base(unitIdx) = FRA.FRACF(setIdx_Base,cIdx)';
    salUnits.minThr_Base(unitIdx) = FRA.FRAThr(setIdx_Base,cIdx)';
    salUnits.FRAFR_Base(unitIdx) = squeeze(mat2cell(FRA.FRASR(:,:,setIdx_Base,cIdx),[16],[22],[1],[ones(1,sum(cIdx))]));

    % salicaylte
    ZeroRate_Sal = squeeze(FRA.FRASR(1,1,setIdx_Sal,cIdx));
    salUnits.ZeroRate_Sal(unitIdx) = ZeroRate_Sal;
    salUnits.CF_Sal(unitIdx) = FRA.FRACF(setIdx_Sal,cIdx)';
    salUnits.minThr_Sal(unitIdx) = FRA.FRAThr(setIdx_Sal,cIdx)';
    salUnits.FRAFR_Sal(unitIdx) = squeeze(mat2cell(FRA.FRASR(:,:,setIdx_Sal,cIdx),[16],[22],[1],[ones(1,sum(cIdx))]));

end

multiWaitbar('Close all');
% clearvars('-except','UFreq','UInt','salUnits');

idx = ~isnan(salUnits.CF_Base) & ~isnan(salUnits.CF_Sal);
salUnits_selected = salUnits(idx,:);


%%
for cc = 1:size(salUnits_selected,1)
    figure;
    subplot(1,2,1)
    imagesc(salUnits_selected.FRAFR_Base{cc})

    subplot(1,2,2);
    imagesc(salUnits_selected.FRAFR_Sal{cc})

end
