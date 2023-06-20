% This script recreates subpanels for van den Berg et al 202X
% "Salicylate on AM"

%% Load analyzed data
% This reads summarized latencies and hit rates of the mice from "results_bX.mat"

clear; close all
addpath('.\functions')
dataPath = '.\Data\Behavior';
DataPathSIM = '.\Data\PCA_SIM\';
load([dataPath '\salicylateResults.mat']); % summary data

if ~exist('SimRes','var') || ~ismember('CDFDiff',fieldnames(SimRes))
    load([DataPathSIM,'SIMResTableEval_Sal.mat'],'SimRes');
end


[~,~]=mkdir('Figures');

FontSize = 14;
saveFigs = 1;

%% false alarm rate
% HR_lat 2 x 6 x 2 x 6 x 3
%        ModFreq x ModDepth x Int x Mouse x condition
%     false alarm rates are distributed across all ModFreq at first
%     ModDepth (= 0). 

% false alarm at 45dB 
MF_idx = 1;
MD_idx = 1;
Int_idx = 1;
FA_45dB = squeeze(HR_lat(MF_idx,MD_idx,Int_idx,:,:));

% false alarm at 60dB
Int_idx = 2;
FA_60dB = squeeze(HR_lat(MF_idx,MD_idx,Int_idx,:,:));

% FA_45dB, FA_60dB: 6 mouse x 3 conditions
% FA_all: 12 mouse*intensity x 3 conditions
FA_all = [FA_45dB;FA_60dB]; 

[p,tbl,stats] = anova2([FA_45dB;FA_60dB],6,off');
% c = multcompare(stats);

%% plot threshold for Poster FENS 2022
close all;
figure('Position',[400,100,500,250]);
uNZMF = uMF(uMF~=0);
AvgThr = mean(Thr_lat,3);
SDThr = std(Thr_lat,0,3);
SEMThr = SDThr./sqrt(length(mice));
axlist = [];
for ii  = 2%1:length(uInt)
    for ff = 1:length(uNZMF)
        index = 1   +   2 *(2-ii) + (ff - 1);
        ax = subplot(1,2,index);
        axlist(index) = ax;
        plot(1:3,squeeze(Thr_lat(ff,ii,:,:)), '-o','Color',[.5,.5,.5]); hold on;
        ebBeh = errorbar(squeeze(AvgThr(ff,ii,:)),squeeze(SEMThr(ff,ii,:)),'-ks','LineWidth',2);
%         title([num2str(uInt(ii)), ' dB  - ', num2str(uNZMF(ff)), ' Hz']);
        title([num2str(uNZMF(ff)), ' Hz']);
        xticks(1:3);xticklabels(Periods);
        set(ax,'XTickLabelRotation',45);
        ylabel('AM detection threshold (dB)');
        
        set(ax,'YDir','reverse')
        ylim(ax,[-20,0]);
        
        yyaxis(ax,'right');
        Md = [0.03,0.06,0.125,0.25,0.5,1];
        dBMd = 20*log10(Md);
        set(ax,'ytick',dBMd,'yticklabel',Md,...
            'ycolor','k');
        ylabel(ax, 'AM detection threshold (m)');
        
        set(ax,'YDir','reverse')
        ylim(ax,[-20,0]);
    end
end
linkaxes(axlist);
xlim(axlist(1),[0.5,3.5]);
xticks(axlist(1),1:3);xticklabels(axlist(1),Periods);
% sgtitle('Effect of Salicylate on AM Detection Threshold');

% add simulation

detLvls = [140]; % best parameters
fIdx = [1,4];
sIdx = ismember(SimRes.thr, detLvls) & abs(SimRes.NDTime - 0.15) < 0.005 & strcmpi(SimRes.condition,'baseline');
sIdx = find(sIdx)';
simColor = [0,.7,.7];

for scnt = 1:length(sIdx)%ss-10:10:ss+50
    sss = sIdx(scnt);
    dP_SIM = SimRes.dP_SIM{sss};
    
    Thr_SIM = SimRes.(['Thr_SIM_',num2str(uInt(ii))]){sss};
    Thr_SIM = Thr_SIM(fIdx,:);
    LineWidth = 3;%scnt;
    for ff = 1:2 % frequency
        ax = subplot(1,2,ff);
        lnSIM = plot(Thr_SIM(ff,:),':','Color',simColor,'MarkerFaceColor',simColor,'LineWidth',LineWidth) ;hold on;%
        legend([ebBeh,lnSIM],{'Behavior','Simulation'},'location','best');
    end
end
%% plot threshold
close all;
figure('Position',[400,100,1000,800]);
uNZMF = uMF(uMF~=0);
AvgThr = mean(Thr_lat,3);
SDThr = std(Thr_lat,0,3);
SEMThr = SDThr./sqrt(length(mice));
axlist = [];
for ii  = 1:length(uInt)
    for ff = 1:length(uNZMF)
        index = 1   +   2 *(2-ii) + (ff - 1);
        ax = subplot(2,2,index);
        axlist(index) = ax;
        plot(1:3,squeeze(Thr_lat(ff,ii,:,:)), '-o'); hold on;
        errorbar(squeeze(AvgThr(ff,ii,:)),squeeze(SEMThr(ff,ii,:)),'-ks','LineWidth',2);
        title([num2str(uInt(ii)), ' dB  - ', num2str(uNZMF(ff)), ' Hz']);
        xticks(1:3);xticklabels(Periods);
        ylabel('Detection Threshold (dB; 20*log(m))');
        
        ylim(ax,[-20,0]);
        
        yyaxis(ax,'right');
        Md = [0.03,0.06,0.125,0.25,0.5,1];
        dBMd = 20*log10(Md);
        set(ax,'ytick',dBMd,'yticklabel',Md,...
            'ycolor','k');
        ylabel(ax, 'Detection Threshold (m)');
         
        ylim(ax,[-20,0]);
    end
end
linkaxes(axlist);
xlim(axlist(1),[0.5,3.5]);
xticks(axlist(1),1:3);xticklabels(axlist(1),Periods);
sgtitle('Effect of Salicylate on AM Detection Threshold');
