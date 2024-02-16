%% Figure 10 Van den Berg*, Wong* et al, 2024, iScience

DataPathPCA = '.\Data\PCA_SIM\';
DataPathSIM = '.\Data\PCA_SIM\';
DataPathNeuro = '.\Data\Neuro\';
DataPathBeh = '.\Data\Behavior\';

%% load results 
% load results from PCA
    load([DataPathPCA,'PCAResTable_Sal.mat'],'PCARes')
    AMKSDs = matfile([DataPathNeuro,'AMKSDs_Sal.mat']);
    

%  load results from Simulation if not present
if ~exist('SimRes','var') || ~ismember('CDFDiff',fieldnames(SimRes))
    load([DataPathSIM,'SIMResTableEval_Sal.mat'],'SimRes');
end
% load results from behavior
results_lat = load([DataPathBeh,'salicylateResults.mat'],'mice', 'Thr_lat','uMF','uMD','uInt');
i = 0;  
uInt = results_lat.uInt;
    nInt = length(uInt);
mice = results_lat.mice;
    nMice = length(mice);
uMF = results_lat.uMF;
uMD = results_lat.uMD;

uMF_Beh = results_lat.uMF;
Thr_lat = results_lat.Thr_lat;

%%
[simPar,ia,~] = unique(SimRes(:,...
    {'input','condition','tempRes_ms',...
    'normalization','demean','useAllMd','bw'}));
    pp = find(... % use which PCA results 
        strcmp(PCARes.input,simPar.input)...
        & strcmp(PCARes.condition,simPar.condition)...
        & PCARes.tempRes_ms == simPar.tempRes_ms ...) < 2*eps(1/64) ...
        & PCARes.normalization == simPar.normalization ...
        & PCARes.demean == simPar.demean ...
        & PCARes.useAllMd == simPar.useAllMd ...
        );
    W = PCARes.V{pp}(:,1);
tempRes_ms = simPar.bw;
tempRes = num2str(tempRes_ms);
KSD = AMKSDs.(['KSD_all_',tempRes]);
tt = AMKSDs.(['tt_',tempRes]);
uAM = AMKSDs.uAM;

PCs = projectPC(KSD,W);
baseTime = [0.8,1];
bIdx = tt > baseTime(1) & tt < baseTime(2); % time index for baseline
PCs = PCs - mean(PCs(:,bIdx),2);

%% Figure 10A: PCA results

Periods = {'Baseline','Salicylate','Washout'};
tOffset = 1; % zero-th point of time axis

% uMd = unique(uAM.Md); uMd = nonzeros(uMd);
uMd = [.125,.25,.5,1];
nMd = length(uMd);
MdLbl = cellfun(@num2str, num2cell(uMd*100), repmat({'%.3g%%'},1,nMd),'UniformOutput', false);
% uMf = unique(uAM.Mf); uMf = nonzeros(uMf); 
uMf = [16,512]; 
nMf = length(uMf);
% uInt = unique(uAM.Intensity); 
uInt = 60;
nInt = length(uInt);
uCond = unique(uAM.condition); 
nCond = length(uCond);

% - define colors for different frequencies
Color = [1,0,0;... 
         ...0,0.5,0;...
         ...0,0,1;...
         0.7,0,0.7;...
         1,0.7,0;...
         0,0.7,0.7;...
         ];
Colors = nan([size(Color),nMd]);
for dd = 1:nMd
    Colors(:,:,dd) = 1 -(1 - Color) * ((dd)/nMd)^2;
end
i = 1;

fig = figure;
fig.Position = [200,200,1000,300];
fig.Units = 'pixels';
fig.Position(3:4) = 96/2.54*[22,5];

% Trial averaged PC1 for different conditions 
for ii = 1:nInt
    ax = [];
 for ff = 1:nMf
     for cc = 1:nCond
%          ax = [ax,subplot(nMf,nCond,cc+nCond*(ff-1))];
         ax = [ax,subplot(1,nMf*nCond,cc+nCond*(ff-1))];
         clear('ll');
         for dd = 1:nMd
             Mf = uMf(ff); Md = uMd(dd); Int = uInt(ii);
             idx = (uAM.Mf == Mf & uAM.Md == Md & uAM.Intensity == Int & strcmp(uAM.condition,uCond(cc)));
             ll(dd) = plot(tt-tOffset,PCs(idx,:),'Color',Colors(ff,:,dd),'LineWidth',1);hold on;
         end
         xlabel('Time re: AM transition (s)');
         xlim([min(tt),max(tt)]-tOffset);
         xline(0,'Color',[0.7,0.7,0.7]);
         xline(1,'Color',[0.7,0.7,0.7]);
         xline(2,'Color',[0.7,0.7,0.7]);
         legend(flip(ll),flip(MdLbl),'location','best','FontSize',8);
         legend('boxoff');
         %     xline(tStart,'Color',[0.3,0.3,0.3],'LineWidth',2);
         %     xline(tEnd,'Color',[0.3,0.3,0.3],'LineWidth',2);
         ylabel('Size of component');
%         title([num2str(Mf,'%d Hz'), ' - ',Periods{cc}]);
     end

 end
 linkaxes(ax);
 xlim(ax,[0.9,1.2]-tOffset)
 ylim(ax,[-50,350])
 set(ax,'FontSize',8)
%  sgtitle(num2str(uInt(ii),'%d dB'));
end

setPDFRes(fig)
saveas(fig,'Figures\Links\Figure10A_PCA_Component.pdf')
% saveas(fig,'Figures\PCA_Component.png')

%% Figure 10B,C,D:

detLvls = [110:10:170];
LineWidth = 1;
MarkerSize = 3;
FontSize = 8;
sIdx = ismember(SimRes.thr, detLvls) & abs(SimRes.NDTime - 0.15) < 0.005 & strcmpi(SimRes.condition,'baseline');
sIdx = find(sIdx)';
sColors = 0.9*hsv(length(sIdx)+1);
sColors(3,:) = []; % manual hack for hard to distinguish green colors
bColors = gray(length(sIdx)+1);

% Simulation results
fig = figure;
fig.Position = [200,200,1000,300];
fig.Units = 'pixels';
fig.Position(3:4) = 96/2.54*[20,10];
clear('ax');

% goodness of fits
[M_all,uThr,uNDTimes] = goodnessHeatmap(SimRes,{'thr','NDTime'},'CDFDiff');
[M_Base60,uThr,uNDTimes] = goodnessHeatmap(SimRes,{'thr','NDTime'},'CDFDiff_Base60');
[M_Sal60,uThr,uNDTimes] = goodnessHeatmap(SimRes,{'thr','NDTime'},'CDFDiff_Sal60');

cRange = [0.10,0.65];
% figure('Position',[10,300,1300,500]);
nCols = 6;
ii = 1;
ax(ii) = subplot(2,nCols,[ii,ii+nCols]);
imagesc(uNDTimes,uThr,M_all,'AlphaData',~isnan(M_all),cRange);hold on;
[~,minIdx] = min(M_all,[],'all','linear');
[minY,minX] = ind2sub(size(M_all),minIdx);
scatter(uNDTimes(minX),uThr(minY),'ow');
colormap(flipud(jet));
title('Overall');%colorbar;
ylabel('Detection level (a.u.)');
% xlabel('Non-decision time (s)')
ii = 2;
ax(ii) = subplot(2,nCols,[ii,ii+nCols]);
imagesc(uNDTimes,uThr,M_Base60,'AlphaData',~isnan(M_Base60),cRange);hold on;
[~,minIdx] = min(M_Base60,[],'all','linear');
[minY,minX] = ind2sub(size(M_Base60),minIdx);
scatter(uNDTimes(minX),uThr(minY),'ow');
colormap(flipud(jet));
title('Baseline');%colorbar;
xlabel('Non-decision time (s)')
ii = 3;
ax(ii) = subplot(2,nCols,[ii,ii+nCols]);
imagesc(uNDTimes,uThr,M_Sal60,'AlphaData',~isnan(M_Sal60),cRange);hold on;
[~,minIdx] = min(M_Sal60,[],'all','linear');
[minY,minX] = ind2sub(size(M_Sal60),minIdx);
scatter(uNDTimes(minX),uThr(minY),'ow');
colormap(flipud(jet));
title('Salicylate');%colorbar;
cb = colorbar(ax(ii),'east');
cb.Position = [0.05,0.15,0.01,0.75];
cb.Label.String = 'CDF difference';cb.Label.FontSize = FontSize;
% xlabel('Non-decision time (s)')



set(ax,'YDir','normal')
for ii = 1:3
scatter(ax(ii),repmat(0.15,1,length(detLvls)),detLvls,20,sColors,'filled','MarkerEdgeColor','k');
end


% thresholds

fIdx = [1,4];

AMThrRange = [-20,0];
% for ii = 2%1:2 % intensity
ii = 2;
intStr = '60';
% figure('Position',[300,200,600,600]);
set(gcf,'defaultAxesColorOrder',[0,0,0; 0,0,0]);
clearvars('ax','ax2');
for scnt = 1:length(sIdx)%ss-10:10:ss+50
    sss = sIdx(scnt);
    dP_SIM = SimRes.dP_SIM{sss};
    
    Thr_SIM = SimRes.(['Thr_SIM_',intStr]){sss};
    Thr_SIM = Thr_SIM(fIdx,:);
    for ff = 1:2 % frequency
        fff = fIdx(ff);
        ax(ff) = subplot(2,nCols,[4,5]+(ff-1)*nCols);
        plot(dP_SIM(fff,:,ii,1),'o--','Color',sColors(scnt,:),'LineWidth',LineWidth,'MarkerSize',MarkerSize); hold on;% 16Hz-allMods-45dB-baseline
        plot(dP_SIM(fff,:,ii,2),'s-','Color',sColors(scnt,:),'MarkerFaceColor',sColors(scnt,:),'LineWidth',LineWidth,'MarkerSize',MarkerSize) % 16Hz-allMods-45dB-salicylate
        
        text(0.05,1-(0.075*(scnt-1)),num2str(SimRes.thr(sss)),'Color',sColors(scnt,:),...
            'Units','normalized','VerticalAlignment','top','FontSize',FontSize)
        
        set(gca,'XTick',1:5,'XTickLabel',uMD);
        
        yline(1,'--','LineWidth',1.5);
        
        title([num2str(uMF(ff),'%dHz')]);%, ' ', num2str(uInt(ii),'%ddB')])
        ylabel('d''')
        ylim([-.7,4])
        xlabel('Mod. depth')
        xlim([0.5,5.5])
        
        ax2(ff) = subplot(2,nCols,[6]+(ff-1)*nCols);
        
        yyaxis(ax2(ff),'left');
        plot(ax2(ff),Thr_SIM(ff,:),'o-','Color',sColors(scnt,:),'MarkerFaceColor',sColors(scnt,:),'LineWidth',LineWidth,'MarkerSize',MarkerSize) ;hold on;%
        xlim([.5,2.5]);xticks([1,2]);xticklabels({'Baseline','Salicylate'});
        title([num2str(uMF(ff),'%dHz')]);%, ' ', num2str(uInt(ii),'%ddB')])
%         if(ff == 1); 
            ylabel('AM threshold (dB)'); 
%         end
% %         if(ff == 2); yticklabels([]);endNonDete
        ylim(AMThrRange);
        set(ax2(ff),'YDir','reverse')
        yyaxis(ax2(ff),'right');
        set(ax2(ff),'YScale','log');%,'YDir','reverse')
        set(ax2(ff),'YDir','reverse')
        ylim(ax2(ff),dB2a(AMThrRange));yticks(uMD);
%         if(ff == 1); yticklabels([]);end
%         if(ff == 2); 
            ylabel('AM threshold (m)'); 
%         end
    end
end
    for ff = 1:2 % frequency
        yyaxis(ax2(ff),'left');
        plot(ax2(ff),squeeze(mean(Thr_lat(ff,ii,:,1:2),3)),'k:s',...
            'LineWidth',2*LineWidth,...
            'MarkerSize',MarkerSize); %'MarkerFaceColor','k'
    end
set([ax,ax2],'FontSize',FontSize);
% end
setPDFRes(fig)
saveas(fig,'Figures\Links\Figure10B-D_PCA_SIMRes.pdf')
% saveas(fig,'Figures\PCA_SIMRes.png')


%% local functions

function a = dB2a(dB)
    a = 10.^(dB./20);
end