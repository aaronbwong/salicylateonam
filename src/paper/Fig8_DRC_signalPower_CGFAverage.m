clearvars;clc;close all;

%% Read data from DRC analysis results
% and compile into a table "data"
% columns include: sp,np,err,PRF,CGF,wtf_reg_params, wtauphi_reg_params,
%                   strfpp, ctxtpp


% STRF parameters
K = 61; J = 8;
freq = logspace(log10(2),log10(64),K);
delay = 5.*(1:J);
% CGF parameters
N = 12; M = 10;
phi = (-N:N)./12;
tau = 5.*(0:M);

addpath('.\functions\');
[data,pairTable] = DRC_load_data();

baseOnly = data.baseOnly;
sigSignal = data.sigSignal;
basePairSel = pairTable.basePairSel;
salPairSel = pairTable.salPairSel;

%% Fig 8A context improvement predictive power for vs STRF


fig = figure; LineWidth = 1;
fig.Position = [150,425,400,400];
sc(1) = scatter(data.strfnpp(baseOnly & sigSignal)*100,data.ctxtnpp(baseOnly & sigSignal)./data.strfnpp(baseOnly & sigSignal),...
    'o','MarkerEdgeColor',[0.7,0.7,0.7],...
    'LineWidth',LineWidth); hold on;
plot([data.strfnpp(basePairSel),data.strfnpp(salPairSel)]'*100,...
    [data.ctxtnpp(basePairSel)./data.strfnpp(basePairSel),...
    data.ctxtnpp(salPairSel)./data.strfnpp(salPairSel)]','-','Color',[0.5,0.5,0]); hold on;
sc(2) = scatter(data.strfnpp(basePairSel)*100,data.ctxtnpp(basePairSel)./data.strfnpp(basePairSel),...
    'ok',...
    'LineWidth',LineWidth); hold on;
sc(3) = scatter(data.strfnpp(salPairSel)*100,data.ctxtnpp(salPairSel)./data.strfnpp(salPairSel),...
    'o','MarkerEdgeColor',[0.7,0.7,0],...
    'LineWidth',LineWidth); 
sc = sc(1:3);
ylabel('pp_{context}/pp_{STRF}');
xlabel('Predictive power STRF (%)');
xlim([-2,70]);
ylim([0,3.5]);
yline(1,'--k');
legend(sc,{'baseline (excl)','baseline (incl)','salicylate'},'Location','northeast');

set(gca,'FontSize',10);
axis square

fig.Renderer = 'painter';
fig.PaperUnits = 'inches';
fig.PaperSize = fig.Position(3:4)./96; %96 dpi
saveas(gcf,'Figures\Links\Figure8A.pdf');

% 
% fig = figure; LineWidth = 1;
% fig.Position = [500,200,350,350];
% sc(1) = scatter(data.ctxtnpp(excSel & sigSignal)*100,data.ctxtnpp(excSel & sigSignal)./data.strfnpp(excSel & sigSignal),...
%     'o','MarkerEdgeColor',[0.7,0.7,0.7],...
%     'LineWidth',LineWidth); hold on;
% plot([data.ctxtnpp(basePairSel),data.ctxtnpp(salPairSel)]'*100,...
%     [data.ctxtnpp(basePairSel)./data.strfnpp(basePairSel),...
%     data.ctxtnpp(salPairSel)./data.strfnpp(salPairSel)]','-','Color',[0.5,0.5,0]); hold on;
% sc(2) = scatter(data.ctxtnpp(baseSel & sigSignal)*100,data.ctxtnpp(baseSel & sigSignal)./data.strfnpp(baseSel & sigSignal),...
%     'ok',...
%     'LineWidth',LineWidth); hold on;
% sc(3) = scatter(data.ctxtnpp(salSel & sigSignal)*100,data.ctxtnpp(salSel & sigSignal)./data.strfnpp(salSel & sigSignal),...
%     'o','MarkerEdgeColor',[0.7,0.7,0],...
%     'LineWidth',LineWidth); 
% sc = sc(1:3);
% ylabel('pp_{context}/pp_{STRF}');
% xlabel('Predictive power context (%)');
% xlim([-3,70]);
% ylim([0,3.5])
% yline(1,'--k');
% legend(sc,{'baseline (excl)','baseline (incl)','salicylate'},'Location','best');
% 
% axis square




%% calculate average CGF

allBaseCGF = cat(3,data.CGF{basePairSel});
avgBaseCGF = mean(allBaseCGF,3);
allBaseCGF_lin = reshape(allBaseCGF,[(M+1)*(2*N+1),size(allBaseCGF,3)]);
allBaseCGF_lin = allBaseCGF_lin - mean(allBaseCGF_lin);
allBaseCGF_lin = allBaseCGF_lin ./ std(allBaseCGF_lin);
[coeffBase, ~, ~, ~, explainedBase] = pca(allBaseCGF_lin');

allSalCGF = cat(3,data.CGF{salPairSel});
avgSalCGF = mean(allSalCGF,3);
allSalCGF_lin = reshape(allSalCGF,[(M+1)*(2*N+1),size(allSalCGF,3)]);
allSalCGF_lin = allSalCGF_lin - mean(allSalCGF_lin);
allSalCGF_lin = allSalCGF_lin ./ std(allSalCGF_lin);
[coeffSal, ~, ~, ~, explainedSal] = pca(allSalCGF_lin');

avgDifCGF = avgSalCGF - avgBaseCGF;
allDifCGF = allSalCGF - allBaseCGF;
allDifCGF_lin = reshape(allDifCGF,[(M+1)*(2*N+1),size(allDifCGF,3)]);
allDifCGF_lin = allDifCGF_lin - mean(allDifCGF_lin);
allDifCGF_lin = allDifCGF_lin ./ std(allDifCGF_lin);
[coeffDif, ~, ~, ~, explainedDif] = pca(allDifCGF_lin');

%% Fig 8B: average CGF figures
fig = figure('Position',[475,520,600,250]); colormap(jet)
cgfmask = ones(2*N+1,M+1);
cgfmask(N+1,1) = 0;
cMaxAll = max([abs(avgBaseCGF),abs(avgSalCGF)],[],'all');
subplot(1,2,1)
imagesc(avgBaseCGF', ...
    'AlphaData',cgfmask, ...
    [-cMaxAll,cMaxAll]);setCGFAxis(gca,'reverse')
colorbar; 
title('Baseline');
set(gca,'Color',[1,1,1]);

subplot(1,2,2)
imagesc(avgSalCGF', ...
    'AlphaData',cgfmask, ...
    [-cMaxAll,cMaxAll]);setCGFAxis(gca,'reverse')
colorbar;
title('Salicylate');
set(gca,'Color',[1,1,1]);

fig.Renderer = 'painter';
fig.PaperUnits = 'inches';
fig.PaperSize = fig.Position(3:4)./96; %96 dpi
saveas(fig,'Figures\Links\Figure8B.pdf');
%
%% Fig 8C-D: collapsed CGF

yrange = [-0.21,0.03];

nUnits = size(pairTable,1);
colRange = 1;
colBaseCGF = squeeze(mean(allBaseCGF(colRange,:,:),1));
colBaseCGF(N+1,:) = nan;
colSalCGF = squeeze(mean(allSalCGF(colRange,:,:),1));
colSalCGF(N+1,:) = nan;

rowRange = N+1; 
rowBaseCGF = (squeeze(mean(allBaseCGF(:,rowRange,:),2)));
rowBaseCGF(1) = NaN;
rowSalCGF = (squeeze(mean(allSalCGF(:,rowRange,:),2)));
rowSalCGF(1) = NaN;

fig = figure('Position',[475,350,600,167]); 
subplot(1,2,1)
errorbar(phi,mean(colBaseCGF'),std(colBaseCGF')./sqrt(nUnits),'Marker','s','Color','k','MarkerFaceColor','k'); hold on;
errorbar(phi,mean(colSalCGF'),std(colSalCGF')./sqrt(nUnits),'Marker','s','Color',[0.7,0.7,0],'MarkerFaceColor',[0.7,0.7,0]);
yline(0,'--k');
ylim(yrange)
xlim([-1.02,1.02])
lgd = legend({'baseline','salicylate'},'Box','off', ...
    'Units','pixels','Position',[175,44,111,32]);
%     'location','best');
ylabel('CGF weight')
xlabel('Phi (oct)')
set(gca,'FontSize',10);

subplot(1,2,2);
errorbar(tau,mean(rowBaseCGF'),std(rowBaseCGF')./sqrt(nUnits),'Marker','s','Color','k','MarkerFaceColor','k'); hold on;
errorbar(tau,mean(rowSalCGF'),std(rowSalCGF')./sqrt(nUnits),'Marker','s','Color',[0.7,0.7,0],'MarkerFaceColor',[0.7,0.7,0]);
yline(0,'--k');
ylim(yrange)
xlim([0,55])
legend({'baseline','salicylate'},'Box','off', ...
    'location','best');
xlabel('Tau (ms)');
set(gca,'FontSize',10);

fig.Renderer = 'painter';
fig.PaperUnits = 'inches';
fig.PaperSize = fig.Position(3:4)./96; %96 dpi
saveas(fig,'Figures\Links\Figure8CD.pdf');
%% REVISION
allBaseSTRF = cat(3,data.STRF{pairTable.basePairSel});
allBasePRF = cat(3,data.PRF{pairTable.basePairSel});
allSalSTRF = cat(3,data.STRF{pairTable.salPairSel});
allSalPRF = cat(3,data.PRF{pairTable.salPairSel});
meanWeightBase = squeeze(mean(allBaseSTRF,[1,2]));
meanWeightSal =  squeeze(mean(allSalSTRF,[1,2]));
meanPRFWeightBase = squeeze(mean(allBasePRF,[1,2]));
meanPRFWeightSal =  squeeze(mean(allSalPRF,[1,2]));

excUnits = meanWeightBase > 0;
inhUnits = meanWeightBase < 0;
excUnits_Sal = meanWeightSal > 0;
inhUnits_Sal = meanWeightSal < 0;

avgExcBaseCGF = mean(allBaseCGF(:,:,excUnits),3);
avgInhBaseCGF = mean(allBaseCGF(:,:,inhUnits),3);
avgExcSalCGF = mean(allSalCGF(:,:,excUnits),3);
avgInhSalCGF = mean(allSalCGF(:,:,inhUnits),3);


fig = figure('Position',[275,220,600,500]); colormap(jet)
cgfmask = ones(2*N+1,M+1);
cgfmask(N+1,1) = 0;
cMaxAll = max([abs(avgExcBaseCGF),abs(avgInhBaseCGF),abs(avgExcSalCGF),abs(avgInhSalCGF)],[],'all');

% Baseline
  % excitatory
subplot(2,2,1)
imagesc(avgExcBaseCGF', ...
    'AlphaData',cgfmask, ...
    [-cMaxAll,cMaxAll]);setCGFAxis(gca,'reverse')
colorbar; 
title('Baseline excitatory');
set(gca,'Color',[1,1,1]);

  % inhibitory
subplot(2,2,3)
imagesc(avgInhBaseCGF', ...
    'AlphaData',cgfmask, ...
    [-cMaxAll,cMaxAll]);setCGFAxis(gca,'reverse')
colorbar;
title('Baseline inhibitory');
set(gca,'Color',[1,1,1]);

% Salicylate
  % excitatory
subplot(2,2,2)
imagesc(avgExcSalCGF', ...
    'AlphaData',cgfmask, ...
    [-cMaxAll,cMaxAll]);setCGFAxis(gca,'reverse')
colorbar; 
title('Salicylate excitatory');
set(gca,'Color',[1,1,1]);
  % inhibitory
subplot(2,2,4)
imagesc(avgInhSalCGF', ...
    'AlphaData',cgfmask, ...
    [-cMaxAll,cMaxAll]);setCGFAxis(gca,'reverse')
colorbar;
title('Salicylate inhibitory');
set(gca,'Color',[1,1,1]);


% ---------------- COLLAPSED CGF ------------------
yrange = [-0.26,0.06];

nExcUnits = sum(excUnits);
nInhUnits = sum(inhUnits);
colRange = 1;
colExcBaseCGF = squeeze(mean(allBaseCGF(colRange,:,excUnits),1));
colExcBaseCGF(N+1,:) = nan;
colExcSalCGF = squeeze(mean(allSalCGF(colRange,:,excUnits),1));
colExcSalCGF(N+1,:) = nan;

colInhBaseCGF = squeeze(mean(allBaseCGF(colRange,:,inhUnits),1));
colInhBaseCGF(N+1,:) = nan;
colInhSalCGF = squeeze(mean(allSalCGF(colRange,:,inhUnits),1));
colInhSalCGF(N+1,:) = nan;

rowRange = N+1; 
rowExcBaseCGF = (squeeze(mean(allBaseCGF(:,rowRange,excUnits),2)));
rowExcBaseCGF(1) = NaN;
rowExcSalCGF = (squeeze(mean(allSalCGF(:,rowRange,excUnits),2)));
rowExcSalCGF(1) = NaN;
rowInhBaseCGF = (squeeze(mean(allBaseCGF(:,rowRange,inhUnits),2)));
rowInhBaseCGF(1) = NaN;
rowInhSalCGF = (squeeze(mean(allSalCGF(:,rowRange,inhUnits),2)));
rowInhSalCGF(1) = NaN;

%%
fig = figure('Position',[475,150,500,667]); 

subplot(2,1,1)
errorbar(phi,mean(colExcBaseCGF'),std(colExcBaseCGF')./sqrt(nExcUnits),'Marker','s','Color','k','MarkerFaceColor','k'); hold on;
errorbar(phi,mean(colExcSalCGF'),std(colExcSalCGF')./sqrt(nExcUnits),'Marker','s','Color',[0.7,0.7,0],'MarkerFaceColor',[0.7,0.7,0]);
errorbar(phi,mean(colInhBaseCGF'),std(colInhBaseCGF')./sqrt(nInhUnits),'Marker','o','Color','k','MarkerFaceColor','none','LineStyle',':'); hold on;
errorbar(phi,mean(colInhSalCGF'),std(colInhSalCGF')./sqrt(nInhUnits),'Marker','o','Color',[0.7,0.7,0],'MarkerFaceColor','none','LineStyle',':');
yline(0,'--k');
ylim(yrange)
xlim([-1.02,1.02])
lgd = legend({['baseline (exc; N = ',num2str(nExcUnits),')'],'salicylate (exc)', ['baseline (inh; N = ',num2str(nInhUnits),')'],'salicylate (inh)'},'Box','off', ...
    ...'Units','pixels',...'Position',[175,44,111,32]);
    'location','best');
ylabel('CGF weight')
xlabel('Phi (oct)')
set(gca,'FontSize',10);


subplot(2,1,2);
errorbar(tau,mean(rowExcBaseCGF'),std(rowExcBaseCGF')./sqrt(nExcUnits),'Marker','s','Color','k','MarkerFaceColor','k'); hold on;
errorbar(tau,mean(rowExcSalCGF'),std(rowExcSalCGF')./sqrt(nExcUnits),'Marker','s','Color',[0.7,0.7,0],'MarkerFaceColor',[0.7,0.7,0]);
errorbar(tau,mean(rowInhBaseCGF'),std(rowInhBaseCGF')./sqrt(nInhUnits),'Marker','o','Color','k','MarkerFaceColor','none','LineStyle',':'); hold on;
errorbar(tau,mean(rowInhSalCGF'),std(rowInhSalCGF')./sqrt(nInhUnits),'Marker','o','Color',[0.7,0.7,0],'MarkerFaceColor','none','LineStyle',':');
yline(0,'--k');
ylim(yrange)
xlim([0,55])
legend({'baseline','salicylate'},'Box','off', ...
    'location','best');
ylabel('CGF weight')
xlabel('Tau (ms)');
set(gca,'FontSize',10);

%%
close all;
% showSTRFs(allBasePRF(:,:,excUnits),num2cell(excUnits_Sal(excUnits)));sgtitle('PRF Base Exc')

showSTRFs(allBaseSTRF(:,:,excUnits))%,num2cell(excUnits_Sal(excUnits)));
sgtitle('Base Exc')
showSTRFs(allBaseSTRF(:,:,inhUnits))%,num2cell(excUnits_Sal(inhUnits)));
sgtitle('Base Inh')
showSTRFs(allSalSTRF(:,:,excUnits))%,num2cell(excUnits_Sal(excUnits)));
sgtitle('Sal Exc')
showSTRFs(allSalSTRF(:,:,inhUnits))%,num2cell(excUnits_Sal(inhUnits)));
sgtitle('Sal Inh')

%%
% required code from DRC_STRF_SalicylateCompare
salColor = [.7,.7,0];
baseColor = [0,0,0];
FontSize = 12;
fig = figure;
fig.Position = [710,100,600,300];
errorbar(mean(Tuning_Base(:,excUnits)./max(Tuning_Base(:,excUnits)),2,'omitnan'),...
    std(Tuning_Base(:,excUnits)./max(Tuning_Base(:,excUnits))./sqrt(nExcUnits),[],2,'omitnan'),...
    'MarkerFaceColor','k','Color',baseColor,'MarkerSize',6,'Marker','s'); hold on;
errorbar(mean(Tuning_Base(:,inhUnits)./max(Tuning_Base(:,inhUnits)),2,'omitnan'),...
    std(Tuning_Base(:,inhUnits)./max(Tuning_Base(:,inhUnits))./sqrt(nInhUnits),[],2,'omitnan'),...
    'MarkerFaceColor','none','Color',baseColor,'MarkerSize',6,'Marker','o','LineStyle',':'); hold on;
errorbar(mean(Tuning_Sal2(:,excUnits)./max(Tuning_Sal2(:,excUnits)),2,'omitnan'),...
    std(Tuning_Sal2(:,excUnits)./max(Tuning_Sal2(:,excUnits))./sqrt(nExcUnits),[],2,'omitnan'),...
    'MarkerFaceColor',salColor,'Color',salColor,'MarkerSize',6,'Marker','s'); hold on;
errorbar(mean(Tuning_Sal2(:,inhUnits)./max(Tuning_Sal2(:,inhUnits)),2,'omitnan'),...
    std(Tuning_Sal2(:,inhUnits)./max(Tuning_Sal2(:,inhUnits))./sqrt(nInhUnits),[],2,'omitnan'),...
    'MarkerFaceColor','none','Color',salColor,'MarkerSize',6,'Marker','o','LineStyle',':'); hold on;
ylabel('Normalized STRF weight');
xticks(1:bandwidth/2:(2*bandwidth+1));
xticklabels([-bandwidth:bandwidth/2:bandwidth]./12);
xlabel('Freq. re. BF (oct)')
set(gca,'FontSize',FontSize);
% legend({'baseline','salicylate'},'location','best')
lgd = legend({['baseline (exc; N = ',num2str(nExcUnits),')'],'salicylate (exc)', ['baseline (inh; N = ',num2str(nInhUnits),')'],'salicylate (inh)'},'Box','off', ...
    'location','best');

%% requires ExcInhAnalysis
nUnits = size(pairTable,1);
for ii = 1:nUnits
    idx = salUnits.mouse == pairTable.mouse(ii) & salUnits.cids == pairTable.cids(ii);
    pairTable.CF_base(ii) = salUnits.CF_Base(idx);
    pairTable.CF_Sal(ii) = salUnits.CF_Sal(idx);
    pairTable.minThr_Base(ii) = salUnits.minThr_Base(idx);
    pairTable.minThr_Sal(ii) = salUnits.minThr_Sal(idx);
    pairTable.ZeroRate_Base(ii) = salUnits.ZeroRate_Base(idx);
    pairTable.ZeroRate_Sal(ii) = salUnits.ZeroRate_Sal(idx);
end

disp('Spont. firing shift')
disp([' excited: ', num2str(mean(pairTable.ZeroRate_Sal(excUnits) - pairTable.ZeroRate_Base(excUnits)),'%.2f'),...
    ' +/- ', num2str(std(pairTable.ZeroRate_Sal(excUnits) - pairTable.ZeroRate_Base(excUnits)),'%.2f')]);
disp([' inhibited: ', num2str(mean(pairTable.ZeroRate_Sal(inhUnits) - pairTable.ZeroRate_Base(inhUnits)),'%.2f'),...
    ' +/- ', num2str(std(pairTable.ZeroRate_Sal(inhUnits) - pairTable.ZeroRate_Base(inhUnits)),'%.2f')]);

disp('Threshold shift')
disp([' excited: ', num2str(mean(pairTable.minThr_Sal(excUnits) - pairTable.minThr_Base(excUnits),'omitnan'),'%.2f'),...
    ' +/- ', num2str(std(pairTable.minThr_Sal(excUnits) - pairTable.minThr_Base(excUnits),'omitnan'),'%.2f')]);
disp([' inhibited: ', num2str(mean(pairTable.minThr_Sal(inhUnits) - pairTable.minThr_Base(inhUnits),'omitnan'),'%.2f'),...
    ' +/- ', num2str(std(pairTable.minThr_Sal(inhUnits) - pairTable.minThr_Base(inhUnits),'omitnan'),'%.2f')]);

%% LOCAL Functions
function setPRFAxis(h,xdir)
    if nargin < 2; xdir = 'normal';end
    K = 61; J = 8;
    freq = logspace(log10(2),log10(64),K);
    delay = 5.*(0:J);
    set(h,'YDir','normal','XDir',xdir)
    set(h,'ytick',1:12:K,'yticklabel',freq(1:12:K));
    xlabel('Delay (ms)');
    set(h,'xtick',1:2:J+1,'xticklabel',delay(1:2:J+1));
    ylabel('Frequency (kHz)');
end

function setCGFAxis(h,xdir)
    if nargin < 2; xdir = 'normal';end
    N = 12; M = 10;
    phi = (-N:N)./12;
    tau = 5.*(0:M);
    set(h,'YDir','normal','XDir',xdir)
    set(h,'ytick',1:N/2:2*N+1,'yticklabel',phi(1:N/2:end));
    set(h,'xtick',1:5:M+1,'xticklabel',tau(1:5:M+1));
    xlabel('Tau (ms)');
    ylabel('Phi (oct)');
%     rectangle('Position',[0.5,N+0.5,1,1],'EdgeColor',[.5 .5 .5],'FaceColor',[.5 .5 .5],'LineWidth',0.05)
%     rectangle('Position',[0.5,N+0.5,1,1],'EdgeColor','none','FaceColor',[.5 .5 .5])
end

function showSTRFs(STRFs,labels)
    nUnits = size(STRFs,3);
    row = floor(sqrt(nUnits));
    col = ceil(nUnits / row);
    rowHeight = 180;
    fig = figure('Position',[100,100,col*rowHeight,row*rowHeight]); 
%     cMax = max(abs(STRFs),[],'all');
    for ii = 1:nUnits
        strf = squeeze(STRFs(:,:,ii))';
        cMax = max(abs(strf),[],'all');
        sp = subplot(row,col,ii);
        imagesc(strf,[-cMax,cMax]);
        colormap(sp,"jet");
        if (nargin > 1)
            if (isstring(labels{ii}))
                title(labels{ii});
            else                       
                title(num2str(labels{ii}));
            end
        end
        setPRFAxis(sp,'reverse')
    end
end