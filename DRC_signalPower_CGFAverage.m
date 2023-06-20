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

% context hyper parameters
context_reg_param = [6,0.5,0.5]';
optstyle = 'ASD';
scale = 'amp';
loadpath = ['Data\Neuro\DRC\',...
    optstyle,'_',...
    replace(mat2str(context_reg_param),{';',',','.','[',']'},{'_','_','x','',''}),...
    '\'];

Mice = [25,26,27,28,29,30,32,38,42,43,49,50,52,57,71,72,73,74];

% read metadata for units
salUnits = readtable('Data\salicylate_units.xlsx');

% initiate progress bar
% progFig = uifigure('Position',[450,400,500,200]);
% d = uiprogressdlg(progFig,'Title',['Analysis started at ',datestr(now,'hh:MM:ss')],...
%         'Message','Please wait...');
    sp = []; np = []; err = [];
    init = 1;

% loop that load data in for each mouse
for m = 1:length(Mice)
    Mouse = num2str(Mice(m),'M%02.0f');
%     d.Value = m/length(Mice);
%     d.Message = 'Loading data';
    load(['Data\Neuro\Spks\', Mouse, '\', Mouse '_SpkS.mat'],'cids','cpos');
    cids = cids';cpos = cpos';
    load(['Data\Neuro\DRC\', Mouse, '_DRC_PSTH.mat'],'DRCSets');
    
    setNums = DRCSets;
    nSets = length(DRCSets);

    for s = 1:nSets
        setNum = setNums(s);   
        load([loadpath,'Results_',Mouse,'_S',num2str(setNum,'%02d'),'_J',num2str(J),'_M',num2str(M),'_',optstyle,'_',scale],'results');   

        nClu = length(results);
        
        salDataSet = nan(nClu,1);
        sp = nan(nClu,1);
        np = nan(nClu,1);
        err = nan(nClu,1);
        c_strf = nan(nClu,1);
        STRF = cell(nClu,1);
        c_ctxt = nan(nClu,1);
        PRF = cell(nClu,1);
        CGF = cell(nClu,1);
        wtf_reg_params = cell(nClu,1);
        wtauphi_reg_params = cell(nClu,1);
        strftpp = nan(nClu,1);
        ctxttpp = nan(nClu,1);
        strfpp = nan(nClu,1);
        ctxtpp = nan(nClu,1);
        for ii = 1:nClu
            salDataSet(ii) = any(salUnits.mouse == Mice(m) & salUnits.cids == cids(ii));              
            sp(ii) = [results(ii).stats.signalpower];
            np(ii) = [results(ii).stats.noisepower];
            err(ii) = [results(ii).stats.error];
            c_strf(ii) = results(ii).strf.ww(1);
            STRF{ii} = flipud(reshape(results(ii).strf.ww(2:end),J,K));
            c_ctxt(ii) = results(ii).full_rank_sparse_rep.c;
            PRF{ii} = results(ii).full_rank_sparse_rep.wtf;
            CGF{ii} = results(ii).full_rank_sparse_rep.wtauphi;
            wtf_reg_params{ii} = results(ii).full_rank_sparse_rep.wtf_reg_params;
            wtauphi_reg_params{ii} = results(ii).full_rank_sparse_rep.wtauphi_reg_params;
            strftpp(ii) = results(ii).strf.tpp;
            ctxttpp(ii) = results(ii).full_rank_sparse_rep.tpp;
            strfpp(ii) = results(ii).strf.pp;
            ctxtpp(ii) = results(ii).full_rank_sparse_rep.pp;
        end
        mouse = repmat(Mice(m),nClu,1);
        set = repmat(setNum,nClu,1);
        data_temp = table(mouse,cids,cpos,set,salDataSet,sp,np,err,c_strf,STRF,c_ctxt,PRF,CGF,wtf_reg_params,wtauphi_reg_params,strftpp,ctxttpp,strfpp,ctxtpp);
        if init == 1
            data = data_temp;init = 0;
        else
            data = [data;data_temp];
        end
    end
    
end
clear('results','ans','cids','cpos','salDataSet','sp','np','err','STRF','PRF','CGF','wtf_reg_params','wtauphi_reg_params','strfpp','ctxtpp');
clear('strftpp','ctxttpp','c_ctxt','c_strf');
clear('set', 'data_temp','mouse','Mouse','nClu','nSets','setNum','setNums','DRCSets');
clear('m','ii', 's','init');
% close(progFig);clear('progFig','d');

% some calculations
data.nnp = data.np./data.sp;
data.strfntpp = data.strftpp./data.sp;
data.ctxtntpp = data.ctxttpp./data.sp;
data.strfnpp = data.strfpp./data.sp;
data.ctxtnpp = data.ctxtpp./data.sp;

nClu = size(data,1);
for ii = 1:nClu
    data.rho(ii) = corr(reshape(data.STRF{ii},[],1),reshape(data.PRF{ii},[],1));
end
%% define parameters & units selection


% units*condition selections
excSel = data.set == 4 & (data.salDataSet == 0 ...
    | data.mouse == 25 | data.mouse == 74); %baseline of units without salicylate
baseSel =  data.set == 4 & data.salDataSet == 1 ...
    & data.mouse ~= 25 & data.mouse ~= 74 ; %baseline of units with salicylate
salSel = data.set == 10& data.salDataSet == 1; %salicylate 
furSel = data.set == 10& data.salDataSet == 0; %furosemide 

% enough signal power
nSigma = 2;
incSel = data.sp > nSigma* data.err;

% select relevant units
baseTable = data(baseSel & incSel,{'mouse','cids'}); % units with acceptable baseline recording
salTable = data(salSel & incSel,{'mouse','cids'}); % units with acceptable salicylate recording
pairTable = intersect(baseTable,salTable); % units with acceptable baseline and salicylate recordings
pairSel = ismember(data(:,{'mouse','cids'}),pairTable); % index of recordings from these units
basePairSel = pairSel & baseSel; % baseline recordings to include
salPairSel = pairSel & salSel;  % salicylate recordings to include
baseTable2 = data(basePairSel,{'mouse','cids'}); % units with acceptable baseline recording
salTable2 = data(salPairSel,{'mouse','cids'}); % units with acceptable salicylate recording
if (all(baseTable2.mouse == salTable2.mouse) && all(baseTable2.cids == salTable2.cids));disp('Order OK. Proceeding with comparison');end

% CGF improvements
pairTable.strfnpp_base = data.strfnpp(basePairSel);
pairTable.ctxtnpp_base = data.ctxtnpp(basePairSel);
pairTable.baseImpr = pairTable.ctxtnpp_base - pairTable.strfnpp_base;
pairTable.strfnpp_sal = data.strfnpp(salPairSel);
pairTable.ctxtnpp_sal = data.ctxtnpp(salPairSel);
pairTable.salImpr = pairTable.ctxtnpp_sal - pairTable.strfnpp_sal;

pairTable.c_base = data.c_ctxt(basePairSel);
pairTable.c_sal = data.c_ctxt(salPairSel);

% Similarity of PRF and CGF upon salicylate
nClu = size(pairTable,1);
for ii = 1:nClu
    sel = data.mouse == pairTable.mouse(ii) & ...
        data.cids == pairTable.cids(ii);
    sel1 = sel & data.set == 4 ;
    sel2 = sel & data.set == 10;
    pairTable.rho_STRFs(ii) = corr(reshape(data.STRF{sel1},[],1), reshape(data.STRF{sel2},[],1));
    pairTable.rho_PRFs(ii) = corr(reshape(data.PRF{sel1},[],1), reshape(data.PRF{sel2},[],1));
    pairTable.rho_CGFs(ii) = corr(reshape(data.CGF{sel1},[],1), reshape(data.CGF{sel2},[],1));
end

%% Fig 10A context improvement predictive power for vs STRF
fig = figure; LineWidth = 1;
fig.Position = [150,425,400,400];
sc(1) = scatter(data.strfnpp(excSel & incSel)*100,data.ctxtnpp(excSel & incSel)./data.strfnpp(excSel & incSel),...
    'o','MarkerEdgeColor',[0.7,0.7,0.7],...
    'LineWidth',LineWidth); hold on;
plot([data.strfnpp(basePairSel),data.strfnpp(salPairSel)]'*100,...
    [data.ctxtnpp(basePairSel)./data.strfnpp(basePairSel),...
    data.ctxtnpp(salPairSel)./data.strfnpp(salPairSel)]','-','Color',[0.5,0.5,0]); hold on;
sc(2) = scatter(data.strfnpp(baseSel & incSel)*100,data.ctxtnpp(baseSel & incSel)./data.strfnpp(baseSel & incSel),...
    'ok',...
    'LineWidth',LineWidth); hold on;
sc(3) = scatter(data.strfnpp(salSel & incSel)*100,data.ctxtnpp(salSel & incSel)./data.strfnpp(salSel & incSel),...
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
saveas(gcf,'Figures\Links\Figure10A.pdf');

% 
% fig = figure; LineWidth = 1;
% fig.Position = [500,200,350,350];
% sc(1) = scatter(data.ctxtnpp(excSel & incSel)*100,data.ctxtnpp(excSel & incSel)./data.strfnpp(excSel & incSel),...
%     'o','MarkerEdgeColor',[0.7,0.7,0.7],...
%     'LineWidth',LineWidth); hold on;
% plot([data.ctxtnpp(basePairSel),data.ctxtnpp(salPairSel)]'*100,...
%     [data.ctxtnpp(basePairSel)./data.strfnpp(basePairSel),...
%     data.ctxtnpp(salPairSel)./data.strfnpp(salPairSel)]','-','Color',[0.5,0.5,0]); hold on;
% sc(2) = scatter(data.ctxtnpp(baseSel & incSel)*100,data.ctxtnpp(baseSel & incSel)./data.strfnpp(baseSel & incSel),...
%     'ok',...
%     'LineWidth',LineWidth); hold on;
% sc(3) = scatter(data.ctxtnpp(salSel & incSel)*100,data.ctxtnpp(salSel & incSel)./data.strfnpp(salSel & incSel),...
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

%% Fig 10B: average CGF figures
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
saveas(fig,'Figures\Links\Figure10B.pdf');
%
%% Fig 10C-D: collapsed CGF

yrange = [-0.21,0.03];

nUnits = sum(basePairSel);
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
saveas(fig,'Figures\Links\Figure10CD.pdf');

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
    ylabel('Freq (kHz)');
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