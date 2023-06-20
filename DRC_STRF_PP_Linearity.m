Mice = [25,26,27,28,29,30,32,38,42,43,49,50,52,57,71,72,73,74];
scales = {'pow','amp','dB40','dB2','dB','bin'};nScales = length(scales);
scaleLbl = {'pow','amp','dB40','dB20','dB0','binary'};
optstyle = 'ASD';

ctxtColor1 = [.5,0.7,0.5];
ctxtColor2 = [0,0.5,0];
ctxtColor3 = [0,0.25,0];

K = 61; J = 8;
freq = logspace(log10(2),log10(64),K);
delay = 5.*(1:J);
% CGF parameters
N = 12; M = 10;
phi = (-N:N)./12;
tau = 5.*(0:M);

context_reg_param = [6,0.5,0.5]';
loadpath = ['Data\Neuro\DRC\',optstyle,'_',...
    replace(mat2str(context_reg_param),{';',',','.','[',']'},{'_','_','x','',''}),...
    '\'];


% read metadata for units
salUnits = readtable('Data\salicylate_units.xlsx');

% initiate progress bar
progFig = uifigure('Position',[450,400,500,200]);
d = uiprogressdlg(progFig,'Title',['Analysis started at ',datestr(now,'hh:MM:ss')],...
        'Message','Please wait...');
    sp = []; np = []; err = [];
    init = 1;
    
% loop that load data in for each mouse
for m = 1:length(Mice)
    Mouse = num2str(Mice(m),'M%02.0f');
    d.Value = m/length(Mice);
    d.Message = 'Loading data';
    load(['Data\Neuro\Spks\', Mouse, '\', Mouse '_SpkS.mat'],'cids','cpos');
    cids = cids';cpos = cpos';
    load(['Data\Neuro\DRC\', Mouse, '_DRC_PSTH.mat'],'DRCSets');
    
    setNums = DRCSets;
    nSets = length(DRCSets);

    for s = 1:nSets
        setNum = setNums(s);   
            load([loadpath,'Results_',Mouse,'_S',num2str(setNum,'%02d'),'_J',num2str(J),'_M',num2str(M),'_',optstyle,'_',scales{1}],'results');   
        nClu = length(results);
        
        salDataSet = nan(nClu,1);
        sp = nan(nClu,1);
        np = nan(nClu,1);
        err = nan(nClu,1);
        strftpp = nan(nClu,nScales);
        ctxttpp = nan(nClu,nScales);
        strfpp = nan(nClu,nScales);
        ctxtpp = nan(nClu,nScales);
        for sc = 1:nScales
            try
                load([loadpath,'Results_',Mouse,'_S',num2str(setNum,'%02d'),'_J',num2str(J),'_M',num2str(M),'_',optstyle,'_',scales{sc}],'results');   
            catch
                continue
            end
%         load(['ASD_0x1_0x5_0x5\Results_',Mouse,'_S',num2str(setNum,'%02d'),'_',optstyle,'_',scales{sc}],'results');   
            for ii = 1:nClu
                salDataSet(ii) = any(salUnits.mouse == Mice(m) & salUnits.cids == cids(ii));              
                sp(ii) = [results(ii).stats.signalpower];
                np(ii) = [results(ii).stats.noisepower];
                err(ii) = [results(ii).stats.error];
                strftpp(ii,sc) = results(ii).strf.tpp;
                ctxttpp(ii,sc) = results(ii).full_rank_sparse_rep.tpp;
                strfpp(ii,sc) = results(ii).strf.pp;
                ctxtpp(ii,sc) = results(ii).full_rank_sparse_rep.pp;
            end
        end
        mouse = repmat(Mice(m),nClu,1);
        set = repmat(setNum,nClu,1);
        data_temp = table(mouse,cids,cpos,set,salDataSet,sp,np,err,strftpp,ctxttpp,strfpp,ctxtpp);
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
close(progFig);clear('progFig','d');

% some calculations
data.nnp = data.np./data.sp;
data.strfntpp = data.strftpp./data.sp;
data.ctxtntpp = data.ctxttpp./data.sp;
data.strfnpp = data.strfpp./data.sp;
data.ctxtnpp = data.ctxtpp./data.sp;

nClu = size(data,1);
% for ii = 1:nClu
%     data.rho(ii) = corr(reshape(data.STRF{ii},[],1),reshape(data.PRF{ii},[],1));
% end

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

% % CGF improvements
pairTable.strfnpp_base = data.strfnpp(basePairSel,:);
pairTable.ctxtnpp_base = data.ctxtnpp(basePairSel,:);
pairTable.baseImpr = pairTable.ctxtnpp_base - pairTable.strfnpp_base;
pairTable.strfntpp_base = data.strfntpp(basePairSel,:);
pairTable.ctxtntpp_base = data.ctxtntpp(basePairSel,:);
pairTable.strfnpp_sal = data.strfnpp(salPairSel,:);
pairTable.ctxtnpp_sal = data.ctxtnpp(salPairSel,:);
pairTable.salImpr = pairTable.ctxtnpp_sal - pairTable.strfnpp_sal;
pairTable.strfntpp_sal = data.strfntpp(salPairSel,:);
pairTable.ctxtntpp_sal = data.ctxtntpp(salPairSel,:);

% 
% pairTable.c_base = data.c_ctxt(basePairSel);
% pairTable.c_sal = data.c_ctxt(salPairSel);

%% Figure 11 A-D 
col = 3;
row = 3;
sc = 3; % 3: 'amp' scale

fig = figure;
fig.Position = [100,100,1000,600];
subplot(row,col,[1,col+1])
nClu = sum(basePairSel);
traces = plot([data.ctxtnpp(basePairSel,sc),data.strfnpp(basePairSel,sc)]'*100,...
    [data.ctxtnpp(salPairSel,sc),data.strfnpp(salPairSel,sc)]'*100,'-',...
    'Color',[.5,.7,.5],'LineWidth',1); hold on;
strfDots = plot(data.strfnpp(basePairSel,sc)*100,data.strfnpp(salPairSel,sc)*100,'o',...
    'Color',[0.5,0.5,0.5],'LineWidth',1.5,'MarkerSize',5);hold on;
strfMean = errorbar(mean(data.strfnpp(basePairSel,sc)*100),... %x
         mean(data.strfnpp(salPairSel,sc)*100),... %y
         std(data.strfnpp(salPairSel,sc)*100/sqrt(nClu)),... %yneg
         std(data.strfnpp(salPairSel,sc)*100/sqrt(nClu)),... %ypos
         std(data.strfnpp(basePairSel,sc)*100/sqrt(nClu)),... %xneg
         std(data.strfnpp(basePairSel,sc)*100/sqrt(nClu)),... %xpos
         's','MarkerFaceColor',[.2,.2,.2],'Color',[.2,.2,.2],'MarkerSize',9,...
         'LineWidth',1);
ctxtDots = plot(data.ctxtnpp(basePairSel,sc)*100,data.ctxtnpp(salPairSel,sc)*100,'o',...
    'Color',[0,0.5,0],'MarkerFaceColor','none','LineWidth',1.5,'MarkerSize',5);hold on;
ctxtMean = errorbar(mean(data.ctxtnpp(basePairSel,sc)*100),... %x
         mean(data.ctxtnpp(salPairSel,sc)*100),... %y
         std(data.ctxtnpp(salPairSel,sc)*100/sqrt(nClu)),... %yneg
         std(data.ctxtnpp(salPairSel,sc)*100/sqrt(nClu)),... %ypos
         std(data.ctxtnpp(basePairSel,sc)*100/sqrt(nClu)),... %xneg
         std(data.ctxtnpp(basePairSel,sc)*100/sqrt(nClu)),... %xpos
         's','MarkerFaceColor',[0,0.25,0],'Color',[0,0.25,0],'MarkerSize',9,...
         'LineWidth',1);
ylabel('Predictive power salicylate (%)');
xlabel('Predictive power baseline (%)');
ylim([0,70]);xlim([0,70]);
plot([-10,100],[-10,100],'k:');
legend([ctxtDots,ctxtMean,strfDots,strfMean],{'Context','(mean+/-sem)','STRF','(mean+/-sem)'},'location','southeast')
axis square;


 
% demonstration of different scales
ax0 = subplot(row,col,2*col+1);
stim = [-Inf,25:5:70];
xx = [0,2:length(stim)];
scales = {'pow','amp','dB40','dB20','dB0','binary'};
Mrkrs = {'s','o','d','^','v','x'};
for ss = 1:length(scales)
   plot(xx,stimScale(stim,scales{ss}),'-k','LineWidth',1,'Marker',Mrkrs{ss});hold on; 
end

set(ax0,'xtick',xx,'xticklabel',num2cell(stim))
% ax0.XTickLabel{1} = 'no sound';
xlabel('Sound intensity (dB)')
ylabel('Stimulus strength');
% legend(scales,'location','northoutside','NumColumns',3)
legend(scales,'location','none','Position',[0.1,0.35,0.25,.05],'NumColumns',3);
legend('boxoff')

% Predictive power of STRF & CTXT models as function of scales
% figure;
maxnpp = 100 * max([pairTable.strfnpp_base,pairTable.strfnpp_sal],[],'all');
maxnpp = max(maxnpp,100 * max([pairTable.ctxtnpp_base,pairTable.ctxtnpp_sal],[],'all'));
[max1,idx1] = max(pairTable.strfnpp_base');
[~,idx1c] = max(pairTable.ctxtnpp_base');

% histogram of best scale (Baseline)
% ax3 = subplot(3,2,5);
ax3 = subplot(row,col,2*col+2);
histogram(idx1,[0.5+(0:nScales)],'FaceColor','k');hold on;
histogram(idx1c,[0.5+(0:nScales)],'FaceColor',ctxtColor2);
xticks(1:nScales);xticklabels(scaleLbl);xlabel('Scale');
ylabel('Best model (# units)');
[max2,idx2] = max(pairTable.strfnpp_sal');
[~,idx2c] = max(pairTable.ctxtnpp_sal');

% histogram of best scale (Salicylate)
% ax4 = subplot(3,2,6);
ax4 = subplot(row,col,2*col+3);
histogram(idx2,[0.5+(0:nScales)],'FaceColor','k');hold on;
histogram(idx2c,[0.5+(0:nScales)],'FaceColor',ctxtColor2);
xticks(1:nScales);xticklabels(scaleLbl);xlabel('Scale');
%ylabel('Best model (# units)');
linkaxes([ax3,ax4]);
ylim(ax3,[0,25])

% plot of Predictive power vs scale (baseline)
% subplot(3,2,[1,3]);
ax1 = subplot(row,col,[2,col+2]);
plot(pairTable.strfnpp_base'*100,'-o','MarkerSize',4,'Color',[.7,.7,.7]);
hold on; plot(pairTable.ctxtnpp_base'*100,'-o','MarkerSize',4,'Color',ctxtColor1);
plot(mean(pairTable.strfnpp_base*100),'-sk','MarkerFaceColor','k','LineWidth',2,'MarkerSize',9);
hold on; plot(mean(pairTable.ctxtnpp_base*100),'-s','Color',ctxtColor2,'MarkerFaceColor',ctxtColor2,'LineWidth',2,'MarkerSize',9);
title('Baseline');
xticks(1:nScales);xticklabels(scaleLbl);xlim(ax3.XAxis.Limits);%xlabel('Scale');
ylim([0,maxnpp+10]);ylabel('STRF predictive power (%)');

% plot of Predictive power vs scale (salicylate)
% subplot(3,2,[2,4]);
ax2 = subplot(row,col,[3,col+3]);
plot(pairTable.strfnpp_sal'*100,'-o','MarkerSize',4,'Color',[.7,.7,.7]);
hold on; plot(pairTable.ctxtnpp_sal'*100,'-o','MarkerSize',4,'Color',ctxtColor1);
hold on; plot(mean(pairTable.strfnpp_sal*100),'-sk','MarkerFaceColor','k','LineWidth',2,'MarkerSize',9);
hold on; plot(mean(pairTable.ctxtnpp_sal*100),'-s','Color',ctxtColor2,'MarkerFaceColor',ctxtColor2,'LineWidth',2,'MarkerSize',9);
title('Salicylate');
xticks(1:nScales);xticklabels(scaleLbl);xlim(ax3.XAxis.Limits);%xlabel('Scale');
ylim([0,maxnpp+10]);%ylabel('Predictive power (%)');

% TODO: subpanel labels ("A","B","C", etc.)
fig.Renderer = 'painter';
fig.PaperUnits = 'inches';
fig.PaperSize = fig.Position(3:4)./96; %96 dpi
saveas(fig,'.\Figures\Links\DRC_PP_Linearity.pdf')



function stimulus = stimScale(stimulus,scale)

    switch scale
        case 'pow'
            stimulus=stimulus-max(max(stimulus));
            stimulus=10.^(stimulus./10);
        case 'amp'
        % --- default (no sound == -Inf)
            stimulus=stimulus-max(max(stimulus));
            stimulus=10.^(stimulus./20);
        case {'dB', 'dB0'}
        % ---  set -inf -> 0;
            idx=isinf(stimulus) & sign(stimulus) < 0;
            stimulus(idx) = 0;
            stimulus=stimulus./max(max(stimulus));
        case {'bin' , 'binary'}
            idx=isinf(stimulus) & sign(stimulus) < 0;
            stimulus(:) = 1;
            stimulus(idx) = 0;
        case {'dB2', 'dB20'}
        % ---  set -inf -> 0;
            refDB = 20;
            stimulus = stimulus - refDB;
            % ---  everything below refDB -> 0;
            idx = sign(stimulus) < 0;
            stimulus(idx) = 0;
            stimulus=stimulus./max(max(stimulus));
        case 'dB40'
            refDB = 40;
            stimulus = stimulus - refDB;
            % ---  everything below refDB -> 0;
            idx = sign(stimulus) < 0;
            stimulus(idx) = 0;
            stimulus=stimulus./max(max(stimulus));    
        case 'dB50'
            refDB = 50;
            stimulus = stimulus - refDB;
            % ---  everything below refDB -> 0;
            idx = sign(stimulus) < 0;
            stimulus(idx) = 0;
            stimulus=stimulus./max(max(stimulus));    
    end
end