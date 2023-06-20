%%

addpath('.\functions\');

salColor = [.7,.7,0];
baseColor = [0,0,0];

DRC_load_data;

%%

allBaseSTRF = cat(3,data.STRF{basePairSel});
allBasePRF = cat(3,data.PRF{basePairSel});
allBaseCGF = cat(3,data.CGF{basePairSel});
allSalSTRF = cat(3,data.STRF{salPairSel});


nClu = size(allBaseSTRF,3);
J = size(allBaseSTRF,1);
K = size(allBaseSTRF,2);
delay = 5*(0:J);
%% Example units
units = [find(pairTable.mouse == 30 & pairTable.cids == 5);...
        find(pairTable.mouse == 30 & pairTable.cids == 262);...
        find(pairTable.mouse == 30 & pairTable.cids == 93);...
        find(pairTable.mouse == 27 & pairTable.cids == 63);...
        ];

fig = figure;
fig.Position = [100,20,600,1000];
FontSize = 16;
for uu = 1:length(units)
    strf = flip(allBaseSTRF(:,:,units(uu)))';
    strf_s = flip(allSalSTRF(:,:,units(uu)))';

    cmax = max(abs([strf(:)]));%max(abs([strf(:);prf(:)]));
    subplot(length(units),2,2*(uu-1)+1)
    imagesc(fliplr(strf),[-cmax,cmax]); colormap(jet);
    colorbar
    setPRFAxis(gca,'reverse')
    set(gca,'FontSize',FontSize)
    
    cmax = max(abs([strf_s(:)]));%max(abs([strf(:);prf(:)]));
    subplot(length(units),2,2*(uu-1)+2)
    imagesc(fliplr(strf_s),[-cmax,cmax]); colormap(jet);
    colorbar
    setPRFAxis(gca,'reverse')
    set(gca,'FontSize',FontSize)
end

%% Extract best frequencies and best delays
BFIdx_base = nan(nClu,1);
DelIdx_base= nan(nClu,1);
BFIdx_sal = nan(nClu,1);
DelIdx_sal= nan(nClu,1);
maxW_base= nan(nClu,1);
maxW_sal= nan(nClu,1);
for u = 1:nClu
%     u = 1;
    strf = squeeze(allBaseSTRF(:,:,u));
    [maxW_base(u),idx] = max(strf,[],'all','linear');
    [DelIdx_base(u),BFIdx_base(u)] = ind2sub([J,K],idx);
    strf = squeeze(allSalSTRF(:,:,u));
    [maxW_sal(u),idx] = max(strf,[],'all','linear');
    [DelIdx_sal(u),BFIdx_sal(u)] = ind2sub([J,K],idx);
end
BFs_base = freq(BFIdx_base)'*1000;
BestDels_base = delay(DelIdx_base)';
BFs_sal = freq(BFIdx_sal)'*1000;
BestDels_sal = delay(DelIdx_sal)';

%% BF - Salicylate vs Baseline
fig = figure;
fig.Position = [710,600,600,420];

[ax1,ax2,sc] = plotScatterChange(BFs_base,BFs_sal,1);
xticks(ax1,freq(1:12:end)*1000);
yticks(ax1,freq(1:12:end)*1000);
xticklabels(ax1,freq(1:12:end));
yticklabels(ax1,freq(1:12:end));
 xlim(ax1,[3e3,64e3]);
xlabel([ax1],'Best frequency - baseline (kHz)')
ylabel([ax1],'Best frequency - salicylate (kHz)')

axis(ax1,'square')
yticks([ax1,ax2],freq(1:12:end)*1000); 
yticklabels([ax1,ax2],freq(1:12:end));
ylim([ax1,ax2],[3e3,64e3]);
ylabel([ax2],'Best frequency (kHz)');
xticks(ax2,1:2);xticklabels({'Base','Sal'});
% set(ax2,'XTickLabelRotation',45);
set([ax1,ax2],'FontSize',FontSize);

%% Descriptive statistics

% best frequency
disp(['Mean change in BF: ',...
    num2str(mean(log2(BFs_sal./BFs_base)),'%.2f'),...
    ' +/- ',...
    num2str(std(log2(BFs_base./BFs_sal)),'%.2f'),...
    ' octaves']);

[~,p,~,~] = ttest(log2(BFs_sal./BFs_base));
disp(['p = ',num2str(p)])

% best delay
disp(['Mean change in best delay: ',...
    num2str(mean(BestDels_sal-BestDels_base),'%.2f'),...
    ' +/- ',...
    num2str(std(BestDels_sal-BestDels_base),'%.2f'),...
    ' ms']);

[~,p,~,~] = ttest(BestDels_sal,BestDels_base);
disp(['p = ',num2str(p)])

%% STRF weight at BF - Salicylate vs Baseline
fig = figure;
fig.Position = [710,150,600,420];

[ax1,ax2,sc] = plotScatterChange(maxW_base,maxW_sal,0);
xlabel([ax1],'max(w_{STRF}) - baseline')
ylabel([ax1],'max(w_{STRF}) - salicylate')

axis(ax1,'square')

ylabel([ax2],'max(w)');
xticks(ax2,1:2);xticklabels({'Base','Sal'});
set([ax1,ax2],'FontSize',FontSize);

disp(['Mean change in w_{STRF}: ',...
    num2str(mean(maxW_sal-maxW_base),'%.2f'),...
    ' +/- ',...
    num2str(std(maxW_sal-maxW_base),'%.2f'),...
    ' ']);

[h,p,ci,stats] = ttest(maxW_sal,maxW_base);
disp(['p = ',num2str(p)])


%% BF centered 
bandwidth = 12; % bins
Tuning_Base = nan(1+2*bandwidth,nClu,1);
Tuning_Sal = nan(1+2*bandwidth,nClu,1);
Tuning_Sal2 = nan(1+2*bandwidth,nClu,1);
Tuning_Base2 = nan(1+2*bandwidth,nClu,1);

for u = 1:nClu
    %%
    BFIdx = BFIdx_base(u);
    if BFIdx > bandwidth && (K - BFIdx) > bandwidth
        STRFRange = BFIdx + (-bandwidth:bandwidth); 
        TunRange = 1:(2*bandwidth+1);
    elseif BFIdx <= bandwidth
        STRFRange = 1:(BFIdx_base(u)+bandwidth);
        TunRange = (bandwidth - BFIdx_base(u) + 2):(2*bandwidth+1); 
    elseif (K - BFIdx) <= bandwidth
        STRFRange = (BFIdx-bandwidth):K;
        TunRange = 1:(bandwidth+1+(K - BFIdx)); 
    end
    
    Tuning_Base(TunRange,u) = allBaseSTRF(DelIdx_base(u),STRFRange,u);
    Tuning_Sal(TunRange,u) = allSalSTRF(DelIdx_sal(u),STRFRange,u);
    
    %% 
    BFIdx = BFIdx_sal(u); %adjusting to the new best frequency
    if BFIdx > bandwidth && (K - BFIdx) > bandwidth
        STRFRange = BFIdx + (-bandwidth:bandwidth); 
        TunRange = 1:(2*bandwidth+1);
    elseif BFIdx <= bandwidth
        STRFRange = 1:(BFIdx_sal(u)+bandwidth);
        TunRange = (bandwidth - BFIdx + 2):(2*bandwidth+1); 
    elseif (K - BFIdx) <= bandwidth
        STRFRange = (BFIdx-bandwidth):K;
        TunRange = 1:(bandwidth+1+(K - BFIdx)); 
    end
    
    Tuning_Sal2(TunRange,u) = allSalSTRF(DelIdx_sal(u),STRFRange,u);
    Tuning_Base2(TunRange,u) = allBaseSTRF(DelIdx_base(u),STRFRange,u);
end

%%
fig = figure;
fig.Position = [710,100,600,200];
% plot(mean(Tuning_Base./max(Tuning_Base),2,'omitnan'),'-sk','MarkerFaceColor','k'); hold on;
errorbar(mean(Tuning_Base./max(Tuning_Base),2,'omitnan'),...
    std(Tuning_Base./max(Tuning_Base)./sqrt(nClu),[],2,'omitnan'),...
    'MarkerFaceColor','k','Color',baseColor,'MarkerSize',6,'Marker','s'); hold on;
errorbar(mean(Tuning_Sal./max(Tuning_Sal),2,'omitnan'),...
    std(Tuning_Sal./max(Tuning_Sal)./sqrt(nClu),[],2,'omitnan'),...
    'MarkerFaceColor',salColor,'Color',salColor,'MarkerSize',6,'Marker','s'); hold on;
ylabel('Norm. w_{STRF}');
xticks(1:bandwidth/2:(2*bandwidth+1));
xticklabels([-bandwidth:bandwidth/2:bandwidth]./12);
xlabel('Freq. re. BF_{base} (oct)')
set(gca,'FontSize',FontSize);
legend({'baseline','salicylate'},'location','best')
%%
fig = figure;
fig.Position = [710,100,600,200];
% plot(mean(Tuning_Base./max(Tuning_Base),2,'omitnan'),'-sk','MarkerFaceColor','k'); hold on;
errorbar(mean(Tuning_Base./max(Tuning_Base),2,'omitnan'),...
    std(Tuning_Base./max(Tuning_Base)./sqrt(nClu),[],2,'omitnan'),...
    'MarkerFaceColor','k','Color',baseColor,'MarkerSize',6,'Marker','s'); hold on;
errorbar(mean(Tuning_Sal2./max(Tuning_Sal2),2,'omitnan'),...
    std(Tuning_Sal2./max(Tuning_Sal2)./sqrt(nClu),[],2,'omitnan'),...
    'MarkerFaceColor',salColor,'Color',salColor,'MarkerSize',6,'Marker','s'); hold on;
ylabel('Norm. w_{STRF}');
xticks(1:bandwidth/2:(2*bandwidth+1));
xticklabels([-bandwidth:bandwidth/2:bandwidth]./12);
xlabel('Freq. re. BF (oct)')
set(gca,'FontSize',FontSize);
legend({'baseline','salicylate'},'location','best')
% %%
% fig = figure;
% fig.Position = [710,100,600,200];
% % plot(mean(Tuning_Base./max(Tuning_Base),2,'omitnan'),'-sk','MarkerFaceColor','k'); hold on;
% errorbar(mean(Tuning_Base2./max(Tuning_Base2),2,'omitnan'),...
%     std(Tuning_Base2./max(Tuning_Base2)./sqrt(nClu),[],2,'omitnan'),...
%     'MarkerFaceColor','k','Color',baseColor,'MarkerSize',6,'Marker','s'); hold on;
% errorbar(mean(Tuning_Sal2./max(Tuning_Sal2),2,'omitnan'),...
%     std(Tuning_Sal2./max(Tuning_Sal2)./sqrt(nClu),[],2,'omitnan'),...
%     'MarkerFaceColor',salColor,'Color',salColor,'MarkerSize',6,'Marker','s'); hold on;
% ylabel('Norm. w_{STRF}');
% xticks(1:bandwidth/2:(2*bandwidth+1));
% xticklabels([-bandwidth:bandwidth/2:bandwidth]./12);
% xlabel('Freq. re. BF_{sal} (oct)')
% set(gca,'FontSize',FontSize);
% legend({'baseline','salicylate'},'location','best')

%% average Pearson correlation between STRF & PRF

dd = data.rho(basePairSel | salPairSel);
disp(['Mean: ',num2str(mean(dd))]);
disp(['SD: ',num2str(std(dd))]);
disp(['Range: ',num2str(min(dd)) ,' - ',num2str(max(dd))]);


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
    rectangle('Position',[0.5,N+0.5,1,1],'EdgeColor','none','FaceColor',[.5 .5 .5])
end