%% get population list
clear; clc; close all;

% PopPath = 'F:\Maurits\Population\';
datapaths;


%% Fig2AB - FRA example plots (threshold elevation)
close all;
% Fig2A example unit 1
mouse1 = 42;
unit1 = 196;
G = figure('Position',[400,200,745,393]);

s1=subplot(1,3,[1 1.3]);
s2=subplot(1,3,[2.7 3]);
load([FRAPath,'\M' num2str(mouse1) '\M' num2str(mouse1) '_FRA_both_data.mat']);
CF = FRA.FRACF;
mThr = FRA.FRAThr;

s1 = plotFRA(s1,FRA, cids==unit1,FRA_Stm.Stm, 2,0,1);
s2 = plotFRA(s2,FRA, cids==unit1,FRA_Stm.Stm, 9,0,1);
title(s1,'Baseline');
title(s2,'Salicylate');
axis(s1,'square');
axis(s2,'square');

saveas(G,[FigPath,filesep,'Fig2A_exampleunit1.pdf']);

%for a colorbar that doesn't squash one of the panels
GCB = figure('Position',[2767,847,745,303]);
s1 = subplot(1,2,2);
s1 = plotFRA(s1,FRA, cids==unit1,FRA_Stm.Stm, 2,1,0);
saveas(GCB,[FigPath,filesep,'Fig2A_exampleunit1_CB.pdf']);

% Fig2B example unit 2
close all;
mouse2 = 30;
unit2 = 104;

F = figure('Position',[2767,811,745,393]);
s1=subplot(1,3,[1 1.3]);
s2=subplot(1,3,[2.7 3]);


load([FRAPath,'\M' num2str(mouse2) '\M' num2str(mouse2) '_FRA_both_data.mat']);
CF = FRA.FRACF;
mThr = FRA.FRAThr;

s1 = plotFRA(s1,FRA, cids==unit2,FRA_Stm.Stm, 2,0,1);
s2 = plotFRA(s2,FRA, cids==unit2,FRA_Stm.Stm, 9,0,1);
title(s1,'Baseline');
title(s2,'Salicylate');

axis(s1,'square');
axis(s2,'square');

saveas(F,[FigPath,filesep,'Fig2B_exampleunit2.pdf']);

%for a colorbar that doesn't squash one of the panels
FCB = figure('Position',[400,200,745,303]);
s1 = subplot(1,2,2);
s1 = plotFRA(s1,FRA, cids==unit2,FRA_Stm.Stm, 2,1,0);
saveas(FCB,[FigPath,filesep,'Fig2B_exampleunit2_CB.pdf']);


%% Fig2C - plot change in population minimum threshold
close all;
mouse1 = 42;
unit1 = 196;
mouse2 = 30;
unit2 = 104;

load([PopPath,'\population_minThr.mat'],'mThr','units');
salicylate_units = units(units(:,4)==1,:);

sel = [isnan(mThr(:,1)) | isinf(mThr(:,1)), isnan(mThr(:,2)) | isinf(mThr(:,2))];
sel = sum(sel,2); sel= sel==0;
UnitsSubSel = salicylate_units(sel,:);
MT = mThr(sel,:);
rng(20240226);
X = rand(size(MT,1),1)-0.5;
F = figure;
s5 = gca;
hold(s5,'on');
plotIdx = [];
for k=1:size(MT,1)
    if UnitsSubSel(k,1) == mouse1 && UnitsSubSel(k,2) == unit1
        plotIdx = [plotIdx; k];
    elseif UnitsSubSel(k,1) == mouse2 && UnitsSubSel(k,2) == unit2
        plotIdx = [plotIdx; k];
    end
    
    plot(s5,[0.5+X(k), 3+X(k)],MT(k,:),'ko-','MarkerSize',4,'MarkerFaceColor','k');
   
end

plot(s5,[0.5+X(plotIdx(1)), 3+X(plotIdx(1))],MT(plotIdx(1),:),'bo-','MarkerSize',4,'MarkerFaceColor','b','LineWidth',2);
plot(s5,[0.5+X(plotIdx(2)), 3+X(plotIdx(2))],MT(plotIdx(2),:),'o-','Color',...
    [0 0.75 0],'MarkerEdgeColor',[0 0.75 0],'MarkerSize',4,'MarkerFaceColor',[0 0.75 0],'LineWidth',2);

%example units means
meanThr = mean(MT,1);
stdThr = std(MT,[],1);

errorbar(s5,[0.5 3],meanThr,stdThr,'ro-','MarkerSize',5,'MarkerFaceColor','r','LineWidth',3);

axprop(s5,16);
set(s5,'XTick',[0.5 3],'XTickLabel',{'Baseline','Salicylate'},'XLim',[-0.5 4],'YLim',[-5 60]);
ylabel(s5,'Minimum threshold (dB SPL)');

saveas(F,[FigPath,'\Fig2C_population_mThr_Salicylate.pdf']);

%% Fig2DE minimum threshold as a function of CF (Supplementary figure 3)
close all;

load([SumPath,filesep,'unitList_all.mat']);
units = units(units(:,4)==1,:);
mice = unique(units(:,1));
nUnits = size(units,1);
nMice = length(mice);

F = figure('Position',[400,200,560,959]);

CF = nan(nUnits,1);
BF = nan(nUnits,1);
minThr = nan(nUnits,2);
cnt = 1;
for k=1:nMice
    load([FRAPath,filesep, 'M' num2str(mice(k)) '\M' num2str(mice(k)) '_FRA_both_data.mat'],'cids','FRA');
    mCF = FRA.FRACF;
    mBF = FRA.FRABF;
    mThr = FRA.FRAThr;
    mouseUnits = units(units(:,1)==mice(k),:);
    nMouseUnits = size(mouseUnits,1);
    for u=1:nMouseUnits
        clustID = mouseUnits(u,2);
        CF(cnt,1) = mCF(1,cids==clustID);
        BF(cnt,1) = mBF(1,cids==clustID);
        if size(mThr,1) <3
            minThr(cnt,:) = [mThr(1,cids==clustID); NaN];
        else
            minThr(cnt,:) = mThr([1 3],cids==clustID);
        end
        cnt = cnt+1;
    end
    
end


%eliminate units with no CF
sel = isnan(CF) | isinf(CF);
CF = CF(~sel);
minThr = minThr(~sel,:);

%eliminate units with non-valid thresholds
sel = isnan(minThr) | isinf(minThr); sel = sum(sel,2); sel = sel>0;
CF = CF(~sel);
rng(20240226);
CFp = CF.*(1+(rand(length(CF),1)-0.5)/14);

minThr = minThr(~sel,:);
dThr = minThr(:,2)-minThr(:,1);


s(1) = subplot(2,1,1);
s(2) = subplot(2,1,2);

% graph 1 
a = s(1);
hold(a,'on');
CF(5)=CF(42);
UCf = unique(CF);
% UCf = UCf(2:end); % to account for the bug in unique producing two 6.72's
nCf = length(UCf);
for u=1:nCf
    sel = CF==UCf(u);
    nUnits = sum(sel);

    xCF = CFp(sel);
    xThr = minThr(sel,:);

%     if nUnits==1
        R = xCF;
%     else
%         R = xCF.*[1+(rand([length(xCF),1])-0.5)/6)];
%     end


    for k=1:nUnits
        if k==1
        plot(a,[R(k) R(k)],[xThr(k,1) xThr(k,2)],'k-');
        z(1) = plot(a,R(k),xThr(k,1),'kd','MarkerFaceColor','k','MarkerSize',4);
        z(2) = plot(a,R(k),xThr(k,2),'o','MarkerFaceColor',[0.7 0.7 0],'MarkerEdgeColor',[0.7 0.7 0],'MarkerSize',4);
        else
        
        plot(a,[R(k) R(k)],[xThr(k,1) xThr(k,2)],'k-');
        plot(a,R(k),xThr(k,1),'kd','MarkerFaceColor','k','MarkerSize',4);
        plot(a,R(k),xThr(k,2),'o','MarkerFaceColor',[0.7 0.7 0],'MarkerEdgeColor',[0.7 0.7 0],'MarkerSize',4);
        end
    end

end

legend(z,'baseline','salicylate');
xlim(a,[5 70]);
ylim(a,[0 70]);
a.XScale = 'log';
xticks(a,[10 20 30 40 50 60]);
xlabel(a,'CF (kHz)');
ylabel(a,'Minimum threshold (dB SPL)');
axprop(a,14);


%graph 2
a = s(2);
scatter(CFp,dThr,30,'ko','MarkerFaceColor','k','MarkerFaceAlpha',0.4,'MarkerEdgeColor','none');

p = polyfit(log10(CF),dThr',2);
lm = fitlm(log10(CF),dThr');
hold(a,'on');
X = 5:70;
Y = (p(1)*(log10(X).^2))  + (p(2)*(log10(X))) + p(3);
plot(X,Y,'k-');

a=gca;
a.XScale = 'log';
xticks(a,[1 10 20 30 40 50 60]);
xlim(a,[5 70]);
ylim(a,[-15 55])
axprop(a,14);
xlabel('CF (kHz)');
ylabel('Change in minimum threshold (dB)');
saveas(F,[FigPath,filesep,'Fig2DE_minThr_CF.pdf']);


%%--Local functions--%

function a = plotFRA(a,FRA, clustNum,Stm, plotSet,CBplot,FACAPlot)


UFreq       =	FRA.UFreq;
NFreq       =	FRA.NFreq;
UInt       =	FRA.UInt;
NInt       =	FRA.NInt;



M   =   FRA.FRASR; % the thing to plot
if size(M,3)<3
    ZSet = 1;
else
    ZSet = [1 3];
end
zMax = max(squeeze(M(:,:,ZSet,clustNum)),[],'all');
zMin = min([ 0, min(squeeze(M(:,:,ZSet,clustNum)),[],'all')]);


% check duration
% sel = zeros(size(Stm));
% for s = 1:NSets
%     tempSel = reshape([Stm.Set] == plotSet(s),size(sel));
%     sel = tempSel | sel;
% end
% 
% xMin = -min([Stm(sel).PreT]) * 1e-3;
% xMax = max([Stm(sel).StimT]+[Stm(sel).PostT]) * 1e-3;


    setNum = plotSet;
    setIdx = find(FRA.FRASetNum == setNum);
    
    %FRA
    h = a;
     % color map of spike rate
    CData = M(:,:,setIdx,clustNum);
    imagesc(h,CData,'AlphaData',~isnan(CData),[zMin,zMax]);
    if (FACAPlot)
    % contour of FACA p-value
    Cont = -log10(FRA.FACApval(:,:,setIdx,clustNum));
    hold(h,'on');contour(h,Cont,[3,3],'w','ShowText','off');hold(h,'off');
%      % contour of max neighbour correlation
%     Cont = (FRA.MaxNeighCorr(:,:,setIdx,clustNum));
%     hold(h,'on');contour(h,Cont,[0.1,0.2,0.5],'r','ShowText','on');hold(h,'off');
    end
    % format and label graph
    set(h,'Xscale','lin','YDir','normal',...
        'FontName','Arial','FontSize',16, ...
        'XTick',2:4:NFreq,'XTickLabel',round(UFreq(2:4:NFreq),1),...
        'YTick',2:4:NInt,'YTicklabel',UInt(2:4:NInt));
    ylim(h,[0.5 14.5]);
    xlabel(h,'Frequency (kHz)')
    ylabel(h,'Intensity (dB SPL)')
    if (CBplot)
        cb = colorbar(h,'eastoutside');
        cb.Label.String = 'Spike rate (spike/s)';
        cb.Label.Rotation = 270;
        P=cb.Label.Position;
        P(1) = 5;
        cb.Label.Position = P;
    end
end   

function axprop(ax,fontsz,tickangle)

if nargin < 3
    tickangle = 0;
end
    
    set(ax,'FontName','Arial','FontSize',fontsz,'XTickLabelRotation',tickangle);

end


