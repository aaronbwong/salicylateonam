%% get population list
clear; clc; close all;
datapaths;

%% FigS2A-B: FRA example plots (threshold elevation)

% Fig S2A: example unit 1
mouse1 = 42;
unit1 = 196;
G = figure('Position',[100,200,745,393]);

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
setPDFRes(G);
saveas(G,[FigPath,'\FigS2A_exampleunit1.pdf']);

%for a colorbar that doesn't squash one of the panels
GCB = figure('Position',[100,200,745,303]);
s1 = subplot(1,2,2);
s1 = plotFRA(s1,FRA, cids==unit1,FRA_Stm.Stm, 2,1,0);
setPDFRes(GCB);
saveas(GCB,[FigPath,'\FigS2A2_exampleunit1_CB.pdf']);

% Fig S2B: example unit 2
close all;
mouse2 = 30;
unit2 = 104;

F = figure('Position',[100,200,745,393]);
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

setPDFRes(F);
saveas(F,[FigPath,'\FigS2B_exampleunit2.pdf']);

%for a colorbar that doesn't squash one of the panels
FCB = figure('Position',[100,200,745,303]);
s1 = subplot(1,2,2);
s1 = plotFRA(s1,FRA, cids==unit2,FRA_Stm.Stm, 2,1,0);
setPDFRes(FCB);
saveas(FCB,[FigPath,'\FigS2B2_exampleunit2_CB.pdf']);


%% Fig S2C: plot change in population minimum threshold
close all;
mouse1 = 42;
unit1 = 196;
mouse2 = 30;
unit2 = 104;

load([PopPath,'\FRA\population_minThr'],'mThr','units');
salicylate_units = units(units(:,4)==1,:);

sel = [isnan(mThr(:,1)) | isinf(mThr(:,1)), isnan(mThr(:,2)) | isinf(mThr(:,2))];
sel = sum(sel,2); sel= sel==0;
UnitsSubSel = salicylate_units(sel,:);
MT = mThr(sel,:);
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

setPDFRes(F);
saveas(F,[FigPath,'\FigS2C_population_mThr_Salicylate.pdf']);


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


