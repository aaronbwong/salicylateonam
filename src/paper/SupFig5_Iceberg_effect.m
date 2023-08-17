%% load data
clear; close all; clc;
datapaths;
load([SumPath,'\AM_results_salicylate.mat']);

MPlotLabels = [0.06,0.93,0.045,0.03; ...
    0.51,0.93,0.045,0.03; ...
    0.06,0.63,0.045,0.03; ...
    0.51,0.63,0.045,0.03; ...
    0.06,0.3,0.045,0.03; ...
    0.51,0.3,0.045,0.03];

% plot VSpp after vs. before for paper (iceberg intensities)

Mf = [16 64 128 512];
Int = [45 60];
Mfcol = mycolors('Mf',4);
col = mycolors('Md',5);
V = VSpp;
nUnits = size(units,1);
mksz = 4;
MMksz = 9;
cnt = 1;
tilecnt = 1;

F = figure; 
F.Units = 'pixels';
F.Position = 96.*[1,1,[21,30]./2.54];
setPDFRes(F)


for f = 1:4

    
    ax(tilecnt) = subplot(3,2,tilecnt); a = ax(tilecnt); hold(a,'on');
    
    for m=1:5
        c = col(m,:);
        for u = 1:nUnits
            v = [V(m,f,Int==45,1,u); V(m,f,Int==60,2,u)]; %select iceberg intensities
            plot(a,v(1),v(2),'o','MarkerFaceColor',c,'MarkerEdgeColor','none','MarkerSize',mksz);
        end
        
    end
    plot(a,[0.01 1],[0.01 1],'k--');
    for m=1:5
        c = col(m,:);
        v1 = squeeze(V(m,f,Int==45,1,:))';
        v2 = squeeze(V(m,f,Int==60,2,:))';
        v = [v1; v2];
        vm = mean(v,2); vs = std(v,[],2); vs = vs/sqrt(nUnits);
        h(m) = errorbar(a,vm(1),vm(2),vs(2),vs(2),vs(1),vs(1),'o','Color',c,'MarkerFaceColor',c,'MarkerEdgeColor','k','MarkerSize',MMksz,'LineWidth',1.5);
    end
    axis(a,'square');
    xlim(a,[0 1]); ylim(a,[0 1]);
    xticks(a,[0 0.5 1]); yticks(a,[0.5 1]);

    if f==4
        leg = legend(h,'0.06','0.125','0.25','0.5','1','Location','southeast','FontSize',12);
        Tleg = title(leg,'Mod. depth','FontName','Arial','FontSize',12);
        leg.Title.NodeChildren.Position = [0.5 0.9 0];

    end

    axprop(a,14);
    cnt = cnt+1; tilecnt = tilecnt +1;
    title([num2str(Mf(f)) ' Hz'],'FontName','Arial','FontSize',16);
end


XL = annotation('textbox',[0.38,0.337,0.297,0.031],'EdgeColor','none','String',"VS_{PP} - baseline - 45dB",'FontSize',14,'FontName','Arial','HorizontalAlignment','center');
YL = annotation('textbox',[0.05,0.533,0.35,0.028],'EdgeColor','none','String',"VS_{PP} - salicylate - 60dB",'FontSize',14,'FontName','Arial','Rotation',90,'HorizontalAlignment','center');

An1 = annotation('textbox',MPlotLabels(1,:),'EdgeColor','none','String',"A",'FontSize',20,'FontName','Arial','FontWeight','bold');
An2 = annotation('textbox',MPlotLabels(2,:),'EdgeColor','none','String',"B",'FontSize',20,'FontName','Arial','FontWeight','bold');
An3 = annotation('textbox',MPlotLabels(3,:),'EdgeColor','none','String',"C",'FontSize',20,'FontName','Arial','FontWeight','bold');
An4 = annotation('textbox',MPlotLabels(4,:),'EdgeColor','none','String',"D",'FontSize',20,'FontName','Arial','FontWeight','bold');

% beeswarm of change in VS thr (only detected --> detected)
load([SumPath,'\AM_results_salicylate.mat']);
nUnits = size(VSthr,4);
UMd = [0.06 0.125 0.25 0.5 1]; NMd = length(UMd);
UMf = [16 64 128 512]; NMf=length(UMf);
UInt = [45 60];

dtsz = 20;
ax(tilecnt) = subplot(3,2,tilecnt); a = ax(tilecnt); hold(a,'on');
a.Position = [0.13,0.07,0.334659090909091,0.225];


cT = [];
dT = [];
m = nan(4,1);
sd = nan(4,1);

for f = 1:NMf
    T1 = squeeze(VSthr(f,Int==45,1,:))';
    T2 = squeeze(VSthr(f,Int==60,2,:))';
    T = [T1; T2];
    tel = ~isnan(T(1,:)) & ~isnan(T(2,:));
    T = T(:,tel);
    
    Td = diff(T,1);
    m(f,1) = mean(Td);
    sd(f,1) = std(Td)/sqrt(length(Td));
    dT = [dT;Td',repmat(f,length(Td),1)];
    cT = [cT; repmat(Mfcol(f,:),length(Td),1)];


end


swarmchart(a,dT(:,2),dT(:,1),dtsz,cT,'o','filled','XJitter','density','XJitterWidth',0.6);
hold(a,'on');
errorbar(a,[1.4 2.4 3.4 4.4],m,sd,'k+','MarkerFaceColor','k');
yline(a,0);
ylim(a,[-25 25]);
xticks(a,[1 2 3 4]);
xticklabels(a,[16 64 128 512]);
title(a,'Temporal');
ylabel('\Delta threshold (dB)');
axprop(a,14); a.YDir = 'reverse';
xlim(a,[0.25 4.9]);
cnt = cnt+1; tilecnt = tilecnt +1;

% beeswarm of change in FR thr (only detected --> detected)

load([SumPath,'\AM_results_salicylate.mat']);
nUnits = size(VSthr,4);
UMd = [0.06 0.125 0.25 0.5 1]; NMd = length(UMd);
UMf = [16 64 128 512]; NMf=length(UMf);
UInt = [45 60];
col = [1 0 0; 0 0.5 0; 0 0 1; 0.75 0 0.75];
ax(tilecnt) = subplot(3,2,tilecnt); a = ax(tilecnt); hold(a,'on');
a.Position=[0.5875,0.07,0.3175,0.225];


m = nan(4,1);
sd = nan(4,1);
cT = [];
dT = [];


for f = 1:NMf
    T1 = squeeze(FRthr(f,Int==45,1,:))';
    T2 = squeeze(FRthr(f,Int==60,2,:))';
    T = [T1; T2];
    tel = ~isnan(T(1,:)) & ~isnan(T(2,:));
    T = T(:,tel);
    Td = diff(T,1);
    m(f,1) = mean(Td);
    sd(f,1) = std(Td)/sqrt(length(Td));
    dT = [dT;Td',repmat(f,length(Td),1)];
    cT = [cT; repmat(Mfcol(f,:),length(Td),1)];
end


swarmchart(a,dT(:,2),dT(:,1),dtsz,cT,'o','filled','XJitter','density','XJitterWidth',0.6);
hold(a,'on');
errorbar(a,[1.4 2.4 3.4 4.4],m,sd,'k+','MarkerFaceColor','k');
xticks(a,[1 2 3 4]);
xticklabels(a,[16 64 128 512]);

yline(a,0);
ylim(a,[-25 25]);
title(a,'Rate');
XLab = xlabel('Modulation frequency (Hz)');
XLab.Position = [-0.733298130280374,31.39322903860981,-1];
% ylabel('Change in threshold (dB)');
axprop(a,14); a.YDir = 'reverse';
xlim(a,[0.25 4.9]);


An5 = annotation('textbox',MPlotLabels(5,:),'EdgeColor','none','String',"E",'FontSize',20,'FontName','Arial','FontWeight','bold');
An6 = annotation('textbox',MPlotLabels(6,:),'EdgeColor','none','String',"F",'FontSize',20,'FontName','Arial','FontWeight','bold');

saveas(F,[FigPath '\FigS5AtoF_population_iceberg_effect.pdf']);

%--Local functions--%
function axprop(ax,fontsz,tickangle)

if nargin < 3
    tickangle = 0;
end
    
    set(ax,'FontName','Arial','FontSize',fontsz,'XTickLabelRotation',tickangle);

end