clearvars; clc;
close all;

datapaths
load('Data\unitList_all.mat');
ResPath = [SpkPath,'\'];
units = units(units(:,4)==1,:);
mice = unique(units(:,1));
nUnits = size(units,1);
nMice = length(mice);

SRate = nan(nUnits,13);
cnt = 1;
for k=1:nMice
    load([ResPath 'M' num2str(mice(k)) '\M' num2str(mice(k)) '_SpkS.mat'],'cids','SponR');
    mouseUnits = units(units(:,1)==mice(k),:);
    nMouseUnits = size(mouseUnits,1);
    for u=1:nMouseUnits
        clustID = mouseUnits(u,2);
        SR = SponR(:,cids==clustID);
        SRate(cnt,1:length(SR)) = SR';
        
        cnt = cnt+1;
    end
    
end

SRate = SRate(:,[1:4 6:13]);


F = figure('Position',[100,200,513,420]);
a= gca;
hold(a,'on');
X = [-130, -100, -60, -20,-10,0,40,80,120,150,170,180,220]; X = X/60;
X = X(1,[1:4 6:13]);
for k=1:nUnits
    S = SRate(k,:);
    S = sqrt(S);
    plot(X,S,'k-o','MarkerFaceColor','k','MarkerSize',4);
end
plot(X,median(sqrt(SRate),'omitnan'),'r-','LineWidth',3);
Yt = a.YTick;
YtL = Yt.^2;
Yt = Yt(1,[1 2 3 5 7 9]);
YtL = string(YtL(1,[1 2 3 5 7 9])');
YtL = cellstr(YtL);
a.YTick = Yt;
a.YTickLabel = YtL;



axprop(a,16);
xlabel(a,'Time re: salicylate (h)');
ylabel(a,'Spontaneous firing rate (spikes/s)');
xlim(a,[-2.3 3.75]);
% F.PaperOrientation = 'landscape';
F.Renderer = 'painter';
F.PaperUnits = 'inches';
F.PaperSize = F.Position(3:4)./96; %96 dpi
saveas(F,[FigPath,'\Fig1B_population_SponR_timecourse.pdf']);

for k=1:nUnits
    S = SRate(k,:);
    Sm(k,1:2) = [mean(S(2:4),'omitnan'), mean(S(6:end),'omitnan')];
end

[p,h,stats] = ranksum(Sm(:,1),Sm(:,2));


idx = 21;
mouse = units(idx,1);
unit = units(idx,2);

load([SpkPath,'\M' num2str(mouse) '\M' num2str(mouse) '_SpkS.mat'],'cids','cpos','SpkS');


% example unit traces

F = figure('Position',[700,200,513,420]);
s1 = subplot(2,6,1:5);
s2 = subplot(2,6,7:11);
w1 = subplot(2,6,6);
w2 = subplot(2,6,12);


load([SumPath,'\Fig1_ExampleTraces_data.mat']);

hold(s1,'on'); hold(s2,'on');
plot(s1,D1f(100:140000),'k-');
plot(s2,D2f(100:140000),'k-');

Spk1 = SpkS{3,cids==unit}; Spk1 = Spk1(Spk1<140000);
Spk2 = SpkS{8,cids==unit}; Spk2 = Spk2(Spk2<140000);

hold(w1,'on');
for w=1:length(Spk1)
    s = Spk1(w,1);
    t = D1f(1,s-30:s+30);
    plot(w1,t,'k-');
end

hold(w2,'on');
for w=1:length(Spk2)
    s = Spk2(w,1);
    t = D2f(1,s-30:s+30);
    plot(w2,t,'k-');
end

linkaxes([s1 s2 w1 w2],'y');

plot(s1,Spk1,450*ones(length(Spk1),1),'r.');
plot(s2,Spk2,450*ones(length(Spk2),1),'r.');

set(s1,'XColor','none','YColor','none','Box','off','XTick',[],'YTick',[]);
set(s2,'XColor','none','YColor','none','Box','off','XTick',[],'YTick',[]);
set(w1,'XColor','none','YColor','none','Box','off','XTick',[],'YTick',[]);
set(w2,'XColor','none','YColor','none','Box','off','XTick',[],'YTick',[]);

yl = s1.YLim;
xlim(w1,[-30 91]);
xlim(w2,[-30 91]);

%plot scale bar for trace
plot(s2,[32000 32000],[-550 -350],'k-');
plot(s2,[32000 32000+15000],[-550 -550],'k-')

%plot scale bar for waveforms
plot(w2,[1 1],[-550 -350],'k-');
plot(w2,[1 31],[-550 -550],'k-')

F.Renderer = 'painter';
F.PaperUnits = 'inches';
F.PaperSize = F.Position(3:4)./96; %96 dpi
saveas(F,[FigPath,'\Fig1A_example_SponR_traces.pdf']);


%--Local functions--%

function axprop(ax,fontsz,tickangle)

if nargin < 3
    tickangle = 0;
end
    
    set(ax,'FontName','Arial','FontSize',fontsz,'XTickLabelRotation',tickangle);

end