%% Fig6A-D: plot HR/FA rates 
clearvars; clc; close all;
saveFigs = 1;
addpath('.\functions\');
addpath(genpath('.\src\'));
datapaths;
load([BehPath,'\SalicylateBehaviorData_HitRates.mat']);
NMf = length(UMf); NMd = length(UMd); UMd = UMd'; UMf = UMf';

F = figure;
F.Units = "pixels";
F.Position = [1,1,13,11.2]./2.54.*96;
setPDFRes(F);

Col = [0 0 0; 1 0 0; 0 0 1];
catch_val = 0.03;
lnwd = 1; 
mksz = 4;

MPlotLabels = [0.07,0.9,0.045,0.0929; ...
    0.52,0.9,0.045,0.0929; ...
    0.07,0.41,0.045,0.0929; ...
    0.52,0.41,0.045,0.0929];

X = [UMd'/1.1;UMd';UMd'*1.1];
Xcatch = [catch_val/1.1; catch_val; catch_val*1.1];
nMice = 6;

tiledlayout(2,2);

cnt = 1; 
for i=[2,1]
    for f=1:NMf
        nexttile;
        a = gca; hold(a,'on');
        for p=1:3
            C = squeeze(FA(i,p,:)); C = C';
            
            
            Cm = mean(C,2);
            Cs = std(C,[],2)/sqrt(nMice);
            errorbar(a,Xcatch(p,1),Cm,Cs,'o','MarkerFaceColor',Col(p,:),'MarkerEdgeColor',Col(p,:),'Color',Col(p,:),'LineWidth',lnwd,'MarkerSize',mksz);
            V = squeeze(HR(:,f,i,p,:));
            Mm = mean(V,2);
            Ms = std(V,[],2)/sqrt(nMice);
            errorbar(a,X(p,:),Mm,Ms,'o-','LineWidth',lnwd,'MarkerFaceColor',Col(p,:),'MarkerEdgeColor',Col(p,:),'Color',Col(p,:),'MarkerSize',mksz);
        end
        
        a.XScale = 'log';
        a.XMinorTick = 'off';
        a.XTick = [catch_val; UMd([1 3 5])];
        a.XTickLabel = [{"C"}; num2cell(UMd([1 3 5]))];%str2cell(string(UMd([1 3 5])))];
        a.XLim = [0.02 1.1];
        a.YLim = [0 1.1];
        a.YTick = [0 0.5 1];
%         xlabel(a,'Modulation depth');
%         ylabel(a,'Mean response delay (s)');
        axprop(a,12);
        title([num2str(UInt(i)) ' dB SPL - ' num2str(UMf(f)) ' Hz']);

        cnt = cnt +1;
    end
end

XL = xlabel(a,'Modulation depth','Position',[0.007157872893006,-0.207338125015334,-1]);
YL = ylabel(a,'Response probability','Position',[0.000023358475611,1.309712754736702]);

An1 = annotation('textbox',MPlotLabels(1,:),'EdgeColor','none','String',"A",'FontSize',20,'FontName','Arial');
An2 = annotation('textbox',MPlotLabels(2,:),'EdgeColor','none','String',"B",'FontSize',20,'FontName','Arial');
An3 = annotation('textbox',MPlotLabels(3,:),'EdgeColor','none','String',"C",'FontSize',20,'FontName','Arial');
An4 = annotation('textbox',MPlotLabels(4,:),'EdgeColor','none','String',"D",'FontSize',20,'FontName','Arial');

if saveFigs; saveas(F,[FigPath,'\Fig6AtoD_hitrate.pdf']);end

%% Fig6E-H: response latencies
clearvars('-except','saveFigs'); clc; close all;
datapaths;
load([BehPath,'\SalicylateBehaviorData.mat']);

UMd = uMD; NMd = length(UMd); UMf = uMF; NMf = length(UMf); UInt = uInt; NInt = length(UInt);
Col = [0 0 0; 1 0 0; 0 0 1];
catch_val = 0.03;
lnwd = 1; 
mksz = 4;

MPlotLabels = [0.07,0.9,0.045,0.0929; ...
    0.52,0.9,0.045,0.0929; ...
    0.07,0.41,0.045,0.0929; ...
    0.52,0.41,0.045,0.0929];

X = [UMd'/1.1;UMd';UMd'*1.1];
Xcatch = [catch_val/1.1; catch_val; catch_val*1.1];
nMice = 6;

F = figure;
F.Position = [1,1,13,11.2]./2.54.*96;
setPDFRes(F);

tiledlayout(2,2);

cnt = 1;
T = RespLatTable;

maxLat = 4;

for i=[2,1]
    for f=1:NMf
        nexttile(cnt);
        a = gca; hold(a,'on');
        for p=1:3
            V = nan(NMd,nMice);
            C = nan(1,nMice);
            for k=1:nMice
                sel  = T.Intensity == UInt(i) & T.Catch == 1 & strcmp(T.Mouse,mice{k,1}) & strcmp(T.condition,Periods{1,p});
                RL = T.RespDelay(sel); RL = RL{1,1}; RL = RL(~isnan(RL));
                C(1,k) = mean(RL);
                T.meanRespDelay(sel) = C(1,k);
                RL = T.RespDelay(sel); RL = RL{1,1}; RL(isnan(RL)) = maxLat;
                T.medianRespDelay(sel) = median(RL);
                for m=1:NMd
                    sel = T.Intensity == UInt(i) & T.ModFreq == UMf(f) & T.ModDepth == UMd(m) & T.Catch == 0 & strcmp(T.Mouse,mice{k,1}) & strcmp(T.condition,Periods{1,p});
                    RL = T.RespDelay(sel); RL = RL{1,1}; RL = RL(~isnan(RL));
                    V(m,k) = mean(RL);
                    T.meanRespDelay(sel) = V(m,k);
                    RL = T.RespDelay(sel); RL = RL{1,1}; RL(isnan(RL)) = maxLat;
                    T.medianRespDelay(sel) = median(RL);
                end

%                 plot(a,UMd,V(:,k),'-','Color',Col(p,:));

            end
            Cm = mean(C,2);
            Cs = std(C,[],2)/sqrt(nMice);
            errorbar(a,Xcatch(p,1),Cm,Cs,'o','MarkerFaceColor',Col(p,:),'MarkerEdgeColor',Col(p,:),'Color',Col(p,:),'LineWidth',lnwd,'MarkerSize',mksz);

            Mm = mean(V,2);
            Ms = std(V,[],2)/sqrt(nMice);
            errorbar(a,X(p,:),Mm,Ms,'o-','LineWidth',lnwd,'MarkerFaceColor',Col(p,:),'MarkerEdgeColor',Col(p,:),'Color',Col(p,:),'MarkerSize',mksz);
        end
        
        a.XScale = 'log';
        a.XMinorTick = 'off';
        a.XTick = [catch_val; UMd([1 3 5])];
        a.XTickLabel = [{"C"}; num2cell(UMd([1 3 5]))];
        a.XLim = [0.02 1.1];
        a.YLim = [0 1.3];
        a.YTick = [0 0.5 1];
        axprop(a,12);
        title([num2str(UInt(i)) ' dB SPL - ' num2str(UMf(f)) ' Hz']);

        cnt = cnt +1;
    end
end

XL = xlabel(a,'Modulation depth','Position',[0.007157872893006,-0.207338125015334,-1]);
YL = ylabel(a,'Response latency (s)','Position',[0.000023358475611,1.6]);

An1 = annotation('textbox',MPlotLabels(1,:),'EdgeColor','none','String',"E",'FontSize',20,'FontName','Arial');
An2 = annotation('textbox',MPlotLabels(2,:),'EdgeColor','none','String',"F",'FontSize',20,'FontName','Arial');
An3 = annotation('textbox',MPlotLabels(3,:),'EdgeColor','none','String',"G",'FontSize',20,'FontName','Arial');
An4 = annotation('textbox',MPlotLabels(4,:),'EdgeColor','none','String',"H",'FontSize',20,'FontName','Arial');

subject = T.Mouse;
batch = nan(size(subject));
batch(ismember(subject,{'16506-02','16506-03'})) = 1;
batch(ismember(subject,{'17966-03','17966-04'})) = 2;
batch(ismember(subject,{'10996-01','10996-02'})) = 3;
T.frequency = T.ModFreq;
T.modDepth = T.ModDepth;
T.intensity = T.Intensity;
T.subject = subject;
T.batch = batch;
writetable(T,'SalicylateMeanLatencyTable.csv')

if saveFigs;saveas(F,[FigPath,'\Fig6EtoH_responsedelay.pdf']);end
%% Fig6legend
close all;
F = figure('Position',[100,100,238,197]);
a= gca;
hold(a,'on');

plot(a,[1 2],[1 2],'ko-','MarkerFaceColor','k','MarkerSize',mksz+2);
plot(a,[1 2],[1 2],'ro-','MarkerFaceColor','r','MarkerSize',mksz+2);
plot(a,[1 2],[1 2],'bo-','MarkerFaceColor','b','MarkerSize',mksz+2);

L = legend(Periods,'location','best','FontSize',14);
ylim(a,[5 10]);
a.XColor = 'none'; a.YColor = 'none'; a.Color = 'none';
setPDFRes(F);
if saveFigs;saveas(F,[FigPath,'\Fig6legend.pdf']);end

%% Fig6I-L plot dPrime 
clearvars('-except','saveFigs'); clc; close all;
datapaths;
load([BehPath,'\SalicylateBehaviorData.mat']);
UMd = uMD; NMd = length(UMd); UMf = uMF; NMf = length(UMf); UInt = uInt; NInt = length(UInt);

F = figure;
F.Position = [1,1,13,11.2]./2.54.*96;
setPDFRes(F);

Col = [0 0 0; 1 0 0; 0 0 1];
catch_val = 0.03;
lnwd = 1; 
mksz = 4;

MPlotLabels = [0.05,0.9,0.045,0.0929; ...
    0.5,0.9,0.045,0.0929; ...
    0.05,0.41,0.045,0.0929; ...
    0.5,0.41,0.045,0.0929];

X = [UMd'/1.1;UMd';UMd'*1.1];
Xcatch = [catch_val/1.1; catch_val; catch_val*1.1];
nMice = 6;

tiledlayout(2,2);

cnt = 1; 
for i=[2 1]
    for f=1:NMf
        nexttile;
        a = gca; hold(a,'on');
        for p=1:3
            M = squeeze(dP_lat(f,:,i,:,p));
            Mm = mean(M,2);
            Ms = std(M,[],2)/sqrt(nMice);
            errorbar(a,UMd,Mm,Ms,'o-','LineWidth',lnwd,'MarkerFaceColor',Col(p,:),'MarkerEdgeColor',Col(p,:),'Color',Col(p,:),'MarkerSize',mksz);
        end
        a.XScale = 'log';
        a.XTick = UMd([1 3 5]);
        a.XLim = [0.05 1.1];
        a.YLim = [-0.25 3];
        a.YTick = [0 1 2 3 4];
      
        axprop(a,12);
        title([num2str(UInt(i)) ' dB SPL - ' num2str(UMf(f)) ' Hz']);
        cnt = cnt +1;


    end
end

XL = xlabel(a,'Modulation depth','Unit','pixels','Position',[-35,-25,-1]);
YL = ylabel(a,'d''_{response latency}','Position',[0.00027895599202,3.7]);

An1 = annotation('textbox',MPlotLabels(1,:),'EdgeColor','none','String',"I",'FontSize',20,'FontName','Arial');
An2 = annotation('textbox',MPlotLabels(2,:),'EdgeColor','none','String',"J",'FontSize',20,'FontName','Arial');
An3 = annotation('textbox',MPlotLabels(3,:),'EdgeColor','none','String',"K",'FontSize',20,'FontName','Arial');
An4 = annotation('textbox',MPlotLabels(4,:),'EdgeColor','none','String',"L",'FontSize',20,'FontName','Arial');

if saveFigs;saveas(F,[FigPath,'\Fig6ItoL_dPrime.pdf']);end
%% Fig6M-P: Thresholds
clearvars('-except','saveFigs'); clc; close all;
datapaths;
load([BehPath,'\SalicylateBehaviorData.mat']);

UMd = uMD; NMd = length(UMd); UMf = uMF; NMf = length(UMf); UInt = uInt; NInt = length(UInt);
Col = [0 0 0; 1 0 0; 0 0 1];
catch_val = 0.03;
lnwd = 1; 
mksz = 4;

P = ["B","S","W"];

MPlotLabels = [0.07,0.9,0.045,0.0929; ...
    0.52,0.9,0.045,0.0929; ...
    0.07,0.41,0.045,0.0929; ...
    0.52,0.41,0.045,0.0929];

X = [UMd'/1.1;UMd';UMd'*1.1];
Xcatch = [catch_val/1.1; catch_val; catch_val*1.1];
nMice = 6;


F = figure;
F.Position = [1,1,13,11.2]./2.54.*96;
setPDFRes(F);

tiledlayout(2,2);

cnt = 1;



for i=[2 1]
    for f=1:NMf
        nexttile;
        a = gca; hold(a,'on');
        for k=1:nMice

            V = squeeze(Thr_lat(f,i,k,:));


            plot(a,[1 2 3],V,'-o','LineWidth',lnwd);

        end
        M = squeeze(Thr_lat(f,i,:,:));
        Mm = mean(M,1,'omitnan');
        Ms = std(M,[],1,'omitnan');
        errorbar(a,[1 2 3],Mm,Ms,'ko-','LineWidth',lnwd+2,'MarkerSize',mksz+2,'MarkerFaceColor','k');
        a.XTick = 1:3;
        a.XTickLabel = P;
        a.XLim = [0 4];
        a.YDir = 'reverse';
        a.YLim = [-18 0];
        a.YTick = [-25 -20 -15 -10 -5 0];
        axprop(a,12);
        title([num2str(UInt(i)) ' dB SPL - ' num2str(UMf(f)) ' Hz']);

        cnt = cnt +1;
    end
end

XL = xlabel(a,'Testing period','Position',[-0.787713921978789,3,-1]);
YL = ylabel(a,'Threshold_{response latency} (dB)','Position',[-6.769207410755154,-21.968345323741012,0]);

An1 = annotation('textbox',MPlotLabels(1,:),'EdgeColor','none','String',"M",'FontSize',20,'FontName','Arial');
An2 = annotation('textbox',MPlotLabels(2,:),'EdgeColor','none','String',"N",'FontSize',20,'FontName','Arial');
An3 = annotation('textbox',MPlotLabels(3,:),'EdgeColor','none','String',"O",'FontSize',20,'FontName','Arial');
An4 = annotation('textbox',MPlotLabels(4,:),'EdgeColor','none','String',"P",'FontSize',20,'FontName','Arial');

if saveFigs;saveas(F,[FigPath,'\Fig6MtoP_thresholds.pdf']);end


%--Local functions--%
function axprop(ax,fontsz,tickangle)

if nargin < 3
    tickangle = 0;
end
    
    set(ax,'FontName','Arial','FontSize',fontsz,'XTickLabelRotation',tickangle);

end
