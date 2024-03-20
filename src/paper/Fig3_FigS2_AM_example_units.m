clear; clc; close all;
NoCol = 1;

datapaths;
FigNums = {'Fig3','FigS2'};

for egunit=1:2

if egunit==1
    mouse1 = 26;
    unit1 = 118;
    IntSel = 60;
    MfSel = 16;
elseif egunit==2
    mouse1 = 42;
    unit1 = 13;
    IntSel = 60;
    MfSel = 64;
end


FigPathUnit = [FigPath,filesep,FigNums{egunit}];%['B:\AMSalicylate Paper\example_units\example_unit' num2str(egunit) '\egunit' num2str(egunit) '_'];
load([AMPath,'\M' num2str(mouse1) '\M' num2str(mouse1) '_AMn_data.mat']);
UInt = AM.UInt;
UMf = AM.UMf;
UMd = AM.UMd;
Sets = AM.SetNum;
mksz = 4;
%----AM raster plot figures up
close all;
%unit1 - baseline 60dB
L = figure('Position',[198,280,545,712]);
a = gca;
plotselAMrasterplot(a,AM_SpkT(:,cids==unit1),AM_Stm.Stm,AM,3,IntSel,MfSel,mksz,NoCol);
[Env,EnvX] = GenAMStim(MfSel,1,-0.2,2.4,1,0.02,2);
% Env = (Env*5)+135;
% plot(a,EnvX,Env,'k-');
plotEnv(a,Env,EnvX,141,2.5);
yline(a,132.5);
ylim(a,[0 148]);
% RYT = get(a,'YTick'); RYT(1,end+1)=141;
% RYTL = get(a,'YTickLabel'); RYTL{end+1,1}='Envelope'; set(a,'YTick',RYT,'YTickLabel',RYTL);
ylabel(a,'Modulation depth');
xlabel(a,'Time re: AM transition (s)');
xticks(a,[0 0.5 1 1.5 2]); xticklabels(a,{-1 -0.5 0 0.5 1});
title(a,'Baseline','FontName','Arial','FontSize',20);
L.Renderer = "painters";
saveas(L,[FigPathUnit 'A' '_egunit' num2str(egunit) '_' 'raster_baseline.pdf']);


% unit1 - salicylate 60dB
R = figure('Position',[199,83,545,712]);
a = gca;
plotselAMrasterplot(a,AM_SpkT(:,cids==unit1),AM_Stm.Stm,AM,8,IntSel,MfSel,mksz,NoCol);
[Env,EnvX] = GenAMStim(MfSel,1,-0.2,2.4,1,0.02,2);
% Env = (Env*5)+135;
% plot(a,EnvX,Env,'k-');
plotEnv(a,Env,EnvX,141,2.5);
yline(a,132.5);
ylim(a,[0 148]);
% RYT = get(a,'YTick'); RYT(1,end+1)=141;
% RYTL = get(a,'YTickLabel'); RYTL{end+1,1}='Envelope'; set(a,'YTick',RYT,'YTickLabel',RYTL);
ylabel(a,'Modulation depth');
xlabel(a,'Time re: AM transition (s)');
xticks(a,[0 0.5 1 1.5 2]); xticklabels(a,{-1 -0.5 0 0.5 1});
title(a,'Salicylate','FontName','Arial','FontSize',20);
R.Renderer = "painters";
saveas(R,[FigPathUnit 'B' '_egunit' num2str(egunit) '_'  'raster_salicylate.pdf']);

% Set up AM Cycle Histograms figure
MdIdx = 6;
col = zeros(4,3);
CT = AM.AMCycT;
bin = [62.5*2,30,14,4]*4;
C = figure('Position',[747,81,273,2*1217]);
C.PaperType = 'a3';

cnt = 1;
%for Md=1
MdSel = 1;
s(cnt) = subplot(5,2,cnt); a = s(cnt);
T = cell2mat(CT{UMd==MdSel,UMf==MfSel,UInt==IntSel,Sets==3,cids==unit1});
plotCycH(a,T,bin(find(UMf==MfSel)-1),col(find(UMf==MfSel)-1,:),MfSel,20);
title(a,'Baseline');
axprop(a,16);
yticks(a,[]);
cnt = cnt +1;
ylabel(a,'1');
xticks(a,[0 0.5 1]);
box(a,'off');
set(a.XAxis,'TickDirection','out');
set(a,'TickLength',[0.05 0.05]);

s(cnt) = subplot(5,2,cnt); a = s(cnt);
T = cell2mat(CT{UMd==MdSel,UMf==MfSel,UInt==IntSel,Sets==8,cids==unit1});
plotCycH(a,T,bin(find(UMf==MfSel)-1),col(find(UMf==MfSel)-1,:),MfSel,20); 
title(a,'Salicylate');
axprop(a,16);
yticks(a,[]);
cnt = cnt +1;
xticks(a,[0 0.5 1]);
box(a,'off');
set(a.XAxis,'TickDirection','out');
set(a,'TickLength',[0.05 0.05]);

%for Md=0.5
MdSel = 0.5;
s(cnt) = subplot(5,2,cnt); a = s(cnt);
T = cell2mat(CT{UMd==MdSel,UMf==MfSel,UInt==IntSel,Sets==3,cids==unit1});
plotCycH(a,T,bin(find(UMf==MfSel)-1),col(find(UMf==MfSel)-1,:),MfSel,20);

axprop(a,16);
yticks(a,[]);
cnt = cnt +1;
ylabel(a,'0.5');
xticks(a,[0 0.5 1]);
box(a,'off');
set(a.XAxis,'TickDirection','out');
set(a,'TickLength',[0.05 0.05]);

s(cnt) = subplot(5,2,cnt); a = s(cnt);
T = cell2mat(CT{UMd==MdSel,UMf==MfSel,UInt==IntSel,Sets==8,cids==unit1});
plotCycH(a,T,bin(find(UMf==MfSel)-1),col(find(UMf==MfSel)-1,:),MfSel,20);

axprop(a,16);
yticks(a,[]);
cnt = cnt +1;
xticks(a,[0 0.5 1]);
box(a,'off');
set(a.XAxis,'TickDirection','out');
set(a,'TickLength',[0.05 0.05]);

%for Md=0.25
MdSel = 0.25;
s(cnt) = subplot(5,2,cnt); a = s(cnt);
T = cell2mat(CT{UMd==MdSel,UMf==MfSel,UInt==IntSel,Sets==3,cids==unit1});
plotCycH(a,T,bin(find(UMf==MfSel)-1),col(find(UMf==MfSel)-1,:),MfSel,20);

axprop(a,16);
yticks(a,[]);
cnt = cnt +1;
ylabel(a,'0.25');
xticks(a,[0 0.5 1]);
box(a,'off');
set(a.XAxis,'TickDirection','out');
set(a,'TickLength',[0.05 0.05]);

s(cnt) = subplot(5,2,cnt); a = s(cnt);
T = cell2mat(CT{UMd==MdSel,UMf==MfSel,UInt==IntSel,Sets==8,cids==unit1});
plotCycH(a,T,bin(find(UMf==MfSel)-1),col(find(UMf==MfSel)-1,:),MfSel,20);
axprop(a,16);
yticks(a,[]);
cnt = cnt +1;
xticks(a,[0 0.5 1]);
box(a,'off');
set(a.XAxis,'TickDirection','out');
set(a,'TickLength',[0.05 0.05]);

%for Md=0.125
MdSel = 0.125;
s(cnt) = subplot(5,2,cnt); a = s(cnt);
T = cell2mat(CT{UMd==MdSel,UMf==MfSel,UInt==IntSel,Sets==3,cids==unit1});
plotCycH(a,T,bin(find(UMf==MfSel)-1),col(find(UMf==MfSel)-1,:),MfSel,20);
axprop(a,16);
yticks(a,[]);
cnt = cnt +1;
ylabel(a,'0.125');
xticks(a,[0 0.5 1]);
box(a,'off');
set(a.XAxis,'TickDirection','out');
set(a,'TickLength',[0.05 0.05]);

s(cnt) = subplot(5,2,cnt); a = s(cnt);
T = cell2mat(CT{UMd==MdSel,UMf==MfSel,UInt==IntSel,Sets==8,cids==unit1});
plotCycH(a,T,bin(find(UMf==MfSel)-1),col(find(UMf==MfSel)-1,:),MfSel,20);
axprop(a,16);
yticks(a,[]);
cnt = cnt +1;
xticks(a,[0 0.5 1]);
box(a,'off');
set(a.XAxis,'TickDirection','out');
set(a,'TickLength',[0.05 0.05]);

%for Md=0.06
MdSel = 0.06;
s(cnt) = subplot(5,2,cnt); a = s(cnt);
T = cell2mat(CT{UMd==MdSel,UMf==MfSel,UInt==IntSel,Sets==3,cids==unit1});
plotCycH(a,T,bin(find(UMf==MfSel)-1),col(find(UMf==MfSel)-1,:),MfSel,20);
axprop(a,16);
yticks(a,[]);
xticks(a,[0 0.5 1]);
box(a,'off');
set(a.XAxis,'TickDirection','out');
set(a,'TickLength',[0.05 0.05]);

cnt = cnt +1;
ylabel(a,'0.06');
xlabel(a,'Phase');

s(cnt) = subplot(5,2,cnt); a = s(cnt);
T = cell2mat(CT{UMd==MdSel,UMf==MfSel,UInt==IntSel,Sets==8,cids==unit1});
plotCycH(a,T,bin(find(UMf==MfSel)-1),col(find(UMf==MfSel)-1,:),MfSel,20);
axprop(a,16);
yticks(a,[]);
xticks(a,[0 0.5 1]);
box(a,'off');
set(a.XAxis,'TickDirection','out');
set(a,'TickLength',[0.05 0.05]);

cnt = cnt +1;
xlabel(a,'Phase');

linkaxes(s,'xy');

saveas(C,[FigPathUnit 'C' '_egunit' num2str(egunit) '_'  'CycHist.pdf']);

% AM VS vs. Md
VS = AM.AMmppVS;

T = figure('Position',[1031,1038,407,259]);


a = gca;
hold(a,'on');
MfIdx = find(UMf==MfSel);
c = col(MfIdx-1,:);
bef = VS(2:end,MfIdx,UInt==IntSel,Sets==3,cids==unit1);
aft = VS(2:end,MfIdx,UInt==IntSel,Sets==8,cids==unit1);
plot(a,UMd(2:end),bef,'o--','MarkerFaceColor',c,'MarkerEdgeColor',c,'Color',c);
plot(a,UMd(2:end),aft,'^-','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',c,'Color',c);
axprop(a,16);
set(a,'XScale','log','YLim',[0 1],'YTick',[0 0.5 1],'XTick',UMd(2:end));
ylabel(a,'VS_{PP}');
if egunit ==1
    legend(a,{'baseline','salicylate'},'location','southeast');
else
    legend(a,{'baseline','salicylate'},'location','northwest');
end
    
saveas(T,[FigPathUnit 'D' '_egunit' num2str(egunit) '_'  'VS.pdf']);


% AM d' VS vs. Md
VS = AM.VSdPrime;

DV = figure('Position',[1030,691,407,259]);


a = gca;
hold(a,'on');
MfIdx = find(UMf==MfSel);
c = col(MfIdx-1,:);
bef = VS(2:end,MfIdx,UInt==IntSel,Sets==3,cids==unit1);
aft = VS(2:end,MfIdx,UInt==IntSel,Sets==8,cids==unit1);
plot(a,UMd(2:end),bef,'o--','MarkerFaceColor',c,'MarkerEdgeColor',c,'Color',c);
plot(a,UMd(2:end),aft,'^-','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',c,'Color',c);
axprop(a,16);
set(a,'XScale','log','YLim',[0 4],'YTick',[0 1 2 3 4],'XTick',UMd);
ylabel(a,'d''_{temporal}');
saveas(DV,[FigPathUnit 'E' '_egunit' num2str(egunit) '_'  'dP_temporal.pdf']);



% AM firing rate vs. Md
X = [0.04; 0.25; 1];
Xl = {'UM','0.25','1'};
F = figure('Position',[1031,342,407,259]);

a = gca;
hold(a,'on');
MfIdx = find(UMf==MfSel);
c = col(MfIdx-1,:);
bef = [FR(1,1,UInt==IntSel,Sets==3,cids==unit1,6); FR(2:end,MfIdx,UInt==IntSel,Sets==3,cids==unit1,6)];
aft = [FR(1,1,UInt==IntSel,Sets==8,cids==unit1,6); FR(2:end,MfIdx,UInt==IntSel,Sets==8,cids==unit1,6)];
plot(a,[0.04; UMd(2:end)],bef,'o--','MarkerFaceColor',c,'MarkerEdgeColor',c,'Color',c);
plot(a,[0.04; UMd(2:end)],aft,'^-','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',c,'Color',c);
axprop(a,16);
ylim(a,[0 inf]);
set(a,'XScale','log','XTick',[0.04, UMd(2:end)'],'XTickLabel',{'UM',0.06, 0.125, 0.25, 0.5, 1});
ylabel(a,'Firing rate (spike/s)');
saveas(F,[FigPathUnit 'F' '_egunit' num2str(egunit) '_'  'firing_rate.pdf']);

% AM d' firing rate vs. Md
DF = figure('Position',[1030,43,407,259]);

dPFR = AM.FRdPrime;
a = gca;
hold(a,'on');
MfIdx = find(UMf==MfSel);
c = col(MfIdx-1,:);
bef = dPFR(2:end,MfIdx,UInt==IntSel,Sets==3,cids==unit1);
aft = dPFR(2:end,MfIdx,UInt==IntSel,Sets==8,cids==unit1);
plot(a,UMd(2:end),bef,'o--','MarkerFaceColor',c,'MarkerEdgeColor',c,'Color',c);
plot(a,UMd(2:end),aft,'^-','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',c,'Color',c);
axprop(a,16);
set(a,'XScale','log','XTick',UMd(2:end));
ylabel(a,'d''_{rate}');
saveas(DF,[FigPathUnit 'G' '_egunit' num2str(egunit) '_'  'dP_rate.pdf']);

end


%--Local functions--%
function axprop(ax,fontsz,tickangle)

if nargin < 3
    tickangle = 0;
end
    
    set(ax,'FontName','Arial','FontSize',fontsz,'XTickLabelRotation',tickangle);

end

function plotselAMrasterplot(a,T,Stm,AM,sSet,sInt,sMf,mksz,noCol)
    
    if (noCol)
        Col = zeros(4,3);
    else
        Col         =   [1 0 0; 0 0.5 0; 0 0 1; 0.75 0 0.75]; 
    end
    Mf          =   Stm.Mf;
    Md          =   Stm.Md;
    Int         =   Stm.Intensity;
    Set         =   Stm.Set;
    USet        =   unique(Set);
    
    ppVS        =   AM.AMmppVS;
    
    UMf         =	AM.UMf;
    NMf         =	AM.NMf;
    UMd         =	AM.UMd;
    NMd         =	AM.NMd;
    UInt        =   AM.UInt;
    NInt        =   AM.NInt;
    
    %modified parameters for VS plotting (excluding UM)
    pUMf        =   UMf(2:end);
    pNMf        =   length(pUMf);
    pUMd        =   UMd(2:end);
    pNMd        =   length(pUMd);

    sepIdx = [];
  
    cnt = 1;
    tick = 1;
    YTicks  =   nan((1+(1*(1)*(1))),1);
    YTickLab=   cell((1+(1*(1)*(1))),1);
    
    plotStore   =   cell((1+(1*(1)*(1))),1);
    MfLabIdx    =   zeros((1+(1*(1)*(1))),1);
    sel = Md == 0 & Int == sInt & Set == sSet;
    nTrls = sum(sel);
    t  =   T(sel,1);
    ax = a;
    for m=1:nTrls
        tt = t{m,1};
        if m == 1
            MfLabIdx(tick,1) = 1;
        end
        
        plotStore{tick,1} = plot(ax,tt,(cnt*ones(length(tt),1)),'k.','MarkerSize',mksz);
        hold(ax,'on');
        
        if m== nTrls/2
            YTicks(tick,1) =  cnt;
            YTickLab{tick,1} = 'UM';
            tick = tick+1;
        end
        
        cnt = cnt+1;
    end

    sepIdx = [sepIdx; cnt+0.5]; %#ok<AGROW>
    
    cnt = cnt +2;
    
    for k=1
        Co = Col(UMf(2:end)==sMf,:);
        
        for p=2:NMd
            sel = Md == UMd(p) & Mf == sMf & Int == sInt & Set == sSet;
            nTrls = sum(sel);
            t   =   T(sel,1);
            
            for z=1:nTrls
                if p==2 && z==1
                    MfLabIdx(tick,1) = 1;
                end
                
                tt = t{z,1};
                %                     plotStore{tick,1} = plot(tt,(cnt*ones(length(tt),1)),'.','MarkerEdgeColor', Col, 'MarkerFaceColor', Col);
                plotStore{tick,1} = plot(ax,tt,(cnt*ones(length(tt),1)),'.','MarkerEdgeColor', Co, 'MarkerFaceColor', Co,'MarkerSize',mksz);
                if z== nTrls/2
                    YTicks(tick,1) =  cnt;
                    YTickLab{tick,1} = num2str(UMd(p));
                    tick = tick+1;
                end
                
                cnt = cnt+1;
            end
            
                sepIdx = [sepIdx; cnt+0.5]; %#ok<AGROW> 
            
            %             keyboard
            cnt = cnt+2;
        end
        cnt = cnt+2;
        
        
    end
    
    sepIdx = sepIdx(1:end-1);

    for k=1:length(sepIdx)
        yline(a,sepIdx(k));
    end
    ylim(ax,[0 cnt+3]);
    yticks(ax,YTicks);
    yticklabels(ax,{'UM',num2str(UMd(2:end))});
    xlim(ax,[-0.1 2.4]);
    axprop(ax,16);
    ylim(ax,[0 135]);
end

function plotCycH(F,CycT,bin,col,Mf,rep)

if (~iscell(CycT))
    C = CycT;
else
    C = [];
    for k=1:size(CycT,1)
        C = [C; CycT{k,1}];
    end
end
edges   =   linspace(0,1,bin);
[N,edges] = histcounts(C,edges);
N = (N/rep)/(Mf*0.5); % multiplied by 0.5 (due to sampling from 0.5s of AM)


histogram(F,'BinEdges',edges,'BinCounts',N,'FaceColor',col,'EdgeColor',col);


y = F.YLim;
yticks(F,[0 y(2)]);
xticks(F,[0 0.5 1]);
xlim([0 1]);

end

function [total_mag,x]=GenAMStim(Mf,Md,startT,endT,transStart,transDur,stimDur)
% startT = -0.2; endT = 2.2;
Fs = 1/5.12e-6;         % sampling rate



% stimulus parameter
%     Md = 1;  % ; modulation depth
%     Mf = 16; % Hz; modulation frequency
%     stimDur = 2; %s

    % rise and fall gates of overall stimulus
    RiseTime = 5e-3; %s
    FallTime = 5e-3; %s
    % transition between UM and AM
%     transStart = 1;%s; transition start
%     transDur = 0.02; %s; transition duration

    % (cosine)phase for the AM transition start
    ModPower    = 1+0.5*Md^2;
    phi = 2*pi-acos(sqrt(ModPower)-1); % phi in radian [0,pi] 

TT = startT:1/Fs:endT;
x = TT;
preIdx = x< transStart ;
transIdx = x>=transStart & x <= transStart+transDur;
postIdx = x > transStart+transDur;
onsetIdx = x >= 0 & x <= RiseTime;
offsetIdx = x >= stimDur - FallTime & x <= stimDur;


% AM envelope
AM_mag = nan(size(x)); AM_mag(transIdx|postIdx) = 1+Md*cos(2*pi*Mf*(x(transIdx|postIdx)-transStart)+phi);

% AM transition "gate"
AM_gate = zeros(size(x)); 
AM_gate(transIdx)= sin((x(transIdx)-1)/transDur *pi/2);
AM_gate(postIdx) = 1;

% UM "envelope"
UM_mag = NaN(size(x)); UM_mag(preIdx|transIdx) = sqrt(ModPower);

% UM transition "gate"
UM_gate = zeros(size(x)); 
UM_gate(preIdx) = 1; 
UM_gate(transIdx) = cos((x(transIdx)-1)/transDur *pi/2); 
UM_gate(postIdx) = 0;

% overall cos2 gating
onOff_gate = ones(size(x));
onOff_gate(x<0) = 0; % before stim
onOff_gate(x>stimDur) = 0; % after stim
onOff_gate(onsetIdx) = 0.5-0.5*cos(x(onsetIdx)/ RiseTime * pi);
onOff_gate(offsetIdx) = 0.5-0.5*cos((x(offsetIdx)-stimDur)/ FallTime * pi);

% combining the results
total_mag = nan(size(x));
total_mag(preIdx) = UM_mag(preIdx);
total_mag(postIdx) = AM_mag(postIdx);
total_mag(transIdx) = sqrt( (AM_mag(transIdx).*AM_gate(transIdx)).^2 + (UM_mag(transIdx).*UM_gate(transIdx)).^2);
total_mag = total_mag .* onOff_gate;

end