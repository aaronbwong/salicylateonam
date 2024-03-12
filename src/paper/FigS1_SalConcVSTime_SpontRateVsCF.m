% Figure S1 Van den Berg*, Wong* et al, 2024, iScience

% FigS1A - Salicylate Concentration vs Time

close all; clear; clc;
datapaths
load([ConcPath,filesep,'salicylate_concentration_data.mat']);
% F = figure; 
% a=gca;
% plot(a,T,C,'ko','MarkerFaceColor','k');

p = polyfit(T,C,1);
% [fitresult, gof] = expfit(T, C);
a=gca;
% axprop(a,16);
rct = area(a,[1.333,2],[500,500],'FaceColor',[1,.7,.7],'EdgeColor','none');
hold(a,'on');
plot(a,T,C,'ko','MarkerFaceColor','k');
X = [0 5];
Y = p(1)*X+p(2);
plot(a,X,Y,'b-','LineWidth',1);
ylim(a,[0 500]);
xticks(a,0:5);
xlim(a,[0 5.5]);
ylabel('Salicylate concentration (mg/L)');
xlabel('Time after salicylate injection (hr)');
% axprop(a,16);
a.FontSize = 16;
F=gcf;

F.Renderer = 'painter';
F.PaperUnits = 'inches';
F.PaperSize = F.Position(3:4)./96; %96 dpi
saveas(F,'.\Figures\Links\FigS1A_salicylate_curve.pdf');

%% FigS1B - Change in SpR against CF (for paper)

clearvars; clc; close all;
datapaths
load([PopPath,filesep,'Population_spontaneous_rate.mat']);
dSFR = diff(SpR);


load([SumPath,filesep,'unitList_all.mat']);
units = units(units(:,4)==1,:);
mice = unique(units(:,1));
nUnits = size(units,1);
nMice = length(mice);

CF = nan(nUnits,1);
BF = nan(nUnits,1);
F = figure;
cnt = 1;
for k=1:nMice
    load([FRAPath filesep 'M' num2str(mice(k)) '\M' num2str(mice(k)) '_FRA_both_data.mat'],'cids','FRA');
    mCF = FRA.FRACF;
    mBF = FRA.FRABF;
    mouseUnits = units(units(:,1)==mice(k),:);
    nMouseUnits = size(mouseUnits,1);
    for u=1:nMouseUnits
        clustID = mouseUnits(u,2);
        CF(cnt,1) = mCF(1,cids==clustID);
        BF(cnt,1) = mBF(1,cids==clustID);
        
        
        cnt = cnt+1;
    end
    
end

sel = isnan(CF) | isinf(CF) | CF==0;
CF = CF(~sel);
SpR = SpR(:,~sel);

S = SpR(2,:) ./ (SpR(1,:)+SpR(2,:));

sel = isnan(S);
CF = CF(~sel);
S = S(~sel);

Sran = (rand([1,length(S)])-0.5)/25;
Sp = S; %+Sran;
CFp = CF.*(1+rand([length(CF),1])/6);
scatter(CFp,Sp,30,'ko','MarkerFaceColor','k','MarkerFaceAlpha',0.4,'MarkerEdgeColor','none');

p = polyfit(log10(CF),S',1);
lm = fitlm(log10(CF),S');
hold on;
X = [5 70];
Y = log10(X)*p(1)+p(2);
plot(X,Y,'k-','LineWidth',1);
a=gca;
a.XScale = 'log';
xticks(a,[1 10 20 30 40 50]);
xlim(a,[5 65]);
% axprop(a,16);
a.FontSize = 16;
xlabel('CF (kHz)');
ylabel('SFR_{salic} / (SFR_{base} + SFR_{salic})');

F.Renderer = 'painter';
F.PaperUnits = 'inches';
F.PaperSize = F.Position(3:4)./96; %96 dpi
saveas(F,[FigPath,filesep,'FigS1B_SponFR_CF.pdf']);