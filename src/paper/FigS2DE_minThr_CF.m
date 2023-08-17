%% minimum threshold as a function of CF (Supplementary figure 3)
clear; clc; close all;
datapaths;

load([PopPath,'\unitList_all.mat']);
units = units(units(:,4)==1,:);
mice = unique(units(:,1));
nUnits = size(units,1);
nMice = length(mice);

F = figure('Position',[576,228,560,959]);

CF = nan(nUnits,1);
BF = nan(nUnits,1);
minThr = nan(nUnits,2);
cnt = 1;
for k=1:nMice
    load([FRAPath, '\M' num2str(mice(k)) '\M' num2str(mice(k)) '_FRA_both_data.mat'],'cids','FRA');
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
rng(2023);%set seed for reproducibility 
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

setPDFREs(F);
saveas(F,[FigPath,'\FigS2DE_minThr_CF.pdf']);

%--Local functions--%
function axprop(ax,fontsz,tickangle)

if nargin < 3
    tickangle = 0;
end
    
    set(ax,'FontName','Arial','FontSize',fontsz,'XTickLabelRotation',tickangle);

end