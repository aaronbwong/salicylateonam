DataPathSIM = '.\Data\PCA_SIM\';
DataPathBeh = '.\Data\Behavior\';

%% load results 

% SimRes_subset = SimRes(sIdx,:);
load([DataPathSIM,'SIMResTableEval_Sal_subset.mat'],'SimRes_subset');
results_lat = load([DataPathBeh,'salicylateResults.mat'],'mice', 'Thr_lat','uMF','uMD','uInt');
i = 0;  
uInt = results_lat.uInt;
    nInt = length(uInt);
mice = results_lat.mice;
    nMice = length(mice);
uMF = results_lat.uMF;
uMD = results_lat.uMD;

uMF_Beh = results_lat.uMF;
Thr_lat = results_lat.Thr_lat;
%%

detLvls = [110:10:170];
LineWidth = 1;
MarkerSize = 3;
FontSize = 8;

sColors = 0.9*hsv(length(detLvls)+1);
sColors(3,:) = []; % manual hack for hard to distinguish green colors
bColors = gray(length(detLvls)+1);

fIdx = [1,4];

AMThrRange = [-20,0];
ii = 2;
intStr = '60';

clear ax

for scnt = 1:length(detLvls)%ss-10:10:ss+50
%     sss = sIdx(scnt);
    dP_SIM = SimRes_subset.dP_SIM{scnt};
    
    Thr_SIM = SimRes_subset.(['Thr_SIM_',intStr]){scnt};
    Thr_SIM = Thr_SIM(fIdx,:);
    for ff = 1:2 % frequency
%         fff = fIdx(ff);
%         ax(ff) = subplot(1,2,ff);
%         plot(dP_SIM(fff,:,ii,1),'o--','Color',sColors(scnt,:),'LineWidth',LineWidth,'MarkerSize',MarkerSize); hold on;% 16Hz-allMods-45dB-baseline
%         plot(dP_SIM(fff,:,ii,2),'s-','Color',sColors(scnt,:),'MarkerFaceColor',sColors(scnt,:),'LineWidth',LineWidth,'MarkerSize',MarkerSize) % 16Hz-allMods-45dB-salicylate
%         
%         text(0.05,1-(0.075*(scnt-1)),num2str(SimRes_subset.thr(scnt)),'Color',sColors(scnt,:),...
%             'Units','normalized','VerticalAlignment','top','FontSize',FontSize)
%         
%         set(gca,'XTick',1:5,'XTickLabel',uMD);
%         
%         yline(1,'--','LineWidth',1.5);
%         
%         title([num2str(uMF(ff),'%dHz')]);%, ' ', num2str(uInt(ii),'%ddB')])
%         ylabel('d''')
%         ylim([-.7,4])
%         xlabel('Mod. depth')
%         xlim([0.5,5.5])
        
        ax2(ff) = subplot(1,2,ff);
        
        yyaxis(ax2(ff),'left');
        plot(ax2(ff),Thr_SIM(ff,:),'o-','Color',sColors(scnt,:),'MarkerFaceColor',sColors(scnt,:),'LineWidth',LineWidth,'MarkerSize',MarkerSize) ;hold on;%
        xlim([.5,2.5]);xticks([1,2]);xticklabels({'Baseline','Salicylate'});
        title([num2str(uMF(ff),'%dHz')]);%, ' ', num2str(uInt(ii),'%ddB')])
%         if(ff == 1); 
            ylabel('AM threshold (dB)'); 
%         end
% %         if(ff == 2); yticklabels([]);endNonDete
        ylim(AMThrRange);
        set(ax2(ff),'YDir','reverse')
        yyaxis(ax2(ff),'right');
        set(ax2(ff),'YScale','log','YDir','reverse')
        ylim(ax2(ff),dB2a(AMThrRange));yticks(uMD);
%         if(ff == 1); yticklabels([]);end
%         if(ff == 2); 
            ylabel('AM threshold (m)'); 
%         end
    end
end
for ff = 1:2 % frequency
    yyaxis(ax2(ff),'left');
    plot(ax2(ff),squeeze(mean(Thr_lat(ff,ii,:,1:2),3)),'k:s',...
        'LineWidth',2*LineWidth,...
        'MarkerSize',MarkerSize); %'MarkerFaceColor','k'
end
set(ax2,'FontSize',FontSize);

%% local functions

function a = dB2a(dB)
    a = 10.^(dB./20);
end