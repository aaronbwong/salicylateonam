
datapaths;

SpontFluo = readtable(fullfile(ImgPath,'SpontFluoSalicylate.csv'),'Delimiter',',');
data = [SpontFluo.GR_pre,SpontFluo.GR_sal,SpontFluo.GR_post];
ss = readtable(fullfile(ImgPath,'SpontFluoSalicylate_CI.csv'),'Delimiter',',');
av = ss.Value(ismember(ss.Group,{'pre','sal','post'}));
cis = [ss.CIs_1(ismember(ss.Group,{'pre','sal','post'}),:),...
        ss.CIs_2(ismember(ss.Group,{'pre','sal','post'}),:)];
sd = std([SpontFluo.GR_pre,SpontFluo.GR_sal,SpontFluo.GR_post;])';
er = cis-av; er(:,2) = -er(:,2);
% er = [sd,sd];

avDiff = ss.Value(ismember(ss.Group,{'sal minus pre','post minus pre'}));
cisDiff = [ss.CIs_1(ismember(ss.Group,{'sal minus pre','post minus pre'}),:),...
            ss.CIs_2(ismember(ss.Group,{'sal minus pre','post minus pre'}),:)];
sdDiff = std([(SpontFluo.GR_sal - SpontFluo.GR_pre) , (SpontFluo.GR_post - SpontFluo.GR_pre)])';
erDiff = cisDiff-avDiff; erDiff(:,2) = -erDiff(:,2);
% erDiff = [sdDiff,sdDiff];


%% plot combined results
close all
fig = figure;
fig.Position = [500,100,330,500];
ax1 = subplot(3,1,[1,2]);
% yyaxis('left')
plot((1:3)-0.1,data','-o','Color',[.7,.7,.7]); hold on;
errorbar((1:3)+0.1,av,er(:,1),er(:,2),'Marker','s', 'MarkerSize', 10, 'MarkerFaceColor', 'auto', 'LineStyle','-','Color','k');
% for ii = 1:length(av)
%     plot([ii,10],[av(ii),av(ii)],':k');
% end
% yyaxis('right')
ylabel('F_{GCaMP} / F_{mRuby}')
ylim([0,5]);
xticks(1:3);xticklabels({'Baseline','Salicylate','Washout'});
set(ax1,'XTickLabelRotation',0);

ax2 = subplot(3,1,3);
errorbar((1:3),[0;avDiff],[0;erDiff(:,1)],[0;erDiff(:,2)],'Marker','o', 'MarkerSize', 6, 'MarkerFaceColor', 'auto', 'LineStyle','none','Color','k');
yline(0,':');
ylim([-1.5, 0.5])
ylabel('Delta (F_{GCaMP} / F_{mRuby})')
xticks(1:3);xticklabels({'Baseline','Salicylate','Washout'});
set(ax2,'XTickLabelRotation',0);


set([ax1,ax2],'FontSize', 12);
pause(0.1);
ax1.Position([1,3]) = ax2.Position([1,3]);
linkaxes([ax1,ax2],'x');
xlim([ax1,ax2],[.5,3.5])

setPDFRes(fig);
saveas(fig,fullfile(FigPath,'Fig1D_SpontFluoFig.pdf'))