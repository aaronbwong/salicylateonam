function plotRespOverview(PSTHs,SortedStim,binEdges,tempRes_ms,W,WLabel,exampleNeurons,exampleNeuronsLabel)

if nargin < 6; WLabel = '|W|';end
if nargin < 7; exampleNeurons = [];end
if nargin < 8; exampleNeuronsLabel = cellfun(@num2str, num2cell(1:length(exampleNeurons)), 'UniformOutput', false);end

%% Response Overview 2 (include onset, exclude average)

    xx = (binEdges(1:end-1)-1 )*1000 + tempRes_ms/2;
    xx_0 = (binEdges(1:end-1) )*1000 + tempRes_ms/2;
    uMf = nonzeros(unique(SortedStim.Mf));nMf = length(uMf);
    uMd = nonzeros(unique(SortedStim.Md));
    uInt = unique(SortedStim.Intensity);
    nCol = nMf + 2;

    nClu = size(PSTHs,3);
    yy = 1:nClu;   

    [~,cidx] = sort(abs(W),'descend');
    unitsNumber = nan(length(exampleNeurons),1);
    for uu = 1:length(exampleNeurons)
        unitsNumber(uu) = find(cidx == exampleNeurons(uu));
    end
    
    % settings
    Y_TICKS = [1,20:20:nClu,nClu];
    Cmap = flipud(gray);
    XLIMS = [-50,150]; X_TICKS = -400:100:400;
    stimIdx = [];
for Int = [45,60]
    fig = figure('Position',[50,150,1400,800]); clear('ax');
    cmax = 0;
    [ax, ~] = tight_subplot(1,nCol,0.005,[.1 .05],[.05 .07]);
    for col = 1:nCol
%         ax(col)=subplot(1,nCol,col);
        axes(ax(col));
        switch col
            case 2 % onset
                stimIdx = SortedStim.Intensity == Int;
            case 1 % PC weights / other param 
                % do nothing
            otherwise  % Mf
                stimIdx = SortedStim.Mf == uMf(col-2) & SortedStim.Md == 1 & SortedStim.Intensity == Int;
        end
        cmax = max(max(mean(PSTHs(stimIdx,:,:),1),[],'all'),cmax);
        switch col
            case 2 % onset
                imagesc(xx_0,yy,squeeze(mean(PSTHs(stimIdx,:,cidx),1))');
                title('Sound onset');
                xlim(XLIMS);xticks(X_TICKS);xlabel('Time (ms)')
%                 ylabel('Units');yticks(Y_TICKS);
                yticks([]);
%                 ax(col).TickLength = [0.01,0];
                xline(0,'b:','LineWidth',1);
%                 addUnitMarker(ax(col),unitsNumber,exampleNeuronsLabel)
            case 1 % PC weights / other param 
                Colors = [0,0,0;1,0,0];
                scatter((W(cidx,1)),1:nClu,25,'k','filled'); hold on;
                scatter(-(W(cidx,1)),1:nClu,25,'k'); hold on;
%                 scatter(abs(W(cidx)),1:nClu,16,Colors(1+(W(cidx) < 0),:),'filled'); hold on;
                ylabel('Units');
                set(ax(col),'YDir','reverse');ylim([1,nClu])
                set(gca,'FontSize',16);yticks(Y_TICKS);
                xlabel(WLabel);xlim([0,Inf])
            otherwise % Mf
                imagesc(xx,yy,squeeze(mean(PSTHs(stimIdx,:,cidx),1))');
                title(num2str(uMf(col-2),'%d Hz'));
                xlim(XLIMS);xticks(X_TICKS);xlabel('Time (ms)')
                yticks([]);
                xline(0,'b:','LineWidth',1);
%                 addUnitMarker(ax(col),unitsNumber,exampleNeuronsLabel)
        end
        if col == nCol; addUnitMarker(ax(col),unitsNumber,exampleNeuronsLabel); end

    end

cb = colorbar(ax(nCol),'Position',[0.95,0.2,0.007,0.6]); cb.Label.String = '# spikes';
linkaxes(ax,'y');

for ii = 2:nCol
	set(ax(ii),'FontSize',16,'CLim',[0,cmax]);colormap(Cmap);
end
% sgtitle([num2str(Int),'dB']);

end
end


function addUnitMarker(curAxes,unitsNumber,unitsMarker)
	xlimits = xlim(curAxes);
    ylimits = ylim(curAxes);
%     x = [xlimits(2),xlimits(2)+0.05*diff(xlimits)];
    for ii = 1:length(unitsNumber)
%         y = repmat(unitsNumber(ii),1,2);
%         annotation('textarrow',x,y,'String',unitsMarker{ii})
        text(curAxes,xlimits(2),unitsNumber(ii),[char(8592)],...
            'VerticalAlignment','middle','HorizontalAlignment','left', 'FontSize', 8, 'Color',[0.75,0,0]);
        text(curAxes,xlimits(2)+0.08*diff(xlimits),unitsNumber(ii),[unitsMarker{ii}],...
            'VerticalAlignment','middle','HorizontalAlignment','left', 'FontSize', 16, 'Color',[0.75,0,0]);
    end
end