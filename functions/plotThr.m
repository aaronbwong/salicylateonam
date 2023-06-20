function [fig, traces] = plotThr(Thr,uMF,mice,Color,LineStyle)
% PLOTTHR plots the threshold of the mice against modulation frequencies
%
if nargin < 5; LineStyle = []; end
if nargin < 4; Color = []; end

    freqLim = [2, 4000];
    thrLim    = [-33, 0];
    
    fig = gcf;
    ax = gca;
    if (~isempty(Color))
        ax.ColorOrder = Color;
    end
    if (~isempty(LineStyle))
        ax.LineStyleOrder = LineStyle;
    end
    hold(ax,'on');
    traces = plot(ax, uMF,Thr','Marker','o','LineWidth',2);
   
    set(ax, 'Ydir', 'reverse');
    ylim(ax, thrLim);
    set(ax, 'xScale','log');
    xlim(ax, freqLim);
    xticklabels = num2cell(uMF);
    set(ax, 'xtick', uMF, 'xticklabel', xticklabels...
            ,'XTickLabelRotation',45, 'XMinorTick','off');
    [traces.DisplayName] = mice{:};
    legend(ax,'Location','northwest'); 

    xlabel(ax, 'Modulation frequency (Hz)');
    ylabel(ax, 'Detection threshold (dB; 20*log(m))');

    yyaxis(ax,'right');
    Md = [0.03,0.06,0.125,0.25,0.5,1];
    dBMd = 20*log10(Md);
    set(ax,'ytick',dBMd,'yticklabel',Md,'Ydir','reverse',...
        'ycolor','k');
    ylim(ax,thrLim);
    ylabel(ax, 'Detection threshold (m)');
    
end
