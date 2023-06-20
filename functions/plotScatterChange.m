function [ax1,ax2,sc] = plotScatterChange(x,y,logscale)
%     fig = figure;
    minAll = min([x(:);y(:)]);
    maxAll = max([x(:);y(:)]);
    ax1 = subplot(1,3,[1,2]);
    sc = scatter(x,y,'MarkerEdgeColor',[.5,.5,.5]);hold on;
    plot([minAll,maxAll],[minAll,maxAll],'--k');
    if (logscale); set(ax1, 'YScale','log','XScale','log'); end
    ax2 = subplot(1,3,3);
    plot([1,2],[x(:)';y(:)'],'-o','Color',[.5,.5,.5]);
    xlim([0.5,2.5]);
    if logscale; set(ax2, 'YScale','log'); end
end