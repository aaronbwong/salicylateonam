function h = plotEnv(ax,env,tt,yOffset,yScale,Color)
    % h returns handle of the plotting 
    % ax: handle of axis in which to plot
    % env: 1D vector the (positive) magnitude of the envelope
    % tt: 1D vector; time axis for env; must be same dimension
    % yOffset: how much to offset the envelope in y-axis
    % yScale: how much to scale the envelope
    % Color: specify color of the envelope (Matlab color code: RGB triplet,
    %        color name, short name or Hexadecimal)
    
    if nargin < 6 || isempty(Color); Color = [.7,.7,.7]; end
%     if nargin < 6 || isempty(Color); Color = 'k'; end
    if nargin < 5 || isempty(yScale); yScale = 1; end
    if nargin < 4 || isempty(yOffset); yOffset = 0; end

    hold(ax,'on');

    % area option
    h = area(tt(:),[yScale*-env(:)+yOffset,2*yScale*env(:)]);
    h(1).EdgeColor = 'none';%Color;
    h(1).FaceColor = 'none';
    h(2).EdgeColor = 'none';%Color;
    h(2).FaceColor = Color;

    % line option
%     h(1) = plot(tt(:),yScale*-env(:)+yOffset,'-','Color',Color);
%     h(2) = plot(tt(:),yScale*env(:)+yOffset,'-','Color',Color);

    % patch option
%     tt_all = [tt(:);flip(tt(:))];
%     env_all = yScale*[env(:);flip(-env(:))] + yOffset;
%     h = fill(tt_all,env_all,Color);
end