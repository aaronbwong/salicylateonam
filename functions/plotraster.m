function [f,YTick,YTickLab] = plotraster(f, EventT,Var, Colour,yinc,ticklevel)
%RASTER creates a raster plot of aligned time-events (spikes, licks etc).

%       Inputs:
%       f           = handle of axes for plotting raster plot
%       EventT      =   1-D cell array of NStim/NTrl x 1, each cell containing
%       event time (in seconds) of one trial/stim
%       Var         =   nTrls x nVar matrix of stimuli variables (1st Var
%       is outermost grouping variable (and grouped by colour)
%       Colour      =   matrix (n x 3) of RGB codes for colouring of the
%       outermost grouping variable
%       Outputs:
%       f           = handle of axes for plotted raster plot
%       YTick       = cell array with NVar + 1 cells
%                     each containing vector with tick mark indices at
%                     different levels (1, 2, 3, etc. variables)
%  written on 31/01/2020 (Maurits van den Berg & Aaron Wong) 

NVar        =   size(Var,2);
NStim       =   size(Var,1);
UVar        =   unique(Var,'rows');
NUVar       =   size(UVar,1);

CVar        =   unique(Var(:,1));
NColor      =   length(CVar);

mksz        =   6;

if nargin < 5; yinc =  flip(10.^[1:NVar]);end
if isempty(yinc); yinc = flip(10.^[1:NVar]); end

%% Creating raster data
[SortVar,idx]     = sortrows(Var);
diffVar = SortVar(1:end-1,:)~=SortVar(2:end,:);
diffVar = [zeros(1,NVar);diffVar];

XX = [];YY = [];CC = [];

Ycnt	=	0;

Ypos    = nan(NStim,1);


for k=1:NStim

    stimNum =   idx(k);
    
    % get spike times and append to XX
    X		=	EventT{stimNum,1};
    XX      =   [XX; X(:)];

    % calculate next position
    Ycnt = Ycnt + 1;
    for v = 1:NVar  % add offset if stimulus is new
        if (diffVar(k,v)); Ycnt = Ycnt + yinc(v); end
    end
    
    % append y points to YY
    Y		=	ones(size(X)) * Ycnt;
    YY      =   [YY; Y(:)];
    Ypos(stimNum) = Ycnt;

    % index for color
    C       =   ones(size(X)) * find(CVar == Var(stimNum,1));
    CC      =   [CC; C(:)];
    
end

%% Plotting 
if (size(Colour,1) > 1) % plotting each color
    hold(f,'off');
    for k=1:NColor
        pColour     =  Colour(k,:);
        sel = (CC == k);
        plot(f,XX(sel),YY(sel),'.','MarkerFaceColor',pColour,'MarkerEdgeColor',pColour);hold(f,'on');
    end
else % plotting single color
    if isempty(Colour)
        pColour     =  'k';
    else
        pColour     =  Colour;
    end
    plot(f,XX,YY,'.','MarkerFaceColor',pColour,'MarkerEdgeColor',pColour,'MarkerSize',mksz)
end

%% Calculate YTicks at each level

YTick       =   cell(NVar+1,1);

for v = 1:NVar
    tempUVar        =   unique(Var(:,1:v),'rows');
    tempNUVar       =   size(tempUVar,1);
    tempYTick       =   nan(tempNUVar,1);
    for k=1:tempNUVar
        sel = ones(NStim,1);
        for p = 1:v
            sel         =   sel & Var(:,p) == tempUVar(k,p);% & Var2 == UVar(k,2);
        end
        tempYTick(k) = mean(Ypos(sel));
    end
    YTick{v} = tempYTick;
    
end
YTickLab = tempUVar(:,ticklevel);
YTick{end} = Ypos;
end

