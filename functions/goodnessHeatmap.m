function [M, uVar1, uVar2] = goodnessHeatmap(SimRes,IndepVar,DepVar,idx)
% SimRes table of simulation resultst
% MatchTable 
% IndepVar cell array of 2 independet variable names
% DepVar string of depedent Variable name

%     nSim = size(SimRes,1);
    
    Var1 = IndepVar{1};
    Var2 = IndepVar{2};
    
    uVar1 = unique(SimRes.(Var1));nVar1 = length(uVar1);
    uVar2 = unique(SimRes.(Var2));nVar2 = length(uVar2);
    
    M = nan(nVar1,nVar2);
    
    useCell = iscell(SimRes.(DepVar));
    
    for i1 = 1:nVar1
        for i2 = 1:nVar2
            sel = SimRes.(Var1) == uVar1(i1) & ...
                SimRes.(Var2) == uVar2(i2);
            if(useCell)
                data = SimRes.(DepVar){sel};
                M(i1,i2) = data(idx);
            else
                M(i1,i2) = SimRes.(DepVar)(sel);
            end
        end
    end
    
end