function [SimRes] = doModel(KSD,tt,W,baseTime,modelTime,thrs,SortedStim,NDTimes)

    % project to PC
    PCs = projectPC(KSD,W);

    % calculate crossing time
    [CrossTime] = calcCrossing(PCs,tt,baseTime,modelTime,thrs);
    
    % crossing time -> respDelay (add non-decision time)
    respWin = 1; offset = 1;
    
    if any(strcmp('condition',SortedStim.Properties.VariableNames))
        uAM = unique(SortedStim(:,{'condition','Intensity','Mf','Md'}),'rows');
    else
        uAM = unique(SortedStim(:,{'Intensity','Mf','Md'}),'rows');
    end
    for nn = 1:length(NDTimes)
        NDTime = NDTimes(nn);
        if nn ==1
            SimRes = simRespLatencies(CrossTime,SortedStim,uAM,thrs,offset,NDTime,respWin);
        else
            SimRes = [SimRes;simRespLatencies(CrossTime,SortedStim,uAM,thrs,offset,NDTime,respWin)];
        end
    end
    
%             [uMF,uMD,uInt,AUCs,dP_SIM,mice,HR_SIM] = calcDPrimeFromLatency(RespLatTable_SIM,minRespDelay,maxRespDelay,respWin);
%         
%         dP_SIM = min(dP_SIM,3.9178);
%         dP_SIM = max(dP_SIM,-3.9178);

end

