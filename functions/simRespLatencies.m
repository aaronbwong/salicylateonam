
function SimRes = simRespLatencies(CrossTime,SortedStim,uAM,thrs,offset,NDTime,maxRespDelay)

    nStm    = size(CrossTime,1);
    nThr    = size(CrossTime,2);
    nUStm   = size(uAM,1);
    
    minRespDelay = 0;
    respWin = 0.75; % for hit rate
    
    RespDelaySim = min(CrossTime - offset + NDTime,maxRespDelay);
    VarNames = uAM.Properties.VariableNames;
    nVars = length(VarNames);
    for tt = 1:nThr
        thr = thrs(tt);
        RespLatTable_SIM = uAM;
        RespLatTable_SIM.Catch = uAM.Md == 0;
        for ss = 1:nUStm
            sel = true(nStm,1);
            for vv = 1:nVars
                sel = sel & ismember(SortedStim.(VarNames{vv}), uAM.(VarNames{vv})(ss));
            end
%             sel =   (SortedStim.Intensity == uAM.Intensity(ss) &...
%                     SortedStim.Mf == uAM.Mf(ss) & ...
%                     SortedStim.Md == uAM.Md(ss) );

            RespLatTable_SIM.RespDelay(ss) = {RespDelaySim(sel,tt)};
            RespLatTable_SIM.Mouse(ss)      = {num2str(tt)};
        end
        RespLatTable_SIM.Properties.VariableNames{'Mf'} = 'ModFreq';
        RespLatTable_SIM.Properties.VariableNames{'Md'} = 'ModDepth';

        if tt == 1
            SimRes = table(thr,NDTime,{RespLatTable_SIM},...%{dP_SIM},...
                    'VariableNames', {'thr','NDTime','RespLatTable_SIM'});%,'dP_SIM'});
        else
            SimRes(tt,:) = table(thr,NDTime,{RespLatTable_SIM});%,{dP_SIM});
        end
    end
    
   
end