function [CrossTime] = calcCrossing(PCs,tt,baseTime,modelTime,thrs)
    nStm    = size(PCs,1);
    bIdx = tt > baseTime(1) & tt < baseTime(2); % time index for baseline
    tIdx = tt > modelTime(1) & tt < modelTime(2);   % time index for model
    nThr = length(thrs);
    CrossSamp = nan(nStm,nThr);
    CrossTime = nan(nStm,nThr);

    PCs = PCs - mean(PCs(:,bIdx),2); % centered at UM noise period
    timeAxis = tt(tIdx);

    for jj = 1:nThr
        thr = thrs(jj);
        if thr == 0; continue; end
        % single direction threshold
        if sign(thr) > 0
            aboveThr = PCs > thr;
        else
            aboveThr = PCs < thr;
        end
        for ii = 1:nStm
            idx = find(aboveThr(ii,tIdx),1);
            if ~isempty(idx); CrossSamp(ii,jj) = idx;end
        end
    end
    CrossTime(~isnan(CrossSamp)) = timeAxis(CrossSamp(~isnan(CrossSamp)));
end