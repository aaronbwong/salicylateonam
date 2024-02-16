function [uMF,uMD,uInt,AUCs,dPrime,mice,HR,NTrl] = calcDPrimeFromLatency(RespLatTable,minRespDelay,maxRespDelay,respWin)
% CALCHRDPRIME calculates dPrime, hit rates for AM sounds
%   Outputs unique values of the three stimulus parameters 
%       uMF (modulation frequency)
%       uMD (modulation depth)
%       uInt (intensity)
%       dPrime: z(HR) - z(FA). False alarm rate (FA) is calculated based on
%       stimulus with zero (0) modulation depth and the same intensity
%       HitRates: fraction of included trials with response latency smaller
%       than RespWin.
%       NTrials: number of trials included in calculation, excluding
%       responses earlier than minRespDelay.
%
    mice = unique(RespLatTable.Mouse,'stable');
    nMice = length(mice);

    uMF     = nonzeros(unique(RespLatTable.ModFreq));
    uMD     =  nonzeros(unique(RespLatTable.ModDepth));
    uInt    =  unique(RespLatTable.Intensity);

    nMF = length(uMF);
    nMD = length(uMD);
    nInt = length(uInt);
    AUCs = nan(nMF,nMD,nInt,nMice);
    dPrime = nan(nMF,nMD,nInt,nMice);
    HR = nan(nMF,nMD+1,nInt,nMice);
    NTrl = nan(nMF,nMD+1,nInt,nMice);

    pd = makedist('Normal','mu',0,'sigma',1);   % creating the standard normal distribution (mu = 0; sigma = 1)
    
    for m = 1:nMice
        Mouse = mice{m};
        for ii = 1:nInt
            Int = uInt(ii);
            catchSel = RespLatTable.Intensity == Int & strcmp(RespLatTable.Mouse,Mouse) ...
                        & RespLatTable.ModFreq == 0 & RespLatTable.ModDepth == 0 ...
                        & RespLatTable.Catch == 1;
            catchDelay = RespLatTable.RespDelay{catchSel};
            catchDelay = catchDelay(~(catchDelay < minRespDelay)); % exclude early response trials (preserves NaNs)
            catchDelay = min(catchDelay,maxRespDelay); % limit value to maximum delay (NaNs will also get maxRespDelay)
           for f = 1:nMF
                MF = uMF(f);
                HR(f,1,ii,m) = mean(catchDelay < respWin);
                NTrl(f,1,ii,m) = length(catchDelay);
                for dd = 1:nMD
                    MD = uMD(dd);
                    testSel = RespLatTable.Intensity == Int & strcmp(RespLatTable.Mouse,Mouse) ...
                        & RespLatTable.ModFreq == MF & RespLatTable.ModDepth == MD ...
                        & RespLatTable.Catch == 0;                    
                    tDelay = RespLatTable.RespDelay{testSel};
                    tDelay = tDelay(~(tDelay < minRespDelay)); % exclude early response trials (preserves NaNs)
                    tDelay = min(tDelay,maxRespDelay); % limit value to maximum delay (NaNs will also get maxRespDelay)
                    
                    HR(f,dd+1,ii,m) = mean(tDelay < respWin);
                    NTrl(f,dd+1,ii,m) = length(tDelay);
                    
                    labels = [zeros(length(catchDelay),1);ones(length(tDelay),1)];
                    scores = [catchDelay(:);tDelay(:)];
                    posclass = 0;
                    [~,~,~,AUCs(f,dd,ii,m)] = perfcurve(labels,scores,posclass);
                end
            end
        end
    end
    
    dPrime =  sqrt(2)*icdf(pd,AUCs);          
end

function tf = matchstrings(findStrings, inStrings)
    if ischar(findStrings);findStrings = {findStrings};end
    
    tf = false(size(inStrings));
    for ii = 1:length(findStrings)
        tf = tf | strcmp(inStrings,findStrings{ii});
    end
end
