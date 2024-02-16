function Thr = thrFromDPrime(dPrime,uMF,uMD,uInt,Int,dThr)
% Inputs
%   dPrime:   4-D matrix with dPrime values. dimensions as follows
%               1. modulation frequency
%               2. modulation depth
%               3. sound intensity
%               4. mouse
%
%   uMF:      unique modulation frequencies (for the 1st dim of dPrime)
%   uMD:      unique modulation depths (for the 2nd dim of dPrime)
%   dThr:     d' value for threshold determination
% 
% Outputs
%   Thr:    Threshold values in 20*log10(m) (dB)
%           does not contain the reference (Mf == 0)

%     nMice = length(mice);
    nMice = size(dPrime,4); 

    uNZMF = uMF(uMF ~= 0);
    NZMFidx = find(uMF~=0);
    nFreq = numel(uMF);
    nNZMF = numel(uNZMF);
    uNZMD = uMD(uMD ~= 0);
    nDepth= numel(uMD);
    nNZMD = numel(uNZMD);
    
    Thr = nan(nNZMF,nMice);
    
    for f = 1:nNZMF
        fidx = NZMFidx(f);
        MF = uMF(fidx);
        for m = 1:nMice
           tempDPrime = dPrime(fidx,uMD~=0,uInt == Int,m);
           tempMD = uNZMD(~isnan(tempDPrime));
           tempDPrime = tempDPrime(~isnan(tempDPrime));
           try
           lowIdx   = find(tempDPrime < dThr, 1,'last');
           if isempty(lowIdx);Thr(f,m) = 20*log10(tempMD(1));continue; end
           highIdx  = find(tempDPrime > dThr & (1:length(tempDPrime)) > lowIdx, 1, 'first');
           if (~isempty(lowIdx) && ~isempty(highIdx))
               Thr(f,m) = interp1(tempDPrime([lowIdx,highIdx]), 20*log10(tempMD([lowIdx,highIdx])),dThr);
           end
           catch ME
               disp(ME)
               keyboard
           end
%            Thr(f,m) = interp1(tempDPrime, 20*log10(tempMD),dThr);
        end
    end
end