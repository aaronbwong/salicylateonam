function Results = getRespLatencies(d,minDate,maxDate,maxBin,RespWin,minRespDelay,sessionFilter)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5, RespWin = Inf; end
if nargin < 6, minRespDelay = 0; end
if nargin < 7 || isempty(sessionFilter); selectSessions = 0; 
elseif(isnumeric(sessionFilter)); selectSessions = sessionFilter; sessionFilter = '';
else; selectSessions = -1;
end
    if (selectSessions)
        selectedSessions = sessionSelect(d,minDate,maxDate,sessionFilter);
    end

    mice = {d.Mouse};
    nMice = length(d);
%     UStim = table;

    vars = {'Intensity','ModFreq','ModDepth','Catch'};
    for m = 1:nMice
        sel = d(m).exptdate' >= minDate & d(m).exptdate' <= maxDate;
        if m == 1
            UStim = unique(d(m).StmMtx(sel,vars),'rows');
        else
            UStim = union(UStim,unique(d(m).StmMtx(sel,vars),'rows'));
        end
    end
    
    NUStim = size(UStim,1);
    
% SoundOnsetDelay = sort(LickDelay(Early==1));
for m = 1:nMice
for ii = 1:NUStim
    sel = true(size(d(m).StmMtx,1),1);
    for v = 1:length(vars)
        sel = sel & d(m).StmMtx.(vars{v}) == UStim{ii,vars(v)}; 
    end
        sel = sel & d(m).exptdate(:) >= minDate & d(m).exptdate(:) <= maxDate;
    if (selectSessions); sel = sel & matchstrings(selectedSessions, d(m).session(:)); end
    tDelay = d(m).RespDelay(sel);
    tDelay = sort(tDelay);
    UStim.RespDelay(ii) = {tDelay};
end
    UStim.Mouse = repmat(mice(m),NUStim,1);
    if m == 1
        Results = UStim;
    else
        Results = [Results;UStim];
    end
end

end

function tf = matchstrings(findStrings, inStrings)
    if ischar(findStrings);findStrings = {findStrings};end
    
    tf = false(size(inStrings));
    for ii = 1:length(findStrings)
        tf = tf | strcmp(inStrings,findStrings{ii});
    end
end
