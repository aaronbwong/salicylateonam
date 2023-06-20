function plotSingleTrialPCs(SortedStim,thr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
uMf = nonzeros(unique(SortedStim.Mf));
nMf = length(uMf);
uMd = nonzeros(unique(SortedStim.Md));
uInt = unique(SortedStim.Intensity);

% [thrIdx,bwIdx,PCIdx,NDIdx] = ind2sub(ModDimMeth,BestModelIdxMeth(2));
   
for mm = 1:4
pp = BestFitParam.Method(mm);
thrIdx = BestFitParam.thrIdx(mm);

for bb = 1:nBws
    tempRes = num2str(bws(bb)*1000);
    PCs = PCs_all(:,:,:,bb);
    PCs = PCs - mean(PCs(:,:,bIdx,1),3); % centered at UM noise period

    disp(['PC ', num2str(pp)]);
    figure('Position',[10,10,1800,1000]); ax = [];
    ii = 1;
    for ff = 1:nMf
        Mf = uMf(ff);
        for dd = length(uMd):-1:1
            Md = uMd(dd);
            ax = [ax; subplot(length(uMd)+1,4,ff+nMf*(length(uMd)-dd))];
            sel = SortedStim.Mf == Mf &SortedStim.Md == Md & SortedStim.Intensity == uInt(ii);
            plot(timeAxis,squeeze(PCs(pp,sel,tIdx))');
            idx = find(sel);
            yline(thrs(thrIdx));
%             yline(-thrs(thrIdx));
            for jj = 1:sum(sel)
%                 xline(RespDelaySim(idx(jj),thrIdx,1,pp,NDIdx) - NDTime + 1,'Color',[0.5,0.5,0.5]);
                time = CrossTime(idx(jj),thrIdx,1,pp);
                if ~isnan(time); xline(time,'Color',[0.5,0.5,0.5]);end
            end
            if dd == length(uMd); title(num2str(Mf));end
        end
        
    end
    % unmodulated
    ax = [ax; subplot(length(uMd)+1,4,1+nMf*length(uMd))];
    sel = SortedStim.Mf == 0 &SortedStim.Md == 0 & SortedStim.Intensity == uInt(ii);
    plot(timeAxis,squeeze(PCs(pp,sel,tIdx))');
            idx = find(sel);
            yline(thrs(thrIdx));
%             yline(-thrs(thrIdx));
            for jj = 1:sum(sel)
%                 xline(RespDelaySim(idx(jj),thrIdx,1,pp,NDIdx) - NDTime + 1,'Color',[0.5,0.5,0.5]);
                time = CrossTime(idx(jj),thrIdx,1,pp);
                if ~isnan(time); xline(time,'Color',[0.5,0.5,0.5]);end
            end


%     linkaxes(ax);
    sgtitle(['Method ', num2str(pp)]);
% end
end
end
end

