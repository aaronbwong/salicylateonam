DataPath = '.\Data\Neuro\';
OutPath = '.\Data\PCA_SIM\';
addpath(genpath('.\functions\'))

%% === PCA ===

AMKSDs = matfile([DataPath,'AMKSDs_Sal.mat']); % trial-averaged KSDs
AMPSTHs = matfile([DataPath,'AMPSTHs_full_Sal.mat']); % single trial PSTHs

bws = AMKSDs.bws;
nBW = length(bws);
bins = AMPSTHs.bins;
nBins = length(bins);
ii = 0;
conditions = {{'baseline','salicylate'},{'baseline'},{'salicylate'}};
condName = {'both','baseline','salicylate'};
nCond = length(conditions);
for normalization = 0:1
for useAllMd = 0:1
for demean = 0:1
    for input = {'PSTH','KSD'}
        switch input{:}
            case 'PSTH'
                uTempRes = bins;
                nTRes = nBins;
            case 'KSD'
                uTempRes = bws;
                nTRes = nBW;
        end
        
        for cc = 1:nCond
            condition = condName{cc};
            for bb = 1:nTRes
                %% --- innermost loop ---
                ii = ii + 1;
                tempRes_ms = uTempRes(bb)*1000;%7.8125;15.625
                tempRes = num2str(tempRes_ms);
                tempResX = replace(tempRes,'.','x');
                fprintf('%d: input=%s; condition=%s; bandwidth=%sms; demean=%d; norm=%d; useAllMod=%d\r',ii,input{:},condition,tempRes,demean,normalization,useAllMd);
                
                switch input{:}
                    case 'PSTH'
                        X = AMPSTHs.(['PSTH_all_' tempResX]);
                        FWHM = tempRes_ms;  % rectangular
                        ERB = tempRes_ms;  % rectangular
                        Stm = AMPSTHs.SortedStim;
                    case 'KSD'
                        X = AMKSDs.(['KSD_all_' tempResX]);
                        FWHM = (2*sqrt(2*log(2))) * tempRes_ms; % gaussian
                        ERB = (sqrt(2*pi)) * tempRes_ms;  % gaussian
                        Stm = AMKSDs.uAM;
                end
                
                % select particular conditions
%                 if cc == 3; keyboard; end
                sel = ismember(Stm.condition,conditions{cc});
                sel = sel & ismember(Stm.Mf,[0,16,512]);
                X = X(sel,:,:);
                Stm = Stm(sel,:);
                
                [NormFactor,D,V,VScore,Vlatent,VVarExp] =...
                    doDiffPCA(X,Stm,normalization,demean,useAllMd);
                
                
                % collect data into table
                if ii == 1
                    PCARes = table(input,{condition},tempRes_ms,FWHM,ERB,normalization,demean,useAllMd,...
                        {V},{VScore},{Vlatent},{VVarExp},...
                        'VariableNames',{'input','condition','tempRes_ms','FWHM','ERB','normalization','demean',...
                        'useAllMd','V','VScore','Vlatent','VVarExp'});
                else
                    PCARes(ii,:) = table(input,{condition},tempRes_ms,FWHM,ERB,normalization,demean,useAllMd,...
                        {V},{VScore},{Vlatent},{VVarExp});
                end
                % ----------------------
            end
        end
end
end
end
end
PCARes = sortrows(PCARes,{'normalization','useAllMd','demean','condition'});
save([OutPath,'PCAResTable_Sal.mat'],'PCARes');disp('Data saved.')
