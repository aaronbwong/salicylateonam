DataPathBeh = '.\Data\Behavior\';
DataPathSIM = '.\Data\PCA_SIM\';
OutPath     = '.\Data\PCA_SIM\';
addpath(genpath('.\functions\'))

%% load results from behavior
behResults = load([DataPathBeh,'salicylateResults.mat'],'mice', 'dP_lat','uMF','uMD','uInt','RespLatTable');

RespLatTable = behResults.RespLatTable;%
RespLatTableUStm = unique(RespLatTable(:,{'Intensity','ModFreq','ModDepth','Catch','condition'}),'rows');
sel = ismember(RespLatTableUStm.ModFreq,[0,16,512]);
RespLatTableUStm = RespLatTableUStm(sel,:);
minRespDelay = 0.075;
maxRespDelay = 1;
tStep = 0.02;
CDF_xx = 0:tStep:maxRespDelay;
nUStm = size(RespLatTableUStm,1);
for ss = 1:nUStm
    sel =   RespLatTableUStm.Intensity(ss) == RespLatTable.Intensity &...
            RespLatTableUStm.ModFreq(ss) == RespLatTable.ModFreq &...
            RespLatTableUStm.ModDepth(ss) == RespLatTable.ModDepth &...
            RespLatTableUStm.Catch(ss) == RespLatTable.Catch &...
            strcmp(RespLatTable.condition,RespLatTableUStm.condition{ss});
    idx = find(sel);
    nAnimals = length(idx);
    for ii = 1:nAnimals
        if ii == 1
            RespDelay = RespLatTable.RespDelay{idx(ii)};
        else
            RespDelay = [RespDelay;RespLatTable.RespDelay{idx(ii)}];
        end
    end
    RespDelay(RespDelay < minRespDelay) = [];
    RespLatTableUStm.RespDelay{ss} = RespDelay;
    CDF = histcounts(RespDelay,CDF_xx,'Normalization','cdf');
    RespLatTableUStm.CDF{ss} = CDF;
end
%% load results from Simulation if not present
if ~exist('SimRes','var')
    load([DataPathSIM,'SIMResTable_Sal.mat'],'SimRes');
end
%% compare results with behavior
nSim = size(SimRes,1);
fprintf('%10d/%10d\r',ss,nSim);
for sim = 1:nSim
    fprintf(repmat('\b',1,22));
    fprintf('%10d/%10d\r',sim,nSim);
    RespLatTable_SIM = SimRes.RespLatTable_SIM{sim};
    sel = ismember(RespLatTable_SIM.ModFreq,[0,16,512]);
    RespLatTable_SIM = RespLatTable_SIM(sel,:);
    nUStm = size(RespLatTable_SIM,1);
    for ss = 1:nUStm
        CDF = histcounts(RespLatTable_SIM.RespDelay{ss},CDF_xx,'Normalization','cdf');
        RespLatTable_SIM.sCDF(ss) = {CDF};
        sel =   RespLatTable_SIM.Intensity(ss)  == RespLatTableUStm.Intensity &...
            RespLatTable_SIM.ModFreq(ss)    == RespLatTableUStm.ModFreq &...
            RespLatTable_SIM.ModDepth(ss)   == RespLatTableUStm.ModDepth &...
            RespLatTable_SIM.Catch(ss)      == RespLatTableUStm.Catch &...
            strcmpi(RespLatTableUStm.condition,RespLatTable_SIM.condition{ss});
        if sum(sel) == 0; keyboard;end
        RespLatTable_SIM.eCDF(ss) =  RespLatTableUStm.CDF(sel);
    end
    selBase = strcmpi(RespLatTable_SIM.condition,'baseline');
    selSal = strcmpi(RespLatTable_SIM.condition,'salicylate');
    sel45 = RespLatTable_SIM.Intensity == 45;
    sel60 = RespLatTable_SIM.Intensity == 60;
    SimRes.CDFDiff(sim)  = rms([RespLatTable_SIM.eCDF{:}] - [RespLatTable_SIM.sCDF{:}]);
    % different conditions
    SimRes.CDFDiff_Base(sim)  = rms([RespLatTable_SIM.eCDF{selBase}] - [RespLatTable_SIM.sCDF{selBase}]);
    SimRes.CDFDiff_Sal(sim)  = rms([RespLatTable_SIM.eCDF{selSal}] - [RespLatTable_SIM.sCDF{selSal}]);
    % different intensities
    SimRes.CDFDiff_45(sim)  = rms([RespLatTable_SIM.eCDF{sel45}] - [RespLatTable_SIM.sCDF{sel45}]);
    SimRes.CDFDiff_60(sim)  = rms([RespLatTable_SIM.eCDF{sel60}] - [RespLatTable_SIM.sCDF{sel60}]);
    % different conditions x intensity combinations
    SimRes.CDFDiff_Base45(sim)  = rms([RespLatTable_SIM.eCDF{selBase&sel45}] - [RespLatTable_SIM.sCDF{selBase&sel45}]);
    SimRes.CDFDiff_Sal45(sim)  = rms([RespLatTable_SIM.eCDF{selSal&sel45}] - [RespLatTable_SIM.sCDF{selSal&sel45}]);
    SimRes.CDFDiff_Base60(sim)  = rms([RespLatTable_SIM.eCDF{selBase&sel60}] - [RespLatTable_SIM.sCDF{selBase&sel60}]);
    SimRes.CDFDiff_Sal60(sim)  = rms([RespLatTable_SIM.eCDF{selSal&sel60}] - [RespLatTable_SIM.sCDF{selSal&sel60}]);
    SimRes.RespLatTable_SIM{sim} = RespLatTable_SIM;
%     if mod(sim,1000)==1; fprintf('%d/%d\r',sim,nSim);end
end
save([OutPath,'SIMResTableEval_Sal.mat'],'SimRes','-v7.3');disp('data saved.');

