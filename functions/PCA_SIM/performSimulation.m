DataPathNeuro = '.\Data\Neuro\';
DataPathPCA = '.\Data\PCA_SIM\';
OutPath = '.\Data\PCA_SIM\';
addpath(genpath('.\functions\'))

%% ########### MODEL ###########

KSD_trl = matfile([DataPathNeuro,'AMKSDs_individual_Sal.mat']); % per-trial KSDs
load([DataPathPCA,'PCAResTable_Sal.mat'],'PCARes')

%% Estimate threshold crossing time

% --- simulation parameters ----
% basic parameters
baseTime = [0.8,1]; modelTime = [1, 2];

% specify parameter ranges to simulate 
thrs = 1:1:200; NDTimes = [0.1:0.01:0.2];
bws = [2]; nBws = length(bws);

% select which PCA result(s) to use
pIdx = [];
pp = find(... % use which PCA results 
    strcmp(PCARes.input,'PSTH')...
    & strcmp(PCARes.condition,'baseline')...
    & PCARes.tempRes_ms == 1/256*1000 ...) < 2*eps(1/64) ...
    & PCARes.normalization == 0 ...
    & PCARes.demean == 0 ...
    & PCARes.useAllMd == 0 ...
    );
pIdx = [pIdx,pp];

% -------------------------------

for pcnt = 1:length(pIdx)
    pp = pIdx(pcnt);
disp(['Using the parameter set #', num2str(pp)]);
disp(PCARes(pp,:));

W = PCARes.V{pp}(:,1);
disp('Simulating resp latency & calculating dprime...');
    for bb = 1:nBws 
        bw = bws(bb); bwX = num2str(bw); % use which bandwith to simulate
        disp(['bw: ' bwX ' ms']);
        tic
        [SimRes] = doModel(KSD_trl.(['KSD_trls_',bwX]),KSD_trl.(['tt_',bwX]),...
                W,baseTime,modelTime,thrs,KSD_trl.SortedStim,NDTimes);
        toc
        nSim = size(SimRes,1);
        SimRes = [table(repmat(bw,nSim,1),'VariableNames',{'bw'}),SimRes];
        SimRes = [repmat(removevars(PCARes(pp,:),{'V','Vlatent','VScore','VVarExp'}),nSim,1),SimRes];
        if bb == 1 && pcnt == 1
            SimRes_all = SimRes;
        else
            SimRes_all = [SimRes_all;SimRes];
        end
    end
end
SimRes = SimRes_all; clear SimRes_all;
nSim = size(SimRes,1);
save([OutPath,'SIMResTable_Sal.mat'],'SimRes','-v7.3');disp('data saved.');

%% Calculate dprime_SIM
minRespDelay = 0.075;
maxRespDelay = 1;
respWin = 0.75;
ss = 0;
disp('Calculating dprime ...');
fprintf('%10d/%10d',ss,nSim);

for ss = 1:nSim
    fprintf(repmat('\b',1,21));
    fprintf('%10d/%10d',ss,nSim);
    
    RespLatTable_SIM = SimRes.RespLatTable_SIM{ss};
    sel = strcmp(RespLatTable_SIM.condition,'baseline');
    [uMF,uMD,uInt,AUCs,dP_SIM_Base,mice,HR_SIM] = calcDPrimeFromLatency(RespLatTable_SIM(sel,:),minRespDelay,maxRespDelay,respWin);
    sel = strcmp(RespLatTable_SIM.condition,'salicylate');
    [uMF,uMD,uInt,AUCs,dP_SIM_Sal,mice,HR_SIM] = calcDPrimeFromLatency(RespLatTable_SIM(sel,:),minRespDelay,maxRespDelay,respWin);
    dP_SIM = cat(4,dP_SIM_Base,dP_SIM_Sal);
    dP_SIM = min(dP_SIM,3.9178);
    dP_SIM = max(dP_SIM,-3.9178);
    SimRes.dP_SIM{ss} = dP_SIM;
end
fprintf('\r');
save([OutPath,'SIMResTable_Sal.mat'],'SimRes','-v7.3');disp('data saved.');

%% Estimate threshold from dprime
disp('Estimating threshold from dprime ...');
    uMF = [16,64,128,512]; uMD = [0.06,0.125,0.25,0.5,1];
    uInt = [45,60]; dThr = 1;
fprintf('%10d/%10d\r',ss,nSim);
for ss = 1:nSim
    fprintf(repmat('\b',1,22));
    fprintf('%10d/%10d\r',ss,nSim);
    dP_SIM = SimRes.dP_SIM{ss};
        dP_SIM = min(dP_SIM,3.9178);
        dP_SIM = max(dP_SIM,-3.9178);
    for ii = 1:2
        Int = uInt(ii);
        SimRes.(['Thr_SIM_',num2str(Int)]){ss} = thrFromDPrime(dP_SIM,uMF,uMD,uInt,Int,dThr);
    end
end
fprintf('\r');

save([OutPath,'SIMResTable_Sal.mat'],'SimRes','-v7.3');disp('data saved.');