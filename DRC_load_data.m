%% Read data from DRC analysis results
% and compile into a table "data"
% columns include: sp,np,err,PRF,CGF,wtf_reg_params, wtauphi_reg_params,
%                   strfpp, ctxtpp

function [data,pairTable] = DRC_load_data()
%
% STRF parameters
K = 61; J = 8;
freq = logspace(log10(2),log10(64),K);
delay = 5.*(1:J);
% CGF parameters
N = 12; M = 10;
phi = (-N:N)./12;
tau = 5.*(0:M);

% context hyper parameters
context_reg_param = [6,0.5,0.5]';
optstyle = 'ASD';
scale = 'amp';
loadpath = ['.\Data\Neuro\DRC\',optstyle,'_',...
    replace(mat2str(context_reg_param),{';',',','.','[',']'},{'_','_','x','',''}),...
    '\'];

Mice = [25,26,27,28,29,30,32,38,42,43,49,50,52,57,71,72,73,74];
% Mice = [26,27,29,30,42,71,72];
% Mice = [30];
% read metadata for units
salUnits = readtable('.\Data\salicylate_units.xlsx');

% initiate progress bar
multiWaitbar('CLOSEALL');
    sp = []; np = []; err = [];
    init = 1;

% loop that load data in for each mouse
for m = 1:length(Mice)
    Mouse = num2str(Mice(m),'M%02.0f');
    multiWaitbar('Mouse', m/length(Mice));
    load(['.\Data\Neuro\Spks\', Mouse, '\', Mouse '_SpkS.mat'],'cids','cpos');
    cids = cids';cpos = cpos';
    load(['.\Data\Neuro\DRC\', Mouse, '_DRC_PSTH.mat'],'DRCSets');
    
    setNums = DRCSets;
    nSets = length(DRCSets);

    for s = 1:nSets
        setNum = setNums(s);   
%         load(['Results_',Mouse,'_S',num2str(setNum,'%02d'),'_J',num2str(J),'_M',num2str(M),'_',optstyle],'results');   
        load([loadpath,'Results_',Mouse,'_S',num2str(setNum,'%02d'),'_J',num2str(J),'_M',num2str(M),'_',optstyle,'_',scale],'results');   

        nClu = length(results);
        
        salDataSet = nan(nClu,1);
        sp = nan(nClu,1);
        np = nan(nClu,1);
        err = nan(nClu,1);
        c_strf = nan(nClu,1);
        STRF = cell(nClu,1);
        c_ctxt = nan(nClu,1);
        PRF = cell(nClu,1);
        CGF = cell(nClu,1);
        wtf_reg_params = cell(nClu,1);
        wtauphi_reg_params = cell(nClu,1);
        strftpp = nan(nClu,1);
        ctxttpp = nan(nClu,1);
        strfpp = nan(nClu,1);
        ctxtpp = nan(nClu,1);
        for ii = 1:nClu
            salDataSet(ii) = any(salUnits.mouse == Mice(m) & salUnits.cids == cids(ii));              
            sp(ii) = [results(ii).stats.signalpower];
            np(ii) = [results(ii).stats.noisepower];
            err(ii) = [results(ii).stats.error];
            c_strf(ii) = results(ii).strf.ww(1);
            STRF{ii} = flipud(reshape(results(ii).strf.ww(2:end),J,K));
            c_ctxt(ii) = results(ii).full_rank_sparse_rep.c;
            PRF{ii} = results(ii).full_rank_sparse_rep.wtf;
            CGF{ii} = results(ii).full_rank_sparse_rep.wtauphi;
            wtf_reg_params{ii} = results(ii).full_rank_sparse_rep.wtf_reg_params;
            wtauphi_reg_params{ii} = results(ii).full_rank_sparse_rep.wtauphi_reg_params;
            strftpp(ii) = results(ii).strf.tpp;
            ctxttpp(ii) = results(ii).full_rank_sparse_rep.tpp;
            strfpp(ii) = results(ii).strf.pp;
            ctxtpp(ii) = results(ii).full_rank_sparse_rep.pp;
        end
        mouse = repmat(Mice(m),nClu,1);
        set = repmat(setNum,nClu,1);
        data_temp = table(mouse,cids,cpos,set,salDataSet,sp,np,err,c_strf,STRF,c_ctxt,PRF,CGF,wtf_reg_params,wtauphi_reg_params,strftpp,ctxttpp,strfpp,ctxtpp);
        if init == 1
            data = data_temp;init = 0;
        else
            data = [data;data_temp];
        end
    end
    
end
clear('results','ans','cids','cpos','salDataSet','sp','np','err','STRF','PRF','CGF','wtf_reg_params','wtauphi_reg_params','strfpp','ctxtpp');
clear('strftpp','ctxttpp','c_ctxt','c_strf');
clear('set', 'data_temp','mouse','Mouse','nClu','nSets','setNum','setNums','DRCSets');
clear('m','ii', 's','init');
multiWaitbar('CLOSEALL');

% some calculations
data.nnp = data.np./data.sp;
data.strfntpp = data.strftpp./data.sp;
data.ctxtntpp = data.ctxttpp./data.sp;
data.strfnpp = data.strfpp./data.sp;
data.ctxtnpp = data.ctxtpp./data.sp;

nClu = size(data,1);
for ii = 1:nClu
    data.rho(ii) = corr(reshape(data.STRF{ii},[],1),reshape(data.PRF{ii},[],1));
end
% %% define parameters & units selection


% units*condition selections
baseOnly = data.set == 4 & (data.salDataSet == 0 ...
    | data.mouse == 25 | data.mouse == 74); %baseline of units without salicylate
baseSel =  data.set == 4 & data.salDataSet == 1 ...
    & data.mouse ~= 25 & data.mouse ~= 74 ; %baseline of units with salicylate
salSel = data.set == 10& data.salDataSet == 1; %salicylate 
furSel = data.set == 10& data.salDataSet == 0; %furosemide 

data.baseOnly = baseOnly;

% enough signal power
nSigma = 2;
incSel = data.sp > nSigma* data.err;
data.sigSignal = incSel;

% select relevant units
baseTable = data(baseSel & incSel,{'mouse','cids'}); % units with acceptable baseline recording
salTable = data(salSel & incSel,{'mouse','cids'}); % units with acceptable salicylate recording
pairTable = intersect(baseTable,salTable); % units with acceptable baseline and salicylate recordings
pairSel = ismember(data(:,{'mouse','cids'}),pairTable); % index of recordings from these units
basePairSel = pairSel & baseSel; % baseline recordings to include
salPairSel = pairSel & salSel;  % salicylate recordings to include
baseTable2 = data(basePairSel,{'mouse','cids'}); % units with acceptable baseline recording
salTable2 = data(salPairSel,{'mouse','cids'}); % units with acceptable salicylate recording
if (all(baseTable2.mouse == salTable2.mouse) && all(baseTable2.cids == salTable2.cids));disp('Order OK. Proceeding with comparison');end

% CGF improvements
pairTable.strfnpp_base = data.strfnpp(basePairSel);
pairTable.ctxtnpp_base = data.ctxtnpp(basePairSel);
pairTable.baseImpr = pairTable.ctxtnpp_base - pairTable.strfnpp_base;
pairTable.strfnpp_sal = data.strfnpp(salPairSel);
pairTable.ctxtnpp_sal = data.ctxtnpp(salPairSel);
pairTable.salImpr = pairTable.ctxtnpp_sal - pairTable.strfnpp_sal;

pairTable.c_base = data.c_ctxt(basePairSel);
pairTable.c_sal = data.c_ctxt(salPairSel);

pairTable.basePairSel = find(basePairSel);
pairTable.salPairSel = find(salPairSel);

% %% Similarity of PRF and CGF upon salicylate
nClu = size(pairTable,1);
for ii = 1:nClu
    sel = data.mouse == pairTable.mouse(ii) & ...
        data.cids == pairTable.cids(ii);
    sel1 = sel & data.set == 4 ;
    sel2 = sel & data.set == 10;
    pairTable.rho_STRFs(ii) = corr(reshape(data.STRF{sel1},[],1), reshape(data.STRF{sel2},[],1));
    pairTable.rho_PRFs(ii) = corr(reshape(data.PRF{sel1},[],1), reshape(data.PRF{sel2},[],1));
    pairTable.rho_CGFs(ii) = corr(reshape(data.CGF{sel1},[],1), reshape(data.CGF{sel2},[],1));
end

end