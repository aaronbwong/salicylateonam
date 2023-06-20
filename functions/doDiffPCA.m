function [NormFactor,D,V,VScore,Vlatent,VVarExp,VPC,VPC_Stm] = doDiffPCA(X,Stm,normalization,demean,useAllMd)

% X: input 
%       nStm x T x nClu;
% Stm: stimulus matrix  
% normalization: 1 = normalize to s.d.; 0 = none;

% VPC: 
%       nClu/nPC x (T*nStm)
% VPC_Stm: 
%       nClu/nPC x T x nStm
    verbose = 0;

    nStm = size(X,1);
    T = size(X,2);
    nClu = size(X,3);

    X = permute(X,[3,2,1]); % nStm * T * nClu -> nClu * T * nStm

    % normalization 
    if normalization == 1
        if verbose; disp('normalize data to s.d.'); end
        X = reshape(X,nClu,T*nStm);
        NormFactor = std(X,0,2);
        X = diag(1./NormFactor) * X;
        X = reshape(X,nClu,T,nStm);
    else
        if verbose; disp('no normalization'); end
        NormFactor = ones(nClu,1);
    end
    
% -- generate difference matrix --
    if verbose; disp('generate difference matrix'); end
    uMd = nonzeros(unique(Stm.Md));
    uMf = nonzeros(unique(Stm.Mf));
    if (useAllMd)
        stmIdx = (Stm.Md ~= 0);
    else
        stmIdx = (Stm.Md == max(uMd));
    end
    ctrlIdx = (Stm.Md == 0);
    D = mean(X(:,:,stmIdx),3) - mean(X(:,:,ctrlIdx),3);

    
    
% -- run PCA on difference matrix --
    if verbose;disp('PCA on difference matrix'); end
    [V,VScore,Vlatent,~,VVarExp,~] = pca(D','Centered',demean); 
    X = reshape(X,nClu,T*nStm);
    
% -- project to PCs ---
    if nargout < 7
        return
    end
    if verbose;disp('projecting data on PCs'); end
    if (demean)
        VPC = V' * (X-mean(X,2));
    else
        VPC = V' * X;
    end
    VPC_Stm = reshape(VPC,nClu,T,nStm);

end
