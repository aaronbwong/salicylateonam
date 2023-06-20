
function PCs = projectPC(KSD,W)
    nStm    = size(KSD,1);
    T       = size(KSD,2);
    nClu    = size(KSD,3);

    X = reshape(KSD,nStm*T,nClu); 
    
    PCs = reshape(W' * X',nStm,T);

end