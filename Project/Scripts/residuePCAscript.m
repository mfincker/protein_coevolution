% Residue PCA script

aaFreq_all = aaSectorDBFreq(membraneSectorDB);

% Normalize by removing the average aa freq
[eigVect_all, eigVal_all, aaEigenBase_all] = residuePCA(aaFreq_all, 2, 1);
%%
aaFreq_under100 = aaSectorDBFreq(membraneSectorDB_under100);
[eigVect_under100, eigVal_under100, aaEigenBase_under100] = residuePCA(aaFreq_under100, 2, 1);
%%
aaFreq_over100 = aaSectorDBFreq(membraneSectorDB_over100);
[eigVect_over100, eigVal_over100, aaEigenBase_over100] = residuePCA(aaFreq_over100, 2, 1);
%%
aaFreq_over10 = aaSectorDBFreq(membraneSectorDB_over10);
[eigVect_over10, eigVal_over10, aaEigenBase_over10] = residuePCA(aaFreq_over10, 2, 1);