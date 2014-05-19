function [clusters,extra] = blast2clust(seqacc) 

[msa,seedseq,seednum] = blast2msa(seqacc,500,1);

% seqacc must be a pdbId in that case !!!!
msa = msa_sort_maeva(seqacc, msa);
% BRAD PRUNING FUNCTION GOES HERE
% e.g. msa2 = phyloprunemsa(msa,cutoff)
[clusters,extra] = miscwrapper(msa,seedseq,seednum);

end