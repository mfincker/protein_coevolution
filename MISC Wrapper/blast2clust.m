function [clusters,extra] = blast2clust(seqacc) 

[msa,seedseq,seednum] = blast2msa(seqacc,500,1);
% BRAD PRUNING FUNCTION GOES HERE
% e.g. msa2 = phyloprunemsa(msa,cutoff)
msa2 = msa; % remove this line when Brad adds function
[clusters,extra] = miscwrapper(msa,seedseq,seednum);

end