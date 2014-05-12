function [refclusts,varargout] = miscwrapper(msa,refseq,refind)
%MISCWRAPPER    wrapper for Mutual Information Sequence Correlation.
%   sectors = MISCWRAPPER(msa,refseqind) for integer MSA (numseq x numpos)
%   and reference sequence residues refseq and position numbers refind will
%   return sectors -- a cell array with lists of residues which co-adapt. 
%
%   [sectors,misc] = MISCWRAPPER(msa,refseqind) also returns a structure 
%   with a numpos x numpos matrix of pair-wise residue co-adaptation scores
%   and the 1xnumseq mapping between the reference sequence and indices in 
%   the matrix, stored as misc.rawmisc and misc.rawats respectively. 
%   Both are re-ordered according to unfiltered sectors. There are also
%   the same versions of misc after bootstrapping, and after being
%   re-ordered according to the sectors.
%
%   Residue numbers are given relative to the reference sequence "refseq", 
%   with indexing according to "refind". MISCWRAPPER uses nwalign to find
%   correspondence between the MSa and the reference sequence.
%
%   MISCWRAPPER prunes fully conserved or highly gappy positions.
%
%   This function requires MeanShiftCluster (c) 2006 Bart Finkston;
%   You can find this function on Mathworks MATLAB Central
%   This function also requires mex-compiled misc.c library
[numseq,numpos] = size(msa);
msa(find(msa==25))=0; % change gaps to zeros
msa(find(msa==23))=0; % change "X" to zeros (unknown aa)

% get indices of reference sequence
ats = msaSeqRef(msa,refseq,refind);

removed_residues = refind(~ismember(refind,ats));

% make sure msa is the right type for c library
msa = uint32(msa);

% Prune away gap-mode positions
removegaps = mode(double(msa),1)==0;

% Prune away conserved positions
removecons = zeros(1,numpos);
for curpos = 1:numpos
    if mode(double(msa(:,curpos)))==0
        removecons(curpos) = 1;
    end
end
msap = msa(:,(find((~removegaps).*(~removecons))));
removed_residues = [removed_residues ats(find(~((~removegaps).*(~removecons))))];
atsp = ats((find((~removegaps).*(~removecons))));

% calculate MSA
n = size(msap,2);
Cp = reshape(1-misc_pc(msap),[n,n]);

% conduct bootstrapping
numboot = round(.01*numseq);
allCp = zeros(n,n,numboot);
for l = 1:numboot
    allCp(:,:,l) = reshape(1-misc_pc(msap(randi(numseq,1,numseq),:)),[n,n]);
end

Cpfinal = median(allCp,3);

% Calculate clusters
[~,~,clusts] = MeanShiftCluster(Cpfinal,1.4,0);

% renumber to reference sequence
% and remove clusters with zero residues
% or with average or below-average scores.
refclusts = {};
avemisc = median(Cpfinal(:));
stdmisc = std(Cpfinal(:));
reorder = 0;%zeros(1,n);
pos = 1;
for l = 1:length(clusts)
    % determine if it's a good cluster, and save it
    if (median(Cpfinal(clusts{l},clusts{l})) - avemisc) > .2*stdmisc
        curclust = atsp(clusts{l});
        curclust(~curclust) = [];
        if length(curclust)>4 && length(curclust) < 100
            refclusts{end+1} = curclust;
            % reorder to show clusters in matrix
%             if ~isempty(clusts{l})
                reorder(pos:(pos+length(clusts{l})-1)) = clusts{l};
                pos = pos + length(clusts{l});
%             end
        end
    end

end

% assign second argument if it was asked for.
if nargout == 2
    misc = struct([]);
    misc(1).rawmisc = Cp;
    misc(1).rawbootmisc = Cpfinal;
    Cpfinal = Cpfinal(reorder,:);
    Cpfinal = Cpfinal(:,reorder);
    Cp = Cp(reorder,:);
    Cp = Cp(:,reorder);
    misc(1).bootmisc = Cpfinal;
    misc(1).misc = Cp;
    misc(1).rawats = atsp;
    misc(1).ats = atsp(reorder);
    misc(1).reorder = reorder;
    misc(1).ignored_residues = removed_residues(removed_residues~=0);
    misc(1).included_residues = atsp(atsp~=0);
    varargout{1} = misc;
end

end