function [combrefclusts,varargout] = miscwrapper(msa,refseq,refind)
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
    if nnz(msa(:,curpos)==msa(1,curpos))==length(msa(:,curpos))
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
numboot = round(.05*numseq);
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
rawrefclusts = {};
avemisc = median(Cpfinal(:));
stdmisc = std(Cpfinal(:));
threshmisc = .3*stdmisc;
reorder = 0;%zeros(1,n);
pos = 1;
for l = 1:length(clusts)
    % determine if it's a good cluster, and save it
    curkeep = Cpfinal(clusts{l},clusts{l});
    if median(curkeep(:)) - avemisc > threshmisc
        curclust = atsp(clusts{l});
        curclust(~curclust) = [];
        if length(curclust)>=3
            refclusts{end+1} = curclust;
            rawrefclusts{end+1} = clusts{l};
            % reorder to show clusters in matrix
%             if ~isempty(clusts{l})
                reorder(pos:(pos+length(clusts{l})-1)) = clusts{l};
                pos = pos + length(clusts{l});
%             end
        end
    end

end

% Now we combine those that have high co-coevolution
taken = zeros(1,length(refclusts));
combrefclusts = {};
for l = 1:length(refclusts)
    if ~taken(l)
        taken(l) = 1;
        combrefclusts{end+1} = refclusts{l};
        for g = (l+1):length(refclusts)
            if ~taken(g)
                curkeep = Cpfinal(rawrefclusts{l},rawrefclusts{g});
                if (median(curkeep(:)) - avemisc) > threshmisc
                    taken(g) = 1;
                    combrefclusts{end} = [combrefclusts{end} refclusts{g}];
                end
            end
        end
    end
end


% assign second argument if it was asked for.
if nargout == 2
    misc = struct([]);
    misc(1).rawmisc = Cp; % co-evolution matrix on all data
    misc(1).rawbootmisc = Cpfinal; % boot-strapped co-evolution matrix
    Cpfinal = Cpfinal(reorder,:);
    Cpfinal = Cpfinal(:,reorder);
    Cp = Cp(reorder,:);
    Cp = Cp(:,reorder);
    misc(1).bootmisc = Cpfinal;  % cluster-reordered bootstrapped matrix
    misc(1).misc = Cp; % cluster-reordered matrix
    misc(1).rawats = atsp; % index-to-residue index mapping
    misc(1).ats = atsp(reorder); % reordered index-to-residue index mappin
    misc(1).reorder = reorder; % how things were reoridered.
    misc(1).ignored_residues = removed_residues(removed_residues~=0); % ignored protein residues numbers
    misc(1).included_residues = atsp(atsp~=0); % included protein residue numbers (not necessarily in a cluster)
    misc(1).rawclusts = refclusts; % clusters recombined if they have cross-coevolution
    misc(1).msap = msap; % pruned msa
    misc(1).msa = msa; % input msa
    varargout{1} = misc;
    
end

end