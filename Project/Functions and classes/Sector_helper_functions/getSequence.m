function [ seq ] = getSequence( data, residueInd,...
							chainIndex, chainID)
%% GETSEQUENCE parses the pdb structure 'data' and look for the
% residues corresponding to the indexes in 'residueInd' in the protein
% sequence characterized by the 'chainIndex' and 'chainID' in the 'data' structure.
	if (sum(strcmp(chainID,'protein complex') + strcmp(chainID,'undefined'))>0)
		seq = 'undefined';
	else
		seqIndex = find(strcmp({data.Sequence.ChainID}, chainID) == 1,1);
        seq = data.Sequence(seqIndex).Sequence(residueInd - ...
              data.DBReferences(chainIndex).seqBegin + 1);
    end
end