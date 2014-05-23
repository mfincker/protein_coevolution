function [ seq ] = getSequence( data, residueInd,...
							chainIndex, chainID)
	if (sum(strcmp(chainID,'protein complex') + strcmp(chainID,'undefined'))>0)
		seq = 'undefined';
	else
		seqIndex = find(strcmp({data.Sequence.ChainID}, chainID) == 1,1);
        seq = data.Sequence(seqIndex).Sequence(residueInd - ...
              data.DBReferences(chainIndex).seqBegin + 1);
    end
end