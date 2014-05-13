function [ proteinLength ] = getProteinLength(data, chainIndex)
	if chainIndex == -1
		proteinLength = -1;
	else
		proteinLength = data.DBReferences(chainIndex).seqEnd - ...
                       		data.DBReferences(chainIndex).seqBegin +1;
    end
end