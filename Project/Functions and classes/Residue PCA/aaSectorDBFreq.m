function [ aaMatrixFreq] = aaSectorDBFreq ( sectorDatabase )
%AASECTORDBCOUNT returns a 20 x n matrix where n is the number
%of sectors in sectorDatabase. Each row corresponds to an amino
%acid - the order of the aa is the same as the order of aa
%from the aacount() output structure. The cell (i,j) is populated
%with the number of occurence of amino acid (i) in sector (j).
	aaMatrixFreq = zeros(20, numel(sectorDatabase));
	for i = 1:numel(sectorDatabase)
		seq = sectorDatabase{i}.Sequence;
		count = aacount(seq);
		count = struct2cell(count);
		count = cell2mat(count);
		aaMatrixFreq(:,i) = count/sectorDatabase{i}.Length;
	end

end