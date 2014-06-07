function [ organismID ] = getOrganismID(data, molID)
%% GETORGANISMID returns the taxonomical id number of the
% organism in which the peptide, characterized by its molID 
% and the pdb record in the form of a pdb matlab structure, is.
	if molID == -1
		organismID = -1;
	else
		src = data.Source;
	    % Reformat the character array to be searchable
	    src = src';
	    src = reshape(src, 1, numel(src));
	    srcSplit = regexp(src, 'MOL_ID: ', 'split');
	    taxonomyId = srcSplit{find( strncmp(srcSplit, [molID ';'], size(molID, 2) + 1) == 1, 1)};
	    taxonomyId = regexp(taxonomyId, 'ORGANISM_TAXID:\s(?<num>\d+)','names');
	    organismID = str2num(taxonomyId.num);
	end
end