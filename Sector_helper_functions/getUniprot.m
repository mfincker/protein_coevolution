function [ uniprotID ] = getUniprot( data )
%GETUNIPROT returns the uniprot id from a pdb file.
% The pdb file must describe a monomer !
%   data is the structure returned by getpdb() or pdbread().
%	Returns an array containing all the uniprot id of the 
%	pdb file.

	uniprot = (regexp([data.DBReferences.database],'UNP   ')-1)/6 +1;

    uniprotID = {};

    for i = 1:size(uniprot,2)
    	id = data.DBReferences(i).dbAccession(1:6);
    	if sum(strcmp(id, uniprotID)) == 0
        	uniprotID = [uniprotID ; {id}]
        end
    end


end