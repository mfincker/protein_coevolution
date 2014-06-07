function [ molID ] = getMolID(data, chainID)
%GETMOLID returns the mol Id letter corresponding to the molID
% of a peptide, characterized by its uniprot accession number, 
%in the given pdb matlab struct.
%   Arguments:
%       - data : pdb matlab struct
%       - chainID : char array of a uniprot number
	if (sum(strcmp(chainID,'protein complex') + strcmp(chainID,'undefined'))>0)
		molID = -1;
	else
		% Find mol_id corresponding to the chainID:
        cmpd = data.Compound;
        cmpd = cmpd';
        cmpd = reshape(cmpd,1,numel(cmpd));
        cmpdFlip = fliplr(cmpd);
        molID = regexp(cmpdFlip,['(' chainID '|' chainID ...
                '[A-Za-a0-9 ,]*?)\s:NIAHC.*?;(?<mol_id>\d*)\s:DI_LOM'],'names');
        % should only return a single mol_id, but might bug if same chainID in two 
        % molecules ...
        molID = fliplr(molID.mol_id);
    end
end