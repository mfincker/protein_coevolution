function [ EC ] = getEC(data, molID)
%% GETEC returns the EC number (if it exists) of the protein
% described by molID in the pdb matlab structure entered as input.
% Input argument:
%   - data: pdb matlab structure
%   - molID: molID letter
% The function returns 'undefined' if the EC number for this peptide
% doesn't exist.

	if molID == -1
		EC = 'undefined';
	else
		cmpd = data.Compound;
        cmpd = cmpd';
        cmpd = reshape(cmpd,1,numel(cmpd));
		cmpdSplit = regexp(cmpd,'MOL_ID: ','split');
        ec = cmpdSplit{find( strncmp(cmpdSplit, [molID ';'], size(molID, 2) + 1) == 1, 1)};
        ec = regexp(ec, 'EC:\s(?<ec>\d+\.\d+\.\d+\.\d+);', 'names');
        if numel(ec) == 1
            EC = ec.ec;
        else
            EC = 'undefined';
        end
    end
end 