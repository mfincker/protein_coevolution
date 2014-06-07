function [ EC ] = getEC(data, molID)
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