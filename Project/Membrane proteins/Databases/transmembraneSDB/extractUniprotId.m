% Extract all uniprot ID from a SDB called DB.
% 06/01/2014 - Maeva

uniprotIDs = {};

for i = 1:numel(DB)
    id = DB{i}.Uniprot;
    % Check if it's really a uniprot ID
    if size(id,2) == 6
    	% Check that the uniprot is not already in the list
        if sum(strcmp(id, uniprotIDs)) == 0
            uniprotIDs = [uniprotIDs ; {id}];
        end
    end
end
    
% Write the list of uniprot id to a file
f = fopen('uniprotId_300514_forCoord.txt','w');
format = '%s\n';
for row = 1:numel(uniprotIDs)
    fprintf(f,format,uniprotIDs{row});
end

