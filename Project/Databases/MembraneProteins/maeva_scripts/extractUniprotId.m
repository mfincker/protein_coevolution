% Extract all uniprot ID:

uniprotIDs = {};
nonUniprotID = {};

for i = 1:numel(DB)
    id = DB{i}.Uniprot;
    if size(id,2) == 6
        if sum(strcmp(id, uniprotIDs)) == 0
            uniprotIDs = [uniprotIDs ; {id}];
        end
    end
end
    
f = fopen('uniprotId_300514_forCoord.txt','w');
format = '%s\n';
for row = 1:numel(uniprotIDs)
    fprintf(f,format,uniprotIDs{row});
end

