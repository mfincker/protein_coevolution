fileName = 'UniqueProteinsDB_Maeva/uniqueMembraneProteinsDB_allMaeva.txt';
uMembraneSectorDB_11 = importdata('./uMembraneSectorDB/uMembraneSectorDB_11GroupSubMPRP.mat');

pdbIdList = cell(1,numel(uMembraneSectorDB_1));
for i = 1:numel(uMembraneSectorDB_1)
    pdbIdList{i} =  uMembraneSectorDB_1{i}.Pdb;
end



file = fopen(fileName);

pdbId = fgetl(file);
while (strcmp(pdbId,'') == 0 && ischar(pdbId))
    group = fgetl(file);
    subgroup = fgetl(file);

    memberProtein = fgetl(file);
    if isempty(memberProtein)
        memberProteins = regexp(memberProtein, '\t', 'split');
    else
        memberProteins = {};
    end

    relatedProtein = fgetl(file);
    if isempty(relatedProtein)
        relatedProteins = regexp(relatedProtein, '\t', 'split');
    else
        relatedProteins = {};
    end

    % Find sectors
    indexes = strcmp(pdbIdList, pdbId);
    indexes = find(indexes);
    for j = indexes
        uMembraneSectorDB_1{j}.Group = group;
        uMembraneSectorDB_1{j}.Subgroup = subgroup;
        uMembraneSectorDB_1{j}.MemberProteins = memberProteins;
        uMembraneSectorDB_1{j}.RelatedProteins = relatedProteins;
    end

    pdbId = fgetl(file);
end
    

save('./uMembraneSectorDB/uMembraneSectorDB_11.mat', 'uMembraneSectorDB_1');