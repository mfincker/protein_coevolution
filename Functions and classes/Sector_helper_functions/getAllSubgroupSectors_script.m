
%% get list of subgroups
redundantsubgroups = {};
for i = 1:length(membraneSectorDB)
    redundantsubgroups{i}= membraneSectorDB{i}.Subgroup;
end
%%
% get unique subgroups
subgroups= unique(redundantsubgroups);



%%

subgroupsectors = {};
for i = 1:length(subgroups)
    subgroupSectors{i} = getSectorsBySubgroup(membraneSectorDB,subgroups{i});
end