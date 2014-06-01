%% Phenotypic enrichment for each sector
%This section is meant to find the number of different phenotypes expressed
%in each sector, and the relative enrichment of each phenotype within each
%sector.
% Sector_phenotypes = cell(2,1,6);
%Create cell arrays for each sector and their phenotypes

%I have to write the second phenotype array first because its size dicates
%the entire size of the array
test = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==2),16);
Sector_phenotypes(:,3,2) = test;

test = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==1),16);
Sector_phenotypes(1:length(test),3,1) = test;

test = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==3),16);
Sector_phenotypes(1:length(test),3,3) = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==3),16);

test = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==4),16);
Sector_phenotypes(1:length(test),3,4)= somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==4),16);

test = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==5),16);
Sector_phenotypes(1:length(test),3,5) = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==5),16);

test = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==6),16);
Sector_phenotypes(1:length(test),3,6) = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==6),16);

comparison = '';
for j = 1:6
    for i = 1:length(Sector_phenotypes)
        if Sector_phenotypes
        temp_name = Sector_phenotypes{i,1,j};
        





    end
end
