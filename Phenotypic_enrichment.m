%% Phenotypic enrichment for each sector
%This section is meant to find the number of different phenotypes expressed
%in each sector, and the relative enrichment of each phenotype within each
%sector.
% Sector_phenotypes = cell(2,1,6);
%Create cell arrays for each sector and their phenotypes

%I have to write the second phenotype array first because its size dicates
%the entire size of the array


test = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==2),16);
Sector_phenotypes(:,4,2) = test;

test = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==1),16);
Sector_phenotypes(1:length(test),4,1) = test;

test = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==3),16);
Sector_phenotypes(1:length(test),4,3) = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==3),16);

test = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==4),16);
Sector_phenotypes(1:length(test),4,4)= somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==4),16);

test = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==5),16);
Sector_phenotypes(1:length(test),4,5) = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==5),16);

test = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==6),16);
Sector_phenotypes(1:length(test),4,6) = somatic_mutations_only_in_sectors(find(cell2mat(somatic_mutations_only_in_sectors(:,17))==6),16);


%%

for j = 1:6
    for i = 1:length(Sector_phenotypes)
        temp_name = Sector_phenotypes{i,4,j};
        if isempty(Sector_phenotypes{1,3,j})
        
        Sector_phenotypes{i,3,j} = temp_name;
        comparison = temp_name;
        Sector_phenotypes{1,1,j} = 1;%length of phenotype vector
        Sector_phenotypes{1,2,j} = 1;%phenotype count
        else
            for k = 1:Sector_phenotypes{1,1,j}
                if strcmp(temp_name,Sector_phenotypes{k,3,j})
                  Sector_phenotypes{k,2,j} = Sector_phenotypes{k,2,j} + 1;
                    break

                elseif ~strcmp(temp_name,Sector_phenotypes{k,3,j}) && k == Sector_phenotypes{1,1,j}
                    Sector_phenotypes{k+1,3,j} = temp_name;
                    Sector_phenotypes{1,1,j} = Sector_phenotypes{1,1,j} +1;
                    
                end
            end
        end
    end
end

%%
    [Sector1_sort,I1] = sort(cell2mat(Sector_phenotypes(:,2,1)),'descend');
    [Sector2_sort,I2] = sort(cell2mat(Sector_phenotypes(:,2,2)),'descend');
    [Sector3_sort,I3] = sort(cell2mat(Sector_phenotypes(:,2,3)),'descend');
    [Sector4_sort,I4] = sort(cell2mat(Sector_phenotypes(:,2,4)),'descend');
    [Sector5_sort,I5] = sort(cell2mat(Sector_phenotypes(:,2,5)),'descend');
    [Sector6_sort,I6] = sort(cell2mat(Sector_phenotypes(:,2,6)),'descend');


%% Create trimmed database of unique mutations for each residue


