function candidate_residues = find_candidate_residues(sector_residues,...
    residues_with_mutations,M,cap,ats)

%% Candidate residues function
%This function takes user inputs to create a list of candidate residue
%sites for compensatory mutations for all mutations in 
%residues_with_mutations.

% sector_residues is a vector of residue numbers of a sector
%residues_with_mutations is a matrix with the residues in the first column
%that have mixed mutations in other species, the second column is the
%wildtype AA at that residue, the third column is the mutation itself, and
%the fourth column is the number of sequences that fix the corresponding
%mutation. M is the covariation matrix
%cap is the number set by the user to recover residues above a certain
%number of mutations that are fixed in other species.
%ats is the position index map for p53

candidate_residues = cell(length(residues_with_mutations(:,1)));

for i = 1:length(residues_with_mutations(:,1))
    %find 
    resind_temp = find(ats==residues_with_mutations{i,1});
    %The value .3 is a relatively high coevolution score
    resinds  = find(M(resind_temp,:)>.3);
    residues = ats(resinds);
    residues = residues(residues >= cap);
    counter = 0;
    trim_residues= ismember(sector_residues,residues);
    new_kept(1) = 0;
    for j = 1:length(trim_residues)
        if(trim_residues(j) == 1 && (residues(j)~=residues_with_mutations{i,1} &&...
                residues(j)~=0))
        counter = counter + 1;
        new_kept(counter) = residues(j);
        end
    end
    candidate_residues{i,1} = residues_with_mutations{i,1};
    candidate_residues{i,2} = residues_with_mutations{i,2};
    candidate_residues{i,3} = residues_with_mutations{i,3};
    candidate_residues{i,4} = new_kept;
    
end
end