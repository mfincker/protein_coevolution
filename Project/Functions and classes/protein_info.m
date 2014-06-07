%% Protein structures
%This script is supposed to take file names of our proteins of interest,
%parse out the protein names, and then create a structure with the fields
%we want. I put arbitrary fields in the body of the for loop. We can decide
%what other fields we want.


files = dir('*-aligned.fasta');
file_length = length(files);
proteins = struct();

for i = 1:file_length
    
    file_names = files(i).name;        
    protein_name = file_names(1:4);

    temp_protein_pdb = getpdb(protein_name);
    proteins(i).Header = temp_protein_pdb.Header;
    proteins(i).Sequence = temp_protein_pdb.Sequence;
    if isfield(temp_protein_pdb,'Helix')
        proteins(i).Helix = temp_protein_pdb.Helix;
    end
    if isfield(temp_protein_pdb,'Sheet')
        proteins(i).Sheet = temp_protein_pdb.Sheet;
    end
    proteins(i).Name = protein_name;
    %And any additional information we may need or want
    
end
