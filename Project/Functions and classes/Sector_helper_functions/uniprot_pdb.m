%This function takes a uniprot ID number, finds the protein associated to
%it, and then searches the file to determine if there is a PDB file for the
%Uniprot protein

function output_file = uniprot_pdb(unitprot_name)

genpept_name = getgenpept(unitprot_name);

%This part of the code searches for any pdb file associated with this
%protein. From experience I know that if there are pdb files associated to
%this protein, the first entry is most relevant to the unitprot protein.
search_string = genpept_name.DBLink;

%I know that the pdb filenames are arranged in this manner PDBsum:filename.
%I will exploit this using the strfind function to find all of the indices
%for which this is true. If the strfind function returns an empty array,
%there is no pdb file, and the function will end with a prompt saying such.
%If there are indices, I will proceed to extract the pdb file name and then
%return the file as an output.
pattern = 'PDBsum:';
patternlength = length(pattern);
indices = strfind(search_string,pattern);
size_ind = size(indices);

if size_ind(1) ~= 0
    startingpoint = indices(1)+patternlength;
    endpoint = startingpoint + 3;
    output_file = search_string(startingpoint:endpoint); 
else
   fprintf('There is no PDB file for this protein \n'); 
   output_file = '';
end
end
