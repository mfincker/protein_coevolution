counter = 0;
for j = 1:length(kept_residues)
    if(kept_residues(j) == 1)
    counter = counter + 1;
    new_kept(counter) = residues(j);
    end
end