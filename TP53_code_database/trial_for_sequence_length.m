%%Trial script
[header sequence] = fastaread('BAC16799-500-aligned.fasta');
breakcount = 0;
for i= 1:length(sequence)
    if i < length(sequence)
        temp1 = length(sequence{1,i});
        temp2 = length(sequence{1,i+1});
        if temp1 ~= temp2
            breakcount = breakcount +1;
            breakage(breakcount) = i;
        end 
    else
        temp1 = length(sequence{1,i});
        temp2 = length(sequence{1,1});
        if temp1 ~= temp2
            breakcount = breakcount +1;
            breakage(breakcount) = i;
        end 
    end
end