function [ membrane ] = getMembrane( data )
%% GETMEMBRANE parses the pdb matlab structure entered as argument
% and look for the word membrane in the keywords
% Returns 1 if it finds it, 0 otherwise.
	keyWords = data.Keywords;
	% Reformat the character array to be searchable
	keyWords = keyWords';
	keyWords = reshape(keyWords,1,numel(keyWords));
	membrane = [strfind(keyWords, 'MEMBRANE') strfind(keyWords, 'membrane')];
	membrane = (sum(membrane) > 0); 
end