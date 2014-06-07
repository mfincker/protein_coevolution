function [ membrane ] = getMembrane( data )
	keyWords = data.Keywords;
	% Reformat the character array to be searchable
	keyWords = keyWords';
	keyWords = reshape(keyWords,1,numel(keyWords));
	membrane = [strfind(keyWords, 'MEMBRANE') strfind(keyWords, 'membrane')];
	membrane = (sum(membrane) > 0); 
end