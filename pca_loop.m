function [ output_args ] = untitled( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
PCAresults = zeros(numel(sectorDatabase), 3) % to keep track of the  lambdas for each sector.
for i = 1:numel(sectorDatabase)
    PCAresults( i, :) = pca( sectorDatabase{i} )
end

% visualize the data in a meaningful way.

end

