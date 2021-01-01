function [ unnested_cell ] = unnestCell( cell_array )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for oIndex = 1:length(cell_array)
   iIndex(oIndex,1) = (length(cell_array{oIndex,1}));
end
sumIndex = sum(iIndex);
unnested_cell = cell(sumIndex,1);
for i = 1:length(cell_array)
    for j = 1:length(cell_array{i,1})
        unnested_cell{i+j,1} = cell_array{i,1}{j,1};
    end
end
end

