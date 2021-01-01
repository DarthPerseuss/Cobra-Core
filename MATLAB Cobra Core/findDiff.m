function [ differences ] = findDiff( solutions, model )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 
    i = size(solutions, 1);
    j = size(solutions, 2);
    differences = cell(i,j);
    FBAsol = optimizeCbModel(model, 'max');
    
    for i = 1:size(solutions, 1)
        for j = 1:size(solutions, 2)
            if isstruct(solutions{i, j}) == 1
                if FBAsol.f - solutions{i,j}.f >= 10^-13
                    differences{i,j} = solutions{i,j};
                else
                    differences{i,j} = NaN;
                end
            end
        end
    end
end

