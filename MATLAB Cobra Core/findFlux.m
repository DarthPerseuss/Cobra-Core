function [core] = findFlux(solution)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


for i = 1:size(solution,1)
    for j = 1: size(solution,2)
        if solution{i, j}.f >= 10^-25;
            core{i,j} = 1;
        else
            core{i,j} = 0;
           continue
        end
    end 
end
end

