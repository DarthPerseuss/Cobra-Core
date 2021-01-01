function [ initCore, fluxes] = findRatio(solutions)
%findRatio Finds ratios of active reactions to find the core
%   The reactions that are active meaning that they have non-zero
%   fluxes under all experimental mediums are assigned to the initial
%   metabolic core of an organism

fluxes = zeros(length(solutions), length(solutions{1}.x));
for i = 1:length(solutions)
    for j = 1:length(solutions{1}.x)
        if solutions{i}.x(j) <= -10^-6 || solutions{i}.x(j)>= 10^6
            fluxes(i,j) = 1;
        end
    end
end

sumActiveFluxes = sum(fluxes);
initCore = [];
for ind = 1:length(sumActiveFluxes)
    if sumActiveFluxes(ind) == length(solutions)
        initCore = [initCore ind];
    end
end
initCore = transpose(initCore);
end
