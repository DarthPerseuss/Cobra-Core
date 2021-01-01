function [FBAsolution] = noUptake( model )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    solution = cell(length(model.lb), 1);
    for i = 1:length(model.lb)
        model = changeRxnBounds(model, model.rxns(i), 0, 'l');
        objRxn = checkObjective(model);
        model = changeObjective(model, objRxn);
        solution{i,1} = optimizeCbModel(model, 'max');
        FBAsolution = solution;
    end
end

