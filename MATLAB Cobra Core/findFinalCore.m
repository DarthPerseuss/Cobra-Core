function finalCore = findFinalCore( model, initialCore)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

model = changeObjective(model, model.rxns(1693));
originalSolution = optimizeCbModel(model);
targetModels = cell([length(initialCore), 1]);
targetSols = cell([length(initialCore), 1]);
finalCore = initialCore;

for i = 1:length(initialCore)
    targetModels{i} = model;
    targetModels{i} = changeRxnBounds(targetModels{i},... 
            targetModels{i}.rxns(initialCore(i)), 0, 'b');
    targetSols{i} = optimizeCbModel(targetModels{i});
    if targetSols{i}.f == originalSolution.f
        finalCore(i) = [];
    end
end

        







