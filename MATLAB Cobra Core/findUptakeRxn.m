function [ uptakeRxns ] = findUptakeRxn( model )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

newModel = cell(length(model.rxns),1);
FBAsolutionArray = cell(length(model.rxns),1);
uptakeRxns = zeros(length(model.rxns),1);

for rxn = 1:length(model.rxns)
    newModel{rxn,1} = changeRxnBounds(model,model.rxns{rxn,1},-0.01,'b');
end
for sol = 1:length(model.rxns)
    FBAsolutionArray{sol,1} = optimizeCbModel(newModel{sol,1});
    if isempty(FBAsolutionArray{sol,1}.x) == true
        uptakeRxns(sol,1) = 0;
    else
        uptakeRxns(sol,1) = 1;
    end
end
end

