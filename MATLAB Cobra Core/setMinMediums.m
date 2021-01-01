function [ minUptMediums, minUptMedSolutions ] = setMinMediums( model, essentialRxns, carbonTrnsRxns, bounds)
%setMinMediums Makes minimal uptake mediums based on a set of
%essential reactions and a set of one-carbon source transport
%reactions.
%   Detailed explanation goes here

minUptMediums = cell([length(carbonTrnsRxns),1]);
minUptMedSolutions = cell([length(carbonTrnsRxns),1]);
model = changeObjective(model, model.rxns(1693));
for i = 1:length(minUptMediums)
    minUptMediums{i} = model;
    minUptMediums{i} = changeRxnBounds(minUptMediums{i},... 
                       essentialRxns, -1000, 'l');
    minUptMediums{i} = changeRxnBounds(minUptMediums{i},... 
                       essentialRxns, 1000, 'u');
    for j = 1:length(minUptMediums)
        if j == i
            minUptMediums{i} = changeRxnBounds(minUptMediums{i},... 
            minUptMediums{i}.rxns(carbonTrnsRxns(j)), bounds, 'l');
            minUptMediums{i} = changeRxnBounds(minUptMediums{i},... 
            minUptMediums{i}.rxns(carbonTrnsRxns(j)), 1000, 'u');
        else
            minUptMediums{i} = changeRxnBounds(minUptMediums{i},... 
            minUptMediums{i}.rxns(carbonTrnsRxns(j)), 0, 'l');
        end
    end
    minUptMedSolutions{i} = optimizeCbModel(minUptMediums{i});
end

