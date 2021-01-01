function [ rxns ] = findRxnsFromNames( model, rxnNames )

ind = cell(length(rxnNames),1);
for i = 1:length(rxnNames)
    ind{i} = find(strcmp(model.rxnNames, rxnNames{i}));
    rxns{i} = model.rxns{ind{i}};
end
rxns = transpose(rxns);
end

