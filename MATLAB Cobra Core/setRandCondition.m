function [solutions, new_models] = setRandCondition( model, inf_gdw, sup_gdw, uptakeRxns, repeats)
%Randomize growth environments for the given model.
%   sets the lower bounds of uptake reactions of a given number models
%   to random values between inf_gdw(negative number) and sup_gdw.

model = changeObjective(model, model.rxns(1693));
new_models = cell([repeats,1]);
solutions = cell([repeats,1]);

for i = 1:repeats
    new_models{i} = model;
end

for i = 1:repeats
    for j = 1:length(uptakeRxns)
        if uptakeRxns(j) == 1
            new_models{i} = changeRxnBounds(new_models{i},... 
            new_models{i}.rxns(j),... 
            ranGenerator(inf_gdw, sup_gdw), 'l');
        else
            continue
        end
    end
end
for i = 1:length(new_models)
    solutions{i} = optimizeCbModel(new_models{i});
end
end
