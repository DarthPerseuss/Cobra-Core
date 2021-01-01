function[opt_solutions] = repeatCalcs(model, inf_gdw, sup_gdw, repeats)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

i = repeats;
j = length(model.rxns);

field1 = 'mets';
field2 = 'metNames';
field3 = 'metFormulas';
field4 = 'genes';
field5 = 'grRules';
field6 = 'rxns';
field7 = 'rxnNames';
field8 = 'subSystems';
field9 = 'csense';
field10 = 'S';
field11 = 'lb';
field12 = 'ub';
field13 = 'b';
field14 = 'c';
field15 = 'rev';
field16 = 'description';

modelStruct(j*i) = struct(field1,[],field2,[],field3,[],field4, [],field5,[],field6,...
    [],field7,[],field8,[],field9,[],field10,[],field11,[],field12,[],field13,[],field14,...
    [],field15,field16,[]);
    

solutions = cell(j, i);
for j = 1:length(model.rxns)
    for i = 1:repeats
        
        combinations = ranGenerator(inf_gdw, sup_gdw, i);
        modelStruct(j*i) = changeRxnBounds(model, model.rxns(j),...
            combinations(i), 'u');  
        %modelStruct(j*i) = changeObjective(model,);
        solutions{j, i} = optimizeCbModel(modelStruct(j*i));
        opt_solutions = solutions;
    end
end
%cell2struct(solutions);
end

