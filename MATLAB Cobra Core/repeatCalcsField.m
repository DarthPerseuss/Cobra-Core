function [ output_args ] = repeatCalcsField( model, inf_gdw, sup_gdw, repeats )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

field1 = 'modelVersion';
field2 = 'rxns';
field3 = 'mets';
field4 = 'S';
field5 = 'rev';
field6 = 'c';
field7 = 'metNames';
field8 = 'lb';
field9 = 'ub';
field10 = 'metCharge';
field11 = 'rules';
field12 = 'genes';
field13 = 'rxnGeneMat';
field14 = 'grRules';
field15 = 'subSystems';
field16 = 'confidenceScores';
field17 = 'rxnReferences';
field18 = 'rxnECNumbers';
field19 = 'rxnNotes';
field20 = 'rxnNames';
field21 = 'metCHEBIID';
field22 = 'metKEGGID';
field23 = 'metPubChemID';
field24 = 'metInChIString';
field25 = 'b';
field26 = 'description';

modelStruct = cell(repeats, 1);

for row = 1:repeats
    modelStruct{row, 1} = struct(field1,[],field2,[],field3,[],field4, [],field5,[],field6,...
    [],field7,[],field8,[],field9,[],field10,[],field11,[],field12,[],field13,[],field14,...
    [],field15,[],field16,[],field17,[],field18,[],field19,[],field20,[],field21,[],...
    field22,[],field23,[],field24,[],field25,[], field26,[]);
end

for c_row = 1:repeats
combinations{c_row, 1} = mat2cell(zeros(length(model.rxns), 1) , length(model.rxns));
end

for model_row = 1: length(model.rxns)
    
    combinations(model_row, 1) = -ranGenerator(inf_gdw, sup_gdw);
    output_args = combinations(model_row);
    
end

