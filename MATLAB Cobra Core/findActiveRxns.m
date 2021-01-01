function [ active_rxns, number ] = findActiveRxns( models, model)
%finActivetxns finds reactions, active under all given environments
%
%   Detailed explanation goes here

FBAsolution = cell(length(models),1);

for i = 1:length(models)
    FBAsolution{i,1} = optimizeCbModel(models{i,1}, 'max');   
end

x_len = length(model.rxns);
array_rxns = zeros(x_len,length(models));

for i = 1:x_len
    for j = 1:length(FBAsolution)
        
         if FBAsolution{j,1}.x(i) ~= 0 % what should this number be?%
                array_rxns(i,j) = 1;
         end 
  
    end
end

sums = sum(array_rxns,2);
q = zeros(x_len,1);
active_rxns = cell(length(x_len),1);

for entry = 1:length(sums)
    
    if sums(entry,1) == length(models) %*******%
       q(entry,1) = 1;
    end
    
    if q(entry,1) == 1
        active_rxns{entry,1} = model.rxnNames{entry,1};
    else
        active_rxns{entry,1} = 0;
    end
    
end
number = sum(q,1); 
end


