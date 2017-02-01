function bnryNodeNames = getBnryNodeNames(transNames,method,param) 
%INPUT:
%transNames is a cell, each component is a char 
%tsType is a cell, each component is a char
%method is char
%param is method-dependent
%Relationship: tsType{i} gives type of index named transNames{i}.

%OUTPUT:
%bnryNodeNames is a cell
%bnryNodeNames{i} is a cell of char giving node names of the index named transNames{i}

bnryNodeNames = cell(length(transNames),1);

switch method
   
    case 'EFI'
        for k =1:length(transNames)
            nodeNameCell = cell(0);
            currentName = transNames{k};
            for j = 1:length(param{k})
                nodeNameCell = [nodeNameCell strcat(currentName,'-',num2str(j))];
            end
            bnryNodeNames{k} = nodeNameCell;
        end
    case 'EWI'
        for k =1:length(transNames)
            nodeNameCell = cell(0);
            currentName = transNames{k};
            for j = 1:param(k)
                nodeNameCell = [nodeNameCell strcat(currentName,'-',num2str(j))];
            end
            bnryNodeNames{k} = nodeNameCell;
        end
    otherwise
       display('This method currently not available!');
       return
end
