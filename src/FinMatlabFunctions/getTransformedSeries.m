function transIndices = getTransformedSeries(indexName,indexType,rawTS)
%INPUT: indexName is char; indexType is char; rawTS is a column matrix

%OUTPUT: transIndices is a cell; transIndices{1} is a collection of
%transformed indices name; transIndices{2} is a matrix, with each column
%being one transformation.
%The i-th component in the name collection corresponds to the i-th column
%of the matrix in transIndices{2}.
%transIndices{3} gives the type of transformed TS under each method.


%The incorporated available transformation/transfromation rules are as of 9 Mar 2015.

transIndices = cell(3,1); 
transNames = {};
transMat = [];
tsType = {};

errorMes = 'Current index type is not decided for transformation yet.';

%The incorporated available transformation/transfromation rules are as of 9 Mar 2015.
%More could be incorporated.


transformedTS = []; %Might be the original TS.
switch indexType
    case {'price','vol'}   %percent change group
        transformedTS = diff(rawTS)./rawTS(1:end-1);
        tsType = {'percent change'};
    case {'ratio','rate'}  %first order difference group
        transformedTS = diff(rawTS);
        tsType = {'fDifference'};
    case {'percent change','fDifference'}
        transformedTS = rawTS(2:end); %no need to transform
        tsType = {indexType}; 
    otherwise
        display(errorMes);
        return
end

longTermScale = {};  %For index with a long term scale, also keep the abs value.
%Decide if the index has long term scale
switch indexType    
    case {'ratio','vol','rate'}
          longTermScale = true;    
    case {'price','percent change','fDifference'}
          longTermScale = false;        
    otherwise
        display(errorMes);
        return
end


if longTermScale  %then keep the original value also
    transNames = {indexName,strcat(indexName,'(Abs)')};
    transMat = [transformedTS rawTS(2:end)]; 
    tsType = {tsType{1},indexType};
else
    transNames = {indexName};   
    transMat = transformedTS;
end




transIndices{1} = transNames;
transIndices{2} = transMat;
transIndices{3} = tsType;
 

        
        
        
        