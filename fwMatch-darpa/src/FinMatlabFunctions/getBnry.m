function bnryCollection = getBnry(transMat,method,param)
%INPUT: 
%transMat is a matrix, each column is for one transformed TS 
%tsType is a cell with each component being a char 
%method is char
% Relationship: tsType{i} gives type of the TS in transMat(:,i)
%param is method-dependent.

%OUTPUT
%bnryCollection is a cell
%bnryCollection{i} is a 0-1 matrix corresponding to the binarization on 
%transMat(:,i) following method with param.


    switch method
        case 'EFI'
             discTS = discEFI(transMat, param);
             sDimVec = [];
             for k = 1:length(param)
                 sDimVec = [sDimVec;length(param{k})];              
             end
             bnryCollection = toBinaryStates(discTS,sDimVec);
        case 'EWI'
             discTS = discEWI(transMat, param);
             bnryCollection = toBinaryStates(discTS,param);
        otherwise
            display('This method currently not available!');
            return        
    end
