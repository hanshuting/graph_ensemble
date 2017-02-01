function bStates = toBinaryStates(statesMat,sDimVec)
%transform from finite states to binary states;
%This is done in the most naive way: if the original state space is of
%dimension 5, then 5 binary states will be generated.
%The returned bStates is a Cell with dimension being the number of columns
%in statesMat; each cell component is the binarized states for one ticker.

T = size(statesMat,1);

bnryStates = cell(size(statesMat,2),1); 

for ticker = 1:length(bnryStates)
    tikStDim = sDimVec(ticker);
    
    tikBStates = zeros(T,tikStDim);  %a matrix
    
    for t = 1:T
        tikBStates(t,statesMat(t,ticker)) = 1;        
    end
    
    bnryStates{ticker} = tikBStates;    
end

bStates = bnryStates;