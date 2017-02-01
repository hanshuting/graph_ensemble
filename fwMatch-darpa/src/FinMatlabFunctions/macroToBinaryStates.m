function macroBStates = macroToBinaryStates(tradeDates,macroDates,macroStatesMat,sDimVec)
%Each non-news-release date has the state "NR" being on.
%sDimVec(k) is number of states for ticker k, including the "NR" state.

T = length(tradeDates);

%number of time series
bnryStates = cell(size(macroStatesMat,2),1);

for ticker = 1:length(bnryStates)
    tikStDim = sDimVec(ticker);
    
    tikBStates = zeros(T,tikStDim);  %a matrix of zeros
    
    for t = 1:T
        
        ind = find(strcmp(macroDates,tradeDates{t}));
        
        if isempty(ind)
            tikBStates(t,1) = 1;  % "NR" is on.   
        else        
        tikBStates(t,macroStatesMat(ind,ticker)+1) = 1;  
        end
    end
    
    bnryStates{ticker} = tikBStates;    
end

macroBStates = bnryStates;