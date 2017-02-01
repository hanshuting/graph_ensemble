function sMatEFI = discEFI(retMat, pMat)
%The naive equal frequency interval discretization.
%pMat is a Cell; each component gives the desired quantile cut-off probabilities
%e.x.  [0.25 0.75 1]
statesMat = [];
for ticker = 1:size(retMat,2)
    tikStatesVec = [];
    
    %tickProbVec = pMat(:,ticker);
    tikQuantiles = quantile(retMat(:,ticker), pMat{ticker}); %quantiles for the current ticker
    
    for k = 1:size(retMat,1)
        %tikState = -1; %unassigned
        
        %Decide which label
        tikState = min(find(tikQuantiles >= retMat(k,ticker)));
        tikStatesVec = [tikStatesVec;tikState];     
    end
    
    statesMat = [statesMat tikStatesVec];
    
end

sMatEFI = statesMat;
