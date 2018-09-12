function sMatEWI = discEWI(retMat, binNum)
%The naive equal length interval discretization.
%binNum is a vector, i-th entry is the # of bins for the i-th time series

%binNum = 10*ones(size(retMat,2),1); %At this moment, naively set at 10 for each
statesMat = [];
%rtpVec = [];
for ticker = 1:size(retMat,2)
    tikStatesVec = [];
    
    [~, cters] = hist(retMat(:,ticker),binNum(ticker));
    
    for k = 1:size(retMat,1)
    tikState = -1; %unassigned
    rtp = min(find(cters >= retMat(k,ticker)));
   
    %Decide which bin
    if rtp == 1
        tikState = 1;
    elseif abs(retMat(k,ticker) - cters(rtp)) < abs(retMat(k,ticker) - cters(rtp-1))
        tikState = rtp;
    elseif abs(retMat(k,ticker) - cters(rtp)) > abs(retMat(k,ticker) - cters(rtp-1))
        tikState = rtp - 1;
    else  %rtp is empty (i.e. current number bigger than all bin centers)
        tikState = binNum(ticker); %Very last bin
    end    
    tikStatesVec = [tikStatesVec;tikState];
    end
    
    statesMat = [statesMat tikStatesVec];
end

sMatEWI = statesMat;