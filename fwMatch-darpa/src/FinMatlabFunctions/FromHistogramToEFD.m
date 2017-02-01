
rVec = retMat(:,1);


%rVec contains daily returns for SP500, with no duplicate values
%300000 bins delibrately make each bin contains at most 2 points.
%hCounts stores number of points in each bin (mostly 0)
%hCenters stores the bin centers
M = 300000;
[hCounts,hCenters] = hist(rVec,M);
binWidth = (max(rVec) - min(rVec))/M;

%recover the boundary points of bins
%Right boundary for i-th bin is hBounds(m + 1)
%pVec contains the bin frequency
hBounds = [hCenters - binWidth/2 max(rVec) ];
pVec = hCounts/sum(hCounts);

%Suppose we want EFD cut-off points with 10% frequency (i.e. 10 bins)
K = 10;
efdCuts = [];
sVec = pVec(1);  %running sum of pVec
for k = 2:M
    sVec = [sVec; sVec(end) + pVec(k)];    
end
%Find the cut-off points
for m = 1:K-1
    efdCuts = [efdCuts;hBounds(min(find(sVec >= m/K)) + 1)]; %add one cut point  
end
efdCuts = [efdCuts;max(rVec)]; 


%%%%%%%%%%%%%%%%%%%Plot the graph%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 30;
[hCounts,hCenters] = hist(rVec,M);
pVec = hCounts/sum(hCounts);

plot(hCenters, pVec)
hold on
plot(efdCuts,zeros(length(efdCuts),1),'rx')






