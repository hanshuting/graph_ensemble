function [bnryStates,nodeNames] = getOriginalBnryStates(rawTSMat,tikName)

%Quantization based on original values.

%INPUT: 
%rawTSMat is a N by 2 matrix containing the time series for one ticker .
%First column is the time stamps, in format of Unix time stamps
%Second column are the index values. The missing values should be
%a blank in the Excel file, which will be read-in as NaN.

%tikName is a string representing ticker's name.

%The default quantization is set as deciles.

%Output:
%bnryStates is the binarized states corresponding to the bucketing scheme used (deciles by default)
%nodeNames is a cell of strings looking like {'tikName-NaN', 'tikName-1','tikName-2'...}

coreTimeStps = [];
coreTS = [];
coreLoc = [];
coreDeltaTS = [];




for t=1:length(rawTSMat(:,1))
    if ~isnan(rawTSMat(t,2))
    coreTimeStps = [coreTimeStps;rawTSMat(t,1)];
    coreTS = [coreTS;rawTSMat(t,2)];
   
    end
end


%sum(isnan(coreDeltaTS))
%max(coreDeltaTS)
%min(coreDeltaTS)

%%%%%%%%Default quantization is based on deciles.%%%%%%
coreEFIStates = discEFI(coreTS, {[0.1:0.1:1]});
bnryStates = zeros(size(rawTSMat,1),11);
nodeNames = cell(1,11);
nodeNames{1} = strcat(tikName,'-NaN');
for k = 2:11
    nodeNames{k} = strcat(tikName,'-',num2str(k-1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t = 1:size(rawTSMat,1)
    currentTime = rawTSMat(t,1);
    
    tLoc = min(find(coreTimeStps == currentTime));
    
    if isempty(tLoc)
        bnryStates(t,1)=1;
    else
        %tLoc
        sLoc = coreEFIStates(tLoc);
        bnryStates(t,1+sLoc) = 1;
    end
end
end