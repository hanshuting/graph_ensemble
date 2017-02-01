function [allBnrySeries,nanInfo] = finBnry(fileName, granularity,method,outputFile)
%fileName is a string.
%granularity is a char specifying the current prevailing granularity:
%'daily','monthly','quarterly','yearly'. 
%Assumption is that the prevailing granularity is the finest among all.
%method is a char specifying method, currently available: "EFI" and "EWI"

%writeToFile is boolean: if true, an Excel file will be generated.

[~,sheets] = xlsfinfo(fileName);
  
[rawNumMat,rawTxt] = xlsread(fileName,sheets{1},'','basic');
indNames = rawTxt(1,2:end);
rawDates = rawNumMat(:,1);
rawIndMat = rawNumMat(:,2:end);
outputDates = rawDates(2:end);  


%[indexType updateFrequency]
%As of 9 Mar 2015.
%indexType: 'price', 'vol', 'percent change','ratio','rate','fDifference'
%updateFrequency: 'daily','monthly','quarterly'
[~,typeInfo] =  xlsread(fileName,sheets{2});

indNum = length(indNames);

nanInfo = cell(0);

for ticker = 1:indNum
    if strcmpi(typeInfo{ticker+1,3},granularity) && sum(isnan(rawIndMat(:,ticker))) > 0
        nanInfo = [nanInfo; {indNames{ticker},sum(isnan(rawIndMat(:,ticker)))}];
    end
end

if ~isempty(nanInfo)
    return
end


allBnrySeries = cell(indNum,1);
for ticker = 1:indNum
    currentFreq = typeInfo{ticker+1,3};
    currentName = indNames{ticker};
    currentType = typeInfo{ticker+1,2};
    currentRawTS = rawIndMat(:,ticker);
    isGranularity = strcmpi(currentFreq, granularity);
    
    %The parameter is decided by frequency.
    param = getMethodParam(method, currentFreq);
    
    switch isGranularity
        case true
            transIndices = getTransformedSeries(currentName,currentType,currentRawTS);
            paramForAll = cell(length(transIndices{1}),1);
            for subInd = 1:length(paramForAll)
                paramForAll{subInd} = param;
            end
            
            switch method
                case 'EWI'
                    paramForAll = cell2mat(paramForAll);
                otherwise
            end
            
             currentNodeNames = getBnryNodeNames(transIndices{1},method,paramForAll);
             currentBnryCollection = getBnry(transIndices{2},method,paramForAll);
             allBnrySeries{ticker} = {currentNodeNames, currentBnryCollection};  %the two components should have the same length.
                         
        case false
            %Need to handle NaNs
            %"core" means non-NaN
            [coreRawDates,coreRawTS] = getCoreTS(rawDates,currentRawTS);
            coreOutputDates = coreRawDates(2:end);
            
            transIndices = getTransformedSeries(currentName,currentType,coreRawTS);
            paramForAll = cell(length(transIndices{1}),1);
            for subInd = 1:length(paramForAll)
                paramForAll{subInd} = param;
            end         
            
            switch method
                case 'EWI'
                    paramForAll = cell2mat(paramForAll);
                otherwise
            end
            
             coreNodeNames = getBnryNodeNames(transIndices{1},method,paramForAll);
             coreBnryCollection = getBnry(transIndices{2},method,paramForAll);
             
             %need to add node "NR" and merge with main dates
             tranNames = transIndices{1};
             currentNodeNames = cell(length(coreNodeNames),1);
            for subTick =1:length(currentNodeNames)
                currentNodeNames{subTick} = [strcat(tranNames{subTick},'-NR') coreNodeNames{subTick}];
            end
            
            currentBnryCollection = mergeToDates(outputDates,coreOutputDates,coreBnryCollection);
            allBnrySeries{ticker} = {currentNodeNames, currentBnryCollection};                   
    end    
end

if nargin > 3
    writeFinBinryToFile(outputDates,allBnrySeries,outputFile);
end

end

function [coreDates,coreTS] = getCoreTS(rawDates,rawTS)
coreDates = [];
coreTS = [];

for t=1:length(rawDates)
    if ~isnan(rawTS(t))
    coreDates = [coreDates;rawDates(t)];
    coreTS = [coreTS;rawTS(t)];
    end
end

end

function mergedBnryCollection = mergeToDates(outputDates,coreOutputDates,coreBnryCollection)
mergedBnryCollection = cell(length(coreBnryCollection),1);

for subInd = 1:length(mergedBnryCollection)
    currentCoreBnryMat = coreBnryCollection{subInd};
    currentBnryMat = zeros(length(outputDates),size(currentCoreBnryMat,2)+1);
    %size(currentBnryMat)
    for t = 1:length(outputDates)
        dateLoc = find(coreOutputDates==outputDates(t));        
        if isempty(dateLoc)
            currentBnryMat(t,1) = 1;  %"NR"
        else       
           stateLoc = find(currentCoreBnryMat(dateLoc,:) ==1 );
           currentBnryMat(t,stateLoc + 1) = 1;
        end
    end
    %size(currentBnryMat)
    mergedBnryCollection{subInd} = currentBnryMat;
end
end



