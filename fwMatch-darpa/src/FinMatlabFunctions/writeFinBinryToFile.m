function  writeFinBinryToFile(outputDates,allBnrySeries,outputFile)
%INPUT: output object from finBnry.m
%OUTPUT: an Excel sheet

fullNodeNameList = cell(0);
fullBnryMatrix = [];

%Turn output dates to 

fullCellList = [{'Time Stamps'};num2cell(outputDates)];


for ticker = 1:length(allBnrySeries)
    tickerInfo = allBnrySeries{ticker};
    subIndNames = tickerInfo{1}; %Each component is a list of nodes
    subIndBnryMat = tickerInfo{2};
    
    for subInd = 1:length(subIndNames)
        fullNodeNameList = [fullNodeNameList subIndNames{subInd}];
        fullBnryMatrix = [fullBnryMatrix subIndBnryMat{subInd}];
    end    
end
fullCellList = [fullCellList [fullNodeNameList;num2cell(fullBnryMatrix)]];
%size(fullCellList)
%size(fullNodeNameList)
%size(num2cell(fullBnryMatrix))
%size([fullNodeNameList;num2cell(fullBnryMatrix)])

%timeCol = [{'Time Stamps'};num2cell(outputDates)];
fid=fopen(outputFile,'w');
numNodes=length(fullNodeNameList);

fprintf(fid,'%s,','Time Stamps');
fprintf(fid,'%s,',fullNodeNameList{1:end-1});
fprintf(fid,'%s\n',fullNodeNameList{end});

for k = 1:size(fullBnryMatrix,1)
    fprintf(fid,'%d,',outputDates(k));
    fprintf(fid,'%d,',fullBnryMatrix(k,1:end-1));
    fprintf(fid,'%d\n',fullBnryMatrix(k,end));
end
fclose(fid);



%csvwrite(outputFile,'fullCellList');