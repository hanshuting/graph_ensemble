%allTimes are the time stamps appearing from either financial and NLP side.
%All in "double"
tic
%finTimes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setu p files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outputFile = '/data/combined_data/finNLPCombo(finProblematic).xlsx';
nlpFile = 'data/combined_data/nlp_hrly_jul16.xlsx';
finFile = 'data/combined_data/FinEconBnryHrly(Problematic).xlsx';

% Read the time strings
[~,finTimesStrAll,~] = xlsread(finFile, 'A:A');
finTimesStr = finTimesStrAll(2:end);
finTimes = datenum(finTimesStr);

[~,nlpTimesStrAll,~] = xlsread(nlpFile, 'A:A');
nlpTimesStr = nlpTimesStrAll(2:end);
nlpTimes = datenum(nlpTimesStr);

allTimes = unique([finTimes;nlpTimes]);

keyboard;




nlpLastRowNum = length(nlpTimes) + 1;
finLastRowNum = length(finTimes) + 1;

rowNum = length(allTimes) + 1;

xlswrite(outputFile,[cell(0) 'time stamps'],'A1:A1');

allTimesStr = cell(length(allTimes),1);

for k = 1:length(allTimesStr)
    allTimesStr{k} = datestr(allTimes(k),'dd-mmm-yyyy HH:MM:SS');
end

timeColumn = strcat('A2:A',num2str(rowNum));
xlswrite(outputFile,allTimesStr,timeColumn);

outputCol = 2;

%handle NLP first.
nlpLastColName = 'GS';
nlpLastColNum = xlscol(nlpLastColName);

nlpCol = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle dates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

finTimes = datenum(finTimesStr);

%nlpTimes
%nlpTimes = datenum(nlpTimesStr);

%allTimes
allTimes = unique([finTimes;nlpTimes]);

tLocFin = []; %assign each time stamp in allTimes a location in the financial time series
tLocNLP = []; %....NLP....

for k = 1:length(allTimes)
    finLoc = min(find(finTimes == allTimes(k)));
    nlpLoc = min(find(nlpTimes == allTimes(k)));
    
    if isempty(finLoc)
        tLocFin = [tLocFin; -1 ];
    else
        tLocFin = [tLocFin; finLoc];
    end
    
    if isempty(nlpLoc)
        tLocNLP = [tLocNLP; -1];
    else
        tLocNLP = [tLocNLP; nlpLoc];
    end   
end


while nlpCol <= nlpLastColNum
    fprintf('nlpCol: %d / %d\n', nlpCol, nlpLastColNum);
    %read-in the current node name
    [~,nlpNode] = xlsread(nlpFile,strcat(xlscol(nlpCol),'1',':',xlscol(nlpCol),'1'));
    nlpNode = nlpNode{1};
    nlpNodeNum = str2num(nlpNode(end));
    
    %binary vector of current node
    nlpVec = xlsread(nlpFile,strcat(xlscol(nlpCol),'2',':',xlscol(nlpCol),num2str(nlpLastRowNum)));
    
    %The first node; an NaN state should be added before it.
    if nlpNodeNum ==1
    
        %display('NLP NaN')
        
        nlpNaNNode = strcat(nlpNode(1:end-2),'_NaN'); %make the NaN node
    
        nlpNodeMat = zeros(length(allTimes),2);
    
        for t = 1:length(allTimes)
            if tLocNLP(t) < 0  % -1, no value for this time stamp
                nlpNodeMat(t,1) = 1; %NaN state is on
            else
                nlpNodeMat(t,2) = nlpVec(tLocNLP(t)); %0 or 1
            end
        end
     
        %write to file
        xlswrite(outputFile,[cell(0) nlpNaNNode],strcat(xlscol(outputCol),'1',':',xlscol(outputCol),'1'));
        xlswrite(outputFile,[cell(0) nlpNode],strcat(xlscol(outputCol + 1),'1',':',xlscol(outputCol + 1),'1'));  
     
        xlswrite(outputFile,nlpNodeMat,strcat(xlscol(outputCol),'2',':',xlscol(outputCol+1),num2str(rowNum)));
     
        outputCol = outputCol + 2;
      
    else
        nlpNodeMat = zeros(length(allTimes),1); %no need to hand NaN
        for t = 1:length(allTimes)
            if tLocNLP(t) >= 0  % there's a value;
                %The NaN has been handled with the first node
                nlpNodeMat(t) = nlpVec(tLocNLP(t)); %0 or 1
            end
        end
        xlswrite(outputFile,[cell(0) nlpNode],strcat(xlscol(outputCol),'1',':',xlscol(outputCol),'1'));
        xlswrite(outputFile,nlpNodeMat,strcat(xlscol(outputCol),'2',':',xlscol(outputCol),num2str(rowNum)));
        
        outputCol = outputCol + 1;
     
    end
    
    
 nlpCol = nlpCol + 1;   
end

%now handle the financial side.
finLastColName = 'AHH';
finLastColNum = xlscol(finLastColName);

finCol = xlscol('ER');
outputCol = xlscol('ND');

while finCol <= finLastColNum
    fprintf('finCol: %d / %d\n', finCol, finLastColNum);
     %read-in the current node name
    [~,finNode] = xlsread(finFile,strcat(xlscol(finCol),'1',':',xlscol(finCol),'1'));
    finNode = finNode{1};
    finNodeNum = str2num(finNode(end)); %if it is 'N' (from NaN), then finNodeNum will be empty.
    
    %binary vector of current node
    finVec = xlsread(finFile,strcat(xlscol(finCol),'2',':',xlscol(finCol),num2str(finLastRowNum)));
    finNodeMat = zeros(length(allTimes),1);
    
    for t = 1:length(allTimes)
        if tLocFin(t) < 0 && strcmpi(finNode(end),'N') %this is an NaN node
            finNodeMat(t) = 1;
        elseif tLocFin(t) >= 0
            finNodeMat(t) = finVec(tLocFin(t));
        end
    end
    xlswrite(outputFile,[cell(0) finNode],strcat(xlscol(outputCol),'1',':',xlscol(outputCol),'1'));
    xlswrite(outputFile,finNodeMat,strcat(xlscol(outputCol),'2',':',xlscol(outputCol),num2str(rowNum)));
    
    outputCol = outputCol + 1;
    
    
    
 finCol = finCol + 1;   
end





