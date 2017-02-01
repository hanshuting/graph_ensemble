function finHrlyBinary(fileName,rowNum,colNum,transformation)


%--------------------------------------------------------------------------
%
%
% PART 1
%
%
%--------------------------------------------------------------------------

%xlsread columns A and B
%testread.xls
%Row 1: Date, Names. Last row: 389

colNum =  xlscol(colNum); %Convert to an integer

excel_col = 1; %Numeric value of column
first_row = '1'; %Permanent for readability
second_row = '2';
last_row = num2str(rowNum); %depends on file
second_last_row = num2str(rowNum-1);
col_name = xlscol(excel_col); %Convert to letter value for input as range
last_col = colNum; %Also depends on file/number of financial indicators

outputFile = 'finHrlyBnry'; 

switch nargin
    
    case 3
         %Each ticker is applied due transformation
         outputFile = strcat(outputFile,'_dueTrans.xlsx');
    case 4

        switch transformation
                case 'original'
                    outputFile = strcat(outputFile,'_original.xlsx');
                case  'difference' 
                    outputFile = strcat(outputFile,'_difference.xlsx');
                case 'percent change' 
                    outputFile = strcat(outputFile,'_perChange.xlsx');
                otherwise
                    display('Please choose transformation from "original","difference" or "percent change"');
                    return
        end
    otherwise
        display('Wrong number of input arguments!')
        return
end



% Reading the first column: times!
temp_range = strcat(col_name,second_row,':',col_name,last_row); %Constructs range from independent strings
%Input should be Unix times
matlab_times = xlsread(fileName,temp_range); %Reads only the dates/times
%matlab_times = cell(0); 
%for k = 2:length(times) %Starting at second row
%    matlab_times=[matlab_times;datenum(times(k),'mm/dd/yyyy HH:MM:SS')]; %Converts Bloomberg times to Matlab times
%end
%Convert back to string values: date_string = datestr(matlab_times)

xlswrite(outputFile,[cell(0) 'time'],'A1:A1');

%abandon first time stamp
xlswrite(outputFile,matlab_times(2:end),strcat('A2:A',second_last_row));

excel_col=excel_col+1; %now excel_col = 2;
excel_col_out = excel_col;  %column number for output file.
ticker_names=cell(0); %Stores the ticker names



%Reading index values of other columns
%and write the binary states to output file.
while (excel_col <= last_col)
    col_name = xlscol(excel_col); %Re-writing variable each iteration
    
    %value range
    temp_range1 = strcat(col_name,second_row,':',col_name,last_row); %Also re-writing variable each iteration
    temp_range2 = strcat(col_name,first_row, ':',col_name,first_row); %ticker name range
    
   
    indicator_values = xlsread(fileName,temp_range1); %Reads indicator values
    indicator_values = [indicator_values;nan(rowNum-1-length(indicator_values),1)];
    [~,indicator_name] = xlsread(fileName,temp_range2); %Stores name of indicator as a string (character array)
    indicator_name = indicator_name{1};
    ticker_names=[ticker_names;indicator_name]; %For some reason, Matlab's having a difficult time storing strings in a cell while keeping the string type
    %numeric_values = cell(0);
    %Converts raw indicator values to numeric values
    %for k = 2:length(indicator_values)
    %    numeric_values=[numeric_values;cell2mat(indicator_values(k))]; %How do you want to handle N/A values? Eliminate cell2mat?
    %end
    
    %size(matlab_times)
    %size(indicator_values)
    function_input=[matlab_times indicator_values];
    %function_output = Liaos_function(indicator_name,function_input)
    
    %function_output  = cell(2,1); % {nodeNames, bnryStates}
    
    switch nargin
        
        case 3
            trans_list = getTransformation(indicator_name);  %The full list is {'original','difference', 'percent change'}
            %trans_list
            for k = 1:length(trans_list)
                
                switch trans_list{k}
                    case 'original'
                        [bnryStates,nodeNames] = getOriginalBnryStates(function_input,indicator_name);
                        bnryStates = bnryStates(2:end,:); %abandon first time stamp
                    case 'difference'
                        [bnryStates,nodeNames] = getDFBnryStates(function_input,indicator_name);
                    case 'percent change'    
                        [bnryStates,nodeNames] =  getDeltaBnryStates(function_input,indicator_name);         
                    otherwise
                        display(trans_list{k});
                        display('Error! No such transformation');
                        return
                end         
                
                %rename the nodes according to transformation.
                for nodeNum = 1:length(nodeNames)
                    switch trans_list{k}
                        case 'percent change'
                            nodeNames{nodeNum} = strcat('percent_change_',nodeNames{nodeNum});
                        otherwise
                            nodeNames{nodeNum} = strcat(trans_list{k},'_',nodeNames{nodeNum});
                    end
                end
                
                for j=1:size(bnryStates,2) %Right now, by default, 11 columns.    
                    col_name_out = xlscol(excel_col_out);
                     %temp_range1 = strcat(col_name_out,second_row,':',col_name_out,last_row); %range for binary states
                     temp_range2_out = strcat(col_name_out,first_row, ':',col_name_out,first_row); %range for nodeNames
        
                    %xlswrite(outputFile,bnryStates(:,j),temp_range1);
                     xlswrite(outputFile,[cell(0) nodeNames{j}],temp_range2_out); %write the node name
                    excel_col_out=excel_col_out+1; %Next node name!
                 end
     
                temp_range1_out = strcat(xlscol(excel_col_out - 11),second_row,':',xlscol(excel_col_out-1),second_last_row); %range for binary states
                 %temp_range1_out
                xlswrite(outputFile,bnryStates,temp_range1_out);
                               
             end
            
        case 4
            switch transformation
                case 'original'
                    [bnryStates,nodeNames] = getOriginalBnryStates(function_input,indicator_name);
                    bnryStates = bnryStates(2:end,:); %abandon first time stamp
                case 'difference'
                    [bnryStates,nodeNames] = getDFBnryStates(function_input,indicator_name);
                case 'percent change'    
                    [bnryStates,nodeNames] =  getDeltaBnryStates(function_input,indicator_name);         
            end
            
                 for j=1:size(bnryStates,2) %Right now, by default, 11 columns.    
                    col_name_out = xlscol(excel_col_out);
                     %temp_range1 = strcat(col_name_out,second_row,':',col_name_out,last_row); %range for binary states
                     temp_range2_out = strcat(col_name_out,first_row, ':',col_name_out,first_row); %range for nodeNames
        
                    %xlswrite(outputFile,bnryStates(:,j),temp_range1);
                     xlswrite(outputFile,[cell(0) nodeNames{j}],temp_range2_out); %write the node name
                    excel_col_out=excel_col_out+1; %Next node name!
                 end
     
                temp_range1_out = strcat(xlscol(excel_col_out - 11),second_row,':',xlscol(excel_col_out-1),second_last_row); %range for binary states
                 %temp_range1_out
                xlswrite(outputFile,bnryStates,temp_range1_out); %write the binary matrix at one time.
        otherwise
            display('Wrong number of input arguments!')
            return
    end
    
    %function_output{1} = nodeNames;
    %function_output{2} = bnryStates;
    
    %write to file

    
    %indicator_values(1:5)

    %save(indicator_name,'function_output');
    %excel_col=excel_col+1;
    excel_col = excel_col + 1;
end

%--------------------------------------------------------------------------
%
%
% PART 2
%
%
%--------------------------------------------------------------------------
% 
% excel_col=1; %Resetting variable to 1
% %last_row = '388'; %Re-initializing, depends on file
% 
% 
% 
% excel_col = excel_col+1;
% 
%         
%         
% for k=1:length(ticker_names)
%     %temp_bin_var=load(char(ticker_names(k))); %ticker_names has cells,
%     %must convert to string in order to pass parameter
%     current_ticker = ticker_names{k};
%     load(current_ticker); %{nodeNames, bnryStates}
%     
%     %{nodeNames, bnryMat} has variable name 'function_output'
%     current_nodeNames = function_output{1};
%     current_bnryMat = function_output{2};
%     %size(current_bnryMat,2)
%     for j=1:size(current_bnryMat,2) %Right now, by default, 11 columns.
%         col_name = xlscol(excel_col); %Re-writing variable each iteration
%         last_row_out = size(current_bnryMat,1) + 1;
%         last_row_out = num2str(last_row_out);
%         temp_range1 = strcat(col_name,second_row,':',col_name,last_row_out); %range for binary states
%         temp_range2 = strcat(col_name,first_row, ':',col_name,first_row); %range for nodeNames
%         
%         xlswrite(outputFile,current_bnryMat(:,j),temp_range1);
%         xlswrite(outputFile,[cell(0) current_nodeNames{j}],temp_range2);
%         excel_col=excel_col+1; %Next column!
%     end
%     
%     delete(strcat(current_ticker,'.mat'));
% end