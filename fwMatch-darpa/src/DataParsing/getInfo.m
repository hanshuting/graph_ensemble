function [data] = getInfo(csvlist)
%% 
% Author: Liang Wu
% Email:  lw2589@columbia.edu
% Date: 2-19-2015
% Read csv file from the provided path
% the function extracts all variables and stores into a structure
% including VariableName, VariableLevels, and indices
% for example, the csv file contains 10 columns, which named as
% var_level1, var_level2, ... var_level10, the function will return
% VariableName: var, VariableLevels: 10, indices: 1:10
%%

% create a struct containing all the variables
field1 = 'varName'; val1 = '';
field2 = 'nLevels'; val2 = 0;
field3 = 'indices'; val3 = zeros(1,val2);
field4 = 'fileName'; val4 = '';
s = struct(field1, val1, field2, val2, field3, val3, field4, val4);

% pre-allocate N structs
N = 1000;
data = repmat(s, N, 1);


k = 1;
count = 0;

for i = 1 : length(csvlist)
    
    % read csv to extract header
    % then convert it into a cell array
    fid = fopen(csvlist{i}, 'r');
    header = fgets(fid);
    fclose(fid);
    tmp = regexp(header,'([^,]*)','tokens');
    out = cat(2,tmp{:});
    
    % parse the file path to retrieve the file name
    str = strsplit(csvlist{i}, {'/','.'});
    filename = str{length(str)-1};

    % parse the variables from header
    % Find the delimiter from the variables
    delimiter = '';
    if any(strfind(out{2},'_'))
        delimiter = '_';
    elseif any(strfind(out{2},'-'))
        delimiter = '-';
    else
        disp('Error: Cannot find the delimiter from the variable name');
    end

    % loop through the header to obtain the meta information
    j = 1;

    while j <= length(out)    
        %count = count + 1;
        if j == 1
            % For some csv files with header starting with 'Date'
            % skip the first string
            if strcmp(out{1}, 'Date')
                j = j + 2;
                c_prev = strsplit(out{2}, delimiter);
            else
                c_prev = strsplit(out{1}, delimiter);
                j = j + 1;
            end
            data(k).varName = c_prev{1};
            data(k).nLevels = 1;
        else
            c = strsplit(out{j}, delimiter);
            if ~(strcmp(c_prev{1}, c{1}))
                k = k + 1;
                data(k).varName = c{1};
                data(k).nLevels = 1;
                c_prev = c;
            else
                data(k).nLevels = data(k).nLevels + 1;
            end
            j = j + 1;
        end
        count = count + 1;
        data(k).indices = (count - data(k).nLevels + 1) : (count);
        data(k).fileName = filename;
    end
    k = k + 1;
end
        
        
        
            
       
   
    
