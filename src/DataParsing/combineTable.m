function [data] = combineTable(csvlist)
%%
% Author: Liang Wu
% Email: lw2589@columbia.edu
% Date: 2-20-2015
%
% This function takes all csv files that are provided
% as a list of full file paths, and then merge them
% into a single cell array.
% 
% It serves as part of the codes for data 
% preprocessing.
%
%%


ftable = cell(1, length(csvlist));

for i = 1 : length(csvlist)
    ftable{i} = dataset('File',csvlist{i},'Delimiter',',');
    ftable{i}.Properties.VarNames{1} = 'Date';
    ftable{i}.Date = datenum(ftable{i}.Date);

%      ftable{i} = readtable(csvlist{i});
%      ftable{i}.Properties.VariableNames{1} = 'Date';
%      ftable{i}.Date = datenum(ftable{i}.Date);
%      ftable{i}.Properties.ObsNames = ftable{i}.Date;
end

data_table = ftable{1};

% % For binary samples
% data_table = cell2mat(data_table);
% variable_count = floor((size(data_table,2)-1)/10);
% sample_count = size(data_table,1);
% binary_samples = zeros(sample_count, variable_count);
% for i = 0:(variable_count-1)
%     binary_samples(:,i+1) = sum(data_table(:,(6+10*i):(10+10*i)),2);
% end


for i = 2 : length(ftable)
    data_table = join(ftable{i},data_table,'type', 'outer','keys', 'Date', 'mergeKeys',true);
end

% convert datenumber back to date
% formatOut = 'mm-dd-yyyy';
% data_table.Date = datestr(data_table.Date, formatOut);
% formatIn = 'MM-dd-yyyy';
% data_table.Date = datetime(data_table.Date, 'InputFormat',formatIn);
data = dataset2cell(data_table);
