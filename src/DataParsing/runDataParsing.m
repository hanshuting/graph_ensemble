% Specail Note: 
% In order to get the code running, first thing to do is to 
% copy all the data into one folder, and remember the folder path.
% The csv files you need are as follows
%
% comenBnryStatesEFI.csv
% comotBnryStatesEFI.csv
% compmBnryStatesEFI.csv
% eq1BnryStatesEFI.csv
% eq2BnryStatesEFI.csv
% fi1BnryStatesEFI.csv
% vol1AbsBnryStatesEFI.csv
% nlp_binary_features.csv
%
% This code actually can access all subfolders and read all csv files
% into the sytem. However, the issue for our data is that there are
% TWO csv files(vol1BnryStatesEFI.csv and vol1AbsBnryStatesEFI.csv) 
% with SAME VARIABLES inside, which causes the merging data failure.
%
% Basically by running this collect_data.m, you will be asked to 
% locate the folder containing your data. After that, the data will
% be saved with name 'All_Data', and the meta-data will be saved as
% 'Meta_Data'.


% call retrieveCSV() to obtain the full paths of all csv files
csvlist = retrieveCSV('data/combined_data/');
% call combineTable() to merge all data into a single data table
data_table = combineTable(csvlist);


data_clean = data_table;
data_clean(1,:) = [];

data_clean = cell2mat(data_clean);
data_clean = data_clean(~any(isnan(data_clean),2),:);

% For binary samples
% data_table = cell2mat(data_table);
% variable_count = floor((size(data_table,2)-1)/10);
% sample_count = size(data_table,1);
% binary_samples = zeros(sample_count, variable_count);
% for i = 0:(variable_count-1)
%     binary_samples(:,2*i+1) = sum(data_table(:,(2+10*i):(6+10*i)),2);
%     binary_samples(:,2*i+2) = sum(data_table(:,(7+10*i):(11+10*i)),2);
% end

% currentFolder = pwd;
% desiredFolder = strrep([currentFolder '/'], '/src/DataParsing', '/cp');
% save it to .mat
% save(fullfile('cp/', 'All_Data'),'data_table');
% 
% save(fullfile('cp/', 'Binary_Sample_Data'),'binary_samples');

save('cp/All_Data.mat', 'data_table', '-v7.3');
% save('cp/Binary_Sample_Data.mat', 'binary_samples');
save('cp/Clean_Data.mat', 'data_clean','-v7.3');

% store all the meta-data into meta-info
meta_info = getInfo(csvlist);
save('cp/Meta_Data.mat', 'meta_info');
% save(fullfile('cp/', 'Meta_Data'), 'meta_info');
