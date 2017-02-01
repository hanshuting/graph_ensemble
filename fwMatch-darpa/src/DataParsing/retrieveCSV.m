function [listOfCSVFiles] = retrieveCSV(topLevelFolder)
%%
% Author: Liang Wu
% Email: lw2589@columbia.edu
% Date: 2-20-2015
%
% This function will loop through all the subfolders
% to search for csv files
%
% Part of the code regarding to the loop through subfolders
% is based on the link below
% http://www.mathworks.com/matlabcentral/answers/uploaded_files/9333/recurse_subfolders.m
%%


% Ask user to confirm or change.
% start_path = '~/';
% topLevelFolder = uigetdir(start_path);
% if topLevelFolder == 0
% 	return;
% end

% topLevelFolder = 'rawData/';
% Get list of all subfolders.
allSubFolders = genpath(topLevelFolder);
% Parse into a cell array.
remain = allSubFolders;
listOfFolderNames = {};
while true
	[singleSubFolder, remain] = strtok(remain, ':');
	if isempty(singleSubFolder)
		break;
	end
	listOfFolderNames = [listOfFolderNames; singleSubFolder];
end
numberOfFolders = length(listOfFolderNames);


% Process all csv files in those folders
% and store them into a cell array
listOfCSVFiles = {};
for k = 1 : numberOfFolders
	% Get this folder
	thisFolder = listOfFolderNames{k};
	
	% Get csv files.
	filePattern = sprintf('%s/*.csv', thisFolder);
	baseFileNames = dir(filePattern);
	numberOfCSVFiles = length(baseFileNames);
	% Now we have a list of all files in this folder.
	
	if numberOfCSVFiles >= 1
		% Go through all those image files.
        for f = 1 : numberOfCSVFiles
            fullFileName = fullfile(thisFolder, baseFileNames(f).name);
            listOfCSVFiles = [listOfCSVFiles; fullFileName];
        end
	end
end


