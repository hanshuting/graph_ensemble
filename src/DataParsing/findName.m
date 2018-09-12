function [ varName, fileName, level ] = findName( index, meta_info )
%%
% Author: Liang WU
% Email: lw2589@columbia.edu
% Date: Ap3. 27, 2015
%
%   This function takes the index of the variable and the meta_info data
%   to find out the name of the variable, as well as which file the
%   variable is from.
%
%%

indArr = {meta_info.indices};

Loc = cellfun(@(x) any(x(:)==index),indArr);

ind_st = find(Loc);

if (isempty(ind_st))
    fprintf('Index Not Found! Please try another index\n');
    return;
else
    varName = meta_info(ind_st).varName;
    fileName = meta_info(ind_st).fileName;
    accLevels = 0;
    for i = 1 : ind_st - 1
        accLevels = accLevels + meta_info(i).nLevels;
    end
    level = index - accLevels;
    if (strfind(lower(fileName), 'fin'))
        fprintf('The searched index belongs to Financial Variables\n');
    elseif (strfind(lower(fileName), 'nlp'))
        fprintf('The searched index belongs to NLP Variables\n');
    else
        fprintf('Error!\n');
        return;
    end
end     

end

