function [ x_train, x_test, variable_names ] = load_columbia_data( train_split )

    % load samples and get only variables that were selected in config file
    load 'cp/columbia_data/Clean_Data.mat';
    load 'cp/columbia_data/Meta_Data.mat';
    
    % iterate through metadata identifying columns and unifying them
    variable_names = {};
    samples = [];
    for i = 1:length(meta_info)
        % check if veriables are over
        info = meta_info(i);
        if info.nLevels == 0
            break;
        end
        variable_names{end+1} = info.varName;
        samples = [samples data_clean(:,info.indices(end)+1)];
    end
    X = logical(samples);
    sample_count = size(X,1);
    
    x_train = X(1:floor(train_split*sample_count),:);
    x_test = X(floor(train_split*sample_count)+1:sample_count,:);
end

