function [x_train, x_test, variable_names] = add_lookback_nodes(base_x_train, base_x_test, time_span, base_variable_names)
%ADD_LOOKBACK_NODES Add time-shift feature node.
%   Adds (time_span - 1) sets of duplicate nodes to x_train and x_test,
%   where the ith set is the activity i timesteps later. For input data of
%   size (T, N) and time_span = K, output data of size (T, KN), where
%   (t, i) = (t + floor((i-1)/N), mod(i-1, N) + 1).
%
%   If base_variable_names is provided, updates variable_names accordingly.
% TODO Clean-up, preallocate
    x_test = base_x_test;
    x_train = base_x_train;
    node_count = size(base_x_train, 2);
    if nargin > 3
        variable_names = base_variable_names;
    else
        variable_names = [];
    end

    for k = 1:(time_span - 1)
        test_block = [base_x_test(1+k:end, :); zeros(k, node_count)];
        x_test = [x_test test_block];
        train_block = [base_x_train(1+k:end, :); zeros(k, node_count)];
        x_train = [x_train train_block];

        if nargin > 3
            make_k_name = @(n) ['(' n ',' int2str(k+1) ')'];
            new_names = cellfun(make_k_name, base_variable_names, 'UniformOutput', 0);
            variable_names = [variable_names new_names];
        end
    end
end
