function [X, variable_names] = add_lookback_nodes(base_X, time_span, base_variable_names)
%ADD_LOOKBACK_NODES Add time-shift feature node.
%   Adds (time_span - 1) sets of duplicate nodes to base_X, where the ith
%   set is the activity i timesteps later. For input data of size (T, N)
%   and time_span = K, output data of size (T, KN), where (t, i) = (t +
%   floor((i-1)/N), mod(i-1, N) + 1).
%
%   If base_variable_names is provided, updates variable_names accordingly.
% TODO Clean-up, preallocate
    X = base_X;
    node_count = size(base_X, 2);
    if nargin > 2
        variable_names = base_variable_names;
    else
        variable_names = [];
    end

    for k = 1:(time_span - 1)
        dup_block = [base_X(1+k:end, :); zeros(k, node_count)];
        X = [X dup_block];

        if nargin > 2
            make_k_name = @(n) ['(' n ',' int2str(k+1) ')'];
            new_names = cellfun(make_k_name, base_variable_names, 'UniformOutput', 0);
            variable_names = [variable_names new_names];
        end
    end
end
