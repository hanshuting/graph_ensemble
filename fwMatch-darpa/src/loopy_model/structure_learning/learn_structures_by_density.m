% LEARN_STRUCTURES_BY_DENSITY
%
% Input
%   samples: logical matrix where each row is a sample and each column a node
%   relative_lambda: relative regularization factor (relative to the lambda
%       that drives all coefficients to zero.
%   densities: vector of densities to be tried.
%   option 'variable_groups': array that marks with integers which variable
%       belongs to each group, so a variable is not predicted by variables of
%       the same group
%
% Return
%   graph_structures: cell array, each one containing a adjacency matrix to
%       the respective threshold

function [graph_structures, numeric_graph_structure] = learn_structures_by_density(samples, relative_lambda, densities, varargin)
    % parse input arguments
    parser = inputParser;
%     is_logical_matrix = @(x) ismatrix(x) && islogical(x);
    is_logical_matrix = @(x) ismatrix(x);
    parser.addRequired('samples', is_logical_matrix);
    parser.addRequired('relative_lambda', @isnumeric);
    parser.addRequired('densities', @isnumeric);
%     parser.addParamValue('variable_groups', [], @isnumeric);
    parser.addParameter('variable_groups', []);     % TODO: validation
    parser.addParameter('lookback', 1, @isscalar);
    parser.parse(samples, relative_lambda, densities, varargin{:});
    node_count = size(samples,2);
    variable_groups = parser.Results.variable_groups;
    lookback = parser.Results.lookback;

    % TODO: real trigger -- currently, iscell(variable_groups) implies
    % loopback_method \in {1,2}, and constraints are imposed on the
    % possible structure
    if iscell(variable_groups)
        coefficients = lasso_node_by_node_group(samples, relative_lambda, 'variable_groups', variable_groups);
    else
        coefficients = lasso_node_by_node(samples, relative_lambda, 'variable_groups', variable_groups);
    end

    % if there are negative values in the symmetric coefficients means that
    % one edge was predicted to be both attractive and repulsive in the two
    % different lassos. We will just zero those variables (at the moment I
    % wrote the code, 10/1820 edges with the Columbia data would have this
    % problem so it is really rare)
    multiplied_coefficients = coefficients .* coefficients';
    negative_values_indexes = find(multiplied_coefficients < 0);
    if negative_values_indexes
        % TODO: This is an odd metric -- first, negative_values_indexes has
        % an entry for each half edge, so its length is twice the number of
        % undirected edges affected.
        % Second, multiplied_coefficients seems inappropriate as some max
        % set of edges to compare against, since if coefficient(i, j) > 0
        % but coefficient(j, i) == 0, edge (i, j) will not be included in
        % length(find(multiplied_coefficients ~= 0)), even though it will
        % show up as a positive edge in numeric_graph_structure and may
        % very well become an actual edge if greater than threshold.
        fprintf('Found %d/%d edges that had contradicting weight signs in both lassos. Zero these coefficients\n',...
            length(negative_values_indexes), length(find(multiplied_coefficients ~= 0)));
        coefficients(negative_values_indexes) = 0;
    end

    % numeric graph structure
    numeric_graph_structure = (coefficients + coefficients');
    % TODO: real trigger -- currently, iscell(variable_groups) implies
    % loopback_method \in {1,2}, and constraints are imposed on the
    % possible structure
    if iscell(variable_groups)
        % DEBUG if lookback_method == 2
        base_node_count = node_count / lookback;
        edge_vector = [vecUT(numeric_graph_structure(1:base_node_count, 1:base_node_count)); ...
            diag(numeric_graph_structure(1:base_node_count, (base_node_count + 1):(2 * base_node_count))); ...
            vecUT(numeric_graph_structure((base_node_count + 1):(2 * base_node_count), (base_node_count + 1):(2 * base_node_count)))];
    else
        edge_vector = vecUT(numeric_graph_structure);
    end
    max_edge_count = numel(edge_vector);

    % for each given density produce an adjacency matrix
    graph_structures = cell(1,numel(densities));
    i = 0;
    for density = densities
        % now let's find what is the threshold we should use to get a certain density
        threshold = quantile(edge_vector, 1 - density);
        if ~threshold
            fprintf('Lambda is too big, such that it is not possible to achieve the desired graph density.\n');
            fprintf('Target was %.f edges, but only able to identify %d non-zero coefficients\n', ...
                max_edge_count * density, numel(find(edge_vector > 0)));
            threshold = 0;
        end
        % now identify the indexes of the edges that will be kept
        i = i + 1;
        graph_structures{i} = logical(numeric_graph_structure > threshold);
        fprintf('Report\n');
        fprintf('Total possible edges: %i\n', max_edge_count);
        fprintf('Density wanted: %f%%\n', density * 100);
        fprintf('Edges wanted: %.f\n', density * max_edge_count);
        fprintf('Non-zero coefficients: %i\n', numel(find(edge_vector > 0)));
        fprintf('Selected edges: %i\n', numel(find(edge_vector > threshold)));
    end
end
