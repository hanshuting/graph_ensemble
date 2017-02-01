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
    parser.addParamValue('variable_groups', [], @isnumeric);
    parser.parse(samples, relative_lambda, densities, varargin{:});
    node_count = size(samples,2);
    variable_groups = parser.Results.variable_groups;
    
    
    % ensures variable groups is a vector of integers with node count size
    if isempty(variable_groups)
        variable_groups = uint8(1:node_count);
    elseif length(variable_groups) ~= node_count
        error('variable_group variable must have one element per node');
    else
        variable_groups = uint8(variable_groups);
    end
    
    coefficients = lasso_node_by_node(samples, relative_lambda, 'variable_groups', variable_groups);
    
    % if there are negative values in the symmetric coefficients means that
    % one edge was predicted to be both attractive and repulsive in the two
    % different lassos. We will just zero those variables (at the moment I
    % wrote the code, 10/1820 edges with the Columbia data would have this 
    % problem so it is really rare)
    multiplied_coefficients = coefficients .* coefficients';
    negative_values_indexes = find(multiplied_coefficients < 0);
    if negative_values_indexes
        fprintf('Found %d/%d edges that had contradicting weight signs in both lassos. Zero these coefficients\n',...
            length(negative_values_indexes), length(find(multiplied_coefficients ~= 0)));
        coefficients(negative_values_indexes) = 0;
    end
    
    % numeric graph structure
    numeric_graph_structure = (coefficients + coefficients');
    
    % for each given density produce an adjacency matrix
    graph_structures = cell(1,numel(densities));
    i = 0;
    for density = densities
        % now let's find what is the threshold we should use to get a certain density
        threshold = quantile(vecUT(numeric_graph_structure), 1 - density);
        if ~threshold
            fprintf('Lambda is too big, such that it is not possible to achieve the desired graph density.\n');
            fprintf('Target was %d edges, but only able to identify %d non-zero coefficients\n', ...
                nchoosek(node_count, 2) * density, length(find(vecUT(numeric_graph_structure) > 0)));
            threshold = 0;
        end
        % now identify the indexes of the edges that will be kept
        i = i + 1;
        graph_structures{i} = logical(numeric_graph_structure > threshold);
        fprintf('Report\n');
        fprintf('Total possible edges: %i\n', nchoosek(node_count, 2));
        fprintf('Density wanted: %f%%\n', density * 100);
        fprintf('Edges wanted: %i\n', density * nchoosek(node_count, 2));
        fprintf('Non-zero coefficients: %i\n', length(find(vecUT(numeric_graph_structure) > 0)));
        fprintf('Selected edges: %i\n', length(find(vecUT(numeric_graph_structure) > threshold)));
    end
end
