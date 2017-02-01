% SIMPLE_INFER_STRUCTURE Infers the structure of MRF of binary variables based on samples 
%  This is a cleaned up version of the INFER_STRUCTURE function, holding
%  only options that matter most.
%
%  Besides that, this model return graph for multiple thresholds without
%  repeating the lasso regression.
%
% Input
%   samples: logical matrix where each row is a sample and each column a node
%   relative_lambda: relative regularization factor (relative to the lambda
%       that drives all coefficients to zero.
%   thresholds: vector of thresholds to be tried.
%
% Return
%   graph_structures: cell array, each one containing a adjacency matrix to
%       the respective threshold

function [graph_structures,numeric_graph_structure] = infer_structures_with_given_densities(samples, relative_lambda, densities)
    % parse input arguments
    parser = inputParser;
    is_logical_matrix = @(x) ismatrix(x) && islogical(x);
    parser.addRequired('samples', is_logical_matrix);
    parser.addRequired('relative_lambda', @isnumeric);
    parser.addRequired('densities', @isnumeric);
    parser.parse(samples, relative_lambda, densities);
    node_count = size(samples,2);
    
    coefficients = lasso_node_by_node(samples,relative_lambda);
    
    % Remove descripancies 
    symmetric_coefficients = coefficients .* coefficients';
    
    % if there are negative values in the symmetric coefficients means that
    % one edge was predicted to be both attractive and repulsive in the two
    % different lassos. We will just zero those variables (at the moment I
    % wrote the code 10/1810 edges would have this problem)
    negative_values_indexes = find(symmetric_coefficients < 0);
    if negative_values_indexes
        fprintf('Found %d/%d edges that had contradicting weight signs in both lassos. Zero these coefficients\n', length(negative_values_indexes), length(find(symmetric_coefficients > 0)));
        symmetric_coefficients(negative_values_indexes) = 0;
    end
    
    % numeric graph structure
    numeric_graph_structure = (coefficients + coefficients');
    
    % for each given density produce an adjacency matrix
    graph_structures = cell(1,numel(densities));
    i = 0;
    for density = densities
        % now let's find what is the threshold we should use to get a certain density
        threshold = quantile(vecUT(symmetric_coefficients), 1 - density);
        if ~threshold
            fprintf('Lambda is too big, such that it is not possible to achieve the desired graph density.\n');
            fprintf('Target was %d edges, but only able to identify %d non-zero coefficients\n', ...
                nchoosek(node_count, 2) * density, length(find(vecUT(symmetric_coefficients) > 0)));
            threshold = 0;
        end
        % now identify the indexes of the edges that will be kept
        i = i + 1;
        graph_structures{i} = logical(symmetric_coefficients > threshold);
        fprintf('Report\n');
        fprintf('Total possible edges: %i\n', nchoosek(node_count, 2));
        fprintf('Density wanted: %f%%\n', density * 100);
        fprintf('Edges wanted: %i\n', density * nchoosek(node_count, 2));
        fprintf('Non-zero coefficients: %i\n', length(find(vecUT(symmetric_coefficients) > 0)));
        fprintf('Selected edges: %i\n', length(find(vecUT(symmetric_coefficients) > threshold)));
    end
end
