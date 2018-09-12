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

function [graph_structures,numeric_graph_structure] = simple_infer_structure(samples, relative_lambda, thresholds)
    % parse input arguments
    parser = inputParser;
    is_logical_matrix = @(x) ismatrix(x) && islogical(x);
    parser.addRequired('samples', is_logical_matrix);
    parser.addRequired('relative_lambda', @isnumeric);
    parser.addRequired('thresholds', @isnumeric);
    parser.parse(samples, relative_lambda, thresholds);
    node_count = size(samples,2);
    
    % each line of the coefficients comes from a regression, there is one
    % regression per node, that generates one coefficient per other node
    % (but diagonal is zero since we cannot use own node to predict itself)
    % We will store these values in the lasso_coefficient square matrix
    lasso_coefficients = zeros(node_count); 
    disp('Starting Lasso Logistic Regression');
    
    % let's try to predict each node based on the others using lasso
    for label_node = 1:node_count
        feature_nodes = setdiff(1:node_count, label_node);
        X = samples(:,feature_nodes);
        Y = samples(:,label_node);
        
        B = lassoglm(X,Y, 'binomial', 'LambdaRatio', relative_lambda, 'NumLambda', 2);
        lasso_coefficients(label_node, feature_nodes) = abs(B(:,1)');
        
        fprintf('.');
    end	
    fprintf('\n');
    
    % make matrix symmetric
    numeric_graph_structure = max(lasso_coefficients, lasso_coefficients');
    
    % for each threshold produce an adjacency matrix
    graph_structures = cell(1,numel(thresholds));
    i = 0;
    for threshold = thresholds
        i = i + 1;
        graph_structures{i} = (numeric_graph_structure >= threshold);
    end
end
