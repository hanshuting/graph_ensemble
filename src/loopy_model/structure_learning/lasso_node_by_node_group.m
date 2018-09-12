function [coefficients] = lasso_node_by_node_group(samples, relative_lambda, varargin)
% LASSO_NODE_BY_NODE_GROUP Tries to predict one variable as a function of
% another and return lasso coefficients in a matrix, where each row refers
% to a lasso logistic regression and each column a coefficient.
%
% Input
%   samples: logical matrix where each row is a sample and each column a node
%   relative_lambda: relative regularization factor (relative to the lambda
%       that drives all coefficients to zero.
%   option 'variable_groups': cell array containing (integer) indices of
%       predictor variables for each node. Variables only predicted by
%       variables indicated. Default predicts using all other variables.
%       Does not require symmetric groupings. An empty group indicates no
%       regression to be solved.
%
% Return
%   coefficients: matrix where each row refers to a lasso logistic
%   regression and each column a coefficient


    % parse input arguments
    parser = inputParser;
    is_logical_matrix = @(x) ismatrix(x);
    is_index_list = @(x) (isempty(x)) || (isvector(x) && all(x > 0));
    is_index_cells = @(x) assert(all(cellfun(is_index_list, x)), ...
        'Must be cell array, all cells index arrays.');

    parser.addRequired('samples', is_logical_matrix);
    parser.addRequired('relative_lambda', @isnumeric);
    parser.addParameter('variable_groups', [], is_index_cells);
    parser.parse(samples, relative_lambda, varargin{:});
    variable_groups = parser.Results.variable_groups;

    node_count = size(samples,2);

    % ensures variable groups is of size node_count, with valid indices
    if isempty(variable_groups)
        variable_groups = repmat(1:node_count, node_count, 1)';
        variable_groups = variable_groups(~eye(size(variable_groups)));
        variable_groups = reshape(variable_groups,node_count - 1,node_count)';
        variable_groups = num2cell(variable_groups, 2)';
    elseif length(variable_groups) ~= node_count
        error('variable_group cell array must have one (possibly empty) element per node');
    end
    valid_indices = @(x) all(x <= node_count);
    % TODO: Validate no self-reference?
    assert(all(cellfun(valid_indices, variable_groups)), ...
        'Invalid indices found in variable_groups.');


    % each line of the coefficients comes from a regression, there is one
    % regression per node, that generates one coefficient per other node
    % that is included in its group. Columns not included in a node's group
    % have coefficient zero.
    % We will store these values in the lasso_coefficients square matrix
    lasso_coefficients = zeros(node_count);
    disp('Starting Lasso Logistic Regression');

    % let's try to predict each node based on the others using lasso
    for label_node = 1:node_count
        % feature nodes are all nodes in the specified group
        feature_nodes = variable_groups{label_node};
        % no feature nodes specified ==> no regression to do
        if isempty(feature_nodes)
            continue;
        end
        X = samples(:,feature_nodes);
        Y = samples(:,label_node);

        B = lassoglm(X,Y, 'binomial', 'LambdaRatio', relative_lambda, 'NumLambda', 2);
        % B = lassoglm(X,Y, 'poisson', 'LambdaRatio', relative_lambda, 'NumLambda', 2); %SH 20160329
        lasso_coefficients(label_node, feature_nodes) = B(:,1)';

        fprintf('.');
    end
    fprintf('\n');

    coefficients = lasso_coefficients;
end
