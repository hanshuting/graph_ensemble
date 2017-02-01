% LASSO_NODE_BY_NODE Tries to predict one variable as a function of another
% and return lasso coefficients in a matrix, where each row refers to a
% lasso logistic regression and each column a coefficient.
%
% Input
%   samples: logical matrix where each row is a sample and each column a node
%   relative_lambda: relative regularization factor (relative to the lambda
%       that drives all coefficients to zero.
%   option 'variable_groups': array that marks with integers which variable
%       belongs to each group, so a variable is not predicted by variables of
%       the same group
%
% Return
%   coefficients: matrix where each row refers to a lasso logistic 
%   regression and each column a coefficient

function [coefficients] = lasso_node_by_node(samples, relative_lambda, varargin)
    % parse input arguments
    parser = inputParser;
    is_logical_matrix = @(x) ismatrix(x);

    parser.addRequired('samples', is_logical_matrix);
    parser.addRequired('relative_lambda', @isnumeric);
    parser.addParamValue('variable_groups', [], @isnumeric);
    parser.parse(samples, relative_lambda, varargin{:});
    variable_groups = parser.Results.variable_groups;
    
    node_count = size(samples,2);
    
    % ensures variable groups is a vector of integers with node count size
    if isempty(variable_groups)
        variable_groups = uint8(1:node_count);
    elseif length(variable_groups) ~= node_count
        error('variable_group variable must have one element per node');
    else
        variable_groups = uint8(variable_groups);
    end
    
    % each line of the coefficients comes from a regression, there is one
    % regression per node, that generates one coefficient per other node
    % (but diagonal is zero since we cannot use own node to predict itself,
    % and for each line columns corresponding to variables in the same group
    % are zero as well).
    % We will store these values in the lasso_coefficients square matrix
    lasso_coefficients = zeros(node_count); 
    disp('Starting Lasso Logistic Regression');
    
    % let's try to predict each node based on the others using lasso
    for label_node = 1:node_count
        % feature nodes are all nodes except the ones in the same group
        % (including the node that is being labeled)
        group_nodes = find(variable_groups == variable_groups(label_node));
        feature_nodes = setdiff(1:node_count, group_nodes);
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
