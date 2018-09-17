% INFER_STRUCTURE Infers the structure of MRF of binary variables based on samples 
%
% Input
%   samples: logical matrix where each row is a sample and each column a node
%   lambda: regularization factor
%   options:
%       'graph_build_method': ['and', 'or', 'average', 'raw']. How to conciliate
%           the fact that sometimes node A is used to predict B, but B is not
%           used to predict A. 'and' methods need both to happen to put an
%           edge between nodes, 'or' only needs one of them, and 'average'
%           checks the lasso coefficient on average is above tolerance.
%           'raw' returns the assymetric matrix which rows are just the
%           output of the lasso regression.
%           Default: 'and'
%       'return_logical': [true, false]. Whether to return a thresholded
%       logical matrix representing the connection or to return the numeric
%       weights predicted for each edge. Default: true
%       'tolerance': Only used if 'return_logical' is true. Threshold for 
%           the lasso coefficients, above which a node is considered 
%           relevant to predict another (i.e. connected). Default: 0.1
%       'use_built_in_lasso': [true, false]. If true uses matlab lasso
%       regression, else uses Boyd binaries. Default: false
%       'use_relative_lambda': interprets lambda not as an absolute value
%       but as a proportions of the maximum lambda. Only works when 
%       'use_built_in_lasso' is set. Default: false.
%       'use_lasso_regression': does regression instead of logistic
%       regression. Option kept so we can make tests. Default: false.


function [ graph_structure,graph_values] = infer_structure( samples, lambda, varargin )
    % parse input arguments
    parser = inputParser;
    is_logical_matrix = @(x) ismatrix(x) && islogical(x);
    parser.addRequired('samples', is_logical_matrix);
    parser.addRequired('lambda', @isnumeric);
    parser.addParamValue('tolerance', 0, @isnumeric);
    parser.addParamValue('graph_build_method', 'and', @ischar);
    parser.addParamValue('return_logical', true, @islogical);
    parser.addParamValue('use_built_in_lasso', false, @islogical);
    parser.addParamValue('use_relative_lambda', false, @islogical);
    parser.addParamValue('use_lasso_regression', false, @islogical);
    parser.addParamValue('preknown_edges',-1,@isscalar);
    parser.parse(samples, lambda, varargin{:});
    
    tolerance = parser.Results.tolerance;
    graph_build_method = parser.Results.graph_build_method;
    return_logical = parser.Results.return_logical;
    use_built_in_lasso = parser.Results.use_built_in_lasso;
    use_relative_lambda = parser.Results.use_relative_lambda;
    use_lasso_regression = parser.Results.use_lasso_regression;
    preknown_edges = parser.Results.preknown_edges;
    
    node_count = size(samples,2);
    
    % each line of the coefficients comes from a regression, there is one
    % regression per node, that generates one coefficient per other node
    % (but diagonal is zero since we cannot use own node to predict itself)
    lasso_coefficients = zeros(node_count); 

    disp('Starting Lasso Logistic Regression');
    
    % let's try to predict each node based on the others using lasso
    for label_node = 1:node_count
        feature_nodes = setdiff(1:node_count,[label_node]);
        X = samples(:,feature_nodes);
        Y = samples(:,label_node);
        
        if use_built_in_lasso
            if use_lasso_regression
                if use_relative_lambda
                    B = lasso(X,Y,'LambdaRatio', lambda, 'NumLambda', 2);
                else
                    B = lasso(X,Y,'Lambda', lambda);
                end
            else
                if use_relative_lambda
                    B = lassoglm(X,Y, 'binomial', 'LambdaRatio', lambda, 'NumLambda', 2);
                else
                    B = lassoglm(X,Y, 'binomial', 'Lambda', lambda);
                end 
            end
            lasso_coefficients(label_node, feature_nodes) = abs(B(:,1)');
        else
            % lasso solver accepts inputs as files, using temp files
            features_file = tempname;
            mmwrite(features_file, samples(:,feature_nodes));
            labels_file = tempname;
            mmwrite(labels_file, samples(:,label_node));

            % solver outputs file, getting temp path for output
            trained_model_file = tempname;

            % call lasso solver
            system(sprintf('thirdparty/l1_logreg-0.8.2/src_c/l1_logreg_train -q -s %s %s %.32f %s', features_file, labels_file, lambda, trained_model_file));

            % get classifier coefficients (except intercept term)
            trained_model = full(mmread(trained_model_file));
            lasso_coefficients(label_node, feature_nodes) = abs(trained_model(2:end));
        end
        
        fprintf('.');
    end	
    
    disp('');

    % now use coefficients to infer graph structure    
    switch graph_build_method
    case 'and'
        graph_structure = min(lasso_coefficients, lasso_coefficients');
    case 'or'
        graph_structure = max(lasso_coefficients, lasso_coefficients');
    case 'average'
        graph_structure = (lasso_coefficients + lasso_coefficients')/2;
    case 'raw'
        graph_structure = lasso_coefficients;
    end
    % if output should be logical, then compare with threshold
    if return_logical
        disp('Returning thresholded estimated structure')
        graph_values = graph_structure;
        graph_structure = triu(graph_structure,1);
        if preknown_edges ~= -1
            [~,edgeInd] = sort(graph_structure(:),'descend');
            edges = zeros(size(graph_structure(:)));
            edges(edgeInd(1:preknown_edges)) = 1;
            graph_structure = reshape(edges,size(graph_structure));
            graph_structure = max(graph_structure,graph_structure');
        else
            graph_structure = graph_structure >= tolerance;
        end
    else
        graph_values = graph_structure;
        disp('Returning non-thresholded estimated structure')
    end
    
    % Make sure diagonal values are 0
    graph_structure = graph_structure .* (1 - diag(ones(1,node_count)));
    
end
