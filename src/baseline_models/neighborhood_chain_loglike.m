function [ avg_test_loglike, avg_train_loglike transition_probabilities] = neighborhood_chain_loglike(x_train, x_test, adjacency_matrix)
    % transition probabilities is a cell array. In each cell we have the
    % conditional probability table for a node that can be accessed like
    % CPT(x_previous_node_i, x_previous_node_j, ...,x_current)
    % where all the x_previous parameters corresponds to neighbors of the
    % node and the node itself, in numeric order.

    node_count = size(x_train, 2);
    train_examples_count = size(x_train, 1);
    
    % initialize tables for transition probabilities
    % each cell will contain a CPT (Conditional Probability Table)
    transition_probabilities = cell(1,node_count);

    % train neighborhood chain
    for i = 1:node_count
        neighborhood = get_ordered_neighborhood_indexes(adjacency_matrix, i);
        neighborhood_size = length(neighborhood);
        
        % create count table with d dimensions, where d is the size of the
        % neighborhood (including node itself)
        count_table = zeros(repmat(2,1,neighborhood_size+1));
        
        % for each sample increment count table
        for j = 2:train_examples_count
            x_previous = x_train(j-1,neighborhood);
            x_current = x_train(j,i);
            count_table_index = num2cell(1 + [x_previous x_current]);
            count_table(count_table_index{:}) = count_table(count_table_index{:}) + 1;
        end
        % laplace smoothing
        count_table = count_table + 1;
        
        sum_tables = zeros(repmat(2,1,neighborhood_size+1));
        indexes = repmat({':'}, 1, neighborhood_size);
        indexes{end+1} = 1;
        sum_tables(indexes{:}) = sum(count_table,neighborhood_size+1);
        indexes{end} = 2;
        sum_tables(indexes{:}) = sum(count_table,neighborhood_size+1);
        transition_probabilities{i} = count_table ./ sum_tables;
    end
    
    % compute avg. test log likelihood
    % append last sample of training to test, so there is an estimate for
    % the first example
    test_examples_count = size(x_test, 1);
    x_test = [x_train(end,:); x_test];
    log_like_acc = 0;
    for i = 2:(test_examples_count+1)
        for j = 1:node_count
           neighborhood_indexes = get_ordered_neighborhood_indexes(adjacency_matrix, j);
           neighborhood_values = num2cell(x_test(i-1,neighborhood_indexes)+1);
           neighborhood_values{end+1} = x_test(i,j) + 1;
           
           prob_table = transition_probabilities{j};
           log_like_acc = log_like_acc + log(prob_table(neighborhood_values{:}));
        end
    end
    avg_test_loglike = log_like_acc / test_examples_count;
    
    % compute avg. train log likelihood
    log_like_acc = 0;
    for i = 2:train_examples_count
        for j = 1:node_count
           neighborhood_indexes = get_ordered_neighborhood_indexes(adjacency_matrix, j);
           neighborhood_values = num2cell(x_train(i-1,neighborhood_indexes)+1);
           neighborhood_values{end+1} = x_train(i,j) + 1;
           
           prob_table = transition_probabilities{j};
           log_like_acc = log_like_acc + log(prob_table(neighborhood_values{:}));
        end
    end
    % notice we divide by train_example_count, which is 1 greater than
    % the number of examples we included in the log_prob. This is
    % because the first sample gets likelihood 100%
    avg_train_loglike = log_like_acc / train_examples_count;
end

