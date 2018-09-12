function [ avg_test_loglike, avg_train_loglike transition_probabilities] = markov_chain_second_order_loglike(x_train, x_test)
    % transition probabilities is a cell array. In each cell we have the
    % conditional probability table for a node that can be accessed like
    % CPT(x_before_previous,x_previous,x_current)

    node_count = size(x_train, 2);
    train_examples_count = size(x_train, 1);
    
    % initialize tables for transition probabilities
    % each cell will contain a 2x2x2 table with conditional probabilities 
    % indexed by (x-2,x-1,x0).
    transition_probabilities = cell(1,node_count);

    % train markov chain
    for i = 1:node_count
        count_table = zeros(2,2,2);
        for j = 3:train_examples_count
            x_before_previous = x_train(j-2,i);
            x_previous = x_train(j-1,i);
            x_current = x_train(j,i);
            count_table(x_before_previous+1,x_previous+1,x_current+1) = count_table(x_before_previous+1,x_previous+1,x_current+1) + 1;
        end
        % laplace smoothing
        count_table = count_table + .001;
        
        sum_tables = zeros(2,2,2);
        sum_tables(:,:,1) = sum(count_table,3);
        sum_tables(:,:,2) = sum(count_table,3);
        transition_probabilities{i} = count_table ./ sum_tables;
    end
    
    % compute avg. test log likelihood
    % append last 2 samples of training to test, so there is an estimate for
    % the first example
    test_examples_count = size(x_test, 1);
    x_test = [x_train(end-1,:); x_train(end,:); x_test];
    log_like_acc = 0;
    for i = 3:(test_examples_count+2)
        for j = 1:node_count
           log_like_acc = log_like_acc + log(transition_probabilities{j}(1+x_test(i-2,j),1+x_test(i-1,j),1+x_test(i,j)));
        end
    end
    avg_test_loglike = log_like_acc / test_examples_count;
    
    % compute avg. train log likelihood
    log_like_acc = 0;
    for i = 3:train_examples_count
        for j = 1:node_count
           log_like_acc = log_like_acc + log(transition_probabilities{j}(1+x_train(i-2,j),1+x_train(i-1,j),1+x_train(i,j)));
        end
    end
    avg_train_loglike = log_like_acc / (train_examples_count-2);
end

