function [ avg_log_likelihood ] = compute_avg_log_likelihood_hidden( node_potentials, edge_potentials, logZ, samples )
%COMPUTE_LIKELIHOOD Summary of this function goes here
%   Here, p(y|x) = 2*(y)^x*(1-y)^(1-x)
%
%   Input:
%       node_potentials: column vector with potentials for all nodes
%       edge_potentials: either symmetric or upper-triangular matrix with
%           edge potentials
%       logZ: log of the partition function
%       samples: one sample per row, one column per node
%

    phi = node_potentials;
    psi = edge_potentials; 

    [sample_count, node_count] = size(samples);
    log_likelihood = 0;
    
    for t=1:sample_count
        % Add constant term (i.e. independent of X's)
        % C^t = -logZ + \sum_{i=1}^K log{(1-{y_i}^t)} + n*log2
        log_likelihood = log_likelihood + sum(log((1 - samples(t,:)') + eps)) - logZ + (node_count * log(2));
        
        % Add the logZ for this sample 
        psi_hat_t = psi; 
        phi_hat_t = log(samples(t,:)' + eps) - log((1-samples(t,:)') + eps) + phi;
        [~,~,~,~,~,~,logZ_hat_t] = run_junction_tree(phi_hat_t, psi_hat_t, 'verbose', true);
        log_likelihood = log_likelihood + logZ_hat_t;
    end
    
    % take the average
    avg_log_likelihood = log_likelihood / sample_count;
end

