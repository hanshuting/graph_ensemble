function [ avg_log_likelihood ] = compute_avg_log_likelihood( node_potentials, edge_potentials, logZ, samples )
%COMPUTE_LIKELIHOOD Summary of this function goes here
%   Input:
%       node_potentials: column vector with potentials for all nodes
%       edge_potentials: either symmetric or upper-triangular matrix with
%           edge potentials
%       logZ: log of the partition function
%       samples: one sample per row, one column per node

    log_likelihood = 0;
    sample_count = size(samples,1);
    
    % add node effects
    log_likelihood = log_likelihood + sum(samples * node_potentials);

    % guarantee edge_potentials has one only in upper triangle
    edge_potentials = edge_potentials - tril(edge_potentials);
    
    % for each sample, add edge effects
    for i = 1:sample_count
        sample = samples(i,:);
        log_likelihood = log_likelihood + sum(sum(edge_potentials .* (double(sample)' * double(sample))));
    end
    
    % subtract partition function
    log_likelihood = log_likelihood - sample_count * logZ;
    
    % take the average
    avg_log_likelihood = log_likelihood / sample_count;
end

