function [ avg_log_likelihood ] = compute_avg_log_likelihood_overcomplete( F, G, logZ, samples, structure )
%COMPUTE_LIKELIHOOD Summary of this function goes here
%   Input:
%       node_potentials: column vector with potentials for all nodes
%       edge_potentials: either symmetric or upper-triangular matrix with
%           edge potentials
%       logZ: log of the partition function
%       samples: one sample per row, one column per node

    sample_count = size(samples,1);
    
    overcomplete_struct = samples_to_overcomplete(samples, structure);
    thetaN = (F * overcomplete_struct.Ut)';
    thetaE = (G * overcomplete_struct.Vt)';
    linearN = vec(thetaN .* overcomplete_struct.YN);
    linearE = vec(thetaE .* overcomplete_struct.YE);
    
    log_likelihood = (sum(linearN) + sum(linearE))/ sample_count;
    log_likelihood = log_likelihood - sum(F(1,:)) - sum(G(1,:));
    
    % take the average
    avg_log_likelihood = log_likelihood  - logZ;
end
