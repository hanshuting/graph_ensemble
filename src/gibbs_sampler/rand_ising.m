% RAND_ISING Creates a random uniform undirected graph in the Ising model.
%
% The graph has random node and edge potentials in ]-1,1[. The graph 
% structure is composed of a spanning tree generated with Wilson's 
% algorithm [1] and random extra edges.
%
% Usage
%   [node_potentials edge_potentials] = randIsing(node_count, extra_edge_density) 
%
% Parameters
%   node_count:          # of nodes in graph
%   extra_edge_density:  in [0,1] represents the probability of adding an 
%                        edge that is not in the spanning tree (optional).
%                        Default is 0.
%
% Options
%   'minimum_edge_weight': this is to avoid generating a graph with edge
%                          weights too low, because then the
%                          relation between those nodes is so weak that
%                          cannot be learned. We are not interested in those
%                          edges in practice. Default: 0.2.
%
% Return values
%   node_potentials: column vector with one potential in ]-1,1[ per node
%   edge_potentials: symmetric sparse matrix with edge weights in ]-1,1[. 
%
% References
% [1] David Bruce Wilson, "Generating Random Spanning Trees More Quickly
%     than the Cover Time" STOC '96 doi>10.1145/237814.237880

function [ node_potentials, edge_potentials ] = rand_ising( node_count, varargin )
    % parsing parameters
    parser = inputParser;
    is_positive = @(x) isnumeric(x) && (x > 0);
    is_between_zero_and_one = @(x) isnumeric(x) && (x >= 0) && (x <= 1);
    parser.addRequired('node_count', is_positive);
    parser.addOptional('extra_edge_density', 0, is_between_zero_and_one);
    parser.addParamValue('minimum_edge_weight', 0.0, is_between_zero_and_one);
    parser.parse(node_count, varargin{:});
    minimum_edge_weight = parser.Results.minimum_edge_weight;
    extra_edge_density = parser.Results.extra_edge_density;

    nodes_in_tree = false(node_count, 1);
    parents = zeros(node_count, 1);
    
    % pick a random root for the spanning tree and add it to the tree
    root = randsample(node_count, 1);
    nodes_in_tree(root) = true;
    
    for i = 1:node_count
        % start from node i and randomly choose a parent for it. If parent
        % is not in the tree, do that recursively until you hit the tree.
        % Note that if this procedure generate a cycle before getting to 
        % the tree, the parent is overwritten and the cycle is broken.
        u = i;
        while ~nodes_in_tree(u)
            % generate array excluding u and sample from it
            candidate_parents = 1:node_count;
            candidate_parents(u) = [];
            parents(u) = randsample(candidate_parents,1);
            u = parents(u);
        end
        
        % go through the path randomly walked in the loop above and add
        % nodes to tree
        u = i;
        while ~nodes_in_tree(u)
            nodes_in_tree(u) = true;
            u = parents(u);
        end
    end
    
    % Now we want to build the conectivity matrix C with random weights. We
    % start by finding the pairs to be connected
    nodes_with_parents = find(parents);
    respective_parents = parents(nodes_with_parents);
    
    % We sort each two paired nodes so that the first node in the pair is 
    % always greater than the second. This makes the conectivity matrix
    % have all its values in the lower-left triangle
    aux = sort([nodes_with_parents respective_parents], 2, 'descend');
    edge_i = aux(:,1); 
    edge_j = aux(:,2);
    
    % Build the matrix with random edge weights 
    edge_weights = (rand(node_count-1, 1) * (1 - minimum_edge_weight) + minimum_edge_weight) .* (2*randi(2, node_count-1, 1) - 3);
    spanning_tree = sparse(edge_i, edge_j, edge_weights, node_count, node_count);    
    
    % Add extra random edges: generate random [-1,1] for all edges not in
    % the spanning tree. Choose only the edges below the density threshold
    % and then scale the remaining range so it goes from 0
    % to 1.
    if extra_edge_density > 0
        extra_edges_prob = sprand(~spanning_tree & tril(ones(node_count,node_count),-1)) / extra_edge_density;
        extra_edges_mask = (extra_edges_prob ~= 0 & extra_edges_prob <= 1) .* ((rand(node_count,node_count) > .5)*2 - 1);
        
        % now compute the values of the extra edges
        extra_edges = (rand(node_count, node_count) * (1 - minimum_edge_weight) + minimum_edge_weight) .* extra_edges_mask;
        
        edge_potentials = spanning_tree + extra_edges;
    else
        edge_potentials = spanning_tree;
    end
    
    % make edge_potentials a symmetric matrix
    edge_potentials = edge_potentials + edge_potentials';
    
    % generate random node weights
    node_potentials = 2*rand(node_count, 1) - 1;
end