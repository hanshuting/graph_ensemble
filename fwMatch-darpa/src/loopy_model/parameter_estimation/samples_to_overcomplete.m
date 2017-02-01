function [ overcomplete_params ] = samples_to_overcomplete(samples, adjacency_matrix)
%CONVERT_SAMPLES_TO_OVERCOMPLETE This functions uses the graph structure
%and the samples in its simple format (one sample per row, one variable per 
%column), and converts the samples to the overcomplete parametrization
%necessary to run the parameter estimation code
            
    sample_count = size(samples,1);
    node_count = size(samples,2);
                
    % creates column array with one element per sample, and the value is
    % always the number of nodes
    Ns = repmat(node_count, sample_count, 1);
                
    % Get the edges from the graph, such that for (row,col), row<col
    [row,col] = find(adjacency_matrix);
    selected_edges_indexes = find(row<col);
    edge_count = numel(selected_edges_indexes);
    edge_list = [row(selected_edges_indexes) col(selected_edges_indexes)];
    edge_list = sortrows(edge_list);
    
    % now create a column cellarray, with one cell per sample. In each cell we
    % will have a matrix with 2 rows and one column per edge. The edges are
    % stored with the smallest index in the first row.
    edges = repmat({edge_list'}, sample_count, 1);
                            
    % YN and YE (overcomplete parametrization of nodes and edges)
    YN = zeros(node_count*sample_count,2);
    YE = zeros(edge_count*sample_count,4);
    for i = 1:sample_count
        % generate the overcomplete parametrization. y_nodes has one row
        % per node with an indicator vector, y_edges has one row per edge
        % with an indicator vector.
        [y_nodes, y_edges] = single_sample_to_overcomplete_hidden(samples(i,:)', 2, edge_list);
        
        % append the overcomplete parametrization of the sample. There is
        % no problem to append since Ns will help separate YN into samples
        % later and the edges variable will help break YE into samples.
        base_index_YN = 1+(i-1)*node_count;
        YN(base_index_YN:(base_index_YN+node_count-1),:) = y_nodes;
        
        base_index_YE = 1+(i-1)*edge_count;
        YE(base_index_YE:(base_index_YE+edge_count-1),:) = y_edges;
    end
    
    % Ut: is an identity matrix for each sample, one concatenated with the 
    % other, horizontally
    Ut = speye(node_count);
    Ut = repmat(Ut, 1, sample_count);
    
    Vt = speye(edge_count);
    Vt = repmat(Vt, 1, sample_count);
%     % Vt: is a matrix that maps a global edge identifier to the edge index
%     % in the edges list. Vt has one row per possible edge in the graph 
%     % (n*(n-1)/2), and in each row one indicator row vector of size
%     % |edges|. This indicator vector indicates if the global edge is
%     % represented in the edge list of the graph.
%     % Finally, we concatenate one of those big matrices per sample,
%     % horizontally
%     global_edge_ids = zeros(edge_count,1);
%     for i = 1:edge_count
%         % compute global edge id. First node has n-1 slots, second node has
%         % n-2 slots and so on. Last node has 0 slots.
%         node_1 = edge_list(i,1);
%         node_2 = edge_list(i,2);
%         global_edge_id = sum((node_count-(node_1-1)):(node_count-1)) + (node_2-node_1);
%         
%         global_edge_ids(i) = global_edge_id;
%     end
%     % creates a sparse matrix for Vt
%     row_count = node_count*(node_count-1)/2;
%     col_count = edge_count;
%     Vt = sparse(global_edge_ids, (1:edge_count), ones(edge_count,1), row_count, col_count);
%     Vt = repmat(Vt, 1, sample_count);
    
    % create structure to be returned
    overcomplete_params = struct();
    overcomplete_params.Ns = Ns;
    overcomplete_params.edges = edges;
    overcomplete_params.YN = YN;
    overcomplete_params.YE = YE;
    overcomplete_params.Ut = Ut;
    overcomplete_params.Vt = Vt;
end

