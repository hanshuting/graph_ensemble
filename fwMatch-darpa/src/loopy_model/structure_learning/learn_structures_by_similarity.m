function [graph_structures, numeric_graph_structure] = learn_structures_by_similarity(samples, densities, variable_groups)
%LEARN_STRUCTURES_BY_SIMILARITY Find edges by similarity metric.
%
% Input
%   samples: logical matrix where each row is a sample and each column a node
%   densities: vector of densities to be tried. Density is taken as a
%       percentage of eligible edges, as defined by variable_groups.
%   variable_groups: cell array containing (integer) indices of nodes eligible
%       to form other side of an edge for each node. Assumes undirected edges,
%       and so symmetric groupings. An empty group indicates no edges should be
%       formed with that node.
%
% Return
%   graph_structures: cell array, each one containing a adjacency matrix to
%       the respective density threshold
%   numeric_graph_structure: Raw similarity scores.
    parser = inputParser;
    is_logical_matrix = @(x) ismatrix(x);
    parser.addRequired('samples', is_logical_matrix);
    parser.addRequired('densities', @isnumeric);
    variable_groups_chk = @(x)validateattributes(x, {'cell'}, {'vector'});
    parser.addOptional('variable_groups', [], variable_groups_chk);
    if nargin < 3
        parser.parse(samples, densities);
        variable_groups = parser.Results.variable_groups;
    else
        parser.parse(samples, densities, variable_groups);
    end
    
    node_count = size(samples,2);
    if isempty(variable_groups)
        % Create default variable_groups (fully connected) if none provided.
        variable_groups = cell(1, node_count);
        for ii = 1:node_count
            variable_groups{ii} = [1:(ii-1) (ii+1):node_count];
        end
    end

    if length(variable_groups) ~= node_count
        error('variable_group cell array must have one (possibly empty) element per node');
    end

    % Convert variable_groups to logical matrix
    eligible_edges = zeros(node_count, node_count);
    for i = 1:node_count
        eligible_edges(i, variable_groups{i}) = 1;
    end
    eligible_edges = logical(eligible_edges);
    assert(all(all(eligible_edges == eligible_edges')), ...
        'variable_groups must be symmetric.');

    similarity = get_correlations(samples);
    coefficients = zeros(size(similarity));
    % Restrict to eligible_edges
    coefficients(eligible_edges) = similarity(eligible_edges);

    % No conflicting edge check necessary: correlations perfectly symmetric

    % numeric graph structure
    numeric_graph_structure = coefficients;

    % Vectorize eligible edge similarity scores
    edge_vector = numeric_graph_structure(triu(eligible_edges, 1));

    max_edge_count = numel(edge_vector);
    num_pos_coefficients = numel(find(edge_vector > 0));

    % for each given density produce an adjacency matrix
    graph_structures = cell(1,numel(densities));
    i = 0;
    for density = densities
        i = i + 1;
        % now let's find what is the threshold we should use to get a certain density
        density_threshold = quantile(edge_vector, 1 - density);
        % Ensure theshold is always positive
        threshold = max(density_threshold, eps);
        % Select inclusive of threshold
        graph_structures{i} = logical(numeric_graph_structure >= threshold);
        edges_found = numel(find(graph_structures{i})) / 2;

        edges_wanted = density * max_edge_count;
        if edges_wanted > edges_found
            fprintf('It is not possible to achieve the desired graph density of %.f%%.\n', density * 100);
            fprintf('Target was %.f edges, but only able to identify %d positive coefficients\n', ...
                edges_wanted, edges_found);
        end
        fprintf('Report\n');
        fprintf('Total possible edges: %i\n', max_edge_count);
        fprintf('Density wanted: %f%%\n', density * 100);
        fprintf('Edges wanted: %.f\n', edges_wanted);
        fprintf('Positive coefficients: %i\n', num_pos_coefficients);
        % Div by 2 to compensate for undirected edge showing up twice in adj matrix
        fprintf('Selected edges: %i\n', edges_found);
    end
end

function [correlations] = get_correlations(samples)
% GET_CORRELATIONS Find correlation between all node poirs.
    correlations = corrcoef(samples);
    % All-zero nodes will produce NaN correlations
    correlations(isnan(correlations)) = 0;
end
