function [yNode, yEdge] = overcompleteLabels(labels, L, edges)
% overcompleteLabels  Convert node labels to overcomplete indicators
%
%   y = overcompleteLabels(labels, L, edges)
%
%   y      : overcomplete binary vector
%   labels : N vector with categorical labels
%   edges  : 2xE matrix of edges

    assert(isvector(labels), 'labels must be vector; did you pass in matrix?');

    N = length(labels);
    E = size(edges, 2);
    
    yNode = zeros(N, L);
    for x = 1:L
        yNode(labels == x, x) = 1;
    end
            
    yEdge = zeros(E, L^2);
    inds = reshape(1:L*L, L, L);
    for e = 1:E
        xi = edges(1,e);
        xj = edges(2,e);
        
        yi = labels(xi);
        yj = labels(xj);
        
        yEdge(e, inds(yi,yj)) = 1;
    end
    
    % Check our work
    assert(all(sum(yNode, 2) == 1), 'yNode not correctly labeled.');
    assert(all(sum(yEdge, 2) == 1), 'yEdge not correctly labeled.');
    
end
