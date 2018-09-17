% Given the edge indices, return its id in Vt
% The edges in Vt are ordered in increasing order of 
% for edges (i,j) such that i < j
% e.g. for graph with 6 nodes (1,2),(1,3), .... ,(3,6),(4,5),(4,6),(5,6)    
function id = getEdgeId(N,edge)
    E = size(edge,1);
    id = zeros(E,1);
    for e=1:E
        i = edge(e,1);
        j = edge(e,2);
        if(i >= j) 
            continue; 
        end
        id(e) = N*(i-1) - sum(1:(i-1)) + (j-i);
    end
    id(find(id == 0)) = [];
end