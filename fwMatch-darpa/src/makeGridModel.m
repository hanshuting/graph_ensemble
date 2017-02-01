function [nNodes, edges] = makeGridModel(R, C)

nNodes = R*C;

idxs = reshape(1:nNodes, R, C);

edges = {};
diffs = [ 1 0 ; 0 1 ; -1 0 ; 0 -1 ];

for r = 1:R
    for c = 1:C
        
        for nd = 1:4
            nr = r + diffs(nd,1);
            nc = c + diffs(nd,2);
                    
            if nr >= 1 && nr <= R && nc >= 1 && nc <= C
                i = idxs(r,c);
                j = idxs(nr,nc);
                
                if i < j                
                    edges{end+1} = [i ; j];
                end
            end
        end
    end
end
   
edges = cell2mat(edges);
nEdges = size(edges, 2);


end

