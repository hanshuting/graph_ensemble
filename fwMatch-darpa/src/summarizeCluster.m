function [edges, withins] = summarizeCluster(A, map)

map = map - min(map) + 1; % for 1-based indexing

for k = 1:max(map);
    mask  = map == k;
    idxs  = find(mask);
    
    nbrs  = zeros(size(idxs));
    
    for ii = 1:length(idxs)
        i = idxs(ii);
        nbrs(ii) = find(A(:,i), 1);
    end
    
    edges{k} = [idxs; nbrs];
    
    withins(k) = length(intersect(idxs, nbrs)) / length(idxs);    
end

end