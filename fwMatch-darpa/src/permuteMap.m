function p = permuteMap(map)

map = map - min(map) + 1;
K = max(map);

p = zeros(size(map));
currIdx = 1;
for k = 1:K    
    mask     = map == k;
    nk       = sum(mask);
    pIdxs    = currIdx:currIdx+nk-1;
    p(pIdxs) = find(mask);
    currIdx = currIdx+nk;
end

end