function p = extractPerm(P)
% extractPerm  Extract from permutation matrix

    [p idxs] = find(P);
    assert(all(vec(idxs) == vec(1:length(idxs))));
    p = p';

end

