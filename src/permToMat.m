function M = permToMat(perm)
% permToMat  Convert permutation vector to column-permutation matrix

    D = length(perm);
    inds = sub2ind([D D], perm, 1:D);
    M = zeros(D, D);
    M(inds) = 1;
            
end

