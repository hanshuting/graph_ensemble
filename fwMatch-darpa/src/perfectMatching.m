function [M, cost] = perfectMatching(A)

    [N, M] = size(A);
    assert(N == M, 'A must be a symmetric weighted adjacency matrix');
    [ri, rj, rw] = findUT(A);
    
    % TODO: possible speed improvements here?
    mVec = blossom_double_mex(N, ri, rj, rw);
    M = sparse(ri, rj, mVec, N, N);
    M = M + M';    
    
    cost = trace(A' * M);

end

