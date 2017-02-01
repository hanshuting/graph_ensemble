function A = biadjacencyToAdjacency(B)

    [M, N] = size(B);
    A = [ zeros(M,M), B ;
          B'        , zeros(N,N) ];    
end

