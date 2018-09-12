function B = normalize(A)
% normalize  Center and normalize weight matrix A to be in [0, 1].

    r = range(vec(A));
    m = min(vec(A));
    
    B = (A - m) / r;

end

