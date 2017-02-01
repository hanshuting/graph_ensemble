function perms = naiveSamplePerm(A, nSamples)

    [D, D2] = size(A);
    assert(D == D2, 'A must be square');
    
    % Compute the partition function
    Z = ryser_small(exp(A));
    
    % Compute everybody's permutations
    M = factorial(D);
    probs = zeros(M, 1);
    for m = 1:M
        P = expandPerm(oneperm(D, m));
        probs(m) = exp(sum(vec(A(P))));        
    end
    
    probs = probs ./ Z;
    assertElementsAlmostEqual(sum(probs), 1);
    
    % Now sample
    idxs = sample(probs, nSamples);
    perms = zeros(nSamples, D);
    for n = 1:nSamples
        perms(n,:) = oneperm(D, idxs(n));
    end
    
end

