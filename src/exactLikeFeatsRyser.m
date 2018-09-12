function logLike = exactLikeFeatsRyser(theta, features, FYInnerProduct)
% exactLikeFeats  Features, exact log-likelihood    

    [M, K] = size(features);
    
    % Hack
    finiteIdxs = ~isinf(theta);
    score = sum(FYInnerProduct(:,finiteIdxs) * theta(finiteIdxs));
    logZ  = 0;    
    
    for m = 1:M           
        logZ = logZ + exactLogZRyser(features(m,:), theta);
    end

    logLike = score - logZ;
    
end

