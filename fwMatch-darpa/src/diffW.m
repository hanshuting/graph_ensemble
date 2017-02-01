function [meanDiff, maxDiff] = diffW(W, X)

    D = size(W, 1);
    a = abs(vec(logSinkhornExp(W) - logSinkhornExp(X)));
    meanDiff = mean(a);
    maxDiff  = max(a);
end

