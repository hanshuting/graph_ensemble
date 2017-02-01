function [m, v] = randomMatchingHammingErr(n)
% randomMatchingError  Compute mean and variance of nNormalized Hamming error of random
%                      matching
%
%   [m, v] = randomMatchingHammingErr(n) returns mean in m and variance in
%   v of the normalized Hamming error of a random matching of n nodes.

% Derivation:
%
% Let X denote the number of incorrect edges in a random matching. First,
% P(X = 0) = 1/n!!, since only one matching out of n!! possibilities will
% yield all of the correct edges. But P(X <= k) = 1/(n-k)!!, since we only
% need (n - k) edges to be matched correctly. Since n and k must both be
% even, we conclude for k > 0,
%
%       P(X = k) = P(X <= k) - P(X <= k - 2)
%                = 1/(n-k)!! - 1/(n-(k-2))!!
%
% Since the terms telescope and k = 0 is the maximum, we see that this
% correctly normalizes.
%
% From this probability mass function, expectation and variance are
% straightforward.

ks  = 2:2:n;
pmf = [1/doubleFact(n), 1 ./ doubleFact(n - ks) - 1 ./ doubleFact(n - (ks - 2))];

assert(sum(pmf) == 1);

% This is the correct normalization. Hamming error should increase by 1
% whenever an edge is predicted wrong. Here, we use the equivalent
% characterization of two nodes being wrong. This double-counts our error,
% so where we would ordinarily divide by n/2 (number of edges), we can
% instead divide by n, and get the same result in the end.

m  = 1/n   * sum((0:2:n)    .* pmf);
m2 = 1/n^2 * sum((0:2:n).^2 .* pmf); 
v  = m2 - m^2;

end

function y = doubleFact(n)
    k = n/2;
    y = 2.^k .* factorial(k);
end
