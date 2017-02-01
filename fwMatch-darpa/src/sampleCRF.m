function samps = sampleCRF(theta, features)
% sampleCRF  Sample a _permutation_ CRF given parameters and features.
%
%   samps = sampleCRF(theta, features) returns an MxN matrix of
%   permutations in samps, row indexing samples. theta is a K-vector while
%   features is an MxK cell array, each entry being an NxN matrix.

[M, K] = size(features);
D = size(features{1}, 1);
samps = zeros(M, D);

for m = 1:M
    A = zeros(D);
    for k = 1:K
        A = A + theta(k) * features{m,k};
    end

    samps(m,:) = sample_perms(exp(A), 1);
end

end

