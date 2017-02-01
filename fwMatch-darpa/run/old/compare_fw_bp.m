%% blaha
D = 5;
B = randn(D);


%% run
[muMin, FBetheHist] = fwBipartite_mex(1/D * ones(D), B);
FBethe = FBetheHist(end);

[muMinFmincon, FBetheFmincon] = minBetheFmincon(1/D * ones(D), B, 1e-6);

A = biadjacencyToAdjacency(B);
[logZ, beliefs, msgs] = bpPerfectMatching(A);
FBetheBP = -logZ;

FBethe
FBetheBP
FBetheFmincon

muMinFmincon - muMin
max(abs(vec(muMin - muMinFmincon)))

% muMin - adjacencyToBiadjacency(beliefs(:,:,2))
