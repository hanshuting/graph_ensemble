%% Load data
load cp/roommatesdata/data.mat;

%% Set params
maxNEdges = 20;
nTrials   = 3;

errs   = zeros(nTrials, maxNEdges);
denoms = zeros(nTrials, maxNEdges);
aucs   = zeros(nTrials, maxNEdges);

%% blah blah
K = size(features, 2);
theta = ones(K, 1);

for ne = 2:2:maxNEdges
    ne
    [errs(:,ne), denoms(:,ne), aucs(:,ne)] = evalThetaUnipartiteSubsample(...
        theta, features(3,:), Ys{3}, ne, nTrials);
end


normalizedErrs = errs ./ denoms;
figure;
boxplot(normalizedErrs); title('Normalized Hamming Errors');
