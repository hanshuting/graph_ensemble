function features = rBiasFeat(X, population)
% rBiasFeat  Bias (all ones) feature

N = length(population);
features{1} = ones(N);
