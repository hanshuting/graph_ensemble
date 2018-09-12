function [W, hist] = learnBetheFakeFeatures(X, varargin)
% learnBetheFakeFeatures

[M, D] = size(X);

K = D * D;
features = cell(M, K);

for m = 1:M
    for k = 1:K
        [i, j] = ind2sub([D D], k);
        features{m,k} = zeros(D);
        features{m,k}(i,j) = 1;
    end
end

[theta, hist] = learnFeats(X, features, @betheLikeFeats, varargin{:});

W = reshape(theta, D, D);


end

