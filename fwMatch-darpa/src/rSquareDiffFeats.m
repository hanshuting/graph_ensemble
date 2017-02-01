function features = rSquareDiffFeats(X, population)
% rDiffFeats  Roommate difference features

[~, K] = size(X);

% Assume that minFeat = 1 for everything.
% NOTE: 0 denotes a MISSING VALUE (see src/roommates.ipynb. Our
% code below handles zero cases specially).

N = length(population);

features = repmat({zeros(N)}, 1, K);

for j = 1:N
    jj = population(j);
    
    for i = 1:(j-1)        
        ii = population(i);        
        
        for k = 1:K
            % If anyone is missing a value, treat the whole pair as
            % "missing" and set the feature matrix to zero.
            if X(ii,k) ~= 0 && X(jj,k) ~= 0
                % Distance feature
%                 dijk = abs(X(ii,k) - X(jj,k));
% CHANGE KT 1/21/14: Use squared differences
                dijk = (X(ii,k) - X(jj,k))^2;                
                features{k}(i,j) = dijk;
            end
        end
    end
    
    fprintf('roommateFeatures: j = %d / %d\n', j, N);
end
