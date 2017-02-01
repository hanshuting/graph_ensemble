function features = roommateFeatures(surveyFile, population)
% allRoommateFeatures  Return pairwise ctors of *all* pairs of students

X = importdata(sprintf('data/%s', surveyFile));
[~, K] = size(X);

% Assume that minFeat = 1 for everything.
% NOTE: 0 denotes a MISSING VALUE (see src/roommates.ipynb. Our
% code below handles zero cases specially).
maxFeat = max(X, [], 1);
medFeat = (1 + maxFeat) / 2;

% Precompute index map
interactionIds = {};
ind = K + 1;
for k = 1:K
    interactionIds{k} = zeros(maxFeat(k));
    
    for vj = 1:maxFeat(k)
        % KT 2 Feb 2015 -- This should really be 1:maxFeat. Symmetry
        % happens between two uses, NOT between the levels of two distinct
        % features.
        for vi = vj:maxFeat(k)
            interactionIds{k}(vi,vj) = ind;
            interactionIds{k}(vj,vi) = ind;
            ind = ind + 1;
        end
    end
end

nTotFeats = ind - 1; % one for distance, one for interactions.

N = length(population);

% Prefill all ctors with zero matrices. In the sequel, we first test for
% nonzero entries (existence of data) before filling in a matrix. Thus, any
% entry with missing data stays zero, and is automatically ignored by the
% algorithm. Neat!

% ctor  = SparseCtor('rows', N, 'cols', N);
% ctors = repmat({ctor}, 1, nTotFeats);

% The loop below construct two kinds of feature matrices:

% Distance ctors (absolute value of difference).
%
% To think about: normalization.
%
% Recall that the questions are answered with sliders. The
% middle value of the slider usually denotes "neutral." But I'm
% not immediately sure how centering, and thus using the
% directionality information, can help. We can't just use the
% raw differences, because then the feature matrix would not be
% symmetric.

% Interaction ctors
%
% Each unordered pair of feature responses gets an entry in
% interactionIds. Given the settings of X(ii,k) and X(jj,k), we
% set the appropriate value to 1.

% JUST WORK WITH DENSE MATRICES AND WASTE MEMORY.
features = repmat({zeros(N)}, 1, nTotFeats);

for j = 1:N
    jj = population(j);
    
    for i = 1:(j-1)        
        ii = population(i);        
        
        for k = 1:K
            % If anyone is missing a value, treat the whole pair as
            % "missing" and set the feature matrix to zero.
            if X(ii,k) ~= 0 && X(jj,k) ~= 0
                % Distance feature
                dijk = abs(X(ii,k) - X(jj,k));
                features{k}(i,j) = dijk;
                
                % Interaction feature: Look up my index
                ik = interactionIds{k}(X(ii,k), X(jj,k));
                features{ik}(i,j) = 1;
            end
        end
    end
    
    fprintf('roommateFeatures: j = %d / %d\n', j, N);
end

% Sparsify the interaction features
for k = K+1:nTotFeats
    features{k} = sparse(features{k});
end

% Add the bias term to the end
features{nTotFeats + 1} = ones(N);

% features = cellfun(@(c) c.makeSparse(), ctors, 'UniformOutput', false);
