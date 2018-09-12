function [obj, dTheta, beliefs] = exactLikeFeats(theta, features, FXBar, varargin)
% exactLikeFeats  Features, exact log-likelihood    

    [M, K] = size(features);
    
    % negative log likelihood
    obj    = 0;
    dTheta = zeros(K, 1);
    
    % Hack
    finiteIdxs = ~isinf(theta);

    %% HACK 2/24 -- USING THIS CODE TO COMPUTE MRFS, ONLY COMPUTE ONE PARTITION FUNCTION.
    % (All the features will be the same here.)
    % This will NOT WORK for true CRFs.
    [logZm, dLogZm] = exactLogZRyser(features(1,:), theta);    
    
    for m = 1:M
        % Batch gradient.        
%         [logZm, dLogZm] = exactLogZRyser(features(m,:), theta);
        
        obj = obj - FXBar(m,finiteIdxs) * theta(finiteIdxs) + logZm;
        
%         obj = obj - FXBar(m,:) * theta + logZm; 
        dTheta = dTheta - FXBar(m,:)' + dLogZm;
    end

    beliefs = [];
    
end

