function [obj, sDTheta, beliefs] = sampleLikeFeats(theta, features, FXBar, varargin)
% sampleLogZ  Features, exact log-likelihood with SAMPLED gradient.   

    p = inputParser;

    p.KeepUnmatched = true; % Enable extraneous options.
    p.addRequired('theta');
    p.addRequired('features', @iscell);    
    p.addRequired('FXBar', @isnumeric);
    
    p.addParamValue('nSamps', 100);    
    p.addParamValue('gradNormDelta', 0.01);
    p.addParamValue('debug', true);
    
    p.parse(theta, features, FXBar, varargin{:});            
    o = p.Results;
    
    % Check to get the obj value
    [obj, dTheta, beliefs] = exactLikeFeats(theta, features, FXBar);        
    [M, K] = size(features);
    D = size(features{1}, 1);
    
    % negative log likelihood        
    
    % Hack
    finiteIdxs = ~isinf(theta);
    
    % Sample dLogZ for sample m
    
    sDTheta = zeros(K, 1);    
    for m = 1:M
        A = zeros(D);
        for k = 1:K
            A = A + theta(k) * features{m,k};
        end
        
        Ys = sample_perms(exp(A), o.nSamps);        
        
        for k = 1:K
            % Gradient of the log-partition function (sample expectation)
            sDLogZm = 0;
            for s = 1:o.nSamps
                Ym = expandPerm(Ys(s,:));                
                sDLogZm = sDLogZm + sum(vec(features{m,k}(Ym)));
            end
            
            % Gradient of score + gradient of log-partition (sample mean)            
            sDTheta(k) = sDTheta(k) - FXBar(m,k) + sDLogZm / o.nSamps;
        end
        fprintf('.');
    end
    fprintf('\n');
    
    beliefs = [];
    
    if o.debug
        fprintf('1-norm diff of dTheta and sDTheta = %g\n', ...
            norm(dTheta - sDTheta, 1));
    end
    
end

