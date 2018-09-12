function [logLikes, exactLogLikes] = permTestLikelihood(Ys, theta, features, rho, lambda, varargin)
    
    % Read stuff...
    p = inputParser;
    p.KeepUnmatched = true; % Enable extraneous options.
    
    p.addRequired('Ys');
    p.addRequired('theta');
    p.addRequired('features');
    p.addRequired('rho');
    
    p.addParamValue('objParams', {});
    p.addParamValue('algo', @BatchFW);
    p.addParamValue('algoParams', {});
    
    p.parse(Ys, theta, features, rho, varargin{:});
    opts = p.Results;
    
    % And compute stuff.
    M = length(Ys);
    D = length(Ys{1});
    
    % First compute the log partition function
    feats = makeFakeFeatures(Ys{1});
    obj   = BipartiteMatchingPredict(theta, feats, 'reweight', rho, opts.objParams{:});
    algo  = opts.algo(obj, opts.algoParams{:});
    
    [objHist, testHist, exitflag] = algo.run();
    
    negLogZ = objHist(end).fval();
    negExactLogZ = -exactLogZRyser(feats, theta);
    
    fprintf('rho = %g, negLogZ = %g, negExactLogZ = %g\n', ...
        rho, negLogZ, negExactLogZ);
    
    % Compute linear terms
    logScores = zeros(M, 1);
    matTheta = reshape(theta, D, D);
    for m = 1:M
        matY         = permToMat(Ys{m});
        logScores(m) = sum(vec(matTheta .* matY));
    end
    
    logLikes      = logScores + negLogZ      - lambda/(2*M) * sum(theta.^2);
    exactLogLikes = logScores + negExactLogZ - lambda/(2*M) * sum(theta.^2);
end