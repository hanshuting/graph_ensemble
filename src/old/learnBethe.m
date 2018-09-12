function [W, hist] = learnBethe(X, varargin)
% learnBethe  Find Bethe-MLE, no features version
%
%   Warning: uses a full DxD intermediate matrix!
%
%   [theta, hist] = learnBethe(X, minBethe) returns parameters theta and an
%   optional history stroy.
%
%   X is an NxD matrix of permutations (bipartite matching) or N cell array
%   of DxD sparse matries (general matching).
%
%   Notation
%    - theta: potentials
%    - mu: marginals/beliefs; muBar: empirical marginals; muMin: argmin
%
%
%   TODO
%    - Features.
%    - Parameterize the FBethe solver (bi/unipartite, Frank-Wolfe/BP)
%    - Regularization.

    p = inputParser;
    p.KeepUnmatched = true; % Enable extraneous options.
    p.addRequired('X', @isnumeric);    
    
    % Initial values
    p.addParamValue('thetaInit', []);    
        
    % Tuning parameters
    % TODO: Scale this with the magnitude of the input data.
    p.addParamValue('initStep', 10);    
    
    % Convergences (this seems a little tight)
    p.addParamValue('TolFun',  1e-6);    
    p.addParamValue('TolX',    1e-6);        
    
    p.addParamValue('rho', []);
    p.addParamValue('stepSizeRule', @(init, it) init / sqrt(it));
    p.addParamValue('normalizeStepLength', false);    
    
    p.addParamValue('exact', false);
    
    % Debug
    p.addParamValue('printIter', true);
    p.addParamValue('thetaTrue', []);
    
    p.parse(X, varargin{:});            
    o = p.Results;    
    
    if iscell(X)
        error('learnBethe:XCell', 'Unipartite X not supported.');
    else
        [N, D] = size(X);
        muBar = meanPerm(X);
    end
    
    if isempty(o.thetaInit)
        W = zeros(D, D);
    else
        W = o.thetaInit;
    end
   
    iter = 1;
    
    fwTolFun = 1e-6;

    % for now, only tol
    while iter <= 2 || abs(hist(end-1).obj - hist(end).obj) > o.TolFun
        if o.exact
            [logZ, muMin] = exactLogZNoFeatRyser(W);
            FBethe = -logZ;            
        else
            %[muMin, FBetheHist] = fwBipartite_mex(1/D * ones(D), W, o.rho, fwTolFun);            
            [muMin, FBetheHist] = fwBipartite_mex_old(1/D * ones(D), W, o.rho, fwTolFun);            
            FBethe = FBetheHist(end);
        end                
        
        stepSz = o.stepSizeRule(o.initStep, iter);
        dW = muMin - muBar;
        if o.normalizeStepLength
            stepSz = stepSz / norm(dW, 2);
        end            
        W = W - stepSz * dW;
                
        obj   = -sum(vec(W .* muBar)) + FBethe;
                
        dMu = norm(muBar - muMin, 1);
        
        % muMin was computed from last iteration's theta.
        % theta, obj are from this iteration.
        hist(iter) = var2struct(iter, muMin, W, obj, dMu, stepSz, FBethe);
        if o.printIter
            disp(hist(iter));
        end            
        
        
        iter = iter + 1;
    end
end

