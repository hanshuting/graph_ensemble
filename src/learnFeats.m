function [theta, hist] = learnFeats(X, features, likeGrad, varargin)
% learnBetheFeats  Find Bethe-MLE, feature version (conditional learning)

    p = inputParser;
    p.KeepUnmatched = true; % Enable extraneous options.
    p.addRequired('X');
    p.addRequired('features', @iscell);    
    p.addRequired('likeGrad', @(f) isa(f, 'function_handle'));    
    
    % Initial values
    p.addParamValue('thetaInit', []);    
    
    % Regularization parameter (no idea how correct)
    p.addParamValue('lambda', 1);

    % 6/3 -- THERE IS A BUG IN WARMSTART. How to fix?
    p.addParamValue('warmStart', false)
    
    p.addParamValue('beliefs', []);
        
    % Tuning parameters
    % TODO: Scale this with the magnitude of the input data.
%     p.addParamValue('initStep', 1e-4);    
    p.addParamValue('initStep', 1);    
    p.addParamValue('normalizeInitStep', true);
    
    % Convergences (this seems a little tight)
    p.addParamValue('TolFun',  1e-4);  
%     p.addParamValue('TolDThetaMax', 1e-2);
    p.addParamValue('maxIter', Inf);
    p.addParamValue('TolDThetaMax', 0.1);    
    
    p.addParamValue('rho', []);
    p.addParamValue('useBackoffRule', true);
    p.addParamValue('stepSizeRule', @(init, it) init / sqrt(it));
    p.addParamValue('normalizeStepLength', false);
    
    % For the likeGrad
    p.addParamValue('likeGradExtraParamValues', {});
    
    % CURRENTLY IGNORED
    p.addParamValue('TolX',    1e-4);        
    
    % Debug
    p.addParamValue('printIter', true);
    p.addParamValue('thetaTrue', []);

    p.parse(X, features, likeGrad, varargin{:});            
    o = p.Results;    

    if o.normalizeInitStep
        o.initStep = o.initStep / o.lambda;
    end
    
    [M, K] = size(features);    
    FXBar = computeFeatInnerProduct(features, X);    

    if isempty(o.thetaTrue)
        thetaTrueNorm = NaN;
    else
        thetaTrueNorm = vec(o.thetaTrue ./ sum(o.thetaTrue));
    end
    
    if isempty(o.thetaInit)
        theta = zeros(K, 1);        
    else
        theta = o.thetaInit;
    end
    stepSz = o.initStep;
    
    % TODO: Perhaps change this?
    
    iter = 1;
    increaseIters = 0;  
    
    beliefs = o.beliefs;
    while iter <= 2 || ...
          abs(hist(end-1).obj - obj) > o.TolFun || ...
          maxDTheta > o.TolDThetaMax 
      
        % TODO: Clean up logic
        if iter > o.maxIter
            break;
        end

        % Batch gradient.
        % TODO: Stochastic gradient. Need to tune step sizes.
        % TODO: Think of some code to make adaptive minibatches.
        
        % Compute partition functions and update objective.
        if ~o.warmStart
            % If no warm-start, erase them.
            beliefs = [];
        end
                
        [obj, dTheta, beliefs] = likeGrad(theta, features, FXBar, 'rho', o.rho, 'X', X, 'beliefs', beliefs, o.likeGradExtraParamValues{:});
        
        % Add the L2 + L1 penalty
        % l1 penalty is unstable...
%        obj    = obj    + o.lambda * (0.5 * sum(theta.^2) + sum(abs(theta)));
%        dTheta = dTheta + o.lambda * (theta + sign(theta));
        obj    = obj    + o.lambda * (0.5 * sum(theta.^2));
        dTheta = dTheta + o.lambda * theta;                       
        
        maxDTheta  = norm(dTheta, Inf); 
        meanDTheta = mean(abs(dTheta));
        
        if o.useBackoffRule
            if iter >= 2 && obj > hist(end).obj
                warning('Objective increased; decreasing step size');
                increaseIters = increaseIters + 1;
                stepSz = o.initStep / (1 + increaseIters);
            end
        else
            stepSz = o.stepSizeRule(o.initStep, iter);
            if o.normalizeStepLength
                % I find that normalizing step lengths is usually a bad idea
                stepSz = stepSz / norm(dTheta, 2);
            end
        end

        % Scale by 1/M for comparability with the no-feature version (and
        % to normalize step size tuning).
        
        theta = theta - stepSz / M * dTheta;                
        
        if o.printIter
            dTheta
            theta
        end

        if isempty(o.thetaTrue)
            thetaDist = NaN;
        else
            thetaDist  = norm(theta - o.thetaTrue, 1) / K;
        end
        
        % muMin was computed from last iteration's theta.
        % theta, obj are from this iteration.
        hist(iter) = var2struct(iter, beliefs, theta, obj, dTheta, thetaDist, maxDTheta, meanDTheta, stepSz);
        if o.printIter
            disp(hist(iter));
        end            
        
        iter = iter + 1;
    end
end


