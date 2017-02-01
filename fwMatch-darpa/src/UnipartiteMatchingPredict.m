classdef UnipartiteMatchingPredict < BCFWObjective
    % UnipartiteMatchingPredict  BCFW interface for bipartite matching predictions.
    %
    % All matrices are _adjacency_ matrices and represented by the 
    % vectorization (column major, of course) of their upper triangular
    % parts.
    %
    %       Essentially, change the constructor (vectorization) and the 
    %       solver call.
    %
    %       Should do by a diff.
       
    properties
        % iterate (actually a vector, as opposed to Ising)
        x
        % constants
        K, N, E
        
        % original problem data        
        theta, features
                
        % (vectorized) problem data
        G, GTimesTheta, rsum, yTrueVec
        % note: GTimesTheta is a column vector.
                
        % options
        opts

        % no bookkeeping, since we have a single sample.
    end
    
    methods
        function th = UnipartiteMatchingPredict(theta, features, varargin)
            % th = UnipartiteMatchingPredict(theta, features)
            %
            %   Initialize the prediction (marginal inference) objective. 
            %   theta    : K x 1 vector of feature weights
            %   features : K cell array of features matrices; must be
            %              symmetric
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Parse options (and set defaults)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
            p = inputParser;
            p.KeepUnmatched = true; % Enable extraneous options.

            p.addRequired('theta', @iscolumn);
            p.addRequired('features', @iscell);
            
            p.addParamValue('dense', false);
            p.addParamValue('beliefs', []);
            p.addParamValue('reweight', 1);
            p.addParamValue('TolGap', 1e-10);
            p.addParamValue('MaxIter', inf);       
            % MaxIter will be the major degree of freedom.
            p.addParamValue('debug', false);            
            p.addParamValue('lineSearch', true);
            p.addParamValue('plotLine', false);
            p.addParamValue('checkStuck', false);
            p.addParamValue('errNegDualityGap', false);
            p.addParamValue('lineSearchOpts', optimset('TolX', 1e-9));
            
            p.addParamValue('YTrue', []);

            p.parse(theta, features, varargin{:});
            th.opts = p.Results;

            th.K = length(features);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Set up problem data and bookkeeping
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            th.features = features;
            th.theta    = theta;
            
            vecUTCurry = @(x) vecUT(x, th.opts.dense);
            
            % Store/compute features.
            th.G = cell2mat(cellfun(vecUTCurry, features, 'UniformOutput', false));
            th.GTimesTheta = th.G * th.theta;

            % Store bookkeeping values.
            F1 = features{1};
            th.N = size(F1, 1);   % number of nodes
            th.E = size(th.G, 1); % number of edges                        

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Set up initial feasible point / warm start.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
            if isempty(th.opts.beliefs)
                th.x = 1/th.N * ones(th.E, 1);
            else
                sz = [th.E, 1];
                assert(all(sz == th.opts.beliefs));
                % We will assume the user supplied us valid beliefs. The
                % typical use-case is passing a saved value of this x in
                % previously.
                th.x = th.opts.beliefs;
            end
            
            % Call this to set up the other dependent temporary vars
            th.moveX(0, 0);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Store and check the ground truth.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                    
            if ~isempty(th.opts.YTrue)
                th.yTrueVec = vecUT(th.opts.YTrue);
                assert(length(th.yTrueVec) == length(th.x));
            end
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% TODO: Suport reweight, particular vector ones (not just 1 param)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
            rho = ones(th.N, 1);
            th.rsum = vecUT(bsxfun(@plus, rho, rho'));
        end        
                
        function xm = getXBlock(th, m)
            error('UnipartiteMatchingPredict:getXBlock', 'Not supported -- UnipartiteMatchingPredict only for one sample!');                        
        end
        
        function err = trainErr(th)
            yHatVec = double(th.x > 0.5);
            % This criterion is wrong.
            err     = full(mean(yHatVec ~= th.yTrueVec));
        end
        
        function [s, dualityGap] = solveLP(th)
            g          = th.grad();
            s          = vec(vecUniMatch(g)); 
            dualityGap = sum(g .* (th.x - s));
        end
        
        function sm = solveLPBlock(th, m)
            % we only have one sample in the prediction part...
            sm = th.solveLP();
        end
        
        function [stepSz, converged] = lineSearchHelper(th, func, maxStep)
            [stepSz, fMin, exitflag] = fminbnd(func, 0, maxStep, th.opts.lineSearchOpts);        
            if exitflag ~= 1
                converged = false;
            end            
            
            if th.opts.checkStuck && func(0) < fMin
                warning('func(0) was smaller than fMin: %g < %g', func(0), fMin);
            end
            converged = true;            
        end        
        
        function [stepSz, converged] = lineSearch(th, dir, maxStep)
            dirConst    = -th.GTimesTheta.' * th.x;
            dirLinCoeff = -th.GTimesTheta.' * dir;
            
            dirHt = @(eta) dirConst + eta*dirLinCoeff + ...
                negBetheEntropy(th.rsum, th.x + eta*dir);
            
            [stepSz, converged] = th.lineSearchHelper(dirHt, maxStep);
        end
        
        function [stepSzm, converged] = lineSearchBlock(th, m, dir, maxStep)     
            error('UnipartiteMatchingPredict:lineSearchBlock', 'Block methods not supported.');
        end                 
                
        function f = fval(th)
            f = -th.GTimesTheta.' * th.x + negBetheEntropy(th.rsum, th.x);
        end                
        
        function theta = computeParams(th)     
            error('UnipartiteMatchingPredict:computeParams', 'Not supported -- UnipartiteMatchingPredict already has aprams!');
        end
        
        function beliefs = reshapeBeliefs(th)
            beliefs = reshape(th.x, th.N, th.N);
        end
        
    end
    
    methods (Access = protected)        
        function g = grad(th)
            g = -th.GTimesTheta + ...
                th.rsum + log(th.x) + (th.rsum - 1).*log(1 - th.x);                
        end                
        
        function gm = gradBlock(th, m)
            error('UnipartiteMatchingPredict:gradBlock', 'Block methods not supported.');
        end
        
        % You should really factor out a class that does NOT support blocks
        function y = ix(th, m)
            error('UnipartiteMatchingPredict:ix', 'Block methods not supported.');
        end
        
        function updateIntermediateValues(th)
            % no-op
        end
    end
    
end

% Method calls are not JITted, so we move out this static function instead.
function H = negBetheEntropy(rsum, x)
    H = sum(x.*log(x) - (rsum - 1).*(1 - x).*log(1 - x)); 
end
