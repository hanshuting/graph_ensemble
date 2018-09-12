classdef BipartiteMatchingPredict < BCFWObjective
    % BipartiteMatching  BCFW interface for bipartite matching problems.
       
    % Hopefully, only a trivial change is needed to go from th to
    % unipartite matchings.
    
    properties
        % iterate (actually a vector, as opposed to Ising)
        x
        % constants
        K
        
        % original problem data
        Ys, features
        
        % (vectorized) problem data
        theta, lambda, G, y, Gt, GtTimesY, rsum
        
        % intermediate data depending on x
        GTimesTheta
        
        % options
        opts
        
        % bookeeping
        N
    end
    
    methods
        function th = BipartiteMatchingPredict(theta, features, varargin)
            % th = BipartiteMatchingPredict(theta, Ys, features)
            %
            %   Initialize the objective.
            %   theta    : K vector of parameters whose likelihood we want to evaluate.
            %   features : M x K cell array. For each k, features{m,k} must
            %   be a square matrix whose rows and columns are the
            %   permutation size.
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Parse options (and set defaults)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
            p = inputParser;
            p.KeepUnmatched = true; % Enable extraneous options.

            p.addRequired('theta', @isnumeric);
            p.addRequired('features', @iscell);
            
            p.addParamValue('YTrue', []);

            p.addParamValue('beliefs', {}, @iscell);
            p.addParamValue('reweight', 1);
            p.addParamValue('TolGap', 1e-10);
            p.addParamValue('MaxIter', inf);            
            p.addParamValue('debug', false);            
            p.addParamValue('lineSearch', true);
            p.addParamValue('plotLine', false);
            p.addParamValue('checkStuck', false);
            p.addParamValue('errNegDualityGap', false);
            p.addParamValue('lineSearchOpts', optimset('TolX', 1e-9));

            p.parse(theta, features, varargin{:});
            th.theta = theta;
            th.opts = p.Results;

            [M2, th.K] = size(features);
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Set up problem data and bookkeeping
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            % Keep unvectorized features for evaluating error... although
            % we could just code the prediction routine to take
            % unfeaturized code... oh well.
            th.features = features;
            th.G = cell2mat(cellfun(@vec, features, 'UniformOutput', false));
            
%             th.y = vertcat(cell2mat(cellfun(@vec, Ymats, 'UniformOutput', false)));
            
            th.GTimesTheta = th.G * th.theta;

            % Store bookkeeping values.
            F1 = features{1};
            th.N = size(F1, 1);   % number of items

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Set up initial feasible point / warm start.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
            if isempty(th.opts.beliefs)
                th.x = 1/th.N * ones(th.N^2, 1);
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
            %% TODO: Suport reweight, particular vector ones (not just 1 param)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
            rho = th.opts.reweight * ones(th.N, 1);
            th.rsum = vec(bsxfun(@plus, rho, rho'));
        end        
                
        function xm = getXBlock(th, m)
            error('BipartiteMatchingPredict:getXBlock', 'Not supported -- UnipartiteMatchingPredict only for one sample!');                        
        end
        
        function err = trainErr(th)
            % Warning: Inaccurate because it is not MPM error.
            % (However, UnipartiteMatchingPredict... It is wrong.)
%             theta = th.computeParams();
            err = evalThetaCellY(th.theta, th.features, th.Ys);
        end
        
        function [s, dualityGap] = solveLP(th)
            g       = th.grad();
            permMat = csaAssignPermMat_mex(reshape(g, th.N, th.N));
            s       = permMat(:);
            
            dualityGap = sum(g .* (th.x - s));
        end
        
        function sm = solveLPBlock(th, m)
            % No need to block stuff
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
            error('BipartiteMatchingPredict:lineSearchBlock', 'Block methods not supported.');
        end                 
                
        function f = fval(th)
            f = -th.GTimesTheta.' * th.x + negBetheEntropy(th.rsum, th.x);
        end                
        
        function theta = computeParams(th) 
            error('BipartiteMatchingPredict:computeParams', 'Not supported -- UnipartiteMatchingPredict already has params!');
        end
        
        function beliefs = reshapeBeliefs(th)
            beliefs = cell(th.M, 1);
            for m = 1:th.M
                beliefs{m} = reshape(th.x(th.ix(m)), ...
                                     th.Ns(m), ...
                                     th.Ns(m));
            end
        end
       
    end
    
    methods (Access = protected)
        function y = ix(th, m)
            y = th.beginRow(m):th.endRow(m);
        end
        
        function updateIntermediateValues(th)
            % no-op
        end               
        
        function g = grad(th)
            g = -th.GTimesTheta + ...
                th.rsum + log(th.x) + (th.rsum - 1).*log(1 - th.x);                
        end
        
        function gm = gradBlock(th, m)
            ix = th.ix(m);
            
            xix = th.x(ix);
            rix = th.rsum(ix);
            
            gm = -(th.G(ix,:) * th.GtTimesYMinusTau)/th.lambda + ...
                rix + log(xix) + (rix - 1).*log(1 - xix);
        end                                                
    end
    
end

% Method calls are not JITted, so we move out this static function instead.
function H = negBetheEntropy(rsum, x)
    H = sum(x.*log(x) - (rsum - 1).*(1 - x).*log(1 - x)); 
end
