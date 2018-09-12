classdef BipartiteMatching < BCFWObjective
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
        lambda, G, y, Gt, GtTimesY, rsum
        
        % intermediate data depending on x
        GtTimesYMinusTau, GtTimesYMinusTauSq
        
        % options
        opts
        
        % bookeeping
        Ns, beginRow, endRow          
    end
    
    methods
        function th = BipartiteMatching(Ys, features, lambda, varargin)
            % th = BipartiteMatching(Ys, features)
            %
            %   Initialize the objective. 
            %   Ys       : M cell vector of permutation vectors.
            %   features : M x K cell array. For each k, features{m,k} must
            %   be a square matrix whose rows and columns are the
            %   permutation size.
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Parse options (and set defaults)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
            p = inputParser;
            p.KeepUnmatched = true; % Enable extraneous options.

            p.addRequired('Ys', @iscell);
            p.addRequired('features', @iscell);
            p.addRequired('lambda', @(x) isnumeric(x) && isscalar(x));
            
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

            p.parse(Ys, features, lambda, varargin{:});
            th.opts = p.Results;

            th.M = length(Ys);
            [M2, th.K] = size(features);
            
            assert(th.M == M2, 'Y and features must have same number of examples');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Set up problem data and bookkeeping
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            % Keep unvectorized features for evaluating error... although
            % we could just code the prediction routine to take
            % unfeaturized code... oh well.
            th.features = features;
            th.Ys       = Ys;
            th.G = cell2mat(cellfun(@vec, features, 'UniformOutput', false));
            
            th.Ns       = zeros(th.M, 1);
            th.beginRow = zeros(th.M, 1);
            th.endRow   = zeros(th.M, 1);                                    
            Ymats = cell(th.M, 1);
            row = 1;
            for m = 1:th.M
                % Vectorize me
                Ymats{m} = expandPerm(Ys{m}, @zeros);                
                
                % Update bookkeeping indices
                th.Ns(m) = length(Ys{m});
                th.beginRow(m) = row;
                th.endRow(m)   = row + numel(Ymats{m}) - 1;
                
                row = th.endRow(m) + 1;
            end

            th.y = vertcat(cell2mat(cellfun(@vec, Ymats, 'UniformOutput', false)));
            th.lambda = lambda;
            th.GtTimesY = th.G'*th.y;           

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Set up initial feasible point.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
            if isempty(th.opts.beliefs)
                th.x = zeros(th.endRow(end), 1);
                for m = 1:th.M
                    th.x(th.ix(m)) = 1/th.Ns(m) * ones(th.Ns(m)^2, 1);
                end
            else
                for m = 1:th.M
                    th.x(th.ix(m)) = th.opts.beliefs{m}(:);                    
                end                
            end
            
            % Call this to set up the other dependent temporary vars
            th.moveX(0, 0);
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% TODO: Suport reweight, particular vector ones (not just 1 param)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
            th.rsum = zeros(th.endRow(end), 1);            
            for m = 1:th.M                
                rho = th.opts.reweight * ones(th.Ns(m));
                th.rsum(th.ix(m)) = vec(bsxfun(@plus, rho, rho'));
            end            
        end        
                
        function xm = getXBlock(th, m)
            xm = th.x(th.ix(m));
        end
        
        function err = trainErr(th)
            theta = th.computeParams();
            err = evalThetaCellY(theta, th.features, th.Ys);
        end
        
        function [s, dualityGap] = solveLP(th)
            g = th.grad();
            s = zeros(size(th.x));
            
            for m = 1:th.M
                gm = reshape(g(th.ix(m)), th.Ns(m), th.Ns(m));
                permMat = csaAssignPermMat_mex(gm);
                s(th.ix(m)) = permMat(:);
            end
            
            dualityGap = sum(g .* (th.x - s));
        end
        
        function sm = solveLPBlock(th, m)
            gm = reshape(th.gradBlock(m), th.Ns(m), th.Ns(m));
            permMat = csaAssignPermMat_mex(gm);
            sm = permMat(:);            
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
            GtTimesDir  = th.G'*dir;
            dirLinCoeff = -2 * sum(th.GtTimesYMinusTau .* GtTimesDir);            
            dirQuadCoeff = sum(GtTimesDir .^ 2);
            
            dirHt = @(eta) 1/(2*th.lambda) * (sum(th.GtTimesYMinusTauSq) + dirLinCoeff.*eta + dirQuadCoeff.*eta.^2) + ...
                negBetheEntropy(th.rsum, th.x + eta*dir);
            
            [stepSz, converged] = th.lineSearchHelper(dirHt, maxStep);
        end
        
        function [stepSzm, converged] = lineSearchBlock(th, m, dir, maxStep)            
            GtTimesDir  = th.G(th.ix(m),:).'*dir;
            
            dirLinCoeff = -2 * sum(th.GtTimesYMinusTau .* GtTimesDir);            
            dirQuadCoeff = sum(GtTimesDir .^ 2);
            
            rsum = th.rsum(th.ix(m));
            x    = th.x(th.ix(m));
            
            % Note that we call negBetheEntropy on a far smaller problem.
            dirHt = @(eta) 1/(2*th.lambda) * (sum(th.GtTimesYMinusTauSq) + dirLinCoeff.*eta + dirQuadCoeff.*eta.^2) + ...
                negBetheEntropy(rsum, x + eta*dir);
            
            [stepSzm, converged] = th.lineSearchHelper(dirHt, maxStep);
            
        end                 
                
        function f = fval(th)
            sq = 1/(2*th.lambda) * sum(th.GtTimesYMinusTauSq);
            f  = sq + negBetheEntropy(th.rsum, th.x);
        end                
        
        function theta = computeParams(th)            
            theta = 1 / th.lambda * th.G' * (th.y - th.x);                        
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
            th.GtTimesYMinusTau = th.G'*(th.y - th.x);
            th.GtTimesYMinusTauSq = th.GtTimesYMinusTau.^2;            
        end               
        
        function g = grad(th)
            g = -(th.G * th.GtTimesYMinusTau)/th.lambda + ...
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
