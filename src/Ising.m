classdef Ising < BCFWObjective
    % Ising  BCFW interface for Ising models.

    properties (Dependent)
        % Flat vector interface for the BatchFW routine. The vectorization
        % format is
        %   x = vertcat(vec(TN'), vec(TE')) = vectorizeNE(TN, TE)
        x
    end
    
    properties
        % vectorized problem data
        lambda, Ut, Vt, YN, YE, UtTimesY, VtTimesY, RN, RE
        % numerical (not 1-hot) label; used to evaluate training error
        YNflat
        
        % intermediate data depending on x; w = mu - y
        TN, TE, UtTimesW, UtTimesWSumSq, VtTimesW, VtTimesWSumSq, edges        
                
        % bookeeping
        Ns, nEdges, ixNode, ixEdge, lastNodeRow, opts
    end
    
    methods
        function th = Ising(YN, YE, Ut, Vt, Ns, edges, lambda, varargin)
            % th = Ising(YN, YE, Ut, Vt, Ns, edges, lambda, varargin)
            %
            %   Initialize the objective. 
            %   YN, YE  : Overcomplete node and edge labels
            %   Ut, Vt  : TRANSPOSED vercomplete node and edge feature matrices
            %             Row index feature and columns index sample/node.
            %             This layout makes our products more efficient.
            %   Ns      : M vector listing number of nodes in each sample            
            %   edges   : M cell vector of 2 x nEdges(m) edge list
            %   lambda  : Regularization constant
            %   varargin : See code for various ParamValue options.
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Parse options (and set defaults)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
            p = inputParser;
            p.KeepUnmatched = true; % Enable extraneous options.

            labelComplete = @(x) ismatrix(x) && all(abs(sum(x, 2)-1) <= eps(1));
            
            p.addRequired('YN', labelComplete);
            p.addRequired('YE', labelComplete);            
            p.addRequired('Ut', @ismatrix);
            p.addRequired('Vt', @ismatrix);
            p.addRequired('Ns', @isvector);            
            p.addRequired('edges', @iscell);
            p.addRequired('lambda', @isscalar);            
            
            p.addParamValue('reweight', 0.5, @isscalar); % 0.5 is TRW.
            % Optionally add test set, and I will compute the Hamming error
            % wrt MAP estimates.
            p.addParamValue('test', {}, @(x) all(isfield(x, {'YN', 'YE', 'Ut', 'Vt', 'edges'})));
            p.addParamValue('TolGap', 1e-10);
            p.addParamValue('debug', true);            
            p.addParamValue('initLabelProp', 0);
            p.addParamValue('errNegDualityGap', false);
            p.addParamValue('lineSearchOpts', optimset('TolX', 1e-7));
            p.addParamValue('checkStuck', true);

            p.parse(YN, YE, Ut, Vt, Ns, edges, lambda, varargin{:});
            th.opts = p.Results;

            th.M = length(Ns);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Set up bookkeeping indices (problem data already vectorized)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            th.Ns                  = Ns;
            th.nEdges              = cellfun(@(x) size(x, 2), edges);            
            [th.ixNode, th.ixEdge] = computeUVIx(Ns, th.nEdges);
                        
            % The bulk move callbacks will put (X_N, X_E) into one long
            % vector. The lastNodeRow field stores where to split.
            th.lastNodeRow = 2*sum(Ns);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Set up bookkeeping for test set, if exists
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            if ~isempty(th.opts.test)
                testNEdges = cellfun(@(x) size(x, 2), th.opts.test.edges);                
                [th.opts.test.ixNode, th.opts.test.ixEdge] = ...
                    computeUVIx(th.opts.test.Ns, testNEdges);
            end            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Set up initial feasible point, mixing uniform and label.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                         
            p = th.opts.initLabelProp;            
            th.TN = (1 - p)*(0.5  * ones(size(YN))) + p*YN;
            th.TE = (1 - p)*(0.25 * ones(size(YE))) + p*YE;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Set up overcomplete reweighting *matrix*.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
            % TODO: Possibly suboptimal; bsxfun may be faster?
            th.RN = zeros(size(YN));
            th.RE = zeros(size(YE));                       
            fprintf('Debug: Ising edgeRhos');
            for p = 1:th.M
                fprintf('+')
                % Expand our constant rhos to a general edge-vector rhos (a
                % constant vector of size nEdges).                
                edgeRhos = th.opts.reweight * ones(th.nEdges(p), 1);
                
                [th.RN(th.ixN(p),:), th.RE(th.ixE(p),:)] = ...
                    makeRhoMatOvercomplete(edgeRhos, edges{p}, Ns(p), 2);
            end
            fprintf('\n');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Copy labels and compute first intermediate values
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
            th.Ut = Ut; th.Vt = Vt; th.YN = YN; th.YE = YE;
            th.edges = edges;
            th.lambda = lambda;
            
            % Compute the flat labels (to compute training error)
            [~, IN] = max(YN,[],2);
            [Ni, Nj] = ind2sub(size(YN),IN);
            th.YNflat = zeros(length(Ni), 1);            
            th.YNflat(Ni) = Nj;
            th.YNflat(find(YN(:,1) == 0.5)) = 1.5;
            
            th.updateIntermediateValues();
        end  
        
        function v = get.x(th)
            v = vertcat(vec(th.TN'), vec(th.TE'));
        end
                        
        function xm = getXBlock(th, m)            
            xmN = th.TN(th.ixN(m),:);
            xmE = th.TE(th.ixE(m),:);
            % TODO: Clean up all this gratuitous use of .'
            xm = vectorizeNE(xmN, xmE);
        end
        
        function moveX(th, stepSz, d)
            [dN, dE] = th.unpack(d);
            
            th.TN = th.TN + stepSz * dN;
            th.TE = th.TE + stepSz * dE;
           
            th.updateIntermediateValues();
        end
        
        function moveXBlock(th, m, stepSz, dm)
            [dmN, dmE] = th.unpackBlock(m, dm);
            
            miN = th.ixN(m);miE = th.ixE(m);
            
            th.TN(miN,:) = th.TN(miN,:) + stepSz * dmN;
            th.TE(miE,:) = th.TE(miE,:) + stepSz * dmE;
                        
            % Update our state variables. This update only accesses the
            % rows of U, V corresponding to sample m. This is because all
            % other inner products with the block-coordinate step dm will
            % be zero. See derivation for details.
            th.UtTimesW = th.UtTimesW + stepSz*th.Ut(:,miN) * dmN;
            th.VtTimesW = th.VtTimesW + stepSz*th.Vt(:,miE) * dmE;
            
            th.updateSumSquares();
        end
        
        function [s, dualityGap] = solveLP(th)
            [GN, GE] = th.grad();
            sN = zeros(size(th.TN));
            sE = zeros(size(th.TE));            
            
            for m = 1:th.M
                % Sorry this is so ugly... It just means work on the GN and
                % GE matrices row-blockwise, storing the results into the
                % corresponding row-blocks of smN and smE.
                [SNt, SEt, eBelow] = ...
                    solveQPBO(GN(th.ixN(m),:).', GE(th.ixE(m),:).', th.edges{m});
                
                % bad code :(
                sN(th.ixN(m),:) = SNt';
                sE(th.ixE(m),:) = SEt';
            end
            
            s = vectorizeNE(sN, sE);

            dgN = frobProd(GN, th.TN - sN);
            dgE = frobProd(GE, th.TE - sE);            
            dualityGap = dgN + dgE;
        end
        
        function sm = solveLPBlock(th, m)
            [gNode, gEdge] = th.gradBlock(m);
            [smNode, smEdge, eBelow] = solveQPBO(gNode.', gEdge.', th.edges{m});
            sm = packOvercomplete(smNode, smEdge);
        end
        
        function [stepSz, converged] = lineSearch(th, dir, maxStep)
            [dN, dE] = th.unpack(dir);
            % Now this part is actually simpler because we don't have to do
            % so much goddamned indexing
            UtTimesD = th.Ut * dN;
            VtTimesD = th.Vt * dE;
            
            const     = th.UtTimesWSumSq + th.VtTimesWSumSq;
            
            linCoeff  = frobProd(th.UtTimesW, UtTimesD) + ...
                        frobProd(th.VtTimesW, VtTimesD);
                   
            quadCoeff = frobProd(UtTimesD, UtTimesD) + ...
                        frobProd(VtTimesD, VtTimesD);
            
            % Hoist property accesses from the inner loop, as these incur
            % overhead.
            RN = th.RN; RE = th.RE; TN = th.TN; TE = th.TE;
            
            dirHt = @(eta) 1/(2*th.lambda) * (const + linCoeff*eta + quadCoeff*eta^2) + ...
                frobProd((1 - RN) .* (TN + eta*dN), log(TN + eta*dN)) + ...                
                frobProd(RE .* (TE + eta*dE),       log(TE + eta*dE));
            
            [stepSz, dirHtMin, exitflag] = fminbnd(dirHt, 0, maxStep, th.opts.lineSearchOpts);        
            % DOUBLE CHECK: Correct sign.
            if exitflag ~= 1
                converged = false;
            end            

            if th.opts.checkStuck && dirHt(0) < dirHtMin
                warning('ht(0) was smaller than htMin: %g < %g', dirHt(0), dirHtMin);
                converged = false;
            end                                    
            converged = true;
        end
        
        % REFACTOR: All the crap actually calling fminbnd should be handled
        % by the parent class. (Nah, parent class would have a fully
        % numeric linesearch.)
        function [stepSzm, converged] = lineSearchBlock(th, m, d, maxStep)            
            [dN, dE] = th.unpackBlock(m, d);

            % Row indices and variable blocks.
            miN = th.ixN(m); TNm = th.TN(miN,:);
            miE = th.ixE(m); TEm = th.TE(miE,:);

            % Tediously precompute all of our quantities.
            UtmTimesD = th.Ut(:,th.ixN(m)) * dN;
            VtmTimesD = th.Vt(:,th.ixE(m)) * dE;
            
            const     = th.UtTimesWSumSq + th.VtTimesWSumSq;
            
            linCoeff  = frobProd(th.UtTimesW, UtmTimesD) + ...
                        frobProd(th.VtTimesW, VtmTimesD);
                    
            quadCoeff = frobProd(UtmTimesD, UtmTimesD) + ...
                        frobProd(VtmTimesD, VtmTimesD);
                        
            oneMinusRhoNm = 1 - th.RN(th.ixN(m),:);
            rhoEm = th.RE(th.ixE(m),:);

            % Fully vectorize everything for speed. This is the tightest
            % bottleneck, and Newton will not help us.
            vecOneMinusRhoNm = oneMinusRhoNm(:).';
            vecRhoEm         = rhoEm(:).';
            
            vecTNm = TNm(:); vecTEm = TEm(:); vecdN = dN(:); vecdE = dE(:);
            
            % Read this as:
            %   Quadratic terms + Singleton Entropy + Pairwise Entropy.            
            dirHt = @(eta) 1/(2*th.lambda) * (const + linCoeff*eta + quadCoeff*eta^2) + ...
                vecOneMinusRhoNm * ((vecTNm + eta*vecdN) .* log(vecTNm + eta*vecdN)) + ...
                vecRhoEm * ((vecTEm + eta*vecdE) .* log(vecTEm + eta*vecdE));
            
            [stepSzm, dirHtMin, exitflag] = fminbnd(dirHt, 0, maxStep, th.opts.lineSearchOpts);        
            if exitflag ~= 1
                converged = false;
            end            
            
%             % TESTING CODE
%             [stepSzmFast, myObj] = th.lineSearchBlockFast(m, d, maxStep);
%             assertElementsAlmostEqual(dirHt(stepSzm), dirHt(stepSzmFast));

%             fprintf('fminbnd got %g; newton got %g; dirtHt diff = %g, newtonObjDiff = %g\n', ...
%                 stepSzm, stepSzmFast, ...
%                 dirHt(stepSzm) - dirHt(stepSzmFast), ...
%                 myObj(stepSzm) - myObj(stepSzmFast));
%             
            if th.opts.checkStuck && dirHt(0) < dirHtMin
                warning('ht(0) was smaller than htMin: %g < %g', dirHt(0), dirHtMin);
                converged = false;
            end                                    

            converged = true;
        end        
        
        function [stepSzm, myObj] = lineSearchBlockFast(th, m, d, maxStep)
            [dN, dE] = th.unpackBlock(m, d);
            
            % Row indices and variable blocks.
            miN = th.ixN(m); TNm = th.TN(miN,:);
            miE = th.ixE(m); TEm = th.TE(miE,:);
            
            % Tediously precompute all of our quantities.
            UtmTimesD = th.Ut(:,th.ixN(m)) * dN;
            VtmTimesD = th.Vt(:,th.ixE(m)) * dE;
            
            % Fully vectorize everything for speed. This is the tightest
            % bottleneck, and Newton will not help us.
            oneMinusRhoNm = 1 - th.RN(th.ixN(m),:);
            rhoEm         = th.RE(th.ixE(m),:);
            
            vecOneMinusRhoNm = oneMinusRhoNm(:).';
            vecRhoEm         = rhoEm(:).';
            
            vecTNm = TNm(:); vecTEm = TEm(:); vecdN = dN(:); vecdE = dE(:);
                        
            % Constant, linear, and quadratic terms of the line search
            % objective.
            c0 = th.UtTimesWSumSq + th.VtTimesWSumSq;
            c1 = frobProd(th.UtTimesW, UtmTimesD) + ...
                frobProd(th.VtTimesW, VtmTimesD);
            c2 = frobProd(UtmTimesD, UtmTimesD) + ...
                frobProd(VtmTimesD, VtmTimesD);
            
            vecOneMinusRhoNm = vec(1 - th.RN(th.ixN(m),:)).';
            vecRhoEm = vec(th.RE(th.ixE(m),:)).';
            
            % Hoist outside of the objective function; property access is slow.
            lam = th.lambda;
            
            function [f, g, h] = obj(eta)
                vecXNm = vecTNm + eta*vecdN;
                vecXEm = vecTEm + eta*vecdE;

                vecLogXNm = log(vecXNm);
                vecLogXEm = log(vecXEm);

                vecXLogXNm = vecXNm .* vecLogXNm;
                vecXLogXEm = vecXEm .* vecLogXEm;                                
                
                f = 1/(2*lam) * (c0 + c1*eta + c2*eta^2) + ...
                    vecOneMinusRhoNm * vecXLogXNm + ...
                    vecRhoEm * vecXLogXEm;

                checkF = 1/(2*th.lambda) * (c0 + c1*eta + c2*eta^2) + ...
                    vecOneMinusRhoNm * ((vecTNm + eta*vecdN) .* log(vecTNm + eta*vecdN)) + ...
                    vecRhoEm * ((vecTEm + eta*vecdE) .* log(vecTEm + eta*vecdE));
                
                assertElementsAlmostEqual(f, checkF);
                
                g = c1/(2*lam) + eta * (c2/lam + ...
                    vecOneMinusRhoNm * (1 + vecLogXNm) + ...
                    vecRhoEm * (1 + vecLogXEm));
                
                h = c2/lam + eta^2 * (vecOneMinusRhoNm * (1./vecXNm + vecLogXNm + 1) + ...
                    vecRhoEm * (1./vecXEm + vecLogXEm + 1));                
            end
            
            % I don't trust Matlab's derivative check...
            f = 0; fOld = Inf;
            eta = 0.5;
            iter = 1;
            % TODO: Need a relative bound.
            while abs(fOld - f) > 1e-6;
                fOld = f;
                [f, g, h] = obj(eta);
                fprintf('Newton LS: iter = %g, eta = %g, f = %g, g = %g, h = %g, g / h = %g\n', iter, eta, f, g, h, g/h);
                eta = eta - g / h;
                if eta <= 0 || eta >= maxStep                    
                    fprintf('Newton went out of bounds to %g; randomly restarting.\n', eta);
                    eta = maxStep * rand;
                end
                iter = iter + 1;
            end
            stepSzm = eta;            
            myObj = @obj;            
        end        
               

	% Training likelihood 
        function f = fval(th)
            sq = 1/(2*th.lambda)*(th.UtTimesWSumSq + th.VtTimesWSumSq);
            f  = sq + frobProd(1 - th.RN, th.TN .* log(th.TN)) + ...
                      frobProd(th.RE,     th.TE .* log(th.TE));
        end  
        
        function params = computeParams(th)                        
            % KT, 10/21 -- Confirmed this is the correct sign (reworked the
            % math on paper.)
            params.F = -(1/th.lambda) * th.UtTimesW.';
            params.G = -(1/th.lambda) * th.VtTimesW.';
            params.lambda = th.lambda;
        end        
	
        function s = score(th) 
            params = th.computeParams();
            s = th.fval() + (th.lambda/2)*(frobProd(params.F,params.F) + frobProd(params.G,params.G));
        end             

        function logZ = partition_function(th,params)
           % compute partition function (linear term minus score divided by number of samples)
           logZ = (th.linearTerm(params) - th.score()) / th.M;
        end
        
        function thetaTx = linearTerm(th,params)
            thetaN = (params.F*th.Ut)';
            thetaE = (params.G*th.Vt)';
            linearN = vec(thetaN .* th.YN);
            linearE = vec(thetaE .* th.YE); 
            thetaTx = (sum(linearN)) + (sum(linearE));
        end
        
        function err = trainErr(th)
            % Compute Hamming error of *maximum posterior marginals.*
            [~, predNflat] = max(th.TN, [], 2);
            Nerr = (floor(th.YNflat) ~= predNflat) & (ceil(th.YNflat) ~= predNflat);
            err = mean(Nerr);
        end
        
        function err = trainMAPErr(th)            
            params = th.computeParams();
            err    = computeMAPErr(params.F, params.G, th.Ut, th.Vt, th.YN, th.ixNode, th.ixEdge, th.edges);
        end
        
        function err = testMAPErr(th)
            params = th.computeParams();
            err    = computeMAPErr(params.F, params.G, th.opts.test.Ut, th.opts.test.Vt, th.opts.test.YN, ...
                                   th.opts.test.ixNode, th.opts.test.ixEdge, th.opts.test.edges);
        end        
    end
    
    methods (Access = protected)
        function y = ix(th, m)
            error('Unsupported.');
        end        
        
        function [vN, vE] = unpack(th, v)
            % Unpack an overcomplete vector for the entire dataset into
            % node and edge matrices conforming to th.TN, th.TE.            
            vN = reshape(v(1:th.lastNodeRow),     fliplr(size(th.TN)))';
            vE = reshape(v(th.lastNodeRow+1:end), fliplr(size(th.TE)))';
        end
        
        function [vN, vE] = unpackBlock(th, m, v)
            % Unpack overcomplete vector for one sample into node and edge
            % matrices conforming to th.TN(th.ixN(m),:), th.TE(th.ixE(m),:).
            %
            % To communicate with the algorithms, we packed the node and
            % edge variables into a single vector. However, we store them
            % in this class in separate variables, to faciliate an easier
            % notation.
            
            splitRow = 2 * th.Ns(m);
            vN = reshape(v(1:splitRow),     2, th.Ns(m)).';
            vE = reshape(v(splitRow+1:end), 4, th.nEdges(m)).';
        end        
        
        function iN = ixN(th, m)
            iN = th.ixNode(m,1):th.ixNode(m,2);
        end
        
        function iM = ixE(th, m)
            iM = th.ixEdge(m,1):th.ixEdge(m,2);
        end
                
        function [GN, GE] = grad(th)
            % Actually just like gradBlock, but with less indexing.
            gHNode = (1 - th.RN) .* (1 + log(th.TN));
            gHEdge = th.RE       .* (1 + log(th.TE));
            
            GN = (1/th.lambda) * (th.Ut.'*th.UtTimesW) + gHNode;
            GE = (1/th.lambda) * (th.Vt.'*th.VtTimesW) + gHEdge;
        end
        
        function [gNode, gEdge] = gradBlock(th, m)            
            gHNode = (1 - th.RN(th.ixN(m),:)) .* (1 + log(th.TN(th.ixN(m),:)));
            gHEdge = th.RE(th.ixE(m),:)       .* (1 + log(th.TE(th.ixE(m),:)));                        
            
            gNode = (1/th.lambda) * (th.Ut(:,th.ixN(m)).'*th.UtTimesW) + gHNode;
            gEdge = (1/th.lambda) * (th.Vt(:,th.ixE(m)).'*th.VtTimesW) + gHEdge;
        end
        
        function updateIntermediateValues(th)
            % Update all state variables. Tricks don't yield any savings
            % here, so we just straightforwardly transcribe the math.
            th.UtTimesW = th.Ut*(th.TN - th.YN);
            th.VtTimesW = th.Vt*(th.TE - th.YE);
            
            th.updateSumSquares();
        end
        
        function updateSumSquares(th)
            th.UtTimesWSumSq = frobProd(th.UtTimesW, th.UtTimesW);
            th.VtTimesWSumSq = frobProd(th.VtTimesW, th.VtTimesW);
        end
        
    end        
end
