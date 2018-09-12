classdef IsingPredict < BCFWObjective
    % IsingPredict  Predictive (marginal inference) objective for Ising model
    
    properties
        % Iterate
        x, TN, TE
        
        % Problem data
        thetaN, thetaE, RN, RE, edges
        
        % Bookkeeping
        N, lastNodeRow, nEdges, opts
    end
    
    methods
        function v = get.x(th)
            v = vectorizeNE(th.TN, th.TE);
        end
        
        function th = IsingPredict(F, G, Ut, Vt, N, edges, varargin)
            p = inputParser;
            p.KeepUnmatched = true; % Enable extraneous options.
            p.addRequired('F', @ismatrix);
            p.addRequired('G', @ismatrix);            
            p.addRequired('Ut', @ismatrix);
            p.addRequired('Vt', @ismatrix);
            p.addRequired('N', @isscalar);            
            p.addRequired('edges', @ismatrix);
            
            p.addParamValue('reweight', 0.5, @isscalar); % 0.5 is TRW.
            % Optionally add test set, and I will compute the Hamming error
            % wrt MAP estimates.
            
            p.addParamValue('TolGap', 1e-10);
            p.addParamValue('errNegDualityGap', false);
            p.addParamValue('lineSearchOpts', optimset('TolX', 1e-7));
            p.addParamValue('checkStuck', false);
            p.addParamValue('YNtrueFlat', []);

            p.parse(F, G, Ut, Vt, N, edges, varargin{:});
            th.opts = p.Results;
            
            th.nEdges = size(edges, 2);
            edgeRhos = th.opts.reweight * ones(th.nEdges, 1);
            [th.RN, th.RE] = makeRhoMatOvercomplete(edgeRhos, edges, N, 2);
            
            th.N           = N;
            th.lastNodeRow = 2*N;
            
            th.thetaN = (F*Ut).';
            th.thetaE = (G*Vt).';
            th.M      = 1;
            
            th.TN = 0.5  * ones(N, 2);
            th.TE = 0.25 * ones(th.nEdges, 4);
            
            th.edges = edges;
        end
        
        function stepSzm = lineSearchBlock(th, m, dir, maxStep)
            error('IsingPredict:lineSearchBlock', 'Not supported -- IsingPredict only for one sample!');
        end
        
        function stepSz = lineSearch(th, dir, maxStep)
            [dN, dE] = th.unpack(dir);
            
            const     = -frobProd(th.TN, th.thetaN) - frobProd(th.TE, th.thetaE);                        
            linCoeff  = -frobProd(dN,    th.thetaN) - frobProd(dE,    th.thetaE);                    
                        
            oneMinusRhoNm = 1 - th.RN;
            rhoEm = th.RE;

            % Fully vectorize everything for speed. This is the tightest
            % bottleneck, and Newton will not help us.
            vecOneMinusRhoNm = oneMinusRhoNm(:).';
            vecRhoEm         = rhoEm(:).';
            
            vecTN = th.TN(:); vecTE = th.TE(:); vecdN = dN(:); vecdE = dE(:);
            
            % Read this as:
            %   Quadratic terms + Singleton Entropy + Pairwise Entropy.            
            dirHt = @(eta) const + eta*linCoeff + ...
                vecOneMinusRhoNm * ((vecTN + eta*vecdN) .* log(vecTN + eta*vecdN)) + ...
                vecRhoEm * ((vecTE + eta*vecdE) .* log(vecTE + eta*vecdE));
            
            [stepSz, dirHtMin, exitflag] = fminbnd(dirHt, 0, maxStep, th.opts.lineSearchOpts);        
            if exitflag ~= 1
                warning('fminbnd did not converge!');
            end            
                        
            if th.opts.checkStuck && dirHt(0) < dirHtMin
                % Exception: Don't move.                
                warning('dirHt(0) was smaller than htMin: %g < %g; not moving.', dirHt(0), dirHtMin);
                stepSz = 0;
%                 htMin = ht(0);
            end                                    
        end
        
        function f = fval(th)
            f = -frobProd(th.TN, th.thetaN) - frobProd(th.TE, th.thetaE) + ...
                frobProd(1 - th.RN, th.TN .* log(th.TN)) + ...
                frobProd(th.RE,     th.TE .* log(th.TE));
        end
        
        function [GN, GE] = grad(th)
            GN = -th.thetaN + (1 - th.RN) .* (1 + log(th.TN));
            GE = -th.thetaE + th.RE .* (1 + log(th.TE));
        end
        
        function [s, dualityGap] = solveLP(th)
            [GN, GE] = th.grad();
            
            [SNt, SEt, eBelow] = solveQPBO(GN.', GE.', th.edges);
            sN = SNt.'; sE = SEt.';
            s = vectorizeNE(sN, sE);
            
            dgN = frobProd(GN, th.TN - sN);
            dgE = frobProd(GE, th.TE - sE);            
            dualityGap = dgN + dgE;            
        end
                
        function err = trainErr(th)
            [~, predYNflat] = max(th.TN, [], 2);
            verr = th.opts.YNtrueFlat ~= predYNflat;
            err = mean(verr);
        end
                   
        function sm = solveLPBlock(th)
            error('IsingPredict:solveLPBlock', 'Not supported -- IsingPredict only for one sample!');            
        end               
        
        function moveX(th, stepSz, d)
            [dN, dE] = th.unpack(d);
            
            th.TN = th.TN + stepSz * dN;
            th.TE = th.TE + stepSz * dE;
        end
        
        function moveXBlock(th, m, stepSz, dir)
            error('IsingPredict:moveXBlock', 'Not supported -- IsingPredict only for one sample!');            
        end
        
        function xm = getXBlock(th, m)
            error('IsingPredict:getXBlock', 'Not supported -- IsingPredict only for one sample!');            
        end
        
        function [xN, xE] = unpack(th, x)
            xN = reshape(x(1:2*th.N),     2, th.N).';
            xE = reshape(x(2*th.N+1:end), 4, th.nEdges).';
        end
    end
        
    methods (Access = protected)
        function y = ix(th, m)
            error('IsingPredict:ix', 'Not implemented -- Not for block use.');
        end
        function updateIntermediateValues(th)
            disp('Nothing to do');
        end
    end    
end

