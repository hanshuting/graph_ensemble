classdef BCFW < AbstractFW
        
    methods
        function th = BCFW(objective, varargin)            
            th = th@AbstractFW(objective, varargin{:});
            th.algoName = 'BCFW';            

%             p = inputParser;
%             p.KeepUnmatched = true;            
%             p.parse(varargin{:});
%             o = p.Results;            
        end
        
        function stepSz = iter(th, t)
            m = randi(th.obj.M);
            sm = th.obj.solveLPBlock(m);
            dm = sm - th.obj.getXBlock(m);

            backupStepSz = 2*th.obj.M / (2*th.obj.M + t);
            stepSz = backupStepSz;

            if th.opts.lineSearchPolicy(th, t)
                [stepSz, converged] = th.obj.lineSearchBlock(m, dm, 1);
                if ~converged
                    warning('BCFW:linesearchDidNotConverge', ...
                        'Linesearch did not converge; backing off to 2N / (2N + t).');
                    stepSz = backupStepSz;
                end
            end
            th.obj.moveXBlock(m, stepSz, dm);
        end
    end
end

