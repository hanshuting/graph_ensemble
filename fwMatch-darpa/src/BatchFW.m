classdef BatchFW < AbstractFW
        
    methods
        function th = BatchFW(objective, varargin)
            th = th@AbstractFW(objective, varargin{:});
            th.algoName = 'BatchFW';
            
%             p = inputParser;
%             p.KeepUnmatched = true;            
%             p.parse(varargin{:});
%             o = p.Results;
        end
        
        function stepSz = iter(th, t)            
            [s, th.dualityGap] = th.obj.solveLP();
            d = s - th.obj.x;
            
            backupStepSz = 2 / (2 + t);
            stepSz = backupStepSz;
            
            if th.opts.lineSearchPolicy
                [stepSz, converged] = th.obj.lineSearch(d, 1);
                if ~converged
                    stepSz = backupStepSz;
                    warning('BatchFW:linesearchDidNotConverge', ...
                        'Linesearch did not converge; backing off to 2 / (2 + t).');                    
                end
            end
                
            th.obj.moveX(stepSz, d);
        end
        
    end
    
end

