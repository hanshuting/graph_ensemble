classdef AbstractDiscLearning < handle
    % AbstractFW  Base class for Frank-Wolfe algorithms. Handles outer loop
    % and convergence checks.
    %
    % For example, the run() loop can handle some callbacks at a given
    % interval.
    
    properties
        algoName
        obj
        %dualityGap
        opts
        totIterTime
    end
    
    methods (Abstract)
        stepSz = iter(th);
    end
    
    methods
        function th = AbstractDiscLearning(objective, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addParamValue('MaxIter', 50000);
            p.addParamValue('printInterval', 15);
            p.addParamValue('testInterval', 15);
            p.addParamValue('outputFcn', []);
            p.addParamValue('testSet', []);
            p.addRequired('objective', @(x) isa(x, 'DiscLearningObjective'));
            p.parse(objective, varargin{:});
            th.obj = objective;
            th.opts = p.Results;
            
            th.totIterTime = 0;
        end
        function [objHist, testHist, exitflag] = run(th)
            t = 1;
            
            %             objHist = struct('iter', {}, 'fval', {}, 'dg', {}, 'relDg', {}, ...
            %                          'totIterTime', {}, ...
            %                          'trainErr', {}, ...
            %                          'dx1', {}, ...
            %                          'stepSz', {});
            
            objHist = [];
            testHist = struct('iter', {}, ...
                'objHistVal', {}, ...
                'MAPErr', {}, ...
                'MPMErr', {});
            
            xOld = th.obj.x;
            while t <= th.opts.MaxIter
                tic;
                stepSz = th.iter(); %%this updates the parameters
                th.totIterTime = th.totIterTime + toc;
                
                if ~isempty(th.opts.outputFcn)
                    th.opts.outputFcn(t, th);
                end
                if mod(t, th.opts.printInterval) == 0
                    fprintf('train err = %f\n',th.obj.trainErr());
                end
                
                if mod(t, th.opts.testInterval) == 0 && ~isempty(th.opts.testSet)
                    MAPErr = th.opts.testSet.MAPErr(th.obj.getTheta());
                    MPMErr = -1;%th.opts.testSet.MPMErr(theta);
                    % TODO: If these intervals are not mod 0 on the same
                    % intervals, objHistVal will be out of date; broken;
                    testHist(end+1) = struct('iter', t, ...
                        'objHistVal', objHist, ...
                        'MAPErr', MAPErr, ...
                        'MPMErr', MPMErr);
                end
                
                t = t + 1;
            end
            exitflag = 0;
        end
        
        
    end
    
end

