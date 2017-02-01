classdef AbstractFW < handle
    % AbstractFW  Base class for Frank-Wolfe algorithms. Handles outer loop
    % and convergence checks.
    %
    % For example, the run() loop can handle some callbacks at a given
    % interval.
    
    properties
        algoName
        obj
        dualityGap
        opts        
        totIterTime
    end
    
    methods (Abstract)
        stepSz = iter(th);        
    end
    
    methods
        function th = AbstractFW(objective, varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            
            p.addRequired('objective', @(x) isa(x, 'BCFWObjective'));
            p.addParamValue('MaxIter', 50000);
            p.addParamValue('MaxTime', inf)
            p.addParamValue('fvalEpsilon', 0); % stopping criteria: if fval decreases less than epsilon in testInterval iterations
            p.addParamValue('printInterval', 10);
            p.addParamValue('testInterval', 100);
            p.addParamValue('trainErrInterval', 1);
            p.addParamValue('printComputedDualityGap', false);
            p.addParamValue('outputFcn', []);
            p.addParamValue('testSet', []);            
            
            p.addParamValue('checkpointFile', []);
            p.addParamValue('checkpointInterval', []);

            % Function which takes myself and iter as an argument and computes
            % whether to do line search or backoff this iteration.
            %
            % Defaults to return true regardless of input.
            
            p.addParamValue('lineSearchPolicy', @(algo, iter) true);
            
            p.parse(objective, varargin{:});
            th.obj = objective;
            th.opts = p.Results;     
                       
            th.totIterTime = 0;            
        end
        
        function [objHist, testHist, exitflag] = run(th, t)
            % Starting iteration is optional.
            if nargin == 1
                t = 1;
            end
            fprintf('DEBUG t = %d, nargin = %d\n', t, nargin);

            objHist = struct('iter', {}, 'fval', {}, 'dg', {}, 'relDg', {}, ...
                         'totIterTime', {}, ...
                         'trainErr', {}, ...
                         'dx1', {}, ...
                         'stepSz', {});
                     
            testHist = struct('iter', {}, ...
                              'objHistVal', {}, ...
                              'MAPErr', {}, ...
                              'MPMErr', {});            
                     
            
            previous_fval = inf;
            xOld = th.obj.x;
            while t <= th.opts.MaxIter && th.totIterTime <= th.opts.MaxTime
                tic;
                stepSz = th.iter(t);
                th.totIterTime = th.totIterTime + toc;

                if ~isempty(th.opts.outputFcn)
                    th.opts.outputFcn(t, th);
                end
                
                % Compute test error online
                testMAPErr = [];
                testMPMErr = [];
                if mod(t, th.opts.testInterval) == 0 && ~isempty(th.opts.testSet)
                    theta  = th.obj.computeParams();
                    MAPErr = th.opts.testSet.MAPErr(theta);
                    MPMErr = th.opts.testSet.MPMErr(theta);
                    
                    % TODO: If these intervals are not mod 0 on the same
                    % intervals, objHistVal will be out of date; broken;
                    testHist(end+1) = struct('iter', t, ...
                        'objHistVal', objHist(end), ...
                        'MAPErr', MAPErr, ...
                        'MPMErr', MPMErr);                    
                end
                                                
                if mod(t, th.opts.printInterval) == 0
                    % Save time by not solving all the LPs if unnecessary
                    % for the printing
                    if th.opts.printComputedDualityGap
                        [~, dg] = th.obj.solveLP();
                    else
                        % If we don't compute a new duality gap, just grab
                        % whatever we've cached. However, it may not always
                        % exist, in which case we'll just print blank.                        
                        dg = th.dualityGap;
                    end
                                   
                    fval  = th.obj.fval();   
                    relDg = abs(dg / fval);
                    
                    objHist(end+1) = struct('iter', t, ...
                        'fval', fval, ...
                        'dg',  dg, ...
                        'relDg', relDg, ...
                        'totIterTime', th.totIterTime, ...
                        'trainErr', th.obj.trainErr(), ...
                        'dx1', norm(th.obj.x - xOld, 1), ...
                        'stepSz', stepSz);
                    
                    % TODO: Factor out
                    xOld = th.obj.x;                  
                                        
                    fprintf('[%s] %s', th.algoName, sprintfStruct(objHist(end)));  
                    
                    if abs(previous_fval - fval) < th.opts.fvalEpsilon
                        fprintf('CONVERGENCE REACHED\n');
                        break;
                    end
                    previous_fval = fval;
                end
                
                t = t + 1;
            end            
            exitflag = 0;            
        end 
        
        function [objHist, paramsHist, exitflag] = fastrun(th)
            % fastrun  Run computations without interleaved error
            % evaluation.
            t = 1;

            objHist = struct('iter', {}, 'fval', {}, 'dg', {}, 'relDg', {}, ...
                         'totIterTime', {}, ...
                         'dx1', {}, ...
                         'stepSz', {});
                     
            xOld = th.obj.x;            
            while t <= th.opts.MaxIter
                tic;
                stepSz = th.iter(t);
                th.totIterTime = th.totIterTime + toc;
                
                paramsHist{t} = th.obj.computeParams();
                
                dg = th.dualityGap;                                        
                fval  = th.obj.fval();
                relDg = abs(dg / fval);
                    
                objHist(t) = struct('iter', t, ...
                    'fval', fval, ...
                    'dg',  dg, ...
                    'relDg', relDg, ...
                    'totIterTime', th.totIterTime, ...
                    'dx1', norm(th.obj.x - xOld, 1), ...
                    'stepSz', stepSz);

                xOld = th.obj.x;                  
                if mod(t, th.opts.printInterval) == 0                                        
                    fprintf('[%s] %s', th.algoName, sprintfStruct(objHist(end)));                                                            
                end
                
                if ~isempty(th.opts.checkpointFile) && ...
                   ~isempty(th.opts.checkpointInterval) && ...
                   mod(t, th.opts.checkpointInterval) == 0               
                   % A full checkpoint is too slow to load. Just save the histories.
                    save('-v7.3', th.opts.checkpointFile, 'objHist', 'paramsHist');
                end
                
                t = t + 1;
            end            
            exitflag = 0;            
        end        
    end
    
end

