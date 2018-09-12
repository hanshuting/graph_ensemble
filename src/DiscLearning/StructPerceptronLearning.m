classdef StructPerceptronLearning < AbstractDiscLearning
    
    properties
        initStepSize
    end
    
    methods
        function th = StructPerceptronLearning(objective, varargin)
            th = th@AbstractDiscLearning(objective, varargin{:});
            th.algoName = 'StructPerceptron';
            
            p = inputParser;
            p.KeepUnmatched = true;
            p.addParamValue('initStepSize', 1.0);
            p.parse(varargin{:});
            o = p.Results;
            th.initStepSize = o.initStepSize;
        end
        
        function stepSz = iter(th)
            %just use a constant step size for now
            stepSz = th.initStepSize;
            m = randi(th.obj.M);
            
            sm = th.obj.solveLPBlock(m);
            lossGradient = 1; %change this when doing better learning. will need to change API, since the loss gradient depends on the objective
            step = stepSz*lossGradient;
            th.obj.updateParams(step,sm,m); %updateParams needs to look at sm, look at y(m) and then update params based on the feats;
            
            %todo: call the LP solver, update the parameters
            %th.obj.
        end
        
    end
    
end

