classdef BCFWObjective < handle
    % BCFWObjective  Interface for block-coordinate Frank Wolfe problems.
    %
    % the Block methods should be thread-safe... global state should
    % somehow be factored to enable distribution... to think about later.
    
    properties
        M, xwavg
    end
    
    methods (Abstract)        

        stepSz = lineSearch(th, dir, maxStep)
        stepSzm = lineSearchBlock(th, m, dir, maxStep)
        
        [s, dualityGap] = solveLP(th)
        sm = solveLPBlock(th);
        % Are you sure you want to implement all of this? To think about.
        
        % Training error, by whatever measure you want.
        err = trainErr(th);
    end
    
    methods (Abstract, Access = protected)
        y = ix(th, m)                 % access indices for block m       
        updateIntermediateValues(th); % called after x changes
    end
    
    methods
        function moveX(th, stepSz, dir)
            th.x = th.x + stepSz*dir;
            th.updateIntermediateValues();            
        end
        
        function moveXBlock(th, m, stepSz, dir)
            th.x(th.ix(m)) = th.x(th.ix(m)) + stepSz*dir;
            th.updateIntermediateValues();
        end
        
        function xm = getXBlock(th, m)
            xm = th.x(th.ix(m));
        end
        
    end
    
end

