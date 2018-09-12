classdef DiscLearningObjective < handle
    % BCFWObjective  Interface for block-coordinate Frank Wolfe problems.
    %
    % the Block methods should be thread-safe... global state should
    % somehow be factored to enable distribution... to think about later.
    
    properties
        M, xwavg, theta, allThetas
    end
    
    methods (Abstract)        
        sm = solveLPBlock(th,m);
        updateParams(step,sm,m);%updateParams needs to look at sm, look at y(m) and then update params based on the feats;
        % Training error, by whatever measure you want.
        err = trainErr(th);
        aT = getTheta(th);
    end
    
end

