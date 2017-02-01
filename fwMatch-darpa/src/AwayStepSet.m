classdef AwayStepSet < handle
    % AwayStepSet  Manage bookkeeping and solve LPs for the away set vertices.
    
    properties
        % Each vertex is stored as a ROW vector.
        verts
        weights                        
        actives
        d
        
        activeRows
        
        BEGIN_VERTS = 2
        nUsed
        size
        capacity
    end       
    
    methods
        function th = AwayStepSet(x0)
            % d:  dimensionality
            % x0: initial point.
            th.capacity = th.BEGIN_VERTS;  
            th.d = length(x0);
            
            th.clearAndReplace(x0);
        end
        
        function [v, maxStep] = solveAwayLP(th, grad)
            vals = th.verts(th.actives,:) * grad;
            [~, i] = min(vals);
            
            ai      = th.activeRows(i);
            v       = th.verts(ai,:).';
            maxStep = th.weights(ai) / (1 - th.weights(vi));
        end                
        
        function updateFWStep(th, fwVert, stepSz)
            if stepSz == 1.0
                % degenerate case; clear out our set.
                th.clearAndReplace(fwVert)
            else
                % Is fwVert already in verts?
                % CAUTION: Your vertices better be integral or
                % half-integral for this to work. Don't compare floating
                % point numbers by equality!
                [tf, row] = ismember(fwVert, th.verts);                
                
                if ~tf
                    % Vertex not already in our set => add it.
                    row = th.addVert(fwVert, stepSz);                    
                end
                
                % Update weights: Downweight everyones and also upweight ours.
                th.weights(th.actives) = (1 - stepSz) * th.weights(th.actives);
                th.weights(row)        = th.weights(row) + stepSz;                
            end
        end
        
        function updateAwayStep(th, awayVert, stepSz)
            [tf, row] = ismember(awayVert, th.verts);
            assert(tf, 'updateAwayStep called with a vertex not in array!');
            
            if stepSz == 1.0
                % degenerate case; drop step
                th.removeVert(row);
            end
            
            th.weights(th.actives) = (1 + stepSz) * th.weights(th.actives);
            th.weights(row)        = th.weights(row) - stepSz;
        end                
    end
    
    methods (Access = private)
        function compact(th)
            th.verts(~th.actives) = [];
            th.weights(~th.actives) = [];
            th.capacity = length(th.weights);            
            th.nUsed = th.capacity;
            assert(th.capacity == th.size);                        
            
            th.actives = true(th.capacity, 1);
            th.activeRows = 1:th.capacity;
        end
        
        function row = addVert(th, v, stepSz)
            % If we took a Frank-Wolfe step to a new vertex, we need to add
            % ourselves. Since we just took the step, our weight will be
            % stepSz.
            if th.nUsed == th.capacity
                th.capacity = th.capacity * 2;
                
                th.verts(th.capacity,:) = 0;
                th.weights(th.capacity) = 0;
                th.actives(th.capacity) = 0;
                
                % Recurse to actually do the now
                row = addVert(th, v, stepSz);
            end
            
            th.nUsed = th.nUsed + 1;
            th.size  = th.size + 1;

            th.verts(th.nUsed,:) = v;
            th.weights(th.nUsed) = stepSz;
            th.actives(th.nUsed) = true;                                                
            
            th.activeRows(th.size) = th.nUsed;
            row = th.nUsed;
        end
        
        function removeVert(th, row)
            th.size = th.size - 1;
            th.actives(row) = false;
            th.activeRows(row) = [];
        end
        
        function clearAndReplace(th, x0)
            th.verts   = zeros(th.capacity, th.d);
            th.weights = zeros(th.capacity, 1);
            th.actives = false(th.capacity, 1);
            th.activeRows = zeros(th.capacity, 1);
            
            th.verts(1,:) = x0.';
            th.weights(1) = 1.0;
            th.actives(1) = true;            
            th.activeRows(1) = 1;     
            
            th.nUsed = 1;            
            th.size  = 1;
        end
            
    end
    
end

