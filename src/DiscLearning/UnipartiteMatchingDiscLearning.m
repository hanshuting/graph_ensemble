classdef UnipartiteMatchingDiscLearning < DiscLearningObjective
    % UnipartiteMatching  BCFW interface for bipartite matching problems.
    %
    % All matrices are _adjacency_ matrices and represented by the 
    % vectorization (column major, of course) of their upper triangular
    % parts.
    %
    % TODO: REFACTOR wrt Bipartite
    %       Essentially, change the constructor (vectorization) and the 
    %       solver call.
    %
    %       Should do by a diff.
       
    properties
        % iterate (actually a vector, as opposed to Ising)
        x
        % constants
        K
        
        % original problem data
        Ys, features
        
        % (vectorized) problem data
         G, y
        
        % options
        opts
        
        % bookeeping
        Ns, Es, beginRow, endRow          
    end
    
    methods
        function th = UnipartiteMatchingDiscLearning(Ys, features, varargin)
            % th = UnipartiteMatching(Ys, features)
            %
            %   Initialize the objective. 
            %   Ys       : M cell vector of sparse matching matrices.
            %   features : M x K cell array. For each k, features{m,k} must
            %              be a symmetric matrix conforming to Ys{m}.
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Parse options (and set defaults)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
            p = inputParser;
            p.KeepUnmatched = true; % Enable extraneous options.

            p.addRequired('Ys', @iscell);
            p.addRequired('features', @iscell);
            p.addParamValue('MaxIter', inf);            
            p.addParamValue('debug', false);            
            p.addParamValue('plotLine', false);
            p.parse(Ys, features, varargin{:});
            th.opts = p.Results;

            th.M = length(Ys);
            [M2, th.K] = size(features);
            
            assert(th.M == M2, 'Y and features must have same number of examples');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Set up problem data and bookkeeping
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            % Keep unvectorized features for evaluating error... although
            % we could just code the prediction routine to take
            % unfeaturized code... oh well.
            th.features = features;
            th.Ys       = Ys;                        
            th.theta = 0.01*randn(size(features,2),1);
            th.G = cell2mat(cellfun(@vecUT, features, 'UniformOutput', false));
            
            th.Ns       = zeros(th.M, 1);
            th.beginRow = zeros(th.M, 1);
            th.endRow   = zeros(th.M, 1);                                    

            % Count number of potential edges in each observation to
            % compute offsets into the vector.
            row = 1;

            % Count nodes and edges of each sample.            
            th.Ns = cellfun(@(X) size(X, 1), Ys);
            % Assume complete graphs.
            th.Es = th.Ns .* (th.Ns - 1) / 2;
            % There is one sample per edge for each sample. This yields the
            % index of the last sample (try it out).
            th.endRow   = cumsum(th.Es);
            % The first row of each sample is one plus the last row of the
            % previous, except for the first one.
            th.beginRow = [1 ; th.endRow(1:end-1) + 1];

            th.y = vertcat(cell2mat(cellfun(@vecUT, Ys, 'UniformOutput', false)));
            %th.GtTimesY = th.G'*th.y;           

            
        end        
      
        function updateParams(th,step,prediction_vector,m)
            inds = ix(th,m);
            feats = th.G(inds,:);
            true_vector = th.y(inds)';
            grad = (true_vector  - prediction_vector)*feats;           
            th.theta = th.theta + step*grad';
            th.allThetas = [th.allThetas th.theta];
        end
        
        
        function err = trainErr(th)
            disp('current iterate')
            n = min(25,length(th.theta));
            th.theta(1:n)'
            disp('avg iterate')
            a = getAvgTheta(th);
            a(1:n)'
            disp('evaluating train error');
%             correctEdges = 0;
%             totalEdges = 0;
%             for m = 1:th.M
%                inds = ix(th,m);
%                pred_vector = solveLPBlock(th, m);
%                true_vector = th.y(inds)';
%                assert(sum(pred_vector) == sum(true_vector));
%                correctEdges = correctEdges + sum(pred_vector & true_vector);
%                totalEdges  = totalEdges + sum(true_vector);
%             end
            %err = 1 - correctEdges/totalEdges;
             err = evalThetaUnipartite(getAvgTheta(th), th.features, th.Ys);
        end
        
        function aT = getTheta(th)
           aT = getAvgTheta(th); 
        end

        function aT = getAvgTheta(th)
            aT = mean(th.allThetas,2);
        end
        function sm = getFeatures(th, m)
            inds = ix(th,m);
            sm =  th.G(inds,:);
        end
        
      
        
        function sm = solveLPBlock(th, m)
            inds = ix(th,m);
            weights = th.G(inds,:)*th.theta;
            sm = vecUniMatch(-weights); 
        end
                
    end
    
    methods (Access = protected)
        function y = ix(th, m)
            y = th.beginRow(m):th.endRow(m);
        end
    end
    
    
end
