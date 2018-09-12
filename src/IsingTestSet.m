classdef IsingTestSet < handle
    properties
        M, st, ixNode, ixEdge, nEdges        
        opts
    end
    
    methods
        function th = IsingTestSet(st, varargin)            
            th.st = st;
            th.nEdges = cellfun(@(x) size(x, 2), st.edges);
            [th.ixNode, th.ixEdge] = computeUVIx(st.Ns, th.nEdges);            
            th.M = size(th.ixNode, 1);
            
            p = inputParser;
            p.addParamValue('maxPredIter', 50);
            p.addParamValue('inferenceOpts', {});
            p.parse(varargin{:});
            th.opts = p.Results;
        end
        
        function err = MAPErr(th, params)    
            nodeCost = -params.F * th.st.Ut;
            edgeCost = -params.G * th.st.Vt;

            mistakes = 0;
            for m = 1:th.M
                ixNm = th.ixNode(m,1):th.ixNode(m,2);
                ixEm = th.ixEdge(m,1):th.ixEdge(m,2);

                [YhatNm, ~, eBelow] = ...
                    solveQPBO(nodeCost(:,ixNm), edgeCost(:,ixEm), th.st.edges{m});

                mistakes = mistakes + sum(vec(YhatNm' ~= th.st.YN(ixNm,:)));
            end

            % th.ixNode(end) = # of nodes.
            err = mistakes / (2*th.ixNode(end));
        end
        
        function [err,likelihood] = MPMErr(th, params)
            mistakes = 0;
            likelihood = 0;
            for m = 1:th.M
                Utm = th.st.Ut(:,th.ixNode(m,1):th.ixNode(m,2));
                Vtm = th.st.Vt(:,th.ixEdge(m,1):th.ixEdge(m,2));

                YNtm = th.st.YN(th.ixNode(m,1):th.ixNode(m,2),:);

                [~, YNtrueFlat] = max(YNtm, [], 2);    

                % trainErr() for IsingPredict *is* MPM.
                tobj  = IsingPredict(params.F, params.G, Utm, Vtm, th.st.Ns(m), th.st.edges{m}, 'YNtrueFlat', YNtrueFlat, th.opts.inferenceOpts{:});
                tprob = BatchFW(tobj, 'MaxIter', th.opts.maxPredIter, 'printInterval', 100, 'linesearch', false, th.opts.inferenceOpts{:});
                tprob.run();
                fprintf('[IsingTest: %d / %d]\n', m, th.M);

                % Adding the negative log partition function
                likelihood = likelihood + tprob.obj.fval();

                mistakes = mistakes + th.st.Ns(m)*tobj.trainErr();
            end
	    
%             likelihood = likelihood / th.M;
            likelihood = (likelihood + th.score(params))/th.M;
            err = mistakes / sum(th.st.Ns);
        end

        function tl = trueLikelihood(th,params)
            tl = th.score(params)/th.M;
%             tl = th.oneValueScore(params)/th.M;
            tl = tl - sum(params.F(1,:)) - sum(params.G(1,:));
            N = size(params.F,2);
            
            % Compute exact log partition function
            [node_weights,edge_weights] = reconstruct_potentials(params);
            [W,~] = construct_graph(edge_weights,N);
            
%             load expt/Parameter_estimation_preknown_edges/Result_20_10000_5_3510.mat model;
%             W = full(triu(model.graph.edge_potentials));
%             node_weights = model.graph.node_potentials;
            
            [~,~,~,~,~,~,logZ] = solveDAI(node_weights,W);
            
            % True likelihood = linear term - logZ
            tl = tl - logZ;
        end
        
        function ll = score(th,params)
            % Adding the linear terms
            thetaN = (params.F*th.st.Ut)';
            thetaE = (params.G*th.st.Vt)';
            linearN = vec(thetaN .* th.st.YN);
            linearE = vec(thetaE .* th.st.YE); 
            
%         ll = (sum(linearN)/length(linearN)) + (sum(linearE)/length(linearE));
            ll = (sum(linearN)) + (sum(linearE));
        
        % Regularization
        %ll = ll - ((th.M * params.lambda)/2)*(frobProd(params.F,params.F) + frobProd(params.G,params.G));
        end
        
       
        function ll = oneValueScore(th,params)
            N = size(params.F,2);
            [node_weights,edge_weights] = reconstruct_potentials(params);
            [W,~] = construct_graph(edge_weights,N);
            
%             load expt/Parameter_estimation_preknown_edges/Result_20_10000_5_3510.mat model;
%             W = full(triu(model.graph.edge_potentials));
%             node_weights = model.graph.node_potentials';
            
            row_mat = repmat([1:N]',1,N);
            col_mat = repmat([1:N],N,1);
            ll = 0;
            for m=1:th.M
                YNtm = th.st.YN(th.ixNode(m,1):th.ixNode(m,2),:);
                YNtm = YNtm(:,2);
                linearE = vec(W.*YNtm(row_mat).*YNtm(col_mat));
                linearN = vec(node_weights'.*YNtm);
                ll = ll + sum(linearE) + sum(linearN);
            end
            
        end
        
    end
    
end

