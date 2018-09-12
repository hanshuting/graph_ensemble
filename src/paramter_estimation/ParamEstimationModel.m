classdef ParamEstimationModel 
    
    properties
        % Original graph if for synthetic data, for real data empty matrix
        graph
        
        % training and Validation samples
        XTrain, XVal
        
        % Struct array for the estimated graphs form structural learning
        % with different parameters
        % Fields : lambda, tolr and graph (NxN adjacency matrix)
        structure
        
        % Struct array of Estimated paramters theta 
        % Fields (learned theta): theta (which has fields F and G ), 
        % Fields (training evaluation) : trainErr, trainMAPErr, trainMPMErr
        % Fields (train likelihood) : trainLikelihoodREgularized (Eq. 6 in
        % Appendix), trainLikelihood (First equation)
        % Fields (test evaluation) : valMAPErr, valMPMErr, valLikelihood
        % Fields (parameters) : slambda, tolr, plambda 
        theta
        
        % parameters
        % slambda - regularizer for structural learning
        % tolr - thereshold for structural learning
        % plambda - regularizer for parameter estimation
        slambda, tolr, plambda
        
        % flags
        add_original_graph
        
    end
    
    methods
        function model = ParamEstimationModel(XTrain,XVal,varargin)
            % model = ParamEstimationModel(XTrain,XVal,varargin)
            % 
            % Inputs (Required):
            % 1. XTrain : #training_samples x #nodes binary matrix of samples
            % 2. XVal : #validation_samples x #nodes binary matrix of samples
            % Inputs (Optional) :
            % 3. graph : an object with two fields, edge_potentials and
            % node_potentials
            % Inputs (Optional and taken as parameters) :
            % 4. numSlambda : number of structural learning regularizers to
            % consider in model selection
            % 5. minSlambda : minimum of log_10(lambda for structural learning)
            % 6. maxSlambda : maximum of log_10(lambda for structural learning)
            % 7. numPlambda, minPlmabda and maxPlambda same as 4,5,6 for
            % regularizer for paramter estimation
            % 8. numTolr, minTolr and maxTolr same as 4,5,6 except we
            % consider linear space instead of log space (no log) for
            % tolerance in structural leanring
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input Parser
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            parser = inputParser;
            is_logical_matrix = @(x) ismatrix(x) && islogical(x) && ~isempty(x);
            ispositive = @(x) x >= 0;
            parser.addRequired('XTrain', is_logical_matrix);
            parser.addRequired('XVal', is_logical_matrix);
            
            parser.addOptional('graph',[],@(x) isfield(x,'edge_potentials') && isfield(x,'node_potentials'));
            
            parser.addParamValue('numSlambda',0,ispositive); 
            parser.addParamValue('minSlambda',1e-04,ispositive);
            parser.addParamValue('maxSlambda',0.9,ispositive);
            
            parser.addParamValue('numPlambda',20,ispositive); 
            parser.addParamValue('minPlambda',1e-03,ispositive);
            parser.addParamValue('maxPlambda',1e09,ispositive);
            
            parser.addParamValue('numTolr',0,ispositive); 
            parser.addParamValue('minTolr',0.05,@isscalar);
            parser.addParamValue('maxTolr',0.5,@isscalar);
            
            parser.addParamValue('add_original_graph',false,@isscalar);
            
            parser.parse(XTrain,XVal,varargin{:});
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Model initialization
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            model.XTrain = parser.Results.XTrain;
            model.XVal = parser.Results.XVal;
            model.graph = parser.Results.graph;
            % lambdas are samples in a logspace
            model.slambda = logspace(log10(parser.Results.minSlambda),log10(parser.Results.maxSlambda),parser.Results.numSlambda);
            model.plambda = logspace(log10(parser.Results.minPlambda),log10(parser.Results.maxPlambda),parser.Results.numPlambda);
            % tolerance is samples in a linear space
            model.tolr = linspace(parser.Results.minTolr,parser.Results.maxTolr,parser.Results.numTolr);
            % compare the original graph among the structures learned
            model.add_original_graph = parser.Results.add_original_graph;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Structural Learning
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function model = StructuralLearning(model)
            % Serial implementation
            % Total numSlambda x numTolr structures leanrned
            % For the other inputs refer to infer_structure.m
            
            for i=1:numel(model.slambda)
                [learned_structures,real_valued_graph] = infer_structures_with_given_densities(model.XTrain, model.slambda(i),model.tolr);
                for j=1:numel(model.tolr)
                    fprintf('Structural Learning with Lambda = %e and Tolerance = %e\n',model.slambda(i),model.tolr(j));
                    indx = sub2ind([numel(model.slambda) numel(model.tolr)],i,j);
                    model.structure(indx).lambda = model.slambda(i);
                    model.structure(indx).tolr = model.tolr(j);
                    model.structure(indx).graph = learned_structures{1};
                    model.structure(indx).realGraph = real_valued_graph;
                    model.structure(indx).max_degree = max(sum(model.structure(indx).graph));
                    if ~isempty(model.graph)
                        % AUC
                        [~,~,~,model.structure(indx).auc] = perfcurve(vecUT((abs(model.graph.edge_potentials) > 0), true),...
                        vecUT(model.structure(indx).realGraph, true), 1);
                        fprintf('AUC : %e\n',model.structure(indx).auc);
                        % TP, FP, FN
                        model.structure(indx).eval = ...
                            eval_graphs((full(model.graph.edge_potentials) ~= 0),model.structure(indx).graph);
                        fprintf('TP : %d, FP : %d, FN : %d Edges Predicted : %d\n',model.structure(indx).eval.TP,...
                            model.structure(indx).eval.FP,model.structure(indx).eval.FN,sum(sum(model.structure(indx).graph ~= 0))/2);
                    end
                end
            end
            % Add the original graph
            if ~isempty(model.graph) && model.add_original_graph
                indx = length(model.structure) + 1;
                model.structure(indx).lambda = Inf;
                model.structure(indx).tolr = Inf;
                model.structure(indx).graph = (abs(full(model.graph.edge_potentials)) > 0);
		model.structure(indx).max_degree = max(sum(model.structure(indx).graph));
                model.structure(indx).realGraph = full(model.graph.edge_potentials);
                model.structure(indx).auc = 1;
                model.structure(indx).eval = eval_graphs((full(model.graph.edge_potentials) ~= 0),...
                    (full(model.graph.edge_potentials) ~= 0));
                fprintf('TP : %d, FP : %d, FN : %d\n',model.structure(indx).eval.TP,...
                            model.structure(indx).eval.FP,model.structure(indx).eval.FN);
            end
            
            % Parallel implementation
            
        end
        
        function model = StructuralLearningPreknownEdges(model,varargin)
            % Serial implementation
            % Total numSlambda x numTolr structures leanrned
            % For the other inputs refer to infer_structure.m
            
            if ~isempty(model.graph)
                % compute number of edges in the original graph
                model.graph.numEdges = full(sum(sum(model.graph.edge_potentials ~= 0))/2);
                fprintf('Number of edges in the original graph : %d\n',model.graph.numEdges);
            end
            
            for indx=1:numel(model.slambda)
                fprintf('Structural Learning with Lambda = %e\n',model.slambda(indx));
               
                model.structure(indx).lambda = model.slambda(indx);
                model.structure(indx).tolr = 0;
                [model.structure(indx).graph,model.structure(indx).realGraph] = ...
                        infer_structure(model.XTrain, model.slambda(indx),varargin{:},'preknown_edges',model.graph.numEdges);
                model.structure(indx).max_degree = max(sum(model.structure(indx).graph));
                
                if ~isempty(model.graph)
                    % AUC
                    [~,~,~,model.structure(indx).auc] = perfcurve(vecUT((abs(model.graph.edge_potentials) > 0), true),...
                        vecUT(model.structure(indx).realGraph, true), 1);
                    fprintf('AUC : %e\n',model.structure(indx).auc);
                    % TP, FP, FN
                    model.structure(indx).eval = ...
                            eval_graphs((full(model.graph.edge_potentials) ~= 0),model.structure(indx).graph);
                    fprintf('TP : %d, FP : %d, FN : %d Edges Predicted : %d\n',model.structure(indx).eval.TP,...
                            model.structure(indx).eval.FP,model.structure(indx).eval.FN,sum(sum(model.structure(indx).graph ~= 0))/2);
                end
            end
            % Add the original graph
            if ~isempty(model.graph) && model.add_original_graph
                indx = length(model.structure) + 1;
                model.structure(indx).lambda = Inf;
                model.structure(indx).tolr = Inf;
                model.structure(indx).graph = (abs(full(model.graph.edge_potentials)) > 0);
                model.structure(indx).max_degree = max(sum(model.structure(indx).graph));
                model.structure(indx).realGraph = full(model.graph.edge_potentials);
                model.structure(indx).auc = 1;
                model.structure(indx).eval = eval_graphs((full(model.graph.edge_potentials) ~= 0),...
                    (full(model.graph.edge_potentials) ~= 0));
                fprintf('TP : %d, FP : %d, FN : %d Edges Predicted : %d\n',model.structure(indx).eval.TP,...
                            model.structure(indx).eval.FP,model.structure(indx).eval.FN,sum(sum(model.structure(indx).graph ~= 0))/2);
            end
            
            % Parallel implementation
            
        end
        
        function model = StructuralLearningChowLiu(model)
            indx = 1;
            fprintf('Structure learning using Chow-Liu\n');
            model.structure(indx).lambda = Inf;
            model.structure(indx).tolr = 0;
            chowliu_struct = treegmFit(model.XTrain);
            model.structure(indx).graph = (abs(chowliu_struct.adjmat) > 0);
            model.structure(indx).max_degree = max(sum(model.structure(indx).graph));
            model.structure(indx).realGraph = chowliu_struct.adjmat;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Parameter estimation training
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function model = ParamterEstimationtraining(model,IsingOpts,BCFWOpts,IsingValOpts)
            % Both parallel and serial implementation possible
            % Serial implementation
            for indx2=1:length(model.structure)
                if indx2 ~= length(model.structure) || model.add_original_graph == false 
                    [i,j] = ind2sub([numel(model.slambda) numel(model.tolr)],indx2);
                    current_slambda = model.slambda(i);
                    current_tolr = model.tolr(j);
                else
                    current_slambda = Inf;
                    current_tolr = 0;
                end
                
                train = DataCRF({model.structure(indx2).graph},{model.XTrain});
                trainObj = IsingTestSet(train, IsingValOpts{:});
                for k=1:numel(model.plambda)
                    if indx2 ~= length(model.structure)
                        indx3 = sub2ind([numel(model.slambda) numel(model.tolr) numel(model.plambda)],i,j,k);
                    else
                        indx3 = numel(model.slambda)*numel(model.tolr)*numel(model.plambda) + k;
                    end
                    
                        
                    fprintf('Parameter Training : Struc. Lambda = %e Tolerance = %e Lambda = %e\n',...
                            current_slambda,current_tolr,model.plambda(k));
                    
                    obj = Ising(train.YN, train.YE, train.Ut, train.Vt, train.Ns, ...
                            train.edges, model.plambda(k),IsingOpts{:},'reweight',2/model.structure(indx2).max_degree);
                    prob = BCFW(obj, BCFWOpts{:});
                    prob.run();
                        
                    model.theta(indx3).theta = prob.obj.computeParams();
                    model.theta(indx3).slambda = current_slambda;
                    model.theta(indx3).tolr = current_tolr;
                    model.theta(indx3).plambda = model.plambda(k);
                    model.theta(indx3).trainMPMErr = prob.obj.trainErr();
                    model.theta(indx3).trainMAPErr = prob.obj.trainMAPErr();
                    model.theta(indx3).trainLikelihood = ...
                            prob.obj.score() / size(model.XTrain,1);
                    model.theta(indx3).structure_id = indx2; 
                    
                   
                end
            end
            
            % Parallel implementation
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Parameter estimation prediction on validation set
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function model = ParamterEstimationPrediction(model,IsingValOpts)
            % Both parallel and serial implementation possible
            % Serial implementation
            for indx2=1:length(model.structure)
                if indx2 ~= length(model.structure) || model.add_original_graph == false 
                    [i,j] = ind2sub([numel(model.slambda) numel(model.tolr)],indx2);
                    current_slambda = model.slambda(i);
                    current_tolr = model.tolr(j);
                else
                    current_slambda = Inf;
                    current_tolr = 0;
                end
                val = DataCRF({model.structure(indx2).graph},{model.XVal});
                valObj = IsingTestSet(val, IsingValOpts{:},'inferenceOpts',...
                    {'printComputedDualityGap', false,'reweight', 2/model.structure(indx2).max_degree});
                for k=1:numel(model.plambda)
                    if indx2 ~= length(model.structure)
                        indx3 = sub2ind([numel(model.slambda) numel(model.tolr) numel(model.plambda)],i,j,k);
                    else
                        indx3 = numel(model.slambda)*numel(model.tolr)*numel(model.plambda) + k;
                    end
                    
                    fprintf('Parameter Prediction on Validation set : Struc. Lambda = %e Tolerance = %e Lambda = %e\n',...
                            current_slambda,current_tolr,model.plambda(k));
                        
                    model.theta(indx3).valMAPErr = valObj.MAPErr(model.theta(indx3).theta);
                    [model.theta(indx3).valMPMErr,model.theta(indx3).valLikelihood]...
                            = valObj.MPMErr(model.theta(indx3).theta);
                    fprintf('Computing True likelihood.....\n');
                    model.theta(indx3).trueValLikelihood = valObj.trueLikelihood(model.theta(indx3).theta);
                    fprintf('True Likelihood is %f\n',model.theta(indx3).trueValLikelihood);
                end
                
            end
            
            % Parallel implementation
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% choose the best model based on validation likelihood
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [opt_theta] = GetOptimalTheta(model)
%             likelihood = [model.theta(1:length(model.theta)-length(model.plambda)).valLikelihood];
            likelihood = [model.theta.trueValLikelihood];
            [ML,MLI] = max(likelihood);
            fprintf('Max Likelihood : %f\n',ML);
            
            opt_theta = model.theta(MLI);
            % Construct the adjacency matrix fir the estimated graph
            % reconstruction via additivity theorem
            [~,edge_weights] = reconstruct_potentials(opt_theta.theta);
            [~,opt_graph] = construct_graph(edge_weights,size(model.XTrain,2));
%             opt_graph(find(opt_graph ~= 0)) = opt_graph(find(opt_graph ~= 0)));
            opt_theta.graph = opt_graph;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Visualization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function visualize_all(model)
            for i=1:length(model.theta)
                model.visualize_theta(model.theta(i));
                pause(0.4);
            end
        end
        
        function visualize_true(model)
            if ~isempty(model.graph)
                tindices = length(model.theta)-length(model.plambda)+1:length(model.theta);
                [~,MLI] = max([model.theta(tindices).valLikelihood]);
                MLI = MLI + length(model.theta)-length(model.plambda);
                model.visualize_theta(model.theta(MLI));
                for i=1:length(tindices)
                    model.visualize_theta(model.theta(tindices(i)));
                    pause(1.4);
                end
            end
        end
        
        function visualize_optimum(model)
            % Get optimal theta
            opt_theta = model.GetOptimalTheta();
            if ~isempty(model.graph)
                model.visualize_theta(opt_theta);
            
                % Visualize the original and estimated graph
                orig_graph = full(abs(full(model.graph.edge_potentials)) > 0);
                plot_graph(orig_graph);        
            end
            plot_graph(opt_theta.graph);
        end
        
        function visualize_theta(model,opt_theta)
            
            ftitle = sprintf('slambda = %e, tolerance = %e, plambda = %e\n',...
                opt_theta.slambda,opt_theta.tolr,opt_theta.plambda);
            
            % True edge potentials
            ThetaTrueE = all_edge_potentials(model.graph.edge_potentials,size(model.XTrain,2));
            
            % Reconstruction of edge and node potentials using additivity
            % theorem
            [node_weights,edge_weights] = reconstruct_potentials(opt_theta.theta);
            
            % Plot scatter plots
            [I,J] = getReverseId(size(model.XTrain,2));
%             plot_scatter_plots(abs(model.graph.node_potentials),exp(node_weights),...
%                     abs(ThetaTrueE),exp(edge_weights),1,ftitle,[I J]);
            plot_scatter_plots(model.graph.node_potentials,node_weights,...
                    ThetaTrueE,edge_weights,1,ftitle,[I J]);
                
%             ThetaHatE = all_edge_potentials(model.structure(opt_theta.structure_id).realGraph...
%                     ,size(model.XTrain,2));
%             plot_scatter_plots(model.graph.node_potentials,node_weights,...
%                     ThetaHatE,edge_weights,2,ftitle);
        end
        
        function parameter_plots(model)
            
            %% Figure 1, slambda Vs AUC and valLikelihood
%             figure(1),
%             opt_theta = model.GetOptimalTheta();
%             
%             lindices = find([model.structure(1:length(model.structure)-1).tolr] == opt_theta.tolr); 
%             lambdas = [model.structure(lindices).lambda];
%             aucs = [model.structure(lindices).auc];
%             
%             tindices = intersect(find([model.theta(1:length(model.theta)).plambda] == opt_theta.plambda),...
%                 find([model.theta(1:length(model.theta)).tolr] == opt_theta.tolr));
%             likelihoods = [model.theta(tindices(1:end-1)).valLikelihood];
%             plot_train_test(2,plambdas,ytrain,ytest,type,legend);

%             print('-depsc','-tiff','-r300','fig/structural_learning_auc_likelihood')

            %% Figure 2, plambda Vs training and test likelihood
           
%             opt_theta = model.GetOptimalTheta();
%             
%             if isempty(model.graph)
%                 tindices = intersect(find([model.theta(1:length(model.theta)).slambda] == opt_theta.slambda),...
%                     find([model.theta(1:length(model.theta)).tolr] == opt_theta.tolr));
%             else
%                 tindices = length(model.theta)-length(model.plambda)+1:length(model.theta);
%             end
%             tindices = 1:length(model.theta);
%             plambdas = [model.theta(tindices).plambda];
%             ytest = [model.theta(tindices).trueValLikelihood];
%             ytrain = [model.theta(tindices).valLikelihood];
%             
%             plot_train_test(2,plambdas,ytrain,ytest,'semilogx','lambda','Log Likelihood','max');
% %             print('-depsc','-tiff','-r300','fig/parameter_estimation_lambda_train_test_likelihhod')
%             
%             %% Figure 3, MAP Error
%             ytrain = [model.theta(tindices).trainMAPErr];
%             ytest = [model.theta(tindices).valMAPErr];
%             plot_train_test(3,plambdas,ytrain,ytest,'semilogx','lambda','MAP Error','min');
% 
%             %% Figure 4, MPM Error
%             
%             ytrain = [model.theta(tindices).trainMPMErr];
%             ytest = [model.theta(tindices).valMPMErr];
%             plot_train_test(4,plambdas,ytrain,ytest,'semilogx','lambda','MPM Error','min');
            
            %% FIGURE 5, surf plot for both the lambdas
            opt_theta = model.GetOptimalTheta();
            opt_plambda = opt_theta.plambda;
            opt_slambda = opt_theta.slambda;
            opt_tolr = opt_theta.tolr;
            
            % surf plot
            Z = size(length(model.plambda),length(model.slambda));
            for i=1:length(model.theta)
                Z(find(model.plambda == model.theta(i).plambda),...
                    find(model.slambda == model.theta(i).slambda)) = model.theta(i).trueValLikelihood;
            end
            figure(5),hold on
            set(gca, 'XScale', 'log', 'YScale', 'log'),grid on
            xlabel('Structural Learning Lambda'),ylabel('Parameter Estimation Lambda')
            surf(model.slambda,model.plambda,Z);
     
            % plambda vs test log likelihood
            figure(6),
            tindices = intersect(find([model.theta.slambda] == opt_slambda),...
                find([model.theta.tolr] == opt_tolr));
            [sp,spi] = sort([model.theta(tindices).plambda]);
            stl = [model.theta(tindices(spi)).trueValLikelihood];
            h6 = semilogx(sp,stl);
            set(h6,'LineWidth',2);set(h6,'Color','b');
            ylabel('Test Log Likelihood');xlabel('\lambda (Parameter estimation)');grid on;
            
            figure(7),
            plot_train_test(7,sp,[model.theta(tindices(spi)).valLikelihood],...
                stl,'semilogx','lambda','Log Likelihood','max');
            
            figure(8),
            tindices = intersect(find([model.theta(1:end).plambda] == opt_plambda),...
                find([model.theta(1:end).tolr] == opt_tolr));
            tindices(3) = [];
            h8 = semilogx([model.theta(tindices).slambda],[model.theta(tindices).trueValLikelihood]);
            set(h8,'LineWidth',2);set(h8,'Color','b');
            ylabel('Test Log Likelihood');xlabel('\lambda (Structure Learning)');grid on;
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Merge Models
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function model = merge_models(model,model2)
            if ~isempty(model2.graph) && model2.add_original_graph 
                model.graph = model2.graph;
                model.add_original_graph = true;
            end
            % change the struct id of the model
            for i=1:length(model2.theta)
                model2.theta(i).structure_id = model2.theta(i).structure_id + length(model.structure);
            end
            model.theta = [model.theta, model2.theta];
            model.structure = [model.structure, model2.structure];
            model.plambda = unique([model.plambda, model2.plambda]);
            model.slambda = unique([model.slambda, model2.slambda]);
            model.tolr = unique([model.tolr, model2.tolr]);
            model.tolr = unique([model.tolr, model2.tolr]);
        end
        
    end
    
end

