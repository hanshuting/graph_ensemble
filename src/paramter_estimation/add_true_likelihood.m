IsingValOpts = { 'maxPredIter', 5000};
for i=1:length(model.theta)
    val = DataCRF({model.structure(1).graph},{model.XVal});
    valObj = IsingTestSet(val, IsingValOpts{:},'inferenceOpts',...
                    {'printComputedDualityGap', false,'reweight', 1});
    model.theta(i).trueValLikelihood = valObj.trueLikelihood(model.theta(i).theta);
    
    params = model.theta(i).theta;
    
end


% Test case with 3 nodes
% N = 3;
% graph = [0 1 1;1 0 1;0 0 1];
% sample = [0 1 1];
% test = DataCRF({graph},{sample});
% 
% % Parameters
% params.F = 5*rand(2,3);
% params.G = 4*rand(4,3);
% 
% % Make them singular
% [node_weights,edge_weights] = reconstruct_potentials(params);
% [W,~] = construct_graph(edge_weights,N);
% 
% % True logZ
% [~,~,~,~,~,~,logZ] = solveDAI(node_weights,W);
% 
% % Linear term 1
% tl1 = params.F(1,1) + params.F(2,2) + params.F(2,3) + ...
%     params.G(2,1) + params.G(2,2) + params.G(4,3);
% 
% % Linear term 2
% tl2 = node_weights(2) + node_weights(3) + edge_weights(3);
% 
% tl1-tl2
% sum(params.F(1,:)) + sum(params.G(1,:))
% 
% % My way
% row_mat = [1 1 1;2 2 2;3 3 3];
% col_mat = [1 2 3;1 2 3;1 2 3];
% YN = [1 0;0 1;0 1];
% YN = YN(:,2);
% linearE = vec(W.*YN(row_mat).*YN(col_mat));
% linearN = vec(node_weights'.*YN);
% tl3 = sum(linearE) + sum(linearN);
% tl1,tl2,tl3
