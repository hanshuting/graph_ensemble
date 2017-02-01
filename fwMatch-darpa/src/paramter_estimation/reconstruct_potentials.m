function [node_weights,edge_weights] = reconstruct_potentials(theta)

% number of nodes
num_nodes = size(theta.F,2);

% Estimated edge potentials 4+1-2-3 (11+00-01-10)
edge_weights = theta.G(4,:) + theta.G(1,:) - theta.G(2,:)  - theta.G(3,:);
            
% Estimated node potentials 
[I,J] = getReverseId(num_nodes);
node_weights = zeros(1,num_nodes);
for i=1:size(theta.G,2)
    node_weights(I(i)) = node_weights(I(i)) + theta.G(2,i) - theta.G(1,i);
    node_weights(J(i)) = node_weights(J(i)) + theta.G(3,i) - theta.G(1,i);
end
node_weights = node_weights + theta.F(2,:) - theta.F(1,:);

end