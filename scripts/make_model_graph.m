
load('./results/model_collection.mat');
best_model = getBestModel(model_collection);
graph = best_model.structure;
node_pot = exp(best_model.theta.node_potentials-best_model.theta.logZ);
edge_pot = exp(best_model.theta.edge_potentials-best_model.theta.logZ);

h = plotGraphWithXY(graph,Coord_active,edge_pot);
%saveas(h,['/vega/brain/users/sh3276/results/luis/best_model.fig']);
