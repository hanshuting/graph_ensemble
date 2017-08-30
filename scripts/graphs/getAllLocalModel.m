function [local_models] = getAllLocalModel(model_collection,coords,ifsave,savename)

savepath = '/vega/brain/users/sh3276/results/luis/fig/m21_sample_d2_d3_local_models/';

model_structs = cellfun(@(x) struct(...
    's_lambda',find(model_collection.s_lambda_sequence==x.s_lambda),...
    'p_lambda',find(model_collection.p_lambda_sequence==x.p_lambda),...
    'test_likelihood',x.test_likelihood,...
    'density',find(model_collection.density_sequence==x.density)...
    ), model_collection.models);

best_by_density = struct([]);
for d = 1:numel(model_collection.density_sequence)
    density_structs = model_structs([model_structs.density] == d);
    dimensions_likelihood_matrix = cell2mat(reshape(struct2cell(density_structs), 4,[])');
    likelihood_grid = full(sparse(dimensions_likelihood_matrix(:,1),...
        dimensions_likelihood_matrix(:,2),...
        dimensions_likelihood_matrix(:,3)));

    best_row = 0;
    best_col = 0;
    best_like = -inf;
    for i = 2:(numel(model_collection.s_lambda_sequence)-1)
        for j = 2:(numel(model_collection.p_lambda_sequence)-1)
            if likelihood_grid(i,j) > likelihood_grid(i-1,j) && ...
                likelihood_grid(i,j) > likelihood_grid(i,j-1) && ...
                likelihood_grid(i,j) > likelihood_grid(i+1,j) && ...
                likelihood_grid(i,j) > likelihood_grid(i,j+1) && ...
                likelihood_grid(i,j) > likelihood_grid(i-1,j-1) && ...
                likelihood_grid(i,j) > likelihood_grid(i-1,j+1) && ...
                likelihood_grid(i,j) > likelihood_grid(i+1,j-1) && ...
                likelihood_grid(i,j) > likelihood_grid(i+1,j+1)
                if likelihood_grid(i,j) > best_like
                    best_row = i;
                    best_col = j;
                    best_like = likelihood_grid(i,j);
                end
            end
        end
     end
     if best_row ~= 0 && best_col ~= 0
         best_by_density(end+1).s_lambda = model_collection.s_lambda_sequence(best_row);
         best_by_density(end).p_lambda = model_collection.p_lambda_sequence(best_col);
         best_by_density(end).like = best_like;
         best_by_density(end).density = model_collection.density_sequence(d);
     end
end

if ~isempty(best_by_density)
    for i = 1:length(best_by_density)
        local_s_lambda = best_by_density(i).s_lambda;
        local_p_lambda = best_by_density(i).p_lambda;
        local_density = best_by_density(i).density;
        local_like = best_by_density(i).like;
        [~, local_model_index] = max(cellfun(@(m) (m.s_lambda==local_s_lambda) ...
            && (m.p_lambda==local_p_lambda) && (m.density==local_density) , model_collection.models));
        local_models{i} = SingleLoopyModel(...
            model_collection.x_train, model_collection.x_test, ...
            model_collection.models{local_model_index}, ...
            model_collection.variable_names);
        
        if ifsave
        % plot model graph
        graph = local_models{i}.structure;
        weight = local_models{i}.theta.edge_potentials;
        h = plotGraphWithXY(graph,coords,weight);
        title(['density ' num2str(local_density) ', likelihood ' num2str(local_like)]);
        saveas(h,[savepath savename '_local_model_' num2str(i) '.fig']);
        end
    end
else
    fprintf('No local extremum found! Searching by likelihood\n');
    model_cell = struct2cell(model_structs);
    lik = squeeze(cell2mat(model_cell(3,1,:)));
    [~,local_model_indx] = sort(lik,'descend');
    for i = 1:round(0.05*length(lik))
        local_models{i} = model_collection.models{local_model_indx(i)};
        if ifsave
        % plot model graph
        graph = local_models{i}.structure;
        weight = local_models{i}.theta.edge_potentials;
        h = plotGraphWithXY(graph,coords,weight);
        title(['density ' num2str(local_models{i}.density) ', likelihood ' num2str(local_models{i}.test_likelihood)]);
        saveas(h,[savepath savename '_local_model_' num2str(i) '.fig']);
        end
    end
end
    
end
