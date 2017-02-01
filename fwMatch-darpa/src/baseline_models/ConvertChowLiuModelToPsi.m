function [CPT,edge_list] = ConvertChowLiuModelToPsi(model)

edge_list = zeros(model.Nedges,2);
CPT = cell(model.Nedges,1);
num_nodes = model.Nnodes;

ne = 1;
for n = 1:num_nodes
    parent = model.pa(n);
%     children = model.edges(find(model.edges(:,2) == n),1);
    if parent == 0 
        continue; 
    end
    if parent == 1
        CPT{ne} = model.CPDs{n} .* repmat(model.CPDs{1}',2,1);
    else
        CPT{ne} = model.CPDs{n};
    end
    edge_list(ne,:) = [n parent];
    ne = ne + 1;
end

end