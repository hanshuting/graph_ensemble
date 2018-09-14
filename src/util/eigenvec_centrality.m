function [node_cent] = eigenvec_centrality(graph)

num_node = size(graph,1);
node_cent = ones(num_node,1);

tol = 1e-3;
count = 0;
while true
    
    count = count+1;
    node_cent_prev = node_cent;
    
    for i = 1:num_node
        node_cent(i) = sum(node_cent_prev.*graph(:,i));
    end
    node_cent = node_cent/max(node_cent);
    
    resd = sum(sqrt((node_cent-node_cent_prev).^2));
%     fprintf('residue is %f at %u\n',resd,count);
    if resd < tol || count>1000
        break;
    end
    
end

end