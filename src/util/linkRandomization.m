function [rand_graph] = linkRandomization(graph)
% Generate random graph using link randomization. In each step two links 
% are selected randomly, and then one of the endpoints of the links are 
% swapped.
% See S. Maslov and K. Sneppen: Science 296, 910 (2002)
% 
% Shuting Han, 2016

graph = graph-tril(graph);
num_edge = sum(sum(logical(graph)));

edge_list = zeros(num_edge,3);
[edge_list(:,1),edge_list(:,2)] = find(graph);

% number of shuffle time for each node
thresh_shuff = 10;

rand_graph = graph;
while true
    
    % randomly select edge pairs
    indx = randperm(num_edge,2);
    edge1 = edge_list(indx(1),1:2);
    edge2 = edge_list(indx(2),1:2);
    new_edge1 = sort([edge1(1),edge2(1)],'ascend');
    new_edge2 = sort([edge1(2),edge2(2)],'ascend');
    
    % check if swapped edges exist
    if ismember(new_edge1,edge_list(:,1:2),'rows') || ...
            ismember(new_edge2,edge_list(:,1:2),'rows')
%         fprintf('edge exists\n');
        continue;
    else
        rand_graph(edge1(1),edge1(2)) = 0;
        rand_graph(edge2(1),edge2(2)) = 0;
        rand_graph(new_edge1(1),new_edge1(2)) = 1;
        rand_graph(new_edge2(1),new_edge2(2)) = 1;
        edge_list(indx(1),1:2) = new_edge1;
        edge_list(indx(2),1:2) = new_edge2;
        edge_list(indx,3) = edge_list(indx,3)+1;
    end
    
    % check if all edges have been shuffled for enough times
%     fprintf('%u\n',sum(rand_graph(:)));
    if all(edge_list(:,3)>thresh_shuff)
        
        % make it a full matrix
        rand_graph = rand_graph+rand_graph';
        
        return;
    end
    
end

end