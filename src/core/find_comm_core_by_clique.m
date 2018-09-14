function [comm_core] = find_comm_core_by_clique(comm,graph,num_stim)

if isempty(comm)
    comm_core = [];
    return;
end

if size(comm{1},find(size(comm)==1,1))~=1
    comm = comm';
end

% find nodes that are connected with stim nodes
num_node = size(graph,1);
stim_indx = num_node-num_stim+1:num_node;
stim_assc_node = cell(num_stim,1);
for n = 1:num_stim
    stim_assc_node{n} = find(graph(stim_indx(n),:)~=0);
end

% find communities
comm_core = cell(num_stim,1);
for n = 1:num_stim

    indx = cellfun(@(x) ~isempty(intersect(x,stim_assc_node{n})),comm);
    indx = indx&cellfun(@(x) isempty(intersect(x,cell2mat(stim_assc_node...
        (setdiff(1:num_stim,n))))),comm);
    
    cr_comm = comm(indx);
    
    node_set = unique(cell2mat(cr_comm));
    num_set = length(node_set);
    mem = histc(cell2mat(cr_comm),1:node_set(end));
    
    % threshold membership
    mem_thresh = zeros(100,1);
    for ii = 1:100
        mem_rand = zeros(num_set,1);
        for jj = 1:length(cr_comm)
            rand_vec = zeros(num_set,1);
            rand_vec(randperm(num_set,length(cr_comm{jj}))) = 1;
            mem_rand = mem_rand+rand_vec;
        end
        mem_thresh(ii) = mean(mem_rand); %quantile(mem_rand,0.95);
    end
%     mem_thresh = 0.1;
    
    comm_core{n} = find(mem>=mean(mem_thresh));

end

end

