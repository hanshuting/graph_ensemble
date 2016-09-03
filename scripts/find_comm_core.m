function [comm_core] = find_comm_core(comm)

if isempty(comm)
    comm_core = [];
    return;
end

if size(comm{1},find(size(comm)==1,1))~=1
    comm = comm';
end

node_set = unique(cell2mat(comm));
num_node = length(node_set);
mem = histc(cell2mat(comm),1:node_set(end));

% threshold membership
mem_thresh = zeros(100,1);
for n = 1:100
    mem_rand = zeros(num_node,1);
    for ii = 1:length(comm)
        rand_vec = zeros(num_node,1);
        rand_vec(randperm(num_node,length(comm{ii}))) = 1;
        mem_rand = mem_rand+rand_vec;
    end
    mem_thresh(n) = mean(mem_rand); %quantile(mem_rand,0.95);
end

comm_core = find(mem>=mean(mem_thresh));

end