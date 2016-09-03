function [core_mem] = find_comm_core_by_rand(comm,comm_rand,num_node,p)

num_shuff = length(comm_rand);
comm_randl = vertcat(comm_rand{:});
if ~isempty(comm_randl)
    comm_sz_rand = cellfun('length',comm_randl);
    bin_range = 1:max(comm_sz_rand);
    hist_comm_sz_rand = hist(comm_sz_rand,bin_range);
    hist_comm_sz_rand = hist_comm_sz_rand/sum(hist_comm_sz_rand);
    cdf_comm_sz_rand = cumsum(hist_comm_sz_rand);
    comm_sz_thresh = find(cdf_comm_sz_rand>1-p,1);
    comm_sz = cellfun('length',comm);
    comm_thresh = comm(comm_sz>=comm_sz_thresh);
else
    comm_thresh = comm;
end

% threshold community membership
mem = histc(cell2mat(comm_thresh),1:num_node);
if ~isempty(comm_randl)
    mem_rand = histc(cell2mat(comm_randl),1:num_node)/num_shuff;
    ccdist = fitdist(mem_rand(:),'normal');
    mem_thresh = icdf(ccdist,0.95);
else
    mem_thresh = 0;
end
core_mem = find(mem>mem_thresh);

end