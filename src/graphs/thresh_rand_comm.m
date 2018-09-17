function [comm_sz_thresh,comm_num,comm_sz] = thresh_rand_comm(comm,p)
% calculate threshold by community size
% comm is a num_shuff-by-1 cell array, each cell contains cells with
% community identities

comml = vertcat(comm{:});
comm_num = length(comml);
comm_sz = cellfun('length',comml);
bin_range = 0:max(comm_sz);
hist_comm_sz = hist(comm_sz,bin_range);
hist_comm_sz = hist_comm_sz/sum(hist_comm_sz);
cdf_comm_sz = cumsum(hist_comm_sz);
comm_sz_thresh = find(cdf_comm_sz>1-p,1);
        
end