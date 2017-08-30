function [core,mc_in_core] = find_core_max_clique(graph,num_stim,mc_minsz)

% find nodes that are connected with stim nodes
num_node = size(graph,1);
stim_indx = num_node-num_stim+1:num_node;
stim_assc_node = cell(num_stim,1);
for n = 1:num_stim
    stim_assc_node{n} = find(graph(stim_indx(n),:)~=0);
end

% find maximal cliques
mc = maximalCliques(graph);
mc_nz = cell(size(mc,2),1);
for n = 1:size(mc,2)
    mc_nz{n} = find(mc(:,n));
end

% threshold by size
mc_sz = cellfun('length',mc_nz);
keepIndx = mc_sz>=3;
mc = mc(:,keepIndx);
mc_nz = mc_nz(mc_sz>=mc_minsz);

% find core cells
% core = cell(num_stim+2,1);
% mc_in_core = cell(num_stim+2,1);
% indx = false(num_stim+2,size(mc,2));
core = cell(num_stim,1);
mc_in_core = cell(num_stim,1);
indx = false(num_stim,size(mc,2));
for n = 1:num_stim

%     indx(n,:) = logical(mc(stim_indx(n),:)~=0&mc(stim_indx(setdiff(1:num_stim,n)),:)==0);
    indx(n,:) = logical(mc(stim_indx(n),:)~=0); % SH 10/17
    
    cmc = mc_nz(indx(n,:));
    node_set = unique(cell2mat(cmc));
    num_set = length(node_set);
    mem = histc(cell2mat(cmc),1:node_set(end));
    
    % threshold membership  
%     mem_thresh = zeros(100,1);
%     for ii = 1:100
%         mem_rand = zeros(num_set,1);
%         for jj = 1:length(cmc)
%             rand_vec = zeros(num_set,1);
%             rand_vec(randperm(num_set,length(cmc{jj}))) = 1;
%             mem_rand = mem_rand+rand_vec;
%         end
%         mem_thresh(ii) = mean(mem_rand); %quantile(mem_rand,0.95);
%     end
    mem_thresh = 0.1;
    
    core{n} = setdiff(find(mem>=mean(mem_thresh)),stim_indx);
%     mc_in_core{n} = cmc;
    mc_in_core{n} = cellfun(@(x) setdiff(x,stim_indx),cmc,'uniformoutput',false);
    
end

% indx(n+1,:) = sum(indx(1:n,:),1)==n;
% indx(n+2,:) = sum(indx(1:n,:),1)==0;
% 
% for n = num_stim+1:num_stim+2
%     
%     cmc = mc_nz(indx(n,:));
%     
%     if isempty(cmc)
%         core{n} = [];
%         mc_in_core{n} = [];
%         continue;
%     end
%     
%     node_set = unique(cell2mat(cmc));
%     num_set = length(node_set);
%     mem = histc(cell2mat(cmc),1:node_set(end));
%     
%     % threshold membership  
%     mem_thresh = zeros(100,1);
%     for ii = 1:100
%         mem_rand = zeros(num_set,1);
%         for jj = 1:length(cmc)
%             rand_vec = zeros(num_set,1);
%             rand_vec(randperm(num_set,length(cmc{jj}))) = 1;
%             mem_rand = mem_rand+rand_vec;
%         end
%         mem_thresh(ii) = mean(mem_rand); %quantile(mem_rand,0.95);
%     end
%     mem_thresh = 0.1;
%     
%     core{n} = setdiff(find(mem>=mean(mem_thresh)),stim_indx);
% %     mc_in_core{n} = cmc;
%     mc_in_core{n} = cellfun(@(x) setdiff(x,stim_indx),cmc,'uniformoutput',false);
%     
% end


end

