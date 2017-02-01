function [LLs] = calc_mc_pot(data,mc,node_pot,edge_pot,logZ,num_stim)
% calculate the clique potentials
% data is num_frame-by_num_neuron
% mc is a cell array of cliques

num_clique = length(mc);
[num_frame,num_node] = size(data);
stim_indx = num_node+1:num_node+num_stim;

LLs = cell(num_clique,1);
for n = 1:num_clique
    
    cc = mc{n};
    LL_stim = zeros(num_stim,num_frame);
    
    for ii = 1:num_frame    
        for jj = 1:num_stim
            
            stim_vec = zeros(1,num_stim);
            stim_vec(jj) = 1;
            
            cc_vec = false(num_node+num_stim,1);
            cc_vec([cc;stim_indx']) = 1;
%             cc_np = node_pot.*cc_vec;
%             cc_ep = edge_pot;
%             cc_ep(~cc_vec,:) = 0;
%             cc_ep(:,~cc_vec) = 0;
            cc_np = node_pot;
            cc_ep = edge_pot;
            cc_data = [data(ii,:).*cc_vec(1:num_node)',stim_vec];
%             cc_np = node_pot(cc_vec);
%             cc_ep = edge_pot(cc_vec,cc_vec);
%             cc_data = [data(ii,cc),stim_vec];
            
            LL_stim(jj,ii) = compute_avg_log_likelihood(cc_np,cc_ep,logZ,cc_data);
            
        end
        
        LLs{n} = LL_stim;
%         LLs(n,ii) = compute_avg_log_likelihood(node_pot(cc),...
%                 edge_pot(cc,cc),logZ,data(ii,cc));

    end
end


end