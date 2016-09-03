function [] = randGraphMC(param)

expt_name = param.expt_name;
ee = param.ee;
num_rand = param.num_rand;
data_path = param.data_path;
result_path_base = param.result_path_base;
cc_type = param.cc_type;
ge_meth = param.ge_meth;

%% cc model
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    
    for e = 1:length(expt_ee)
        
        fprintf('processing %s_%s...\n',expt_name{n},expt_ee{e});
        model_path = [result_path_base expt_name{n} '\models\'];
        cc_path = [result_path_base expt_name{n} '\cc_' cc_type '\']; 
        save_path = [result_path_base expt_name{n} '\rand_graph_' ge_meth '\'];
        
        crf_graph = load([model_path expt_name{n} '_' expt_ee{e} '_loopy_best_model.mat']);
        crf_graph = crf_graph.graph;
        load([cc_path expt_name{n} '_' expt_ee{e} '_cc_graph.mat']);
        
        %% CC
        % generate random graphs
        num_node = size(cc_graph,1);
        num_edge = sum(cc_graph(:))/2;
        cc_rand_graph = zeros(num_node,num_node,num_rand);
        
        % make random graphs
        for ii = 1:num_rand
            cc_rand_graph(:,:,ii) = mkRandomGraph(num_node,num_edge);
        end
        
        % find mc
        mc_cc_rand = cell(num_rand,1);
        for ii = 1:num_rand
            mc = maximalCliques(cc_rand_graph(:,:,ii));
            mc_cc_rand_nz = cell(size(mc,2),1);
            for jj = 1:size(mc,2)
                mc_cc_rand_nz{jj} = find(mc(:,jj));
            end
            mc_cc_rand{ii} = mc_cc_rand_nz;
        end
        
        %% CRF
        % in case I want to switch to 'on' edges only
        if strcmp(param.ge_meth,'on')
            num_node = size(crf_graph,1);
            num_edge = sum(sum(logical(tril(crf_graph))));
            edge_list = zeros(num_edge,2);
            [edge_list(:,2),edge_list(:,1)] = find(tril(crf_graph));
            G_crf = zeros(num_node,num_node);
            for ii = 1:size(G_crf,2)
                node_1 = edge_list(ii,1);
                node_2 = edge_list(ii,2);
                G_crf(node_1,node_2) = G_crf(4,ii);%-G(3,i)-G(2,i);
                G_crf(node_2,node_1) = G_crf(4,ii);%-G(3,i)-G(2,i);
            end
            G_crf = G_crf>0;
        elseif strcmp(param.ge_meth,'full')
            G_crf = crf_graph;
        else
            error('Undefined graph edge option ''%s''\n',param.ge_meth);
        end
        
        % generate random graphs
        num_node = size(G_crf,1);
        num_edge = sum(G_crf(:))/2;
        crf_rand_graph = zeros(num_node,num_node,num_rand);
        
        % make random graphs
        for ii = 1:num_rand
            crf_rand_graph(:,:,ii) = mkRandomGraph(num_node,num_edge);
        end
        
        % find mc
        mc_crf_rand = cell(num_rand,1);
        for ii = 1:num_rand
            mc = maximalCliques(crf_rand_graph(:,:,ii));
            mc_crf_rand_nz = cell(size(mc,2),1);
            for jj = 1:size(mc,2)
                mc_crf_rand_nz{jj} = find(mc(:,jj));
            end
            mc_crf_rand{ii} = mc_crf_rand_nz;
        end
        
        %% save results
        save([save_path expt_name{n} '_' expt_ee{e} '_mc_cc_rand.mat'],...
            'mc_cc_rand','cc_rand_graph','-v7.3');
        save([save_path expt_name{n} '_' expt_ee{e} '_mc_crf_rand.mat'],...
            'mc_crf_rand','crf_rand_graph','-v7.3'); 
        
    end
end

end
