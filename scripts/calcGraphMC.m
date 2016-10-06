function [] = calcGraphMC(param)

expt_name = param.expt_name;
ee = param.ee;
data_path = param.data_path;
result_path_base = param.result_path_base;
ge_type = param.ge_type;

% make xcorr graph
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    
    for e = 1:length(expt_ee)
        
        fprintf('processing %s_%s...\n',expt_name{n},expt_ee{e});
        
        % set path
        model_path = [result_path_base '\' expt_name{n} '\models\'];
        cc_path = [result_path_base '\' expt_name{n} '\cc\'];
        save_path = [result_path_base '\' expt_name{n} '\mc\'];
        
        % load graphs
        crf_graph = load([model_path expt_name{n} '_' expt_ee{e} ...
            '_loopy_best_model_' ge_type '.mat']);
        load([cc_path expt_name{n} '_' expt_ee{e} '_cc_graph.mat']);
        crf_graph = crf_graph.graph;
        
        % find maximal cliques
        mc = maximalCliques(crf_graph);
        mc_crf = cell(size(mc,2),1);
        for ii = 1:size(mc,2)
            mc_crf{ii} = find(mc(:,ii));
        end
        
        mc = maximalCliques(cc_graph);
        mc_cc = cell(size(mc,2),1);
        for ii = 1:size(mc,2)
            mc_cc{ii} = find(mc(:,ii));
        end
        
        % save results
        save([save_path expt_name{n} '_' expt_ee{e} '_mc_crf_' ge_type '.mat'],...
            'mc_crf','-v7.3');
        save([save_path expt_name{n} '_' expt_ee{e} '_mc_cc.mat'],'mc_cc','-v7.3');
        
    end
end

end