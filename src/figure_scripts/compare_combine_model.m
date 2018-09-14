function [] = compare_combine_model(param)

% parameters
expt_name = param.expt_name;
ee = param.ee;
indv_ee = param.indv_ee;
num_shuff = param.num_shuff;
k = param.k;
p = param.p;
cc_type = param.cc_type;
comm_type = param.comm_type;
data_path = param.data_path;
fig_path = param.fig_path.graph_prop;
save_path = param.result_path.stats_path;
result_path_base = param.result_path_base;
savestr = param.savestr;
ccode_path = param.ccode_path;

load(ccode_path);

for n = 1:length(expt_name)
    
    expt_ee = ee{n}{1};
    expt_indv_ee = indv_ee{n};
    
    model_path = [result_path_base expt_name{n} '\models\']; 
    comm_path = [result_path_base expt_name{n} '\comm_' comm_type '\'];
    
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    load([comm_path expt_name{n} '_' expt_ee '_loopy_comm_k_' num2str(k) '.mat']);
    load([comm_path expt_name{n} '_' expt_ee '_loopy_rand_graph_comm_k_' num2str(k) '.mat']);
    best_model = load([model_path expt_name{n} '_' expt_ee '_loopy_best_model.mat']);
    num_node = size(best_model.graph,1);
    core_mem = find_comm_core(comm_loopy,comm_loopy_rand,num_node,p);
    
    core_mem_indv = cell(length(expt_indv_ee),1);
    for e = 1:length(expt_indv_ee)
        
        load([comm_path expt_name{n} '_' expt_indv_ee{e} '_loopy_comm_k_' num2str(k) '.mat']);
        load([comm_path expt_name{n} '_' expt_indv_ee{e} '_loopy_rand_graph_comm_k_' num2str(k) '.mat']);
        best_model = load([model_path expt_name{n} '_' expt_indv_ee{e} '_loopy_best_model.mat']);
        num_node = size(best_model.graph,1);
        core_mem_indv{e} = find_comm_core(comm_loopy,comm_loopy_rand,num_node,p);
        
    end
    
    % plot
    rr = 0.5;
    figure;set(gcf,'color','w')
    subplot(1,2,1)
    plotGraphHighlight([Coord_active;0,0;0,max(Coord_active(:,1))],core_mem,'k',rr)
    subplot(1,2,2)
    hold on
    plotGraphHighlight([Coord_active;0,0],core_mem_indv{1},mycc.red,rr)
    plotGraphHighlight([Coord_active;0,max(Coord_active(:,1))],core_mem_indv{2},mycc.blue,rr)
    
end

end