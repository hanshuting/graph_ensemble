function [] = plot_ens_circ_graph(param)

% parameters
expt_name = param.expt_name;
ee = param.ee;
ge_type = param.ge_type;
data_path = param.data_path;
fig_path = param.fig_path.core;
result_path_base = param.result_path_base;
savestr = param.savestr;
ccode_path = param.ccode_path;

load(ccode_path);

for n = 1:length(expt_name)
    
    expt_ee = ee{n}{1};

    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    best_model = load([model_path expt_name{n} '_' expt_ee ...
        '_loopy_best_model_' ge_type '.mat']);
    num_stim = length(setdiff(vis_stim,0));
    num_node = size(best_model.graph,1);
    
    %% find ensembles
    % load results: 'core_crf','core_svd'
    load([result_path_base '\' expt_name{n} '\core\' expt_ee '_crf_svd_core.mat']);
    
    % currenly only works for two stimuli
    num_ens = cellfun('length',core_crf);
    noncore = setdiff(1:num_node,cell2mat(core_crf(1:num_stim)));
    sort_indx = [reshape(core_crf{1},1,[]),reshape(core_crf{2},1,[]),noncore];
    edgec = cell(num_node,3);
    edgec(:,1) = num2cell(1:num_node);
    edgec(:,2) = repmat({'all'},1,num_node);
    edgec(1:num_ens(1),3) = mat2cell(ones(num_ens(1),1)*mycc.red,...
        ones(1,num_ens(1)),3);
    edgec(num_ens(1)+1:sum(num_ens(1:2)),3) = mat2cell(ones(num_ens(2),1)*mycc.blue,...
        ones(1,num_ens(2)),3);
    edgec(sum(num_ens(1:2))+1:num_node,3) = mat2cell(ones(length(noncore),1)*mycc.gray_light,...
        ones(1,length(noncore)),3);
    nodec = zeros(num_node,3);
    nodec(1:num_ens(1),:) = repmat(mycc.red,num_ens(1),1);
    nodec(num_ens(1)+1:sum(num_ens(1:2)),:) = repmat(mycc.blue,num_ens(2),1);
    nodec(sum(num_ens(1:2))+1:num_node,:) = repmat(mycc.gray_light,length(noncore),1);
    
    % plot
    figure; set(gcf,'color','w','position',[2022 313 894 372]);
    visGraphCirc(best_model.graph(sort_indx,sort_indx),'edgeColor',edgec,...
        'nodeColor',nodec);
    title(expt_name{n},'interpreter','tex')
    

end