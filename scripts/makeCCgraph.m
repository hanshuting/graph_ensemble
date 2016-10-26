
function [] = makeCCgraph(param)

expt_name = param.expt_name;
ee = param.ee;
num_shuff = param.num_shuff;
data_path = param.data_path;
shuff_path_base = param.shuff_path_base;
p = param.p;

% make xcorr graph
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    num_node = size(Spikes,1);
    
    for e = 1:length(expt_ee)
        
        fprintf('processing %s_%s...\n',expt_name{n},expt_ee{e});
        
        % load real and random models
        shuff_path = [shuff_path_base 'shuffled_' ...
            expt_name{n} '_' expt_ee{e} '_loopy\'];
        
        % load data
        load([data_path expt_name{n} '\' expt_name{n} '_' expt_ee{e} '.mat']);
        
        % calculate correlation
        cc = corr(data);
        
        % threshold correlation with shuffled data
        shuff_cc = zeros(size(cc,1),size(cc,2),num_shuff);
        for i = 1:num_shuff
            load([shuff_path 'shuffled_' expt_name{n} '_' expt_ee{e} '_'...
                num2str(i) '.mat']);
            shuff_cc(:,:,i) = corr(data);
        end
        cc_thresh = zeros(size(cc));
        for i = 1:size(cc,1)
            for j = 1:size(cc,2)
                if all(isnan(shuff_cc(i,j,:)))
                    continue;
                end
                ccdist = fitdist(squeeze(shuff_cc(i,j,:)),'normal');
                cc_thresh(i,j) = icdf(ccdist,1-p);
            end
        end

        % generate graph
        cc_graph = cc>cc_thresh;
        cc_graph = cc_graph-diag(diag(cc_graph));
        cc_weight = cc.*cc_graph;
        
        % save results
        ccpath = [param.result_path_base '\' expt_name{n} '\cc\'];
        save([ccpath expt_name{n} '_' expt_ee{e} '_cc_graph.mat'],...
            'cc_graph','cc_weight','-v7.3');
        
    end
end

end