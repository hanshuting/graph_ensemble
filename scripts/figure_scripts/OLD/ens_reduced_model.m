function [] = ens_reduced_model(param)

% parameters
expt_name = param.expt_name;
ee = param.ee;
ge_type = param.ge_type;
data_path = param.data_path;
fig_path = param.fig_path.core;
save_path = param.result_path.stats_path;
result_path_base = param.result_path_base;
ccode_path = param.ccode_path;
linew = param.linew;
wsz = param.winsz;

load(ccode_path);

%% initialize
stats_all = [];
msim_all = [];

%% go through experiments
expt_count = 0;
num_win = ceil(param.maxf/wsz)+1;
for n = 1:length(expt_name)
    
    expt_ee = ee{n}{1};
    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    load([data_path expt_name{n} '\Pks_Frames.mat']);
    
    num_stim = length(unique(vis_stim))-1;
    data_high = Spikes(:,Pks_Frame)';
    vis_stim_high = vis_stim(Pks_Frame);
    
    num_frame = length(Pks_Frame);
    trunc_vec = wsz:wsz:min([param.maxf,floor(num_frame/wsz)*wsz]);
    
    % find ensembles
    core_crf = cell(num_win,1);
    for ii = 1:length(trunc_vec)
        best_model = load([model_path expt_name{n} '_' expt_ee '_' ...
            num2str(trunc_vec(ii)) '_loopy_best_model_' ge_type '.mat']);
        core_crf{ii} = find_core_max_clique(best_model.graph,num_stim);
    end
    % last one is the full model
    best_model = load([model_path expt_name{n} '_' expt_ee ...
        '_loopy_best_model_' ge_type '.mat']);
    core_crf{num_win} = find_core_max_clique(best_model.graph,num_stim);
    
    % predict
    for s = 1:num_stim
        expt_count = expt_count+1;
        true_label = double(vis_stim_high==s)';
        for ii = 1:length(trunc_vec)
            [~,~,~,sim_avg,acc,prc,rec] = core_cos_sim(core_crf{ii}{s},...
                data_high,true_label);
            msim_all(expt_count,ii,:) = sim_avg;
            stats_all(expt_count,ii,:) = [acc,prc,rec];
        end
        if length(trunc_vec)<num_win-1
            for ii = length(trunc_vec)+1:num_win-1
                msim_all(expt_count,ii,:) = NaN;
                stats_all(expt_count,ii,:) = [NaN,NaN,NaN];
            end
        end
        % last one is full data
        [~,~,~,sim_avg,acc,prc,rec] = core_cos_sim(core_crf{num_win}{s},...
            data_high,true_label);
        msim_all(expt_count,num_win,:) = sim_avg;
        stats_all(expt_count,num_win,:) = [acc,prc,rec];
    end
    
end

save([save_path 'reduced_model_pred.mat'],'core_crf','msim_all','stats_all','-v7.3');

%% plot
binsz = 20;
wr = 0.5;
trunc_vec = wsz:wsz:param.maxf+wsz;
sample_step = 50;

figure
set(gcf,'color','w')
set(gcf,'position',[2041 532 728 441]);
set(gcf,'paperpositionmode','auto')

% similarity
subplot(2,2,1);
gcapos = get(gca,'position');
dims = size(msim_all);
plot_box_seq_pair(reshape(msim_all,[dims(1),dims(2),1,dims(3)]),trunc_vec,...
    binsz,wr,sample_step,linew,mycc)
ylim([0 0.5])
xlabel('number of frames');ylabel('similarity');box off
set(gca,'position',gcapos);

% acr
dims = size(stats_all);
subplot(2,2,2);
gcapos = get(gca,'position');
plot_box_seq_single(reshape(stats_all(:,:,1),[dims(1),dims(2),1]),...
    trunc_vec,wr,sample_step,linew,mycc)
ylim([0 1])
xlabel('number of frames');ylabel('accuracy');box off
set(gca,'position',gcapos);

% prc
subplot(2,2,3);
gcapos = get(gca,'position');
plot_box_seq_single(reshape(stats_all(:,:,2),[dims(1),dims(2),1]),...
    trunc_vec,wr,sample_step,linew,mycc)
ylim([0 1])
xlabel('number of frames');ylabel('precision');box off
set(gca,'position',gcapos);

% rec
subplot(2,2,4);
gcapos = get(gca,'position');
plot_box_seq_single(reshape(stats_all(:,:,3),[dims(1),dims(2),1]),...
    trunc_vec,wr,sample_step,linew,mycc)
ylim([0 1])
xlabel('number of frames');ylabel('recall');box off
set(gca,'position',gcapos);

saveas(gcf,[fig_path 'reduced_data_model_pred_stats_' ge_type '.pdf'])


end