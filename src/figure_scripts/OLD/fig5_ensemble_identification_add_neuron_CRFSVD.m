function [] = fig5_ensemble_identification_add_neuron_CRFSVD(param)

% parameters
expt_name = param.expt_name;
ee = param.ee;
num_shuff = param.num_shuff;
p = param.p;
ge_type = param.ge_type;
data_path = param.data_path;
fig_path = param.fig_path.core;
save_path = param.result_path.stats_path;
result_path_base = param.result_path_base;
savestr = param.savestr;
ccode_path = param.ccode_path;
OSI_thresh = param.OSI_thresh;
mc_minsz = param.mc_minsz;
linew = param.linew;

load(ccode_path);

%% initialize
num_expt = length(expt_name);
num_crf_svd = zeros(num_expt,1);
num_crf_osi = zeros(num_expt,1);
num_crf = zeros(num_expt,1);
num_svd = zeros(num_expt,1);
num_osi = zeros(num_expt,1);
num_cell = zeros(num_expt,1);

num_mc_in_core = [];

expt_count = 0;

%% graph properties
for n = 1:length(expt_name)
    
    expt_ee = ee{n}{1};

    model_path = [result_path_base '\' expt_name{n} '\models\']; 
    
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    svd_data = load([data_path 'ensembles\' expt_name{n} '_core_svd.mat']);
    best_model = load([model_path expt_name{n} '_' expt_ee ...
        '_loopy_best_model_' ge_type '.mat']);
    num_stim = length(unique(vis_stim))-1;
    
    %% find ensembles
    % this code now only works with two stimuli
    % SVD
    core_svd = cell(num_stim,1);
    for ii = 1:num_stim
        for jj = 1:length(svd_data.svd_state)
            if strcmp(num2str(ii),svd_data.svd_state{jj})
                core_svd{ii} = svd_data.core_svd{jj};
                break;
            end
        end
    end
    
    % CRF max clique
    [core_crf,mc_in_core] = find_core_max_clique(best_model.graph,num_stim,mc_minsz);
    num_in_core = cellfun('length',mc_in_core);
    num_mc_in_core(end+1:end+length(num_in_core)) = num_in_core;
    
    %% plot ensemble
    rr = 1;
    figure;
    set(gcf,'color','w','position',[2041 430 543 338]);
    set(gcf,'paperpositionmode','auto')
    for ii = 1:num_stim
        
        expt_count = expt_count+1;
        
        crf_svd = intersect(core_crf{ii},core_svd{ii});
        
        num_cell(expt_count) = size(Spikes,1);
        num_crf(expt_count) = length(core_crf{ii});
        num_svd(expt_count) = length(core_svd{ii});
        num_crf_svd(expt_count) = length(crf_svd);
        
        subplot(1,num_stim,ii);
        plotCoreOverlay(Coord_active,core_crf{ii},core_svd{ii},mycc.orange,...
            mycc.green,rr)
        
    end
    
    print(gcf,'-dpdf','-painters',[fig_path expt_name{n} '_' expt_ee '_mc_svd_core.pdf'])
    
    %% save result
    save([result_path_base '\' expt_name{n} '\core\' expt_ee '_mc_svd_core_' ...
        ge_type '.mat'],'core_svd','core_crf','-v7.3');
    
end

%% plot number of maximal cliques in core
figure;set(gcf,'color','w','position',[2008 525 577 219])
boxwd = 0.2;hold on;
h = boxplot(num_mc_in_core,'positions',0.5,'width',...
    boxwd,'colors',mycc.black);
setBoxStyle(h,linew);
xlim([0 1])
gcapos = get(gca,'position');
ylabel('# maximal cliques in core')
set(gca,'linewidth',linew)
set(gca,'position',gcapos);
box off

saveas(gcf,[fig_path 'num_mc_in_core.pdf']);

%% plot stats
stepsz = 0.5;
ww = 0.3;

figure
set(gcf,'color','w');
set(gcf,'position',[2096 452 447 233])
set(gcf,'paperpositionmode','auto')

% ensemble number
subplot(1,2,1); hold on
h = boxplot(num_svd./num_cell,'positions',stepsz,'width',ww,'colors',mycc.green);
setBoxStyle(h,linew);
h = boxplot(num_crf./num_cell,'positions',2*stepsz,'width',ww,'colors',mycc.orange);
setBoxStyle(h,linew);
xlim([0 3*stepsz])
ylim([0 max([num_svd./num_cell;num_crf./num_cell])])
set(gca,'xtick',[1,2]*stepsz,'xticklabel',{'SVD','CRF'})
ylabel('Cells %')
box off

% CRF+SVD
subplot(1,2,2); hold on
h = boxplot(num_crf_svd./(num_crf+num_svd),'positions',stepsz,'width',ww,'colors',mycc.black);
setBoxStyle(h,linew);
xlim([0 3*stepsz])
ylim([0 1])
ylabel('Nshared%')
set(gca,'xcolor','w')
box off

print(gcf,'-dpdf','-painters',[fig_path expt_ee '_mc_svd_core_nums_' ...
    ge_type '.pdf'])

end
