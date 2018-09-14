function [] = fig3_random_graph_comm(param)

% parameters
expt_name = param.expt_name;
ee = param.ee;
num_rand = param.num_rand;
k_seq = param.k_seq;
k = param.k;
p = param.p;
ge_meth = param.ge_meth;
fig_path = param.fig_path.graph_prop;
result_path_base = param.result_path_base;
ccode_path = param.ccode_path;
savestr = param.savestr;

load(ccode_path);

%% 
nexpt = sum(cellfun('length',ee));
num_clique_pair = zeros(nexpt,length(k_seq),2);
num_comm_pair = zeros(nexpt,length(k_seq),2);
sz_comm_pair = cell(nexpt,length(k_seq),2);

% load data
expt_count = 0;
for n = 1:length(expt_name)
    expt_ee = ee{n};
    for e = 1:length(expt_ee)

        fprintf('processing %s_%s...\n',expt_name{n},expt_ee{e});
        expt_count = expt_count+1;

        % load real and random models
        model_path = [result_path_base expt_name{n} '\models\'];
        comm_path = [result_path_base expt_name{n} '\comm_' ge_meth '\']; 
        rand_path = [result_path_base expt_name{n} '\rand_graph_' ge_meth '\'];
        load([model_path expt_name{n} '_' expt_ee{e} '_loopy_best_model.mat']);
        
        % store results
        for ii = 1:length(k_seq)
            
            % load clique_loopy_rand, comm_loopy_rand, loopy_rand_graph
            load([comm_path expt_name{n} '_' expt_ee{e} ...
                '_loopy_rand_graph_comm_k_' num2str(k_seq(ii)) '.mat']);
            % load clique_loopy, comm_loopy
            load([comm_path expt_name{n} '_' expt_ee{e} '_loopy_comm_k_' ...
                num2str(k_seq(ii)) '.mat']);
            % load num_rand_clique, num_rand_comm, sz_rand_comm
            load([rand_path expt_name{n} '_' expt_ee{e} '_loopy_rand_graph_stats_k_' ...
                num2str(k_seq(ii)) '.mat']);
            
            num_node = size(loopy_rand_graph,1);
            
            % clique number
            sz_clique = cellfun('length',clique_loopy);
            sz_clique = sz_clique(sz_clique>k_seq(ii));
            if ~isempty(sz_clique)
                sz_clique = num2cell(sz_clique);
                num_clique = cellfun(@(x) nchoosek(x,k_seq(ii)),sz_clique,'uniformoutput',false);
                num_clique_pair(expt_count,ii,1) = sum(cell2mat(num_clique)); %/num_node;
            else
                num_clique_pair(expt_count,ii,1) = 0;
            end
            num_clique_pair(expt_count,ii,2) = mean(num_rand_clique); %/num_node;

            % community number and size
            % first column, real data
            if ~isempty(comm_loopy)
                sz_comm_pair{expt_count,ii,1} = cellfun('length',comm_loopy);
                num_comm_pair(expt_count,ii,1) = length(sz_comm_pair{expt_count,ii,1}); %/num_node;
            else
                sz_comm_pair{expt_count,ii,1} = 0;
                num_comm_pair(expt_count,ii,1) = 0;
            end
            % second column, random graph data
            cr_rand_comm = sz_rand_comm;
            for loopvar = 1:num_rand
                if cr_rand_comm{loopvar}==0
                    cr_rand_comm{loopvar} = [];
                end
            end
            sz_comm = cell2mat(cr_rand_comm');
            if isempty(sz_comm)
                sz_comm = 0;
            end
            sz_comm_pair{expt_count,ii,2} = sz_comm; %/num_node;
            num_comm_pair(expt_count,ii,2) = mean(cellfun('length',cr_rand_comm)); %/num_node;

        end
    end

end

% save([save_path 'comm_clique_stats.mat'],'sz_comm_pair','num_comm_pair',...
%     'num_clique_pair','-v7.3');


%% rank sum test
% clique number
p_num_clique = zeros(length(k_seq),1);
for ii = 1:length(k_seq)
    [~,p_num_clique(ii)] = ranksum(num_clique_pair(:,ii,1),num_clique_pair(:,ii,2),...
        'tail','right');
end

% community number
p_num_comm = zeros(length(k_seq),1);
for ii = 1:length(k_seq)
    [~,p_num_comm(ii)] = ranksum(num_comm_pair(:,ii,1),num_comm_pair(:,ii,2),...
        'tail','right');
end

% community size
sz_comm_inq = sz_comm_pair(:,k_seq==k,:);
sz_comm_real = cell2mat(squeeze(sz_comm_inq(:,:,1)));
sz_comm_rand = cell2mat(squeeze(sz_comm_inq(:,:,2)));

% convert to distribution
bin_range = 0:1:60;
hist_comm_real = hist(sz_comm_real,bin_range);
hist_comm_real = hist_comm_real/sum(hist_comm_real);
hist_comm_rand = hist(sz_comm_rand,bin_range);
hist_comm_rand = hist_comm_rand/sum(hist_comm_rand);
cdf_comm_rand = cumsum(hist_comm_rand);

[~,p_sz_comm] = ranksum(hist_comm_rand,hist_comm_real,'tail','right');


%% plot
linew = 1;
scsz = 30;

figure;
set(gcf,'color','w','position',[1969 602 683 178],'PaperPositionMode','auto');
gridsz = 0.2;

% clique number
subplot(1,3,1)
hold on;
yma = 1.1*max(num_clique_pair(:));
for ii = 1:length(k_seq)
    scatter((k_seq(ii)-gridsz)*ones(size(num_clique_pair(:,ii,2))),...
        num_clique_pair(:,ii,2),scsz,mycc.black,'+','linewidth',linew);
    scatter((k_seq(ii)+gridsz)*ones(size(num_clique_pair(:,ii,1))),...
        num_clique_pair(:,ii,1),scsz,mycc.orange,'+','linewidth',linew);
    plot([(k_seq(ii)-gridsz)*ones(size(num_clique_pair(:,ii,2))),...
        (k_seq(ii)+gridsz)*ones(size(num_clique_pair(:,ii,1)))]',...
        [num_clique_pair(:,ii,2),num_clique_pair(:,ii,1)]','color',mycc.gray);
    if p_num_clique(ii)<p
        scatter(k_seq(ii),yma,20,'k*');
    end
end
xlim([k_seq(1)-1,k_seq(end)+1]);
ylim([0 yma])
xlabel('k');ylabel('N_c_l_i_q_u_e');
set(gca,'xtick',k_seq,'linewidth',linew)
% legend('real data','random graph')
box off

% community number
subplot(1,3,2)
hold on;
yma = 1.1*max(num_comm_pair(:));
for ii = 1:length(k_seq)
    scatter((k_seq(ii)-gridsz)*ones(size(num_comm_pair(:,ii,2))),...
        num_comm_pair(:,ii,2),scsz,mycc.black,'+','linewidth',linew);
    scatter((k_seq(ii)+gridsz)*ones(size(num_comm_pair(:,ii,1))),...
        num_comm_pair(:,ii,1),scsz,mycc.orange,'+','linewidth',linew);
    plot([(k_seq(ii)-gridsz)*ones(size(num_comm_pair(:,ii,2))),...
        (k_seq(ii)+gridsz)*ones(size(num_comm_pair(:,ii,1)))]',...
        [num_comm_pair(:,ii,2),num_comm_pair(:,ii,1)]','color',mycc.gray);
    if p_num_comm(ii)<p
        scatter(k_seq(ii),yma,20,'k*');
    end
end
xlim([k_seq(1)-1,k_seq(end)+1]);
ylim([0 yma])
xlabel('k');ylabel('N_c_o_m');
set(gca,'xtick',k_seq,'linewidth',linew)
% legend('real data','random graph')
box off

% community size
subplot(1,3,3)
histogram(sz_comm_rand,bin_range,'normalization','probability','facecolor',mycc.black)
hold on;
histogram(sz_comm_real,bin_range,'normalization','probability','facecolor',mycc.orange)
sz_thresh = find(cdf_comm_rand>1-p,1);
plot([sz_thresh,sz_thresh],[0,max(hist_comm_rand)],'k--');
xlabel('s^c^o^m');ylabel('percentage');
xlim([bin_range(1) bin_range(end)])
ylim([0,max(hist_comm_rand)])
set(gca,'linewidth',linew)
legend('CRF_r_a_n_d','CRF')
box off

saveas(gcf,[fig_path 'rand_graph_comm_' savestr '.fig'])
saveas(gcf,[fig_path 'rand_graph_comm_' savestr '.pdf'])

%% probability of observing a k-clique in random graphs

end
