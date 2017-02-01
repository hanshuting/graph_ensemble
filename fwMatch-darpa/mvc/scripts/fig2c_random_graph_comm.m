%% --------- random graph properties --------------

% parameters
expt_name = 'm24_vis_opto';%'m21_sample_d2_d3';
ee = {'vis_01','vis_02'};
data_path = 'C:\Shuting\fwMatch\data\';
model_path = ['C:\Shuting\fwMatch\results\' expt_name '\models\'];
save_path = ['C:\Shuting\fwMatch\results\' expt_name '\rand_graph\'];

num_rand = 100;
k = 2:6;

%% load data
load([data_path expt_name '.mat']);

for n = 1:length(ee)
    
    load([model_path expt_name '_' ee{n} '_loopy_best_model.mat']);
    
    %% generate random graphs
    num_node = size(graph,1);
    num_edge = sum(graph(:))/2;
    rand_graph = zeros(num_node,num_node,num_rand);

    for i = 1:num_rand
        rand_graph(:,:,i) = mkRandomGraph(num_node,num_edge);
    end

    %% clique and community properties
    num_rand_clique = zeros(length(k),num_rand);
    num_rand_comm = zeros(length(k),num_rand);
    sz_rand_comm = cell(length(k),num_rand);
    for i = 1:length(k)
        fprintf('k=%u\n;',k(i));
%         % number of node that can form a complete graph
%         Nc = floor(sqrt(2*num_edge));
% 
%         % number of extra node that can form k cliques
%         Ne = floor((num_edge-Nc)/k(i));
%         
%         % maximum number of k cliques
%         Nk = nchoosek(Nc,k(i))+nchoosek(Ne*k(i),k(i)-1);
    
        for j = 1:num_rand
        
            fprintf('random graph #%u\n',j);
            [comm,clique,~] = k_clique(k(i),rand_graph(:,:,j));
        
            % number of k cliques
            sz_clique = cellfun('length',clique);
            sz_clique = sz_clique(sz_clique>k(i));
            if ~isempty(sz_clique)
                sz_clique = num2cell(sz_clique);
                num_clique = cellfun(@(x) nchoosek(x,k(i)),sz_clique,'uniformoutput',false);
                num_rand_clique(i,j) = sum(cell2mat(num_clique));
            else
                num_rand_clique(i,j) = 0;
            end
        
            % number and size of communities
            if ~isempty(comm)
                sz_rand_comm{i,j} = cellfun('length',comm);
                num_rand_comm(i,j) = length(sz_rand_comm{i,j});
            else
                sz_rand_comm{i,j} = 0;
                num_rand_comm(i,j) = 0;
            end
    
        end
    end
        
    save([save_path ee{n} '_rand_graph_stats.mat'],'num_rand_clique',...
        'num_rand_comm','sz_rand_comm','-v7.3');

end

%% compare with real data

% expt_name = {'m21_d2_vis',...
%              'm21_sample_d2_d3',...
%              'm24_vis_opto',...
%              'mf2_v1'};
% ee = {{'01','02','on','off'},...
% {'vis_01','vis_02','po_vis_01','po_vis_02','all_vis_01','all_vis_02','all'},...
% {'vis_01','vis_02','vis_on','vis_off'},...
% {'vis_01','vis_02','vis_03','vis_04','vis_01_all','vis_02_all','vis_03_all','vis_04_all','vis_all','vis_all_all'}};
expt_name = {'m21_d2_vis','mf2_v1'};
ee = {{'01','02'},{'vis_01_all','vis_02_all','vis_04_all'}};

num_rand = 100;
k = 2:6;
save_path = 'C:\Shuting\fwMatch\results\rand_graph\';

% initialize
nexpt = sum(cellfun('length',ee));
num_clique_pair = zeros(nexpt,length(k),2);
num_comm_pair = zeros(nexpt,length(k),2);
sz_comm_pair = cell(nexpt,length(k),2);

% load data
expt_num = 0;
for n = 1:length(expt_name)
    expt_ee = ee{n};
    for e = 1:length(expt_ee)
        
        fprintf('processing %s_%s...\n',expt_name{n},expt_ee{e});
        expt_num = expt_num+1;
        
        % load real and random models
        model_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\models\'];
        rand_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\rand_graph\'];
        load([model_path expt_name{n} '_' expt_ee{e} '_loopy_best_model.mat']);
        load([rand_path expt_ee{e} '_rand_graph_stats.mat']);
        
        % store results
        for i = 1:length(k)
            
            % clique number
            [comm,clique,~] = k_clique(k(i),graph);
            sz_clique = cellfun('length',clique);
            sz_clique = sz_clique(sz_clique>k(i));
            if ~isempty(sz_clique)
                sz_clique = num2cell(sz_clique);
                num_clique = cellfun(@(x) nchoosek(x,k(i)),sz_clique,'uniformoutput',false);
                num_clique_pair(expt_num,i,1) = sum(cell2mat(num_clique));
            else
                num_clique_pair(expt_num,i,1) = 0;
            end
            num_clique_pair(expt_num,i,2) = mean(num_rand_clique(i,:));
            
            % community number and size
            % first column, real data
            if ~isempty(comm)
                sz_comm_pair{expt_num,i,1} = cellfun('length',comm);
                num_comm_pair(expt_num,i,1) = length(sz_comm_pair{expt_num,i,1});
            else
                sz_comm_pair{expt_num,i,1} = 0;
                num_comm_pair(expt_num,i,1) = 0;
            end
            % second column, random graph data
            cr_rand_comm = sz_rand_comm(i,:);
            for loopvar = 1:num_rand
                if cr_rand_comm{loopvar}==0
                    cr_rand_comm{loopvar} = [];
                end
            end
            sz_comm = cell2mat(cr_rand_comm');
            if isempty(sz_comm)
                sz_comm = 0;
            end
            sz_comm_pair{expt_num,i,2} = sz_comm;
            num_comm_pair(expt_num,i,2) = mean(cellfun('length',cr_rand_comm));
            
        end
        
        
    end
end

% save([save_path 'comm_clique_stats.mat'],'sz_comm_pair','num_comm_pair',...
%     'num_clique_pair','-v7.3');

%% plot
p = 0.05;

% ----- clique number ----- %
% paired t-test
p_num_clique = zeros(length(k),1);
for i = 1:length(k)
    [~,p_num_clique(i)] = ttest(num_clique_pair(:,i,1),num_clique_pair(:,i,2),...
        'tail','right');
end

% plot
h = figure;set(h,'color','w');
gridsz = 0.2;
hold on;
for i = 1:length(k)
    scatter((k(i)-gridsz)*ones(size(num_clique_pair(:,i,2))),...
        num_clique_pair(:,i,2),'k','linewidth',1);
    scatter((k(i)+gridsz)*ones(size(num_clique_pair(:,i,1))),...
        num_clique_pair(:,i,1),'k','filled','linewidth',1);
    plot([(k(i)-gridsz)*ones(size(num_clique_pair(:,i,2))),...
        (k(i)+gridsz)*ones(size(num_clique_pair(:,i,1)))]',...
        [num_clique_pair(:,i,2),num_clique_pair(:,i,1)]','k-');
    if p_num_clique(i)<p
        scatter(k(i),1.2*max(max(num_clique_pair(:,i,:))),20,'k*');
    end
end
xlim([k(1)-1,k(end)+1]);
xlabel('k');ylabel('number of cliques');
legend('real data','random graph')
box off

% ----- community number ----- %
% paired t-test
p_num_comm = zeros(length(k),1);
for i = 1:length(k)
    [~,p_num_comm(i)] = ttest(num_comm_pair(:,i,1),num_comm_pair(:,i,2),...
        'tail','right');
end
    
% plot
h = figure;set(h,'color','w');
gridsz = 0.2;
hold on;
for i = 1:length(k)
    scatter((k(i)-gridsz)*ones(size(num_comm_pair(:,i,2))),...
        num_comm_pair(:,i,2),'k','linewidth',1);
    scatter((k(i)+gridsz)*ones(size(num_comm_pair(:,i,1))),...
        num_comm_pair(:,i,1),'k','filled','linewidth',1);
    plot([(k(i)-gridsz)*ones(size(num_comm_pair(:,i,2))),...
        (k(i)+gridsz)*ones(size(num_comm_pair(:,i,1)))]',...
        [num_comm_pair(:,i,2),num_comm_pair(:,i,1)]','k-');
    if p_num_comm(i)<p
        scatter(k(i),1.2*max(max(num_comm_pair(:,i,:))),20,'k*');
    end
end
xlim([k(1)-1,k(end)+1]);
xlabel('k');ylabel('number of community');
legend('real data','random graph')
box off

% ----- community size ----- %
k_inq = 3;
sz_comm_inq = sz_comm_pair(:,k==k_inq,:);
sz_comm_real = cell2mat(squeeze(sz_comm_inq(:,:,1)));
sz_comm_rand = cell2mat(squeeze(sz_comm_inq(:,:,2)));

% convert to distribution
bin_range = 0:1:60;
hist_comm_real = hist(sz_comm_real,bin_range);
hist_comm_real = hist_comm_real/sum(hist_comm_real);
hist_comm_rand = hist(sz_comm_rand,bin_range);
hist_comm_rand = hist_comm_rand/sum(hist_comm_rand);
cdf_comm_rand = cumsum(hist_comm_rand);

% t-test
[~,p_sz_comm] = ttest2(hist_comm_rand,hist_comm_real,'tail','right');
    
% plot
h = figure;set(h,'color','w');
histogram(sz_comm_rand,bin_range,'normalization','probability','facecolor',[0.5 0.5 0.5])
hold on;
histogram(sz_comm_real,bin_range,'normalization','probability','facecolor',[1 0.2 0.2])
sz_thresh = find(cdf_comm_rand>1-p,1);
plot([sz_thresh,sz_thresh],[0,max(hist_comm_rand)],'k--');
xlabel('community size');ylabel('percentage');
legend('random graph','real data')
box off


%% probability of observing a k-clique in random graphs


