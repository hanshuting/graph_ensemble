%% --------- random graph properties --------------
% parameters
expt_name = 'mf2_v1';%'m21_d2_vis';%'m21_sample_d2_d3';
%ee = {'01','02','on','off'};
%ee = {'vis_01','vis_02','po_vis_01','po_vis_02','all_vis_01','all_vis_02','all'};
ee = {'vis_01','vis_02','vis_03','vis_04','vis_01_all','vis_02_all','vis_03_all','vis_04_all','vis_all','vis_all_all'};
data_path = '/home/sh3276/work/fwMatch/data/';
model_path = ['/home/sh3276/work/fwMatch/results/' expt_name '/models/'];
save_path = ['/home/sh3276/work/fwMatch/results/' expt_name '/rand_graph/'];

num_rand = 100;
k = 2:6;

%% load data
load([data_path expt_name '.mat']);

for n = 1:length(ee)
    
    fprintf('experiment %s_%s...\n',expt_name,ee{n});
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
%         fprintf('k=%u\n;,k(i)');
%         % number of node that can form a complete graph
%         Nc = floor(sqrt(2*num_edge));
% 
%         % number of extra node that can form k cliques
%         Ne = floor((num_edge-Nc)/k(i));
%         
%         % maximum number of k cliques
%         Nk = nchoosek(Nc,k(i))+nchoosek(Ne*k(i),k(i)-1);
    
        parfor j = 1:num_rand
        
%             fprintf('random graph #%u\n',j);
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