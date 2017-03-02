
dpath = 'C:\Shuting\graph_ensemble\data\pa_raw\';
spathbase = 'C:\Shuting\graph_ensemble\data\';
ecid = {'511509529','511507650','511510650','511510670','511510718','511510855'};
tf_seq = 1;
stim_seq = [1,3];
expt1 = 20000;
expt2 = 40000;

for n = 1:length(ecid)

    load([dpath 'ecid_' ecid{n} '_fr.mat']);
    load([dpath 'ecid_' ecid{n} '_coords.mat']);
    
    vis_stim(vis_stim==5) = 1;
    vis_stim(vis_stim==6) = 2;
    vis_stim(vis_stim==7) = 3;
    vis_stim(vis_stim==8) = 4;
    vis_stim_all = vis_stim;
    Spikes_all = Spikes;
    high_vec = zeros(1,size(Spikes,2));
    high_vec(Pks_Frame)=1;

    %% extract TF=1 for all expt
    spath = [spathbase 'pa_' ecid{n} '_TF1_2stim\'];
    sname = ['pa_' ecid{n}];
    if exist(spath,'dir')~=7
        fprintf('created directory: %s\n',spath);
        mkdir(spath);
    end

    % expt data
    keep_indx = vis_stim==stim_seq(1)|vis_stim==stim_seq(2);
    Spikes = Spikes_all(:,tf==1&keep_indx);
    vis_stim = vis_stim_all(tf==1&keep_indx);
    pkframe = find(high_vec(tf==1&keep_indx));
    Pks_Frame = pkframe;
    Coord_active = coords;
    save([spath sname '_TF1.mat'],'Spikes','vis_stim','Coord_active');
    save([spath 'Pks_Frames.mat'],'Pks_Frame');

    % add neuron data
    data = Spikes_all(:,tf==1 & keep_indx & high_vec)';
    vis_stim = vis_stim_all(tf==1 & keep_indx & high_vec)';
    for ii = 1:2
        data(:,end+1) = vis_stim==stim_seq(ii);
    end
    save([spath sname '_TF1_2stim_vis_high_add_neuron.mat'],'data');

    
end

