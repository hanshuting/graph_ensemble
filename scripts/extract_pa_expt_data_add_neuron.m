
dpath = 'C:\Shuting\graph_ensemble\data\pa_raw\';
spathbase = 'C:\Shuting\graph_ensemble\data\';
ecid = {'511509529','511507650','511510650','511510670','511510718','511510855'};
tf_seq = [1,2,4,8,15];
expt1 = 20000;
expt2 = 40000;

for n = 1:length(ecid)

    load([dpath 'ecid_' ecid{n} '_fr.mat']);
    load([dpath 'ecid_' ecid{n} '_coords.mat']);

    vis_stim(vis_stim==5)=1;
    vis_stim(vis_stim==6)=2;
    vis_stim(vis_stim==7)=3;
    vis_stim(vis_stim==8)=4;
    vis_stim_all = vis_stim;
    Spikes_all = Spikes;
    high_vec = zeros(1,size(Spikes,2));
    high_vec(Pks_Frame)=1;

    %% extract TF=1 for all expt
    spath = [spathbase 'pa_' ecid{n} '_TF1\'];
    sname = ['pa_' ecid{n}];
    if exist(spath,'dir')~=7
        mkdir(spath);
    end

    % expt data
    Spikes = Spikes_all(:,tf==1);
    vis_stim = vis_stim_all(tf==1);
    pkframe = find(high_vec(tf==1));
    Pks_Frame = pkframe;
    Coord_active = coords;
    save([spath sname '_TF1.mat'],'Spikes','vis_stim','pkframe','Coord_active');
    save([spath 'coords.mat'],'Coord_active');
    %     save([spath 'Pks_Frame.mat'],'Pks_Frame');

    % add neuron data
    data = Spikes_all(:,tf==1&high_vec)';
    vis_stim = vis_stim_all(tf==1&high_vec)';
    for ii = 1:4
        data(:,end+1) = vis_stim==ii;
    end
%     save([spath sname '_TF1_vis_high_add_neuron.mat'],'data');

    % test data with other TFs
    for ii = tf_seq
        data = Spikes_all(:,tf==ii&high_vec)';
        vis_stim = vis_stim_all(tf==ii&high_vec)';
        save([spath sname '_TF' num2str(ii) '_vis_high.mat'],'data','vis_stim');
    end
    
    %% extract TF=1 for the first experiment
    spath = [spathbase 'pa_' ecid{n} '_1_TF1\'];
    sname = ['pa_' ecid{n}];
    if exist(spath,'dir')~=7
        mkdir(spath);
    end

    % expt data
    Spikes = Spikes_all(:,1:expt1);
    Spikes = Spikes(:,tf(1:expt1)==1);
    vis_stim = vis_stim_all(1:expt1);
    vis_stim = vis_stim(tf(1:expt1)==1);
    pkframe = high_vec(1:expt1);
    Pks_Frame = pkframe;
    pkframe = find(pkframe(tf(1:expt1)==1));
    save([spath sname '_1_TF1.mat'],'Spikes','vis_stim','pkframe','Coord_active');
    save([spath 'Pks_Frame.mat'],'Pks_Frame');
    save([spath 'coords.mat'],'Coord_active');

    % add neuron data
    data = Spikes_all(:,1:expt1);
    data = data(:,tf(1:expt1)==1&high_vec(1:expt1))';
    vis_stim = vis_stim_all(1:expt1);
    vis_stim = vis_stim(tf(1:expt1)==1&high_vec(1:expt1))';
    for ii = 1:4
        data(:,end+1) = vis_stim==ii;
    end
%     save([spath sname '_1_TF1_vis_high_add_neuron.mat'],'data');

    % test data with other two expts
    data = Spikes_all(:,expt2:end);
    data = data(:,tf(expt2:end)==1&high_vec(expt2:end))';
    vis_stim = vis_stim_all(expt2:end);
    vis_stim = vis_stim(tf(expt2:end)==1&high_vec(expt2:end))';
%     save([spath sname '_23_TF1_vis_high.mat'],'data','vis_stim');
    
    % test data with other TFs
    for ii = tf_seq
        data = Spikes_all(:,1:expt1);
        data = data(:,tf(1:expt1)==ii&high_vec(1:expt1))';
        vis_stim = vis_stim_all(1:expt1);
        vis_stim = vis_stim(tf(1:expt1)==ii&high_vec(1:expt1))';
        save([spath sname '_1_TF' num2str(ii) '_vis_high.mat'],'data','vis_stim');
    end
    
    for ii = tf_seq
        data = Spikes_all(:,expt2:end);
        data = data(:,tf(expt2:end)==ii&high_vec(expt2:end))';
        vis_stim = vis_stim_all(expt2:end);
        vis_stim = vis_stim(tf(expt2:end)==ii&high_vec(expt2:end))';
        save([spath sname '_23_TF' num2str(ii) '_vis_high.mat'],'data','vis_stim');
    end
    
end

