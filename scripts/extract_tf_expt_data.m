
dpath = 'C:\Shuting\graph_ensemble\data\pa_raw\';
ecid = {'ecid_511507650','ecid_511509529','ecid_511510650','ecid_511510718','ecid_511510855'};
spath = 'C:\Shuting\graph_ensemble\data\';
target_tf = [1,15];

for n = 1:length(ecid)
    
    load([dpath ecid{n} '_spikes.mat']);
    vis_stim_full = vis_stim;
    
    for ii = 1:length(target_tf)
        
        data = Spikes(:,tf==target_tf(ii))';
        vis_stim = vis_stim_full(tf==target_tf(ii));
        
        expt_spath = [spath ecid{n} '_tf_' num2str(target_tf(ii))];
        if exist(expt_spath,'dir')~=7
            mkdir(expt_spath);
        end
        save([expt_spath '\' ecid{n} '_tf_' num2str(target_tf(ii)) '.mat'],...
            'data','vis_stim','tf','-v7.3');
        
    end
end