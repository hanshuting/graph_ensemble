function [] = truncDataAddNeuron(param,wsz)
% generate truncated data by desired window size

% parameters
expt_name = param.expt_name;
data_path = param.data_path;

for n = 1:length(expt_name)
    
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    load([data_path expt_name{n} '\Pks_Frames.mat']);
    data_high = Spikes(:,Pks_Frame)';
    vis_stim_high = vis_stim(Pks_Frame);
    
    num_node = size(data_high,2);
    num_frame = length(Pks_Frame);
    trunc_vec = wsz:wsz:floor(num_frame/wsz)*wsz;
    
    % add stim associated nodes
    num_stim = length(unique(vis_stim))-1;
    for ii = 1:num_stim
        data_high(:,num_node+ii) = vis_stim_high==ii';
    end
    
    % truncate data
    for ii = 1:length(trunc_vec)
        data = data_high(1:trunc_vec(ii),:);
        save([data_path expt_name{n} '\' expt_name{n} '_all_high_add_neuron_'...
            num2str(trunc_vec(ii)) '.mat'],'data');
    end

end

end