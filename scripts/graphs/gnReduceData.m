function [] = gnReduceData(param)

data_path = param.data_path;
rand_perc = param.rand_perc;
expt_name = param.expt_name;
ee = param.ee;
num_rep = 30;

num_expt = length(expt_name);
for n = 1:num_expt
    
    expt_ee = ee{n};
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    load([data_path expt_name{n} '\Pks_Frames.mat']);
    data_raw = Spikes(:,Pks_Frame)';
    num_node = size(data_raw,2);
    
    for ii = 1:length(rand_perc)
        num_nn = floor(num_node*rand_perc(ii));
        for jj = 1:num_rep
            cell_indx = randperm(num_node,num_nn);
            data = data_raw(:,cell_indx);
            save([data_path expt_name{n} '\' expt_name{n} '_' expt_ee{1} '_' ...
                num2str(rand_perc(ii)) '_' num2str(jj) '.mat'],'data','cell_indx');
        end
    end

end