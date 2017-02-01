% generate high activity frames

% parameters
rng(1000);
expt_name = {'m21_d2_vis'};
ee = {{'01_high','02_high'}};
test_ee = {{'01','02','off','spon'}};
num_shuff = 100;
k = 3;
data_path = 'C:\Shuting\fwMatch\data\';

% make shuffled data
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    expt_test = test_ee{n};
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    num_node = size(Spikes,1);
    
    for e = 1:length(expt_ee)
        
        fprintf('processing %s_%s...\n',expt_name{n},expt_ee{e});
        
        for m = 1:length(expt_test)
        
            % load data
            core_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\core\'];
            load([core_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k) ...
                '_core.mat']);
            load([data_path expt_name{n} '\' expt_name{n} '_' expt_test{m} '.mat']);
            full_data = data;
            num_frame = size(data,1);

            %---------- svd ----------%
            cell_indx = ee_core_svd;
            num_spike = zeros(num_frame,num_shuff);
            for i = 1:num_shuff
                data = shuffle(full_data(:,cell_indx)','time')';
                num_spike(:,i) = sum(data,2);
            end
        
        end
        
    end
end