% extract ensemble data for training

rng(1000);
expt_name = {'m21_d2_vis'};
ee = {{'01_high','02_high'}};
test_ee = {{'01_high','02_high','off_high','spon_high'}};
% ee = {{'01_high','02_high'},{'vis_01_all','vis_02_all','vis_04_all'}};
num_shuff = 100;
k = 3;
p = 0.05;
perc_train = 0.8;

%% extract training and test data
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    
    for e = 1:length(expt_ee)
         
        % set path
        data_path = ['C:\Shuting\fwMatch\data\' expt_name{n} '\']; 
        core_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\core\'];
        load([core_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k) ...
            '_core.mat']);
        load([data_path expt_name{n} '_' expt_ee{e} '.mat']);
        
        num_frame = size(data,1);
        num_node = size(data,2);
        full_data = data;
        num_train = round(perc_train*num_frame);
        indxTrain = randperm(num_frame);
        indxTest = indxTrain(num_train+1:end);
        indxTrain = indxTrain(1:num_train);
        
        %---------- svd ----------%
        % core training
        cell_indx = ee_core_svd;
        data = full_data(indxTrain,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_svd_core_train_k_'...
            num2str(k) '.mat'],'data');
        % core test
        data = full_data(indxTest,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_svd_core_test_k_'...
            num2str(k) '.mat'],'data');
        % noncore training
        cell_indx = setdiff(1:num_node,cell_indx);
        data = full_data(indxTrain,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_svd_noncore_train_k_'...
            num2str(k) '.mat'],'data');
        % noncore test
        data = full_data(indxTest,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_svd_noncore_test_k_'...
            num2str(k) '.mat'],'data');
        
        %---------- cc mem ----------%
        % core training
        cell_indx = core_mem_cc;
        data = full_data(indxTrain,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_cc_mem_core_train_k_'...
            num2str(k) '.mat'],'data');
        % core test
        data = full_data(indxTest,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_cc_mem_core_test_k_'...
            num2str(k) '.mat'],'data');
        % noncore training
        cell_indx = setdiff(1:num_node,cell_indx);
        data = full_data(indxTrain,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_cc_mem_noncore_train_k_'...
            num2str(k) '.mat'],'data');
        % noncore test
        data = full_data(indxTest,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_cc_mem_noncore_test_k_'...
            num2str(k) '.mat'],'data');
        
        %---------- loopy mem ----------%
        % core training
        cell_indx = core_mem_loopy;
        data = full_data(indxTrain,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_loopy_mem_core_train_k_'...
            num2str(k) '.mat'],'data');
        % core test
        data = full_data(indxTest,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_loopy_mem_core_test_k_'...
            num2str(k) '.mat'],'data');
        % noncore training
        cell_indx = setdiff(1:num_node,cell_indx);
        data = full_data(indxTrain,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_loopy_mem_noncore_train_k_'...
            num2str(k) '.mat'],'data');
        % noncore test
        data = full_data(indxTest,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_loopy_mem_noncore_test_k_'...
            num2str(k) '.mat'],'data');
        
        %---------- cc centrality ----------%
        % core training
        cell_indx = core_cent_cc;
        data = full_data(indxTrain,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_cc_cent_core_train_k_'...
            num2str(k) '.mat'],'data');
        % core test
        data = full_data(indxTest,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_cc_cent_core_test_k_'...
            num2str(k) '.mat'],'data');
        % noncore training
        cell_indx = setdiff(1:num_node,cell_indx);
        data = full_data(indxTrain,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_cc_cent_noncore_train_k_'...
            num2str(k) '.mat'],'data');
        % noncore test
        data = full_data(indxTest,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_cc_cent_noncore_test_k_'...
            num2str(k) '.mat'],'data');

        
        %---------- loopy centrality ----------%
        % core training
        cell_indx = core_cent_loopy;
        data = full_data(indxTrain,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_loopy_cent_core_train_k_'...
            num2str(k) '.mat'],'data');
        % core test
        data = full_data(indxTest,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_loopy_cent_core_test_k_'...
            num2str(k) '.mat'],'data');
        % noncore training
        cell_indx = setdiff(1:num_node,cell_indx);
        data = full_data(indxTrain,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_loopy_cent_noncore_train_k_'...
            num2str(k) '.mat'],'data');
        % noncore test
        data = full_data(indxTest,cell_indx);
        save([data_path expt_name{n} '_' expt_ee{e} '_loopy_cent_noncore_test_k_'...
            num2str(k) '.mat'],'data');

    end
    
end

%% extract data on test conditions
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    expt_test = test_ee{n};
    
    for e = 1:length(expt_ee)
         
        % set path
        data_path = ['C:\Shuting\fwMatch\data\' expt_name{n} '\']; 
        core_path = ['C:\Shuting\fwMatch\results\' expt_name{n} '\core\'];
        load([core_path expt_name{n} '_' expt_ee{e} '_k_' num2str(k) ...
            '_core.mat']);
            
        for m = 1:length(expt_test)    
                        
            load([data_path expt_name{n} '_' expt_test{m} '.mat']);

            num_frame = size(data,1);
            num_node = size(data,2);
            full_data = data;

            %---------- svd ----------%
            % core training
            cell_indx = ee_core_svd;
            data = full_data(:,cell_indx);
            save([data_path expt_name{n} '_' expt_ee{e} '_svd_core_' ...
                expt_test{m} '_k_' num2str(k) '.mat'],'data');
            % noncore training
            cell_indx = setdiff(1:num_node,cell_indx);
            data = full_data(:,cell_indx);
            save([data_path expt_name{n} '_' expt_ee{e} '_svd_noncore_' ...
                expt_test{m} '_k_' num2str(k) '.mat'],'data');

            %---------- cc mem ----------%
            % core training
            cell_indx = core_mem_cc;
            data = full_data(:,cell_indx);
            save([data_path expt_name{n} '_' expt_ee{e} '_cc_mem_core_' ...
                expt_test{m} '_k_' num2str(k) '.mat'],'data');
            % noncore training
            cell_indx = setdiff(1:num_node,cell_indx);
            data = full_data(:,cell_indx);
            save([data_path expt_name{n} '_' expt_ee{e} '_cc_mem_noncore_'...
                expt_test{m} '_k_' num2str(k) '.mat'],'data');

            %---------- loopy mem ----------%
            % core training
            cell_indx = core_mem_loopy;
            data = full_data(:,cell_indx);
            save([data_path expt_name{n} '_' expt_ee{e} '_loopy_mem_core_'...
                expt_test{m} '_k_' num2str(k) '.mat'],'data');
            % noncore training
            cell_indx = setdiff(1:num_node,cell_indx);
            data = full_data(:,cell_indx);
            save([data_path expt_name{n} '_' expt_ee{e} '_loopy_mem_noncore_'...
                expt_test{m} '_k_' num2str(k) '.mat'],'data');

            %---------- cc centrality ----------%
            % core training
            cell_indx = core_cent_cc;
            data = full_data(:,cell_indx);
            save([data_path expt_name{n} '_' expt_ee{e} '_cc_cent_core_'...
                expt_test{m} '_k_' num2str(k) '.mat'],'data');
            % noncore training
            cell_indx = setdiff(1:num_node,cell_indx);
            data = full_data(:,cell_indx);
            save([data_path expt_name{n} '_' expt_ee{e} '_cc_cent_noncore_'...
                expt_test{m} '_k_' num2str(k) '.mat'],'data');

            %---------- loopy centrality ----------%
            % core training
            cell_indx = core_cent_loopy;
            data = full_data(:,cell_indx);
            save([data_path expt_name{n} '_' expt_ee{e} '_loopy_cent_core_'...
                expt_test{m} '_k_' num2str(k) '.mat'],'data');
            % noncore training
            cell_indx = setdiff(1:num_node,cell_indx);
            data = full_data(:,cell_indx);
            save([data_path expt_name{n} '_' expt_ee{e} '_loopy_cent_noncore_'...
                expt_test{m} '_k_' num2str(k) '.mat'],'data');
        
        end
        
    end
    
end
