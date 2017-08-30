
function [] = makeShuffledData(param)

expt_name = param.expt_name;
ee = param.ee;
num_shuff = param.num_shuff;
data_path = param.data_path;
shuff_path_base = param.shuff_path_base;

% make shuffled data
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    load([data_path expt_name{n} '\' expt_name{n} '.mat']);
    num_node = size(Spikes,1);
    
    for e = 1:length(expt_ee)
        
        fprintf('processing %s_%s...\n',expt_name{n},expt_ee{e});
        
        % set save path
        shuff_path = [shuff_path_base 'shuffled_' ...
            expt_name{n} '_' expt_ee{e} '_loopy\'];
        if ~exist(shuff_path)
            mkdir(shuff_path);
        end
        
        % load data
        load([data_path expt_name{n} '\' expt_name{n} '_' expt_ee{e} '.mat']);
        data_raw = data;
        
        % shuffle
        for i = 1:num_shuff
            data = shuffle(data_raw','exchange')';
            save([shuff_path 'shuffled_' expt_name{n} '_' expt_ee{e} '_'...
                num2str(i) '.mat'],'data');
        end
        
    end
end

