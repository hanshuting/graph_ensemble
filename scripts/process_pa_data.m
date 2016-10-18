% script for calculating binary spike matrix from PA dataset

dpath = 'C:\Shuting\graph_ensemble\data\pa_raw\';
ecid = {'ecid_511507650','ecid_511509529','ecid_511510650','ecid_511510718',...
    'ecid_511510855','ecid_511510670'};

qnoise = 0.05;
num_id = length(ecid);
for n = 1:num_id
    
    load([dpath ecid{n} '.mat']); % raw, corrected, dff, stim_info
    
    % get spikes
    dff_filt = lowpass_filter_ct2(dff);
    grad = gradient(dff_filt);
    th1 = quantile(grad(:),qnoise);
    th2 = quantile(grad(:),1-qnoise);
    grad_th = grad;
    grad_th(grad<=th1) = NaN;
    grad_th(grad>=th2) = NaN;
    baseline = 5*nanstd(grad_th(:))+nanmean(grad_th(:));
    Spikes = grad>baseline;
    
    ii = 10;
    subplot(2,1,1)
    hold off;plot(corrected(ii,1:10000));
    subplot(2,1,2)
    hold off;plot(grad(ii,1:10000));
    hold on;plot([1,10000],baseline*[1,1],'r');
    
    % extract stim vector
    vis_stim = zeros(1,size(grad,2));
    tf = zeros(1,size(grad,2));
    for ii = 1:length(stim_info.orientation)
        if ~isnan(stim_info.orientation(ii))
            vis_stim(stim_info.start(ii):stim_info.end(ii)) = stim_info.orientation(ii)/45+1;
            tf(stim_info.start(ii):stim_info.end(ii)) = stim_info.temporal_frequency(ii);
        else
            vis_stim(stim_info.start(ii):stim_info.end(ii)) = -1;
        end
    end
    
    % save result
    save([dpath ecid{n} '_spikes.mat'],'Spikes','vis_stim','tf','-v7.3');
    
end
