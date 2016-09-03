
% load data
load('C:\Shuting\fwMatch\data\m21_d2_vis.mat');
load('C:\Shuting\fwMatch\data\ensembles\Core_m21_d2_vis.mat');
load('C:\Shuting\fwMatch\data\ensembles\m21_d2_vis_edos.mat');

% sort by ensemble index
num_frame = size(Spikes,2);
num_ensemble = size(Pools_coords,3);
eindx = cell(num_ensemble,1);
spikes_sorted = [];
acm = zeros(num_ensemble+1,1);
for i = 1:num_ensemble
    single_eindx = squeeze(Pools_coords(:,3,i));
    single_eindx = single_eindx(single_eindx~=0);
    eindx{i} = single_eindx;
    acm(i+1) = acm(i)+length(single_eindx);
    spikes_sorted(acm(i)+1:acm(i+1),:) = Spikes(single_eindx,:);
end
eindx_all = cell2mat(eindx);
oindx = setdiff(1:size(Spikes,1),eindx_all);
acm(end+1) = acm(end)+length(oindx);
spikes_sorted(end+1:end+length(oindx),:) = Spikes(oindx,:);

% plot
h = figure;set(h,'color','w');
imagesc(~spikes_sorted);colormap(gray);
hold on;
for i = 1:num_ensemble+1
    plot([0 num_frame],[acm(i),acm(i)],'b--');
end
set(gca,'tickdir','out');
xlabel('frame number')
ylabel('cell number')

% save figure
saveas(h,'C:\Shuting\fwMatch\paper\figures\matlab_fig\1e');
saveas(h,'C:\Shuting\fwMatch\paper\figures\matlab_fig\1e.pdf');

% plot stimuli
figure;set(gcf,'color','w')
imagesc(vis_stim'==0);colormap(gray);
set(gca,'xtick',[],'ytick',[]);
saveas(gcf,'C:\Shuting\fwMatch\paper\figures\matlab_fig\1e_stim');
saveas(gcf,'C:\Shuting\fwMatch\paper\figures\matlab_fig\1e_stim.pdf');


