% extract svd
fname = 'pa_511510855_TF1';
load(['C:\Shuting\graph_ensemble\data\pa_svd\SVD_' fname '.mat']);

seq = [7,3,2,1];
svd_state = {'1','2','3','4'};
core_svd = cell(4,1);
for n = 1:length(seq)
    core_svd{n} = Pools_coords(:,3,seq(n));
    core_svd{n} = core_svd{n}(core_svd{n}~=0);
end
disp(core_svd);
save(['C:\Shuting\graph_ensemble\data\ensembles\' fname '_core_svd.mat'],...
    'core_svd','svd_state');