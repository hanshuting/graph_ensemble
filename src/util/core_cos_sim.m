function [pred,sim_core,sim_thresh,sim_avg,acc,prc,rec] = core_cos_sim(core,data,true_label)
% calculate cosine similarity using ensembles
% data is num_frame-by-num_neuron

% noise quantile
qnoise = 0.6;

% coefficient of S.D.
c = 3;
% if nargin<3 || ifspont==1 % if the dataset contains spont activity
%     c = 3;
% else % if vis stim only
%     c = 5;
% end

num_node = size(data,2);
core_vec = zeros(num_node,1);
core_vec(core) = 1;
sim_core = 1-pdist2(data,core_vec','cosine')';
sim_avg = [mean(sim_core(true_label==0)),mean(sim_core(true_label==1))];
% sim_thresh = 3*std(sim_core(sim_core<=quantile(sim_core,qnoise)))+...
%     mean(sim_core(sim_core<=quantile(sim_core,qnoise)));
th = quantile(sim_core(:),qnoise);
sim_core_th = sim_core;
sim_core_th(sim_core>=th) = NaN;
sim_thresh = c*nanstd(sim_core_th(:))+nanmean(sim_core_th(:));
if isnan(sim_thresh)
    sim_thresh = 0;
end
pred = sim_core>sim_thresh;

% make sure the dimensions are the same
pred = reshape(pred,[],1);
true_label = reshape(true_label,[],1);

% prediction statistics
TP = sum(pred==1&true_label==1);
TN = sum(pred==0&true_label==0);
FP = sum(pred==1&true_label==0);
FN = sum(pred==0&true_label==1);
acc = (TP+TN)/(TP+TN+FN+FP);
prc = TP/(TP+FP);
rec = TP/(TP+FN);

end
