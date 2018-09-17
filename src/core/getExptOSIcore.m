function [core_osi] = getExptOSIcore(Spikes,vis_stim)
% Spikes is num_neuron-by-num_frame

OSI_thresh = 0.4;

% calculate OSI
[OSI,OSIstim] = calcOSI(Spikes,vis_stim);

for i = 1:length(unique(vis_stim))-1
    core_osi = find((OSI>OSI_thresh)&(OSIstim==i));
%         save([savepath expt_name{n} '_core_OSI_' num2str(i) '.mat'],'core_osi');
end
    
end    