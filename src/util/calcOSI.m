function [OSI,OSIstim] = calcOSI(data,stim)
% NOTE: this code doesn't really calculate R_ortho correctly at this moment
% this function calculates the orientation selectivity index of neurons
% INPUT:
%     data: N-by-T spike matrix
%     stim: T-by-1 stimulation vector, 0 should represent no stim
% OUTPUT:
%     OSI: N-by-1 vector
%     OSIstim: N-by-1 vector, preferred stimulation
% 
% Shuting Han, 2016

num_cell = size(data,1);

stim_seq = unique(stim);
stim_seq = setdiff(stim_seq,0);
num_stim = length(stim_seq);

% get averaged response for each stimulation
rstim = zeros(num_cell,num_stim);
for i = 1:num_stim
    rstim(:,i) = sum(data(:,stim==stim_seq(i)),2)/sum(stim==stim_seq(i));
end

[~,OSIstim] = max(rstim,[],2);

OSI = zeros(num_cell,1);
for i = 1:num_cell
    r_ortho = sum(rstim(i,setdiff(stim_seq,OSIstim(i))));
    OSI(i) = (rstim(i,OSIstim(i))-r_ortho)/(rstim(i,OSIstim(i))+r_ortho);
end

end