function [cum_hist] = calc_cum_dist(data,bins)
% calculate cumulative distribution with given data and bins
% currently only returns the mean across samples

num_sample = size(data,2);
cum_hist_all = zeros(length(bins),num_sample);

for i = 1:num_sample
    cum_hist_all(:,i) = histc(data(~isnan(data(:,i)),i),bins);
    cum_hist_all(:,i) = cum_hist_all(:,i)/sum(cum_hist_all(:,i));
    cum_hist_all(:,i) = cumsum(cum_hist_all(:,i),'reverse');
end

cum_hist = mean(cum_hist_all,2);

end