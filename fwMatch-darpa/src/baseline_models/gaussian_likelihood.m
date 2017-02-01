function log_likelihood = gaussian_likelihood(x_train, x_test)

covariance = cov(x_train);
det_covariance = det(covariance);

[num_samples, num_var] = size(x_test);
mu = repmat(mean(x_train),[num_samples,1]);

log_likelihood = mean(-0.5*diag(((x_test - mu)/covariance) * (x_test - mu)') ...
    - (num_var/2)*log(2*pi) - 0.5*det_covariance);

end