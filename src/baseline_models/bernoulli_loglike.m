function [ avg_test_loglike, avg_train_loglike ] = bernoulli_loglike(x_train, x_test)
    mu = (sum(x_train) / size(x_train, 1));
    log_mu = log(mu + eps);
    log_one_minus_mu = log(1 - mu + eps);
    
    rep_train_log_mu = repmat(log_mu, size(x_train,1), 1);
    rep_train_log_one_minus_mu = repmat(log_one_minus_mu, size(x_train,1), 1);
    avg_train_loglike = sum(sum(rep_train_log_mu .* x_train + rep_train_log_one_minus_mu .* (1-x_train))) / size(x_train, 1);
    
    rep_test_log_mu = repmat(log_mu, size(x_test,1), 1);
    rep_test_log_one_minus_mu = repmat(log_one_minus_mu, size(x_test,1), 1);
    avg_test_loglike = sum(sum(rep_test_log_mu .* x_test + rep_test_log_one_minus_mu .* (1-x_test))) / size(x_test, 1);
end

