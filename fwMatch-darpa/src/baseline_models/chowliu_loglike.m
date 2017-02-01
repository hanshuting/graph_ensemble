function [ avg_test_loglike, avg_train_loglike ] = chowliu_loglike( x_train, x_test )
    model = treegmFit(x_train);
    
    avg_train_loglike = sum(treegmLogprob(model, x_train)) / size(x_train,1);
    avg_test_loglike = sum(treegmLogprob(model, x_test)) / size(x_test,1);
end

