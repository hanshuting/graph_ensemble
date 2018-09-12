function [who_train who_test] = kfold_sets(N,K,k)

who_test = K-k+1:K:N;

who = 1:N;
who(who_test) = 0;

who_train = who(who>0);