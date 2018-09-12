function A = A_bruteforce(model,theta_ij,theta_i,x0)

% computes exact (non-variational) log-partition function
X = tf_table(model.nnodes,model.nvals);
v = zeros(size(X,2),1);
where = 0;
for i=1:size(X,2)
    x = X(:,i);
    % if this is actually a point that obeys constraints
    if prod(double(x(find(x0))==x0(find(x0))))
        where = where + 1;
        for i=1:model.nnodes
            v(where) = v(where) + theta_i(x(i),i);
        end
        for n=1:model.ncliques
            i = model.pairs(n,1);
            j = model.pairs(n,2);
            who = x(i)+(x(j)-1)*model.nvals;
            v(where) = v(where) + theta_ij(who,n);
        end
    end
end
A = log_sum_exp(v(1:where));

end