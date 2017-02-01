function x = tf_table(N,nvals)

x = [];
for n=1:N
    mine = [];
    for val=1:nvals
        mine = [mine; val*ones(nvals^(n-1),1)];
    end
    x = [repmat(x,nvals,1) mine];
    %x = [repmat(x,2,1) [zeros(2^n,1); ones(2^n,1)]];
end
x = x';
end