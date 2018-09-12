function loss = testCRF(testSamps, theta, features)

lapjvTol = 1e-6;
loss = 0;

for m = 1:M
    A = zeros(D);
    for k = 1:K
        A = A - theta(k) * features{m,k};
    end
    
    Yhat = lapjv(A', lapjvTol);

    loss = loss + sum(Yhat ~= testSamps(m,:));    
end

loss = loss / (D * M);

end

