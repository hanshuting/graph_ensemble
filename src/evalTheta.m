function err = evalTheta(theta, features, Ytrue)
%evalTheta  Evaluate Hamming error of theta on features and Ytrue.

[M, N] = size(Ytrue);
K = size(features, 2);
D = size(features{1,1}, 1);

theta

assert(size(Ytrue, 1) == M, 'features and Ytrue do not have the same # of samples.');

err = 0;

for m = 1:M    
    % Construct specific A matrix for sample m        
    A = zeros(D, D);    
    for k = 1:K
        A = A + theta(k) .* features{m,k};
    end

    A
    
    % Our model requires the MAP assignment (maximization), while the
    % subrountine performs minimization.
    Yhat = csaAssignPerm(-A);

    Yhat
    Ytrue(m,:)
    
    err = err + sum(Yhat ~= Ytrue(m,:));
end

% Normalized Hamming error
err = err / (D * M);

end

