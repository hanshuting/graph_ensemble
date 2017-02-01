function err = evalThetaCellY2(theta, features, Ytrue)
%evalTheta  Evaluate Hamming error of theta on features and Ytrue, cell Ytrue.

M = length(Ytrue);
K = size(features, 2);
D = size(features{1,1}, 1);

err = 0;

for m = 1:M    
    % Construct specific A matrix for sample m        
    A = zeros(D, D);    
    for k = 1:K
        A = A + theta(k) .* features{m,k};
    end   
    
    % Our model requires the MAP assignment (maximization), while the
    % subrountine performs minimization.
    %Yhat = csaAssignPerm(-A);
    [assignment,cost] = munkres(-A);
    Yhat = assignment;

    err = err + sum(Yhat ~= Ytrue{m});
end

% Normalized Hamming error
err = err / (D * M);

end

