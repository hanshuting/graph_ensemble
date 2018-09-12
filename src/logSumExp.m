function y = logSumExp(x, dim)
% y = logSumExp(x, [dim]) is equivalent to y = log(sum(x, [dim])) with extra precision.
% See http://machineintelligence.tumblr.com/post/4998477107/the-log-sum-exp-trick

    A = max(real(x(:)));
    if nargin == 1
        y = A + log(sum(exp(x - A)));
    else
        y = A + log(sum(exp(x - A), dim));
    end

end

