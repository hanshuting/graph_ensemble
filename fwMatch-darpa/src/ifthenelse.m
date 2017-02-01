function y = ifthenelse(pred, iftrue, iffalse)
% ifthenelse  Functional if, for lambdas.
    if pred
        y = iftrue();
    else
        y = iffalse();
    end
end

