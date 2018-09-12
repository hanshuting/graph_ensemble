function dEnt = modOneDEntropyDerivative(beliefs, reweight)
riPlusRj = bsxfun(@plus, reweight, reweight');
dEnt = riPlusRj + log(beliefs) + (riPlusRj - 1) .* log(1 - beliefs);
% Technically this isn't right but whatever... maybe David Sontag's new
% paper has some way to deal with this?
if any(vec(beliefs == 0 | beliefs == 1))
    error('modOneDEntropyDerivative:boundary', 'Hit the boundary; unbounded gradient');
end
end
