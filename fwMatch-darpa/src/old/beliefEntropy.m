function h = beliefEntropy(B)

    h = sum(B(:) .* log(B(:)));

end

