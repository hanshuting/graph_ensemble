%% Blah
D = 6;
M = factorial(D);

for m = 1:M
    p = oneperm(D, m);
    m2 = rankperm(p);
    assert(m == m2, 'Inversion fail');
end
