T = 100;

for t = 1:T
    %% Make data
    N = 100;
    n = 20;
    A = rand(N);

    [YTrue, CTrue] = blossom(A);

    %% blah blah blah
    [Asub, Ysub] = sampleSubmatching(A, YTrue, n);

    Asub
    Ysub

    %% Rematch
    % If this was done correctly, we should recover YsubTrue == Ysub
    [YsubTrue, ~] = blossom(Asub);
    assert(all(vec(YsubTrue == Ysub)));
end