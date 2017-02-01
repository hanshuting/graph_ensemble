%% Test it
D = 5;
A = rand(D);
B = rand(D);

a1 = sum(vec(A .* B))
a2 = trace(A' * B)
