%% init with random matrix
D = 5;
M = 5000;
A = rand(5,5);

%% init with systematic matrix
D = 5;
M = 5000;
A = eye(5);


%% Code
hubers = sample_perms(exp(A), M);
naives = naiveSamplePerm(A, M);

rankhubers = rankperms(hubers);
ranknaives = rankperms(naives);

Nperm = factorial(D);

figure; hist(rankhubers, Nperm); title('Huber-Law Samples');
figure; hist(ranknaives, Nperm); title('Naive Samples');

[h, p, stat] = kstest2(rankhubers, ranknaives)
