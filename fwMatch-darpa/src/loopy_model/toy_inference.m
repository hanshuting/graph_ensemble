
F = [1 4 3;2 5 4];
G = [1 5 3;3 7 5;2 6 4;4 8 6];
edge_list = [1 2;1 3; 2 3]';

% Observed variables - 1,2
sample = [1 0 1];
Ut = eye(3);
Vt = eye(3);
N = 3;
YN = [0 1; 1 0;0 1];
YNtrueFlat = max(YN,[],2);

% Change F and G
sub = -1000;
F = [sub 4     3
     2   sub 4];

G = [sub 3
     sub 5
     6   sub    
     8   sub];
edge_list = [1 3; 2 3]';
Ut = eye(3);
Vt = eye(2);
N = 3;

obj = IsingPredict(F, G, Ut, Vt, N, edge_list,'YNtrueFlat', YNtrueFlat);
prob = BatchFW(obj, 'MaxIter', 100,'printInterval', 10, 'linesearch', false);
prob.run();

logZ = -1*prob.obj.fval();

