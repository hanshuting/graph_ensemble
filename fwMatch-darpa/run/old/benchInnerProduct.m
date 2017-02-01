%% Parameter setup
N = 1e6;
M = 100;
T = 100;

chunkSz = 1e4;

beginRow = 1:chunkSz:N;
endRow   = chunkSz:chunkSz:N;
nChunks = N / chunkSz;

times = zeros(2,2,T);
% times(1,1) = no transpose, cached idxs
% times(1,2) = no transpose, computed idxs
% times(2,1) = transpose, cached idxs
% times(2,2) = transpose, computed idxs

%% Run them all
for t = 1:T
    %% One iter
    U = rand(N, M);
    v = rand(chunkSz, 4);
    i = randi(nChunks);
    
    tic;
    y11 = U(beginRow(i):endRow(i),:).' * v;
    times(1,1,t) = toc;
    
    tic;
    y12 = U(chunkSz*(i-1)+1:chunkSz*i,:).' * v;
    times(1,2,t) = toc;
    
    tic;
    y21 = v' * U(beginRow(i):endRow(i),:);
    times(2,1,t) = toc;
    
    tic;
    y22 = v' * U(chunkSz*(i-1)+1:chunkSz*i,:);
    times(2,2,t) = toc;
    
end
