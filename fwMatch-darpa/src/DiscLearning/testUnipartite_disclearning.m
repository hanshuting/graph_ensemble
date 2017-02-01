clear;

addpath(genpath('../BCFWstruct/'));

%% Make fake data
M = 99;
D = 4;
K = 5;

% (Only the UT matters.)
Y1 = {sparse([1 3], [2 4], [1 1], 4, 4)};
Y2 = {sparse([1 2], [2 4], [1 1], 4, 4)};

Y = [Y1;Y2;Y1];
Ys = repmat(Y, M/3, 1);
Yss = cell(1,M);
Y1s = Y1{1,1};
Y2s = Y2{1,1};

for i = 1:M
   if(rand > 0.5); Y = Y1s; else Y = Y2s; end
   Yss{i} = vecUT(Y); 
end

%[ri, rj, rv] = findUT(rand(D));
feats = cell(M, K);

for i = 1:M
    for k = 1:K
        feats{i,k} = rand(D);
    end
end


obj = UnipartiteMatchingDiscLearning(Ys, feats);
patterns = cell(1,M);
for i = 1:M
   patterns{i} = i; 
end
   



 lossCB = @(param, y, ybar) sum(abs(y(:) - ybar(:)));

    
if(false)
    parm.patterns = patterns;
    parm.labels = Yss ;
    parm.lossFn = lossCB ;
    parm.constraintFn  = @constraintCB ;
    parm.featureFn = @featureCB ;
    parm.dimension = K ;
    parm.obj = obj;
    args = ' -c 1.0 -o 1 -v 1 ' ;
    model = svm_struct_learn(args, parm) ;
else
    
    param = [];
    param.patterns = patterns;
    param.labels = Yss;
    param.lossFn = lossCB;
    param.oracleFn = @constraintCB;
    param.featureFn = @featureCB;
    param.obj = obj;
    
    options = [];
    options.lambda = 1e-2;
    options.gap_threshold = 0.1; % duality gap stopping criterion
    options.num_passes = 100; % max number of passes through data
    options.do_line_search = 1;
    options.debug = 1; % for displaying more info (makes code about 3x slower)

   [model, progress] = solverBCFW(param, options);
%[model, progress] = solverFW(param, options);
%[model, progress] = solverSSG(param, options);


 
end

% trainer = StructPerceptronLearning(obj, 'printInterval', 5);
% trainer.run()