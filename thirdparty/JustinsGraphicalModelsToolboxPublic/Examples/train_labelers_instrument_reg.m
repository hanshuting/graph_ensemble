function [params times fvals histMinFunc testResult] = train_labelers_instrument_reg(loss_spec,params0,testmode,...
    imsdir,labelsdir,nvals,feat_params,edge_params,maxims,reg,checkpoint_interval,checkpoint_file,maxtime)
% *KT* note: callers experiment_binarydenoising
% *KT* note: times measures cumulative time per objective function
%            evaluation. This is NOT the same as an outer loop (LBFGS)
%            iteration, since the function may be called several times per
%            iteration. (But it usually not too often, thanks to being
%            Newton.)
% *KT* note: I added histMinFunc which is updated once per minFunc iteration.
%
%            QUESTION: Can linesearch cause params to change?
%            NOTE 1: I store the *structured* parameters.
%            NOTE 2: eval loops through Aim -- all images. (Though it's all
%            somehow cached so I won't penalize data access?) So once
%            *evaluation* of eval should be considered an effective pass.
%               ^^ confused mumbo-jumbo; ignore. histMinFunc is one entry
%               per ITERATION of minFunc. That's right. ^^
%
% *KT* note: testResult.L = final log likelihood
%                      .E = MPM error.
%
% loss spec could be
% trunc_em_5
% pert_ul_1e-5
% ...
% whatever

if nargin < 3
    testmode = 0;
end

rescale   = .5;%1;
rescale   = 1;
rescale   = 1;
cachedata = 1; % load all ims, precompute feats, store in memory

Aims    = [dir([imsdir    '/*.png']);dir([imsdir    '/*.jpg']);dir([imsdir    '/*.bmp'])];
Alabels = [dir([labelsdir '/*.png']);dir([labelsdir '/*.jpg']);dir([labelsdir '/*.bmp']);];

if exist('maxims','var')
    Aims    = Aims   (1:maxims);
    Alabels = Alabels(1:maxims);
end

if ~exist('maxtime','var')
    maxtime = inf;
end

M = 1:length(Aims);

%reg = 1e-3;
%if strcmp(loss_spec(1:6),'pseudo')
%    reg = 1e-4;
%elseif strcmp(loss_spec(1:6),'piecewise')
%    reg = 1e-4;
%elseif strcmp(loss_spec(end-1:end),'_0')
%    reg = 1e-4;
%end
%reg = 1e-5;

% store all models
modelAR  = cell(1000,1000);
modelZ   = cell(length(Aims),1);
for n=1:length(Aims)
    im    = imread([imsdir    '/' Aims(n).name]);
    im    = double(im)/255;
    if size(im,3)==1, im=repmat(im,[1 1 3]); end;
    label = double(imread([labelsdir '/' Alabels(n).name]));
        
    
%     'get rid of'
    label(label==0) = 1;
    
    if rescale ~= 1
        im    = imresize(im   ,rescale,'bilinear');
        label = imresize(label,rescale,'nearest');
    end
    
    [ly lx] = size(label);
    if isempty(modelAR{ly,lx}) % make models at reduced resolution
        %fprintf('generating a model of size %d x %d ...\n',ly,lx);
        modelAR{ly,lx} = gridmodel(ly,lx,nvals);
    end
    modelZ{n} = modelAR{ly,lx};
    
    % compute features on reduced resolution image
    if n==1
        feats  = featurize_im(im,feat_params);
        nfeat  = size(feats,3);
        [efeats enames] = edgeify_im(im,edge_params,modelAR{ly,lx}.pairs,modelAR{ly,lx}.pairtype);
        nefeat = size(efeats,2);
    end
    
end

% parfor wants these pre-defined.
IM     = cell(length(Aims),1);
LABEL  = cell(length(Aims),1);
FEATS  = cell(length(Aims),1);
EFEATS = cell(length(Aims),1);
if cachedata
    %fprintf('precomputing and storing all features, images...\n');
    parfor n=1:length(Aims)
        imm    = imread([imsdir    '/' Aims(n).name]);
        imm    = double(imm)/255;
        if size(imm,3)==1, imm=repmat(imm,[1 1 3]); end;
        labell  = double(imread([labelsdir '/' Alabels(n).name]));
%         'get rid of'
        labell(labell==0) = 1;
        if rescale ~= 1
            imm     = imresize(imm    ,rescale,'bilinear');
            labell  = imresize(labell ,rescale,'nearest');
            %featss  = imresize(featss ,rescale,'bilinear');
            %efeatss = imresize(efeatss,rescale,'bilinear');
        end
        [ly lx lz] = size(imm);
        %model   = modelAR{ly,lx};
        model    = modelZ{n};
        featss  = featurize_im(imm,feat_params);
        efeatss = edgeify_im(imm,edge_params,model.pairs,model.pairtype);
        
        IM{n}     = imm;
        LABEL{n}  = labell;
        FEATS{n}  = featss;
        EFEATS{n} = efeatss;
    end
end

%Gadj = 1;


% *KT* QUESTION: When are the values ever used? Whatever... no matter
%make parameters
p.F = 0*.1*randn(nvals,nfeat);
%p.G = .0001*randn(nvals*nvals,2);
p.G = 0*.0001*randn(nvals*nvals,nefeat);

% *KT* Instrument function (called once per minFunc iter, to be fair.)

histMinFunc = struct('iter', {}, ...
                     'x', {}, ...
                     'b_i', {}, ...
                     'b_ij', {}, ...
                     'F', {}, ...
                     'G', {}, ...
                     'fval', {}, ...
                     'firstorderopt', {}, ...
                     'totTime', {});
                 
instrumentTime = timetic;
tic(instrumentTime);

instrument_b_i = [];
instrument_b_ij = []; 
function outputFcn(x,optimValues,state,varargin)
    pause(instrumentTime);
    if strcmp(state, 'iter')
        x_unpacked = unpackstruct(x, p);
        
        histMinFunc(end+1) = struct('iter', optimValues.iteration, ...
            'x', x, ...
            'b_i', instrument_b_i, ...
            'b_ij', instrument_b_ij, ...
            'F', x_unpacked.F, ...
            'G', x_unpacked.G, ...
            'fval', optimValues.fval, ...
            'firstorderopt', optimValues.firstorderopt, ...
            'totTime', toc(instrumentTime));

        if mod(optimValues.iteration, checkpoint_interval) == 0
            fprintf('Saving checkpoint to %s.\n', checkpoint_file);
            save('-v7.3', checkpoint_file);
        end
    end
    start(instrumentTime);
end

times = [];
fvals = [];

time0 = tic;
function [L grad E] = eval(params,N,nograd)
    if toc(time0)>maxtime
        L = 0;
        grad = zeros(size(params));
        E = 0;
        return
    end
    
    if nargin < 3
        nograd = 0;
    end
    
    p  = unpackstruct(params  ,p);
    
    grad = 0*params;
        
    if nargin < 2
        N = 1:length(Aims);
    end
    
    T = zeros(length(N),1);
    L = zeros(length(N),1);
    E = zeros(length(N),1);
    
    % parfor wants these pre-declared for some reason
    b_out     = cell(length(N));
    im_out    = cell(length(N));
    label_out = cell(length(N));
    
%     parfor n=1:length(Aims)
    % *KT* removed parfor for accurate single core timings.
    for n=1:length(Aims)    
    %for i=1:length(N)
        
        g = unpackstruct(params*0,p);
        
        if cachedata==1
            im     = IM{n};
            label  = LABEL{n};
            feats  = FEATS{n};
            efeats = EFEATS{n};
        else
            im    = imread([imsdir    '/' Aims(n).name]);
            im    = double(im)/255;
            if size(im,3)==1, im=repmat(im,[1 1 3]); end;
            label  = double(imread([labelsdir '/' Alabels(n).name]));
            if rescale ~= 1
                im     = imresize(im    ,rescale,'bilinear');
                label  = imresize(label ,rescale,'nearest');
                %feats  = imresize(feats ,rescale,'bilinear');
                %efeats = imresize(efeats,rescale,'bilinear');
            end
            [ly lx] = size(label)
            %model   = modelAR{ly,lx};
            model   = modelZ{n};
            feats   = featurize_im(im,feat_params);
            efeats  = edgeify_im(im,edge_params,model.pairs,model.pairtype);
        end
        
        [ly lx] = size(label);
        %model = modelAR{ly,lx};
        model = modelZ{n};
        rho   = .5;
        
%         efeats2 = edgeify_im(im,edge_params,model.pairs,model.pairtype);
%         norm(efeats-efeats2)
        
        %x = featurize_im(im,feat_params);
        x = feats;
        x = reshape(x,size(x,1)*size(x,2),size(x,3));
        y = label(:);
        z = efeats;
                
        %[myL b_ij b_i g.F g.G] = crf_linear_independent_pairtypes(model,p.F,p.G,x,y,rho,loss_spec);
        if nograd
            %tic
            [myL b_ij b_i] = crf_linear_linear(model,p.F,p.G,x,z,y,rho,loss_spec);
            %toc
        else
            %tic
            [myL b_ij b_i g.F g.G] = crf_linear_linear(model,p.F,p.G,x,z,y,rho,loss_spec);
            %toc
        end
        
        % copy to the outside function so that our instrument function
        % outputFcn the next time it is called
        instrument_b_ij = b_ij;
        instrument_b_i  = b_i;
        
        %g.G = g.G/Gadj;
        
        if 1%ischar(testmode) || mod(length(fvals)+1,10)==0
            if (strcmp(loss_spec,'pseudo') || strcmp(loss_spec,'piecewise') || strcmp(loss_spec,'piecewise_shared'))
            % use TRW for visualization purposes
            [~, ~, b_i] = crf_linear_linear(model,p.F,p.G,x,z,y,rho,'em_trw_5');
            %[~, ~, b_i] = crf_linear_linear(model,p.F,p.G,x,z,y,rho,'em_mnf_1e-5');
            end
        end
        
        L(n)      = myL          ;
        grad(:,n) = packstruct(g);
        
        [~,ypred] = max(b_i,[],1);
        E(n) = sum(ypred(y>0)' ~= y(y>0));
        T(n) = sum(double(y>0));
        
        if n==1 || ischar(testmode)
            b_out{n} = zeros(ly,lx,nvals);
            for val=1:nvals
                b_out{n}(:,:,val) = reshape(b_i(val,:),ly,lx);
            end
            im_out{n} = im;
            label_out{n} = label;
        end
        
        
        % *KT*: I don't care about visuals!
%         if testmode
%             testmode_stuff = regexp(testmode,'_','split');
%             lab1 = testmode_stuff{1};
%             warning off all
%             colormap gray
%             imwrite(im_out{n}                    ,sprintf('rez/im_%d_%s.png' ,n,lab1));
%             imwrite(miximshow(label_out{n},nvals),sprintf('rez/lab_%d_%s.png',n,lab1));
%             imwrite(miximshow(b_out{n},nvals)    ,sprintf('rez/b_%d_%s.png'  ,n,testmode));
%             warning on all
%         end
    end
    
    %numel(label)
    E    = sum(E)/sum(T);
    L    = sum(L)/sum(T);
    grad = sum(grad,2)/sum(T);
    
    % numerical hardening-- if you get a nan, return a huge loss, hopefully
    % force backtracking
    if isbad(L) || isbad(grad)
        L    = 2^100;
        grad = 0*params;
    end
    
    scaling = p;
    scaling.F = 1 + 0*p.F;
    %scaling.G = 1./Gadj + 0*p.G;
    scaling.G = 1 + 0*p.G;
    scaling = packstruct(scaling);
    
    %reg  = 0;%1e-3;
    L    = L    + reg*sum((params.*scaling).^2);
    grad = grad + 2*reg*(params.*scaling.^2);

    fprintf('sum(T) = %g\n', sum(T));
    
%     L   = L + reg*sum(p.G(:).^2);
%     g   = unpackstruct(params*0,p);
%     g.G = g.G + 2*reg*p.G;
%     grad = grad + packstruct(g);

    times = [times; toc(time0)];
    fvals = [fvals; L];
    
    if mod(length(fvals),10)==0
        % *KT*: TODO: Hook in your iterate extraction tracker. (Or even
        % outside.)
        
        % comment out these lines to run on a machine with no display...
%         figure(1), imshow(im_out{1})
%         figure(2), imshow(miximshow(label_out{1},nvals))
%         figure(3), imshow(miximshow(b_out{1},nvals)); title(E)
%         figure(4), plot(params)
%         figure(5),
%         fparams = p.G(1,:);
%         D = 0:.01:2;
%         %r=.001;
%         %diff = r./(r+D);
%         diff = D;
%         R = 0*diff;
%         maxk = length(fparams)/2-1;
%         where = 1;
%         for n=0:maxk
%             R = R + fparams(where)*cos(pi*n*diff);
%             where = where+1;
%             R = R + fparams(where)*sin(pi*n*diff);
%             where = where+1;
%         end    
%         plot(D,R);
%         drawnow
    end

    
%     g = unpackstruct(grad,p);
%     g.F = 0*g.F;
%     grad = packstruct(g);
end

if nargin>=2 && ~isempty(params0)
    params = params0;
    %p0     = unpackstruct(params,p);
    %p0.G   = p0.G/50;
    %params = packstruct(p0);
else
    params = packstruct(p);
end

% this function fits ignoring G completely
    function [L grad E] = evalF(params)
        p = unpackstruct(params,p);
        p.G = 0*p.G;
        [L grad E] = eval(packstruct(p));
        g = unpackstruct(grad,p);
        g.G = 0*g.G;
        grad = packstruct(g);
    end

%reg = 1e-5;
if ~testmode
    p2 = p;
    p2.F =   1 + 0*p.F;
    p2.G =   1 + 0*p.G;
    scaling = packstruct(p2);
    
    options = [];
    options.Method = 'lbfgs';
    % backtracking line search-- used in all experiments except when explicitly changed for
    % perturbation vs. backpropagation vs. truncation experiments
    options.LS      = 2;
    %options.LS      = 4;
    %options.LS_init = 1;
    options.TolX    = 1e-7;
    options.MaxIter = 1e6/1000;
    options.MaxFunEvals = 1e6/1000;
    options.outputFcn = @outputFcn;
    
    % *KT* NOTE: This minimizes *evalF*, **not** eval! Fitting independent
    % nodes is easy. We use this value to initialize minimization over
    % eval.
    if isempty(params0)
        params      = minFunc(@evalF,params,options);
        p  = unpackstruct(params,p);
        p0 = p;
        save p0 p0
    end
    
%     p0 = [];
%     load p0;
%     p.F = p0.F;
%     params = packstruct(p);

    % *KT* NOTE: What the hell, why do we restart? No matter... our
    % instrumenting function will still log everything. Although the
    % iteration counter at the end may be off. BEWARE!
    %
    % Remember: The instrument function (outputFcn) writes everything
    % itself. You don't have to worry about minFunc being called twice...
    % you already forgot this point!
    
    % restart LBFGS 
    params_min = params;
    for retrys=1:1
        if toc(time0)<maxtime
            params_min  = minFunc(@eval ,params_min,options);
        end
    end
    
    params = params_min;
end

if testmode
% No regularization! So L *is* actually plain likelihood.
reg = 0;
[testResult.L grad testResult.E] = eval(params,1:length(Aims),1); % 1 turns off gradient computation
fprintf('loss_spec: %s  L: %f  E: %f\n', loss_spec, testResult.L, testResult.E); %E/length(Aims)
beep; pause(1); beep; pause(1); beep
end

if testmode
    params = [testResult.L testResult.E];
end

end
