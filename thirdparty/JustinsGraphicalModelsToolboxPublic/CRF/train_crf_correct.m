function [p, history] = train_crf_correct(feats,efeats,labels,models,loss_spec,crf_type,options)

% required inputs:
%   features  - a cell array of features (should be nvars x nfeats)
%   efeats    - a cell array of edge features
%   labels    - a cell array of true labels
%   models    - a model structure defining the graph (or a cell array if desired)
%   (the above three can be a function taking an integer input instead of a cell array)
%   loss_spec - a string specifying what loss to use
%   crf_type  - a string specifying the type of crf to train

% optional inputs: (put in a structure array as options.nvals=5, etc.)
%   nvals   - an integer value specifying the number of output values  (default from labels)
%   maxtime - the maximum number of time to spend training (in seconds) (default inf)
%   maxiter - the maximum number of training iterations to use (passed to minFunc) (default 1000)
%   reg     - the amount of ridge regularization to apply (vs. loss per pixel) (default 0)
%   viz     - a function that takes a cell array of univariate marginals and displays them (default none)
%   rho     - edge appearance probability for TRW (one for all edges) default 1
%   print_times - binary if you want to print stats on speed
%   gradual     - binary if you want to gradually increase data and decrease iterations
%   derivative_check - turn derivative check on

% outputs
%   p - a structure array with the given parameters

% author: justin domke
% email (my first name).(my last name)@rit.edu

if isa(feats,'function_handle')
    if ~isfield(options,'N')
        error('if giving a function handle for labels, must provide N as optional input');
    end
    N = options.N;
else
    N = length(feats);
end

% could do lots more input checking here...
% model sizes fit with inputs
% labels are integers
% labels obey nvals

% default values for optional parameters
if ~isfield(options,'nvals')
nvals = 0; for i=1:N, label = fun_or_cell(labels,i); nvals = max(nvals,round(max(label(:)))); end
end
maxtime = inf;
maxiter = 1000;
reg     = 0;
viz     = @(marg) [];
rho     = 1;
gradual = 0;
%N       = nan; % no meaningful default
print_times = 0;
opt_display = 0;
opt_TolX    = 1e-7;
opt_TolFun  = 1e-7;
derivative_check = 0;
checkpoint_interval = false;
checkpoint_file     = false;
save_beliefs        = false;
p = [];

new_scale = 1;

% overwrite with user inputs
parseoptions(options);
derivative_check;

% if only a single model, repeat it
if ~isa(models,'function_handle') && ~iscell(models)
    models = repmat({models},N,1);
end

% if no edge features, put in constants
if isempty(efeats)
    for n=1:N
        npairs = size(models{n}.pairs,1);
        efeats{n} = ones(npairs,1);
    end
end

feat   = fun_or_cell(feats ,1);
efeat  = fun_or_cell(efeats,1);
nfeat  = size( feat,2);
nefeat = size(efeat,2);

if strcmp(crf_type,'linear_linear')
    if isempty(p)
        p.F = zeros(nvals  ,nfeat);
        p.G = zeros(nvals^2,nefeat);
    end
else
    error('unsupported crf_type')
end

time0  = tic;
ncalls = 0;

% *KT* Instrument function (called once per minFunc iter, to be fair.)
history = struct('iter', {}, ...
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

instrument_b_i = {};
instrument_b_ij = {}; 
function outputFcn(x,optimValues,state,varargin)
    pause(instrumentTime);
    if strcmp(state, 'iter')
        x_unpacked = unpackstruct(x, p);

        if save_beliefs
            b_i  = cell2mat(instrument_b_i);
            b_ij = cell2mat(instrument_b_ij);
        else
            b_i  = [];
            b_ij = [];
        end
        
        history(end+1) = struct('iter', optimValues.iteration, ...
            'x', x, ...
            'b_i', b_i, ...
            'b_ij', b_ij, ...
            'F', x_unpacked.F, ...
            'G', x_unpacked.G, ...
            'fval', optimValues.fval, ...
            'firstorderopt', optimValues.firstorderopt, ...
            'totTime', toc(instrumentTime));

        if checkpoint_interval > 0 && ...
            mod(optimValues.iteration, checkpoint_interval) == 0
       
            fprintf('Saving checkpoint to %s.\n', checkpoint_file);
            save('-v7.3', checkpoint_file);            
        end
    end
    start(instrumentTime);
end

function [L grad E] = eval(params)
    if toc(time0)>maxtime
        L = 0;
        grad = zeros(size(params));
        E = 0;
        return
    end
        
    p  = unpackstruct(params  ,p);
    
    grad = 0*params;
    
    T = zeros(N,1);
    L = zeros(N,1);
    E = zeros(N,1);
    
    % parfor wants this pre-declared for some reason
    b_i     = cell(N,1);
    
    % Clear the instrumentation variables (to save the new beliefs)
    instrument_b_i  = cell(N,1);
    instrument_b_ij = cell(N,1);
    
%     parfor n=1:N
    for n=1:N
        g = unpackstruct(params*0,p);
        
        model = models{n};

        % streamlined -- has to be structs, not functions
        y = labels{n};
        x = feats{n};
        z = efeats{n};
                
        if strcmp(crf_type,'linear_linear')
            [L(n) b_ij b_i{n} g.F g.G] = crf_linear_linear(model,p.F,p.G,x,z,y,rho,loss_spec);            
            instrument_b_i{n}  = b_i{n}.';
            instrument_b_ij{n} = b_ij.';
        else
            error('unsupported crf_type')
        end
        
        grad(:,n) = packstruct(g);

        [~,ypred] = max(b_i{n},[],1);
        E(n) = sum(ypred(y>0)' ~= y(y>0));
        T(n) = sum(double(y(:)>0));
    end
    
    ncalls = ncalls + 1;
    
    viz(b_i);
        
    E    = sum(E);
    L    = sum(L);
    grad = sum(grad,2);
    
    % numerical hardening-- if you get a nan, return a huge loss, hopefully
    % force backtracking
    if isbad(L) || isbad(grad)
        L    = 2^100;
        grad = 0*params;
    end
    
    % regularization
    L    = new_scale * (L    + 0.5*reg*sum(params.^2));
    grad = new_scale * (grad + reg*params);

    if print_times
        fprintf('time: %f / %f  (%f per call)  loss: %f  error: %f \n', ...
            toc(time0), maxtime, toc(time0)/ncalls, L, E);
        disp('First 10 columns of F and its gradient');
        ncol = min(10, size(p.F, 2));
        p.F(:,1:ncol)
        g = unpackstruct(grad, p);
        g.F(:,1:ncol)
    end
    
end

opt_options = [];
opt_options.LargeScale = 'off';
% opt_options.HessUpdate = 'bfgs';
opt_options.HessUpdate = 'sd';
opt_options.MaxIter = maxiter;
opt_options.DerivativeCheck = derivative_check;
opt_options.Display = opt_display;
opt_options.TolX    = opt_TolX;
opt_options.TolFun  = opt_TolFun;
opt_options.outputFcn = @outputFcn;

params = packstruct(p);

if gradual
    N0 = N;
    maxiter0 = maxiter;
    N = 16;
    while 1
        N = 2*N;
        if N>N0, N = N0; end
        maxiter = max(5,ceil(maxiter/2));
        
        fprintf('N: %d  maxiter: %d  time: %f \n', N, maxiter, toc(time0));
        opt_options.MaxIter = maxiter;
        params = minFunc(@eval,params,opt_options);
        if N==N0
            break;
        end
    end
else
    params = minFunc(@eval,params,opt_options);
end
    
p = unpackstruct(params,p);

end

function parseoptions(options)
if exist('options','var')
    names = fieldnames(options);
    for i=1:length(names)
        name = names{i};
        val  = options.(name);
%         if ~(strcmp(name,'maxtime') ...
%             || strcmp(name,'maxiter') ...
%             || strcmp(name,'nvals') ...
%             || strcmp(name,'reg') ...
%             || strcmp(name,'rho') ...
%             || strcmp(name,'print_times') ...
%             || strcmp(name,'N') ...
%             || strcmp(name,'gradual') ...
%             || strcmp(name,'derivative_check') ...
%             || strcmp(name,'opt_display') ...
%             || strcmp(name,'opt_TolX') ...
%             || strcmp(name,'viz')) ...
%             || strcmp(name,'checkpoint_file') ...
%             || strcmp(name,'checkpoint_interval') ...
%             || strcmp(name,'add_score') ...
%             error(['unrecognized option: ' name]);
%         end

        % Why bother checking?
        assignin('caller',name,val);
    end
end
end

function out = fun_or_cell(fun,n)
if iscell(fun)
    out = fun{n};
elseif isa(fun,'function_handle')
    out = fun(n);
else
    error('fun must be a cell or a function')
end
end
