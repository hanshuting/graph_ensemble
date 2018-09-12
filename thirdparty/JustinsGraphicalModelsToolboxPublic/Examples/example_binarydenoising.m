function example_binarydenoising

N     = 1;
siz   = 50;
rho   = .5;
nvals = 2;

% make a graph
model = gridmodel(siz,siz,nvals);

% make a bunch of data
for n=1:N
    x{n} = round(imfilter(rand(siz),fspecial('gaussian',50,7),'same','symmetric'));
    % extremely difficult noise pattern -- from perturbation paper
    t = rand(size(x{n}));
    noiselevel = 1.25; % in pert paper 1.25
    y{n} = x{n}.*(1-t.^noiselevel) + (1-x{n}).*t.^noiselevel; 
end

% make features and labels
for n=1:N
    feats{n}  = [y{n}(:) 1+0*x{n}(:)];
    labels{n} = x{n}+1;
end

% no edge features here
efeats = []; % none

    % visualization function
    function viz(b_i)
        % here, b_i is a cell array of size nvals X nvars
        for n=1:N
            subplot(3,N,n    ); imshow(reshape(b_i{n}(2,:),siz,siz));
            subplot(3,N,n+  N); imshow(reshape(feats{n}(:,1),siz,siz));
            subplot(3,N,n+2*N); imshow(reshape(labels{n}-1,siz,siz));
            
        end
        xlabel('top: marginals  middle: input  bottom: labels')
        drawnow
    end

loss_spec = 'trunc_cl_trw_5';

crf_type  = 'linear_linear';
options.derivative_check = 'off';
options.viz         = @viz;
options.rho         = rho;
options.print_times = 1;
options.nvals       = nvals;

p = train_crf(feats,efeats,labels,model,loss_spec,crf_type,options);

% make a test image
x = round(imfilter(rand(siz),fspecial('gaussian',50,7),'same','symmetric'));
% extremely difficult noise pattern -- from perturbation paper
t = rand(size(x));
y = x.*(1-t.^noiselevel) + (1-x).*t.^noiselevel; 
feats  = [y(:) 1+0*x(:)];
labels = x+1;

[b_i b_ij] = eval_crf(p,feats,efeats,model,loss_spec,crf_type,rho);

b_i = reshape(b_i',[siz siz nvals]);

[~,label_pred] = max(b_i,[],3);
error = mean(label_pred(:)~=labels(:))

end