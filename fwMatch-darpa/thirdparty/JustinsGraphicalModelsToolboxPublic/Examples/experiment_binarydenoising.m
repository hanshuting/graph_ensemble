imsdir      = 'data/JustinBinary/BerkeleyBinaryDenoisingTrain5/';
labelsdir   = 'data/JustinBinary/BerkeleyBinaryDenoisingTrainLabels/';
imsdirT     = 'data/JustinBinary/BerkeleyBinaryDenoisingTest5/';
labelsdirT  = 'data/JustinBinary/BerkeleyBinaryDenoisingTestLabels/';
nvals       = 2;
maxims      = 32;

%feat_params = {{'patches',0},{'hog',8},{'position'},{'const'}}; % {'labpatches',0},
feat_params = {{'patches',0},{'const'}};

edge_params = {{'const'},{'pairtypes'}};

% This is from Justin
% loss_spec ={'trunc_ul_trw_0',...
%             'em_trw_1e-4',...
%             'bbp_ul_trw_1e-4',...
%             'bbp_cl_trw_1e-4',...
%             'bbp_uquad_trw_1e-4',...
%             'pseudo',...
%             'piecewise',...
%             'bbp_acc5_trw_1e-4',...
%             'bbp_acc15_trw_1e-4',...
%             'bbp_acc50_trw_1e-4'};

% trunch_ul_trw_0 *is* surrogate likelihood!
loss_spec ={'trunc_ul_trw_0'};

errs = {'1.25','1.5','5'};

testmode = 0;
for j=1:length(errs);
    err = errs{j};
    imsdir      = sprintf('data/JustinBinary/BerkeleyBinaryDenoisingTrain%s/',err);
    labelsdir   = 'data/JustinBinary/BerkeleyBinaryDenoisingTrainLabels/';
    imsdirT     = sprintf('data/JustinBinary/BerkeleyBinaryDenoisingTest%s/',err);
    labelsdirT  = 'data/JustinBinary/BerkeleyBinaryDenoisingTestLabels/';
    
    params_univ = train_labelers(...
          'trunc_ul_trw_0',[],testmode,...
          imsdir,labelsdir,nvals,feat_params,edge_params,maxims);
       
    for i=1:length(loss_spec)
        %params0  = [];
        % initialize smoothed loss to surr. likelihood solution.
        if i <= 7
            params0 = params_univ;
        else
            params0 = params{j,2};
        end
        
        [params{j,i} times{j,i} fvals{j,i}] = train_labelers(...
           loss_spec{i},params0,testmode,...
           imsdir,labelsdir,nvals,feat_params,edge_params,maxims);
    end
end

save experiment_binarydenoising

for j=1:length(errs);
    err = errs{j};
    imsdir      = sprintf('data/JustinBinary/BerkeleyBinaryDenoisingTrain%s/',err);
    labelsdir   = 'data/JustinBinary/BerkeleyBinaryDenoisingTrainLabels/';
    imsdirT     = sprintf('data/JustinBinary/BerkeleyBinaryDenoisingTest%s/',err);
    labelsdirT  = 'data/JustinBinary/BerkeleyBinaryDenoisingTestLabels/';
    
    for i=1:length(loss_spec)
        params0  = params{j,i};
        testmode = ['train_' loss_spec{i} '_' err];
        fprintf('training error for %s ...\n', loss_spec{i});
        colormap gray
        train_labelers(loss_spec{i},params0,testmode,...
            imsdir,labelsdir,nvals,feat_params,edge_params,maxims);
        testmode = ['test_' loss_spec{i} '_' err];
        fprintf('test    error for %s ...\n', loss_spec{i});
        colormap gray
        train_labelers(loss_spec{i},params0,testmode,...
            imsdirT,labelsdirT,nvals,feat_params,edge_params);
    end
   
    
    % change filenames for all images
    % this presumably won't work on windows
    for i=1:32
        system(sprintf('mv rez/im_%d_train.png rez/im_%d_train_%s.png',i,i,errs{j}));
    end
    for i=1:100
        system(sprintf('mv rez/im_%d_test.png rez/im_%d_test_%s.png',i,i,errs{j}));
    end
end
