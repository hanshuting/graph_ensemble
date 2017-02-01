function [feats names] = edgeify_im(im,edge_params,pairs,pairtype)

% F = edgeify_im(im,{{'patches',2},{'hog',32}});
% take a ly x lx image and return a lz x npairs image of features
%
% pairtypes (which must be last) takes all the previously existing
% features, and multiplies them by an indicator function for each possible
% pairtype (must be provided as third input).  This is useful, for example,
% to have separate parameters for vertical and horizontal links.

nthresh = 10;

[ly lx lz] = size(im);

nfeat = 0;
for i=1:length(edge_params)
    fp = edge_params{i};
    if strcmp(fp{1},'const')
        nfeat = nfeat + 1;
    elseif strcmp(fp{1},'diffthresh')
        %nfeat = nfeat + length(fp{2});
        nfeat = nfeat + nthresh;
    elseif strcmp(fp{1},'difffourier')
        maxk = fp{2};
        nfeat = nfeat + (maxk+1)*2;
    elseif strcmp(fp{1},'sobelthresh')
        %nfeat = nfeat + length(fp{2});
        nfeat = nfeat + nthresh;
    elseif strcmp(fp{1},'sobelfourier')
        maxk = fp{2};
        nfeat = nfeat + (maxk+1)*2;
    elseif strcmp(fp{1},'edges')
        nfeat = nfeat + 5;
    elseif strcmp(fp{1},'canny')
        nfeat = nfeat + 10;
    elseif strcmp(fp{1},'pairtypes')
        if i~=length(edge_params)
            error('pairtypes must be last!');
        end
        if nargin < 4
            error('must provide pairtype as 4th input');
        end
        nfeat = nfeat * max(pairtype);
    else
        error('unsupported feature type: %s', fp{1})
    end
end
npairs = size(pairs,1);
feats = zeros(npairs,nfeat);

where = 1;
for i=1:length(edge_params)
    fp = edge_params{i};
    if strcmp(fp{1},'const')
        %feats = cat(3,feats,ones(ly,lx));
        feats(:,where) = ones(1,npairs);
        names{where} = 'const';
        where=where+1;
    elseif strcmp(fp{1},'diffthresh')
        diff = 0;
        for z=1:lz
            diff = diff + (im(pairs(:,1)+ly*lx*(z-1))-im(pairs(:,2)+ly*lx*(z-1))).^2;
        end
        diff = sqrt(diff);
        diff = diff(:);
        %diff2 = sort(diff);
        for n=1:nthresh
            %thresh = diff2(round(length(diff)*n/11));
            thresh = .5*n/nthresh; % similar to above strategy in practice
            feats(:,where) = double(diff>thresh);
            names{where} = ['diff-' num2str(n)];
            where = where+1;
        end
    elseif strcmp(fp{1},'difffourier')
        diff = 0;
        for z=1:lz
            diff = diff + (im(pairs(:,1)+ly*lx*(z-1))-im(pairs(:,2)+ly*lx*(z-1))).^2;
        end
        diff = sqrt(diff);
        diff = diff(:)';
        %diff = min(1,diff*3);
        %r=.001; diff = r./(r+diff);
        
        maxk = fp{2};
        for n=0:maxk
            feats(:,where) = cos(pi*n*diff);
            names{where} = ['difffourier-' num2str(n)];
            where = where+1;
            feats(:,where) = sin(pi*n*diff);
            names{where} = ['difffourier-' num2str(n)];
            where = where+1;
        end    
    elseif strcmp(fp{1},'sobelthresh')
%         for n=1:nthresh
%             thresh = n*.2/nthresh;
%             E = edge(rgb2gray(im),'sobel',thresh);
%             feats(:,where) = double(max( E(pairs(:,1)), E(pairs(:,2)))==0);
%             names{where} = ['sobel-' num2str(n)];
%             where = where+1;
%         end
        [E,thresh,gv,gh] = edge(rgb2gray(im),'sobel');
        g = sqrt(gv.^2+gh.^2);
        diff = max(g(pairs(:,1)), g(pairs(:,2)));
        for n=1:nthresh
            thresh = n*.2/nthresh;
            %feats(:,where) = double((diff>thresh) & (diff <= thresh_next));
            feats(:,where) = double((diff>thresh));
            names{where} = ['sobel-' num2str(n)];
            where = where+1;
        end
    elseif strcmp(fp{1},'sobelfourier')
        [E,thresh,gv,gh] = edge(rgb2gray(im),'sobel');
        g = sqrt(gv.^2+gh.^2);
        diff = max(g(pairs(:,1)), g(pairs(:,2)));
        
        maxk = fp{2};
        for n=0:maxk
            feats(:,where) = cos(pi*n*diff);
            names{where} = ['sobelfourier-' num2str(n)];
            where = where+1;
            feats(:,where) = sin(pi*n*diff);
            names{where} = ['sobelfourier-' num2str(n)];
            where = where+1;
        end
        
    elseif strcmp(fp{1},'canny')
        for thresh=[.5]
            for sig=[1 5 15]
            E = edge(rgb2gray(im),'canny',thresh,sig);
            feats(:,where) = double(max( E(pairs(:,1)), E(pairs(:,2))));
            names{where} = 'canny';
            where = where+1;
            end
        end
    elseif strcmp(fp{1},'edges')
        for n=1:5
            if n==1
                type = 'sobel';
            elseif n==2
                type = 'prewitt';
            elseif n==3
                type = 'roberts';
            elseif n==4
                type = 'log';
            else
                type = 'canny';
            end
            E = edge(rgb2gray(im),type);
            feats(:,where) = double(max( E(pairs(:,1)), E(pairs(:,2))));
            names{where} = 'edges';
            where = where+1;
        end
    elseif strcmp(fp{1},'pairtypes')
        oldfeats = feats;
        oldnames = names;
        oldwhere = where;
        where = 1;
        for ii=1:oldwhere-1
            for jj=1:max(pairtype)
                feats(:,where) = oldfeats(:,ii).*double(pairtype(:)==jj);
                names{where} = [oldnames{ii} ' type ' num2str(jj)];
                where = where+1;
            end
        end
        %[oldwhere-1 nfeat]
    else
        error('unsupported feature type: %s',fp{1})
    end
end