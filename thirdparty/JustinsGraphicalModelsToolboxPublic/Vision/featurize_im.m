function [feats names] = featurize_im(im,feat_params)

% F = featurize_im(im,{{'patches',2},{'hog',32}});
% take a ly x lx image and return a ly x lx x lz image of features
%
% each feature takes some parameters:
% patches(k)       : the RBG  intensities in the surrounding 2k+1 x 2k+1 patch
% hsvpatches(k)    : the HSV  intensities in the surrounding 2k+1 x 2k+1 patch
% labpatches(k)    : the LAB  intensities in the surrounding 2k+1 x 2k+1 patch
% graybpatches(k)  : the gray intensities in the surrounding 2k+1 x 2k+1 patch
% const            : a constant of one
% position         : the y and x position
% hog(k)           : Histogram of Oriented Gradient features
% colhist(nbins,k) : compute a color histogram with resolution nbins in
%                  : each dimension about each pixel in a 2k1+ x 2k+1 patch
% lbp              : compute local binary patterns

% 'hog' features require piotr dollar's toolbox
% from http://vision.ucsd.edu/~pdollar/toolbox/doc/

[ly lx lz] = size(im);

%feats = [];

nfeat = 0;
for i=1:length(feat_params)
    fp = feat_params{i};
    if strcmp(fp{1},'patches')
        buf_size = fp{2};
        nfeat = nfeat + (1+2*buf_size)^2*lz;
    elseif strcmp(fp{1},'hsvpatches')
        buf_size = fp{2};
        nfeat = nfeat + (1+2*buf_size)^2*lz;
    elseif strcmp(fp{1},'labpatches')
        buf_size = fp{2};
        nfeat = nfeat + (1+2*buf_size)^2*lz;
        if ~exist('colorspace')
            error(['to use lab features you must install Pascal Getreuers colorspace.m from' ...
                   'http://www.mathworks.com/matlabcentral/fileexchange/28790-colorspace-transformations'])
        end
    elseif strcmp(fp{1},'fouriercolor')
        maxk  = fp{2};
        nfeat = nfeat + 2*(maxk+1)^lz;
    elseif strcmp(fp{1},'graypatches')
        buf_size = fp{2};
        nfeat = nfeat + (1+2*buf_size)^2;
    elseif strcmp(fp{1},'const')
        nfeat = nfeat + 1;
    elseif strcmp(fp{1},'position')
        %order = fp{2};
        %nfeat = nfeat + 2*order;
        nfeat = nfeat + 2;
    elseif strcmp(fp{1},'hog')
        if ~exist('hog')
            error(['to use hog features you must install piotr dollars toolbox' ...
                   'from http://vision.ucsd.edu/~pdollar/toolbox/doc/'])
        end
        nfeat = nfeat + 36;
    elseif strcmp(fp{1},'colhist')
        rez = fp{2};
        nfeat = nfeat + rez^lz;
    elseif strcmp(fp{1},'lbp')
        nfeat = nfeat + 16;
        if ~exist('lbp')
            error(['to use lbp features you must get code from' ...
                   'from http://www.cse.oulu.fi/MVG/Downloads/LBPMatlab'])            
        end
    elseif strcmp(fp{1},'texton')
        C = fp{2};
        nfeat = nfeat + size(C,2);
    elseif strcmp(fp{1},'daisy')
        nfeat = nfeat + 200;
    elseif strcmp(fp{1},'fourier')
        %if i~=length(feat_params)
        %    error('fourier must be last!');
        %end
        maxk = fp{2};
        nfeat = 2*(maxk+1)^nfeat;
    else
        error('unsupported feature type: %s', fp{1})
    end
end
feats = zeros(ly,lx,nfeat);
%fprintf('nfeat: %d\n',nfeat);

where = 1;
for i=1:length(feat_params)
    fp = feat_params{i};
    if strcmp(fp{1},'patches')
        buf_size = fp{2};
        im2    = reflectim(im   ,buf_size);
        for dy=-buf_size:buf_size
            for dx=-buf_size:buf_size
                for z=1:lz
                    stuff = im2((1:ly)+buf_size+dy,(1:lx)+buf_size+dx,z);
                    %feats = cat(3,feats,stuff);
                    feats(:,:,where) = stuff;
                    names{where} = 'patches';
                    where=where+1;
                end
            end
        end
    elseif strcmp(fp{1},'hsvpatches')
        buf_size = fp{2};
        im2    = reflectim(rgb2hsv(im),buf_size);
        for dy=-buf_size:buf_size
            for dx=-buf_size:buf_size
                for z=1:lz
                    stuff = im2((1:ly)+buf_size+dy,(1:lx)+buf_size+dx,z);
                    %feats = cat(3,feats,stuff);
                    feats(:,:,where) = stuff;
                    names{where} = 'hsvpatches';
                    where=where+1;
                end
            end
        end
    elseif strcmp(fp{1},'labpatches')
        buf_size = fp{2};
        im2    = reflectim(colorspace('rgb->lab',im),buf_size);
        for dy=-buf_size:buf_size
            for dx=-buf_size:buf_size
                for z=1:lz
                    stuff = im2((1:ly)+buf_size+dy,(1:lx)+buf_size+dx,z);
                    %feats = cat(3,feats,stuff);
                    feats(:,:,where) = stuff;
                    names{where} = 'labpatches';
                    where=where+1;
                end
            end
        end
    elseif strcmp(fp{1},'fouriercolor')
        maxk = fp{2};
        
        C = tf_table(lz,maxk+1)-1;
        
        %X2 = zeros(2*size(C,2),size(X,2));
        for j=1:size(C,2)
            c = C(:,j);
            % compute c'*X;
            rez = zeros(ly,lx);
            for z=1:lz
                rez = rez + c(z)*im(:,:,z);
            end
            feats(:,:,where) = cos(pi*rez);
            names{where} = 'fouriercolor';
            where = where+1;
            feats(:,:,where) = sin(pi*rez);
            names{where} = 'fouriercolor';
            where = where+1;
        end        
    elseif strcmp(fp{1},'graypatches')
        buf_size = fp{2};
        im2    = reflectim(rgb2gray(im),buf_size);
        for dy=-buf_size:buf_size
            for dx=-buf_size:buf_size
                for z=1
                    stuff = im2((1:ly)+buf_size+dy,(1:lx)+buf_size+dx,z);
                    %feats = cat(3,feats,stuff);
                    feats(:,:,where) = stuff;
                    names{where} = 'graypatches';
                    where=where+1;
                end
            end
        end
    elseif strcmp(fp{1},'const')
        %feats = cat(3,feats,ones(ly,lx));
        feats(:,:,where) = ones(ly,lx);
        names{where} = 'const';
        where=where+1;
    elseif strcmp(fp{1},'position')
        maxorder = fp{2};
        [xpos ypos] = meshgrid(1:lx,1:ly);
        xpos = xpos/lx;
        ypos = ypos/ly;
        for xorder=0:maxorder
            for yorder=0:maxorder
                if (xorder ~=0 || yorder~=0) && xorder+yorder <= maxorder
                    feats(:,:,where) = xpos.^xorder .* ypos.^yorder;
                    names{where} = 'position-x';
                    where=where+1;
                end
            end
        end
    elseif strcmp(fp{1},'hog')
        % rez - train   / test
        % 2   - 0.194550 / 0.2056700
        % 4   - 0.185340 / 0.193420
        % 8   - 0.176620 / 
        % 16  - 0.162500 / 0.1766800 - 0.168660 / 0.182520
        % 32  - 0.150070 / 0.1752900 - 0.157080 / 0.180830
        % 8+32- 0.142270 / 0.1637600
        
        rez = fp{2};        
        %H = hog(im,rez,9,10);
%         H = hog(single(reflectim(im,rez)),rez,9,10);
        H = hog(reflectim(im,rez),rez,9,10);
        
        lzH = size(H,3);
        H2 = imresize(H,[ly lx],'bilinear');
        feats(:,:,where:where+lzH-1) = H2;
        for where2=where:where+lzH-1
            names{where2} = 'hog';
        end
        where=where+lzH;

%         rez = fp{2};
%         H = myhog(im);
%         lzH = size(H,3);
%         feats(:,:,where:where+lzH-1) = H;
%         where=where+lzH;

%         r = fp{2};
%         H = hog(im,r,9,10);
%         H2 = zeros(size(H,1)+2,size(H,2)+2,size(H,3));
%         H2(2:end-1,2:end-1,:) = H;
%         
%         H2 = imresize(H2,[size(im,1),size(im,2)],'bilinear');        
%         lzH = size(H,3);
%         feats(:,:,where:where+lzH-1) = H2;
%         where=where+lzH;
        
        
        %figure(1), imshow(im)
        %figure(2), imshow(hogdraw(H,25))
    elseif strcmp(fp{1},'colhist')
        rez = fp{2};
        buf = fp{3};
        F = histify(im,rez,buf);
        feats(:,:,where:where+size(F,3)-1) = F;
        for where2=where:where+size(F,3)
            names{where2} = 'colhist';
        end
        where=where+size(F,3);
    elseif strcmp(fp{1},'lbp')
        SP=[-1 0; 0 -1; -0 1; 1 0];
        I2 = lbp(rgb2gray(im),SP,0,'i')+1;
        % lbp is smaller by one pixel...
        I2 = [ones(1,lx-2); I2; ones(1,lx-2)];
        I2 = [ones(ly,1) I2 ones(ly,1)];
        
        F = zeros(ly,lx,16);
        for y=1:ly
            for x=1:lx
                F(y,x,I2(y,x))=1;
            end
        end
        feats(:,:,where:where+15) = F;
        for where2=where:where+15
            names{where2} = 'lbp';
        end
        where=where+16;
    elseif strcmp(fp{1},'texton')
        % these seem to be utterly useless...
        C = fp{2};
        T = MR8textonify(im,C);
        for z=1:size(C,2)
            F(:,:,where) = double(T==z);
            names{where} = 'texton';
            where=where+1;
        end
    elseif strcmp(fp{1},'daisy')
        D = compute_daisy(im);
        myF = reshape(D.descs,[lx ly 200]); myF = permute(myF,[2 1 3]);
        F(:,:,where:where+199) = myF;
        where = where+200;
    elseif strcmp(fp{1},'fourier')
        % throw away all the old features and replace them with fourier
        % features!
        maxk = fp{2};
        maxz = where-1;
        oldfeats = feats(:,:,1:maxz);
        C = tf_table(maxz,maxk+1)-1;
        
        nfeat = 2*(maxk+1)^maxz;
        feats = zeros(ly,lx,nfeat);
        
        where = 1;
        for j=1:size(C,2)
            c = C(:,j);
            % compute c'*X;
            rez = zeros(ly,lx);
            for z=1:size(oldfeats,3)
                rez = rez + c(z)*oldfeats(:,:,z);
            end
            feats(:,:,where) = cos(pi*rez);
            names{where} = 'fourier';
            where = where+1;
            feats(:,:,where) = sin(pi*rez);
            names{where} = 'fourier';
            where = where+1;
        end
        %fprintf('[%d %d]\n',nfeat,where);
    else
        error('unsupported feature type: %s',fp{1})
    end
end
