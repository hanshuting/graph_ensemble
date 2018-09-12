function C = makeclusters(imsdir,nclusters,method,maxims)

% go through imsdir, run filters, collect data, do cultering

Aims = [dir([imsdir    '/*.png']);dir([imsdir    '/*.jpg']);dir([imsdir    '/*.bmp'])];
Aims = Aims(1:5:end); 
'changing # ims!' 
nims = length(Aims);
if nargin==4
    if maxims > nims
        error('maxims is greater than # ims in director');
    end
    nims = maxims;
end

rescale   = 1;
ndata_per_im = 100;

% get the data (100 per im)
if strcmp(method,'MR8')
    data = zeros(8,nims*ndata_per_im);
elseif strcmp(method,'canny')
    data = zeros(11*11*1,nims*ndata_per_im);
else
    error('undefined cluster feature type');
end
where = 1;
for m=1:nims
    im    = double(imread([imsdir    '/' Aims(m).name]))/255;
    if ndims(im)==3
        im    = rgb2gray(im);
    end
    im = imresize(im,rescale,'bilinear');
    [ly lx] = size(im);
    
    if strcmp(method,'MR8')
        f = MR8fast(im);
        a = ceil(rand(1,ndata_per_im)*size(f,2));
        data(:,where:where+ndata_per_im-1) = f(:,a);
        where = where+1000;
    elseif strcmp(method,'canny')
        %E = edge(im,'canny',.1)/3 + edge(im,'canny',.2)/3 + edge(im,'canny',.5)/3;
        %E   = double(edge(im,'canny'));
        
        %im2 = imresize(im,.5,'bilinear');
        %E2  = double(edge(im2,'canny'));
        %E2  = imresize(E2,size(E),'bilinear');
        
        im3 = imresize(im,.25,'bilinear');
        E3  = double(edge(im3,'canny',.1));
        E3  = imresize(E3,[ly lx],'nearest');
        
        %E  = reflectim(E ,50);
        %E2 = reflectim(E2,50);
        E3 = reflectim(E3,50);
        
        for reps=1:ndata_per_im
            y = 50+ceil(rand*(ly));
            x = 50+ceil(rand*(lx));
            %e  = E (y-5  :1:y+5  ,x-5  :1:x+5  );
            %e2 = E2(y-5*2:2:y+5*2,x-5*2:2:x+5*2);
            e3 = E3(y-5*4:4:y+5*4,x-5*4:4:x+5*4);
            
            %data(:,where) = [e(:); e2(:); e3(:)];
            data(:,where) = e3(:);
            where = where+1;
        end
        warning off all
        %imshow([E E2 E3])
        imshow([E3])
        warning on all
        drawnow
    end
end

[IDX,C] = kmeans(data',nclusters);
C=C';

for i=1:nclusters
    who = find(IDX==i,1);
    C(:,i) = data(:,who);
end