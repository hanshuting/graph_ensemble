function [bethe,iters,mytime,margs] = estper(W,path_to_binary,tol)

num = round(10000000*rand);
if(nargin == 2)
    tol = 1e-6;
end

fname = ['tmp_bpweights' num2str(num)];

eval(['save ' fname ' -ascii W']);

tries = 1;

while(~exist([fname 'out.txt'], 'file') && tries < 2)
    if (tries>1)
        disp(['Something went wrong, retrying. Try number ' num2str(tries)]);
        drawnow;
    end
    
    command = [path_to_binary ' ' fname ' ' num2str(size(W,1)) ' ' fname 'out.txt ' num2str(tol) ' ' fname 'out.txt.marginals'];
    %if (size(W,1)>50)
    %    command = ['nice -n 19 ' command];
    %    disp('Nicing process since N>50');
    %end
    tic; 
    system(command);
    mytime = toc;

    tries = tries+1;
end

tmp=load([fname 'out.txt']);
bethe=tmp(1);
iters = tmp(2);
margs = load([fname 'out.txt.marginals']);
system(['rm ' fname '*']);

