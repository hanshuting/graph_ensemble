function [bethe,iters] = estper(W);

num = round(10000000*rand);

fname = ['tmp_bpweights' num2str(num)];

eval(['save ' fname ' -ascii W']);

tries = 1;

global estper_exe;

while(~exist([fname 'out.txt'], 'file'))
    if (tries>1)
        disp(['Something went wrong, retrying. Try number ' num2str(tries)]);
        drawnow;
    end
    
    command = [estper_exe ' ' fname ' ' num2str(size(W,1)) ' ' fname 'out.txt'];

    %if (size(W,1)>50)
    %    command = ['nice -n 19 ' command];
    %    disp('Nicing process since N>50');
    %end

    system(command);

    tries = tries+1;
end

tmp=load([fname 'out.txt']);
bethe=tmp(1);
iters = tmp(2);
system(['rm ' fname '*']);

