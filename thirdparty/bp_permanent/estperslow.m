function [per,G] = estperslow(W);

THRESHOLD = 1e-20;

N = size(W,1);

G = makebmatchexpanded(W,1);
%G = makenewmatch(W);

go = 1;

iters = 1;

oldb = G.b;

while(go)

    G=bp(G);

    %P = zeros(N);
    P = G.b;

    change = sum(abs(G.b(:)-oldb(:)));

    if (change<THRESHOLD)
        go = 0;
    end

    oldb = G.b;

    if (iters>500*N)
        go=0;
        disp('reached maximum iterations');
    end
    
    %fb(iters) = bethe(G);
    %plot(fb); drawnow;
    
    iters = iters+1;
end


per = exp(-bethe(G));