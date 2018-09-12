function [Z, msgs, B] = bpPerfectMatching(G, varargin)
tol = 10^(-9);
maxIters = 100;

[n, m] = size(G);
p = inputParser;
%G(i,j) should be the weight on edge (i,j) in the graph defined by the
%non-zero entries of G

%'reweight' is an optional symmetric matrix specifying the reweighting parameters
defaultC = ones(n,1);
addOptional(p,'reweight',defaultC,@isnumeric);

%'msgs' is an optional initialization of the msgs
defmitoj = rand(n, n, 2);
for i = 1:n
    for j = 1:n
        if G(i,j) == 0
            defmitoj(i,j,1) = 0;
            defmitoj(i,j,2) = 0;
        end
    end
end
addOptional(p,'msgs',defmitoj,@isnumeric);

%control damping parameter
addOptional(p,'damp',.95,@isnumeric);

parse(p,varargin{:});

C = p.Results.reweight;
mitoj = p.Results.msgs;
damp = p.Results.damp;

conv = 0;

%used to test for convergence
temp = mitoj(1:n, 1:n, 1:2);

%iterate until messages converge
iter = 1;

while conv == 0
    
    for i = 1:n
        
        %compute the msgs from i to j
        for j = find(G(i,:))
            if G(i,j) ~= 0
                %compute the product of all incoming msgs to i evaluated at
                %0 except j
                intoi0 = 1;
                for k =  find(G(i,:))
                    if G(i,k) ~= 0  && k~=j
                        intoi0 = intoi0 *mitoj(k, i, 1)^C(k)*mitoj(i, k, 1)^(C(i)-1);
                    end
                end
                
                %compute the sum of the product of all incoming msgs to i with
                %exactly one edge equal to one except j
                intoi1 = 0;
                
                for k =  1:n
                    if G(i,k) ~= 0  && k ~= j
                        tmp = 1;
                        for l = 1:n
                            if G(i,l) ~= 0 && l ~= k  && l ~= j
                                tmp = tmp *mitoj(l, i, 1)^C(l)*mitoj(i, l, 1)^(C(i) - 1);
                            end
                        end
                        intoi1 = intoi1 + tmp*G(i,k)*mitoj(k, i, 2)^C(l)*mitoj(i, k, 1)^(C(i) - 1);
                    end
                end
                
                zero = intoi1;
                one = intoi0;
                
                %normalize if the messages are not too small
                nm = zero + one;
                if nm > 10^(-8)
                    zero = zero/nm;
                    one = one/nm;
                else
                    if(zero > one)
                        zero = 1; one = 0;
                    else
                        one = 1;zero = 0;
                    end
                end
                %
                mitoj(i, j, 1) = (1-damp)*mitoj(i, j, 1) + (damp)*zero;
                mitoj(i, j, 2) = (1-damp)*mitoj(i, j, 2) + (damp)*one;
                
            end
        end
    end
    
    %check for convergence
    if max(max(max(abs(temp - mitoj(1:n, 1:n, 1:2))))) < tol || iter > maxIters
        conv = 1;
    else
        temp = mitoj(1:n, 1:n, 1:2);
    end
    iter = iter +1;
end

%fprintf('converged %d\n',iter)
%calculate beliefs at convergence
msgs = mitoj(1:n, 1:n, 1:2);
Z = 0;
b = zeros(n,n,2);

%edge beliefs
for i = 1:n
    for j = 1:n
        if G(i,j) ~= 0
            b(i,j,1) = mitoj(i,j,1)^C(i)*mitoj(j,i,1)^C(j);
            b(i,j,2) = mitoj(i,j,2)^C(i)*mitoj(j,i,2)^C(j)*G(i,j);
            
            %            assert(nm > 10^(-10),num2str(nm));
            nm = b(i,j,1) + b(i,j,2);
            
            if(nm > 10^(-10) && ~(isinf(b(i,j,1)) || isinf(b(i,j,2))))
                b(i,j,1) = b(i,j,1)/nm;
                b(i,j,2) = b(i,j,2)/nm;
            else
                if( b(i,j,1) > b(i,j,2))
                    b(i,j,1) = 1;  b(i,j,2) = 0;
                else
                    b(i,j,2) = 1;  b(i,j,1) = 0;
                end
            end
            
            
%             b(j,i,1) = b(i,j,1);
%             b(j,i,2) = b(i,j,2);
%             
            %disp([b(i,j,1), b(i,j,2)]);
            assert(~isnan(xlogx(b(i,j,1))));
            assert(~isnan(xlogx(b(i,j,2))));

            Z = Z + G(i,j)*b(i,j,2) + (C(i) + C(j) - 1)*xlogx(b(i,j,1)) - xlogx(b(i,j,2));
        end
    end
end
% disp(Z);
% end

B = b(1:n,1:n,2);
%Z = exp(Z);

