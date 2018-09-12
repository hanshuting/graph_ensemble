
if (3>2)

export LD_LIBRARY_PATH=/home/jebara/slick/src/mosek64/mosek/5/tools/platform/linux64x86/bin/:$LD_LIBRARY_PATH

export MOSEKLM_LICENSE_FILE=/home/jebara/slick/src/mosek64/mosek/5/licenses/mosek.lic

end


addpath  /home/jebara/slick/src/mosek64/mosek/5/toolbox/r2007a
mosekopt
help quadprog
H=eye(3,3);f=[1;2;3];A=eye(3,3);b=[4;7;8];quadprog(H,f,A,b)
load data.mat
Z=zeros(10,10);
Zor = zeros(10,10);
Zand=zeros(10,10);
for i=1:10
for j=1:10
C=exp(i);
alpha=exp(j-1)-1;
lcnst=0.000001*exp((i-1)*2+(j-1)/5);
tolrcnst=exp(-j);
[probsucc,fracdisagree,auctable,counts,graphret,truegraph] = evaluateAlgorithms(graphs,Xs,60,lcnst,tolrcnst,C,alpha);
Z(i,j)=auctable(3); 
Zor(i,j)=auctable(1); 
Zand(i,j)=auctable(2);
end
end





Z=zeros(10,10);
Zor = zeros(10,10);
Zand=zeros(10,10);
for i=1:10
for j=1:10
C=exp(i);
alpha=0.25*exp(j-1)-1;
lcnst=0.000001*exp((i-1)*2+(j-1)/5);
tolrcnst=exp(-j);
[probsucc,fracdisagree,auctable,counts,graphret,truegraph] = evaluateAlgorithms(graphs,Xs,60,lcnst,tolrcnst,C,alpha);
Z(i,j)=auctable(3); 
Zor(i,j)=auctable(1); 
Zand(i,j)=auctable(2);
end
end

CC=zeros(10,10);
AA=zeros(10,10);
LL=zeros(10,10);
for i=1:10
for j=1:10
C=exp(i);
alpha=0.25*exp(j-1)-1;
lcnst=0.000001*exp((i-1)*2+(j-1)/5);
CC(i,j)=C;
AA(i,j)=alpha;
LL(i,j)=lcnst;
end
end
