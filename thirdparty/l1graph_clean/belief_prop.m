%assumes binary variables; and no node potentials
function [mun,mue] = belief_prop(graph)

n = size(graph,1);
k = 2;
theta_n = zeros(n*k,1);
theta_e = zeros(n*k);
for s = 1:n
	for t = 1:n
		theta_e((s-1)*k + 1,(t-1)*k + 1) = graph(s,t);
		theta_e((s-1)*k + 2,(t-1)*k + 2) = graph(s,t);
		
		theta_e((s-1)*k + 1,(t-1)*k + 2) = -1 * graph(s,t);
		theta_e((s-1)*k + 2,(t-1)*k + 1) = -1 * graph(s,t);
	end
end
theta_n;
theta_e;
bb = (graph~=0) * 1;
[mun,mueret] = sumprod_mex(theta_n,theta_e,bb,n,k,10);

mue = zeros(n,n,k,k);
for s = 1:n
	for t = 1:n
		for i = 1:k
			for j = 1:k
				mue(s,t,i,j) = mueret((s-1)*k + i,(t-1)*k + j);
			end
		end
	end
end

%size(mun)
%'is mue symm'
%sum(sum(mueret ~= mueret'))	
%pause
%'checking mun'
%for s= 1:n
%	sum(mun(s,:)) == 1
%end
%pause
em1 = zeros(2);
em2 = zeros(2);
for s = 1:n
	Nbs = find(graph(s,:));
	for idt = 1:length(Nbs)
		t = Nbs(idt);
		em1(:,:) = mue(s,t,:,:);
		em2(:,:) = mue(t,s,:,:);
		if(sum(sum(em1 ~= em2')) > 0)
			em1
			em2
			%pause
		end
	end
end