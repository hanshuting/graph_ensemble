%ASSUMES the samples X are binary with values {-1,1}
%OUTPUT: graphret_gen: chow liu with known degree
%	 graphret_hub: chow liu with known degree and known hub

function [graphret_gen,graphret_hub] = chowliu(X,d)

n = size(X,1);
p = size(X,2);

probn = zeros(p,2);
probn(:,1) = sum((X == -1),1)/n;
probn(:,2) = 1 - probn(:,1);

probe = zeros(p,p,2,2);

%CALCULATE EDGE PROBS
for s = 1:p
	for t = (s+1):p
		for l = 1:n
			for i = 0:1
				for j = 0:1
					if((X(l,s) == (2*i - 1)) && (X(l,t) == (2*j - 1)))
						probe(s,t,i+1,j+1) = probe(s,t,i+1,j+1) + 1; 
					end
				end
			end
		end
		probe(s,t,:,:) = probe(s,t,:,:)/sum(sum(probe(s,t,:,:)));
	
	end
end

%CALCULATE MI WEIGHTS
W = zeros(p);
for s = 1:p	
	for t = (s+1):p
		for i = 0:1
			for j = 0:1 
				prbst = probe(s,t,i+1,j+1);
				if(prbst == 0)
					continue;
				end
				W(s,t) = W(s,t) + prbst * log(prbst/(probn(s,i+1) * probn(t,j+1)));
			end
		end
	end
end

W_gen = W;
graphret_gen = zeros(p);
%KRUSKAL MWST
nodeid = 1:p;
nume = 0;
itr = 0;
while(1)
	itr = itr + 1;
	
	[er,ri] = max(W_gen);
	[maxw,rc] = max(er);
	if(maxw == 0)
		break;
	end
	s = ri(rc);
	t = rc;
	W_gen(s,t) = 0;
	W_gen(t,s) = 0;
	
	if(nodeid(s) == nodeid(t))	
		continue;
	else
		nodeid(nodeid == nodeid(t)) = nodeid(s);
		graphret_gen(s,t) = 1;
		graphret_gen(t,s) = 1;
		nume = nume + 1;
	end
	if(nume == d)
		break;
	end
end

%GREEDY WITH HUB KNOWN
graphret_hub = zeros(p);
if(sum(W(1,:) > 0) < d)
	graphret_hub(1,:) = (W(1,:) > 0) * 1;
	graphret_hub(:,1) = graphret_hub(1,:);
else
	[sval,sid] = sort(W(1,:),'descend');
	graphret_hub(1,sid(1:d)) = 1;
	graphret_hub(sid(1:d),1) = 1;
end
