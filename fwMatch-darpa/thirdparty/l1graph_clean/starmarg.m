function [margn,marge] = starmarg(graph)

p = size(graph,1);
margn = zeros(p,2);
marge = zeros(p,p,2,2);


mesgin = zeros(p,2);
mesgout = zeros(p,2);

k = 2;

Nbh = find(graph(1,:));
for ids = 1:length(Nbh)
	s = Nbh(ids);
	for x_h = 0:1
		mesgin(s,x_h+1) = exp(graph(1,s) * (2* x_h - 1)) + exp(graph(1,s) * -1 * (2 * x_h - 1));
	end
end
for ids = 1:length(Nbh)
	s = Nbh(ids);
	for x_s = 0:1
		for x_h = 0:1
			tmps = exp(graph(1,s) * (2 * x_s - 1) * (2 * x_h - 1));
			for idu = 1:length(Nbh)
				u = Nbh(idu);
				if(u == s)
					continue;
				end
				tmps = tmps * mesgin(u,x_h+1);
			end
			mesgout(s,x_s+1) = mesgout(s,x_s+1) + tmps;
		end
	end
end

for x_h = 0:1
	tmph = 1;
	for idu = 1:length(Nbh)
		u = Nbh(idu);
		tmph = tmph * mesgin(u,x_h+1);
	end
	margn(1,x_h+1) = tmph;
end
margn(1,:) = margn(1,:)/sum(margn(1,:));

for ids = 1:length(Nbh)
	s = Nbh(ids);
	for x_s = 0:1
		margn(s,x_s + 1) = mesgout(s,x_s + 1);
		for x_h = 0:1
			tmpsh = exp(graph(1,s) * (2 * x_s - 1) * (2 * x_h - 1));
			for idu = 1:length(Nbh)
				u = Nbh(idu);
				if(u == s)
					continue;
				end
				tmpsh = tmpsh * mesgin(u,x_h+1);
			end
			marge(s,1,x_s+1,x_h+1) = tmpsh;
			marge(1,s,x_s+1,x_h+1) = tmpsh;
		end
	end
	msum = sum(sum(marge(s,1,:,:)));
	marge(s,1,:,:) = marge(s,1,:,:)/msum;
	marge(1,s,:,:) = marge(1,s,:,:)/msum;
	
	margn(s,:) = margn(s,:)/sum(margn(s,:));
end
	