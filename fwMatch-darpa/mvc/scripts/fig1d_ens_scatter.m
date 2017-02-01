
load('C:\Shuting\fwMatch\data\m21_d2_vis.mat')
load('C:\Shuting\fwMatch\data\ensembles\Core_m21_d2_vis.mat')

ens_coord1 = Pools_coords(:,1:2,1);
ens_coord1 = ens_coord1(sum(ens_coord1,2)~=0,:);
ens_coord2 = Pools_coords(:,1:2,2);
ens_coord2 = ens_coord2(sum(ens_coord2,2)~=0,:);

circ_sz = 100;

h = figure;set(gcf,'color','w');
hold on;
scatter(ens_coord1(:,1),-ens_coord1(:,2),circ_sz,'r','filled');
scatter(Coord_active(:,1),-Coord_active(:,2),circ_sz,'k')
axis off equal
saveas(h,'C:\Shuting\fwMatch\paper\figures\matlab_fig\1d_scatter_ens1');
saveas(h,'C:\Shuting\fwMatch\paper\figures\matlab_fig\1d_scatter_ens1.pdf');

h = figure;set(gcf,'color','w');
hold on;
scatter(ens_coord2(:,1),-ens_coord2(:,2),circ_sz,'r','filled');
scatter(Coord_active(:,1),-Coord_active(:,2),circ_sz,'k')
axis off equal
saveas(h,'C:\Shuting\fwMatch\paper\figures\matlab_fig\1d_scatter_ens2');
saveas(h,'C:\Shuting\fwMatch\paper\figures\matlab_fig\1d_scatter_ens2.pdf');