clear all; clc ;
addpath ../../Noor/GPS_D17pt_D246/bin_util/
load ../figure5/samples_allpriors_again_paper_mag0037_AS004.mat
samples = change_okada(unburnedsamples);
mommagsam = momentmag(samples);
cents = centroids(samples);

% authors in the order of the table:
% 1. GCMT
% 2. NEIC
% 3. NIED F net
% 4. Ito et al.
% 5. Nishimura et al.
% 6. Asano and Iwata
% 7. Takenaka et al.
% 8. Horikawa et al.
% 9. Kobayashi et al.
% 10. Sekiguchi et al.
% 11. Ozawa et al. 
% Locations
latlon = [130.15 130.131 130.17 130.18 130.297 130.1918 nan 130.2011 130.2486 130.2482 130.2016; ...
33.72 33.807 33.74 33.74 33.683 33.7363 nan 33.721 33.7074 33.6985 33.7264];

utmlonlat = FO_CrdTrans(latlon,1);
close
% strikes
strikes = [122 125 122 123 118 122 124 122 123 126 122];
% moment magnitudes
mommags = [6.612 6.639 6.565 6.586 6.596 6.677 6.565 6.474 6.637 6.68 nan];
% dip
dips = [88 117 93 94 101 93 87 91 92.3 93 90];
%color specs
colsplot = [1 0 0; 0 1 1; 1 .5 0;0 .5 1; 1 1 0; 0 1 0;0 0 1; .5 0 1; 0 .4 .3; 1 0  1; .5 .5 .5 ];
depths = [12 10 11 10 7.4 8.14 nan 12 6 4 5]; 

subplot(2,6,[1 2]);
%addpath ~/Desktop/softwares/gkde2/
%addpath ../../Noor/GPS_D17pt_D246/bin_util/
%load slipcolor
%p = gkde2([samples(1:100:end,6) samples(1:100:end,7)]);
%contour(p.x,p.y,p.pdf,1000)
%colormap(slipcolor)
hold on 
lw = 30; 

fp = confellipse2([samples(1:100:end,6) samples(1:100:end,7)],.68);
set(fp,'color','k','Linewidth',4)
fp = confellipse2([samples(1:100:end,6) samples(1:100:end,7)],.95);
set(fp,'color','k','Linewidth',2)

plot(utmlonlat(1,1),utmlonlat(2,1),'.','MarkerSize',lw,'color',colsplot(1,:))
plot(utmlonlat(1,2),utmlonlat(2,2),'.','MarkerSize',lw,'color',colsplot(2,:))
plot(utmlonlat(1,3),utmlonlat(2,3),'.','MarkerSize',lw,'color',colsplot(3,:))
plot(utmlonlat(1,4),utmlonlat(2,4),'.','MarkerSize',lw,'color',colsplot(4,:))
plot(utmlonlat(1,5),utmlonlat(2,5),'.','MarkerSize',lw,'color',colsplot(5,:))
plot(utmlonlat(1,6),utmlonlat(2,6),'.','MarkerSize',lw,'color',colsplot(6,:))
plot(utmlonlat(1,7),utmlonlat(2,8),'.','MarkerSize',lw,'color',colsplot(7,:))
plot(utmlonlat(1,8),utmlonlat(2,8),'.','MarkerSize',60,'color',colsplot(8,:))
plot(utmlonlat(1,9),utmlonlat(2,9),'.','MarkerSize',lw,'color',colsplot(9,:))
plot(utmlonlat(1,10),utmlonlat(2,8),'.','MarkerSize',lw,'color',colsplot(10,:))
plot(utmlonlat(1,11),utmlonlat(2,8),'.','MarkerSize',lw,'color',colsplot(11,:))

grid on
axis equal; 
%axis([605 612 3731 3735]); 
box('on')
% plot lat lon
% subplot(2,6,[1 2]);
% [nx,ny]= hist(samples(:,6),100);
% plot(ny,nx,'k','Linewidth',4)
% set(gca,'YTickLabel',[]);
% ca = [min([utmlonlat(1,[1:3 5:end]) ny]) max([utmlonlat(1,[1:3 5:end]) ny]) min(nx) max(nx)];
% hold on
% plot([utmlonlat(1,1) utmlonlat(1,1)],[ca(3) ca(4)],'Color',colsplot(1,:),'Linewidth',lw)
% plot([utmlonlat(1,2) utmlonlat(1,2)],[ca(3) ca(4)],'Color',colsplot(2,:),'Linewidth',lw)
% plot([utmlonlat(1,3) utmlonlat(1,3)],[ca(3) ca(4)],'Color',colsplot(3,:),'Linewidth',lw)
% plot([utmlonlat(1,5) utmlonlat(1,5)],[ca(3) ca(4)/3],'Color',colsplot(5,:),'Linewidth',lw)
% plot([utmlonlat(1,6) utmlonlat(1,6)],[ca(4)/3 2*ca(4)/3],'Color',colsplot(6,:),'Linewidth',lw)
% plot([utmlonlat(1,8) utmlonlat(1,8)],[2*ca(4)/3 ca(4)],'Color',colsplot(8,:),'Linewidth',lw)
% plot([utmlonlat(1,9) utmlonlat(1,9)],[ca(3) ca(4)],'Color',colsplot(9,:),'Linewidth',lw)
% xlabel('X center [km]')
% % 606.5373  608.3650  609.2913  620.2079  607.4263  607.4388       NaN  607.4263  608.5503
% axis(ca)
% subplot(2,6,[3 4])
% [nx,ny]= hist(samples(:,7),100);
% plot(ny,nx,'k','Linewidth',4)
% set(gca,'YTickLabel',[]);
% ca = [min([utmlonlat(2,[1:3 5:end]) ny]) max([utmlonlat(2,[1:3 5:end]) ny]) min(nx) max(nx)];
% hold on
% plot([utmlonlat(2,1)+.1 utmlonlat(2,1)+.1],[ca(3) ca(4)],'Color',colsplot(1,:),'Linewidth',lw)
% plot([utmlonlat(2,2) utmlonlat(2,2)],[ca(3) ca(4)/4],'Color',colsplot(2,:),'Linewidth',lw)
% plot([utmlonlat(2,3) utmlonlat(2,3)],[ca(4)/4 2*ca(4)/4],'Color',colsplot(3,:),'Linewidth',lw)
% plot([utmlonlat(2,4) utmlonlat(2,4)],[ca(3) ca(4)],'Color',colsplot(4,:),'Linewidth',lw)
% plot([utmlonlat(2,5) utmlonlat(2,5)],[ca(3) ca(4)/2],'Color',colsplot(5,:),'Linewidth',lw)
% plot([utmlonlat(2,6) utmlonlat(2,6)],[2*ca(4)/4 3*ca(4)/4],'Color',colsplot(6,:),'Linewidth',lw)
% plot([utmlonlat(2,8) utmlonlat(2,8)],[ca(3)/2 ca(4)],'Color',colsplot(8,:),'Linewidth',lw)
% plot([utmlonlat(2,9) utmlonlat(2,9)],[3*ca(4)/4 ca(4)],'Color',colsplot(9,:),'Linewidth',lw)
% xlabel('Y center [km]')
% % 3.7313    3.7336    3.7336    3.7274    3.7347    3.7336       NaN    3.7347    3.7336
% axis(ca)
lw = 3; 
subplot(2,6,[3 4])
[nx,ny]= hist(-cents(:,3),100);
plot(ny,nx,'k','Linewidth',4)
set(gca,'YTickLabel',[]);
ca = [min([-cents(:,3)' ny]) max([-cents(:,3)' ny]) min(nx) max(nx)];
hold on
plot([depths(1) depths(1)],[ca(3) ca(4)/2],'Color',colsplot(1,:),'Linewidth',lw)
plot([depths(2) depths(2)],[ca(3) ca(4)/2],'Color',colsplot(2,:),'Linewidth',lw)
plot([depths(3) depths(3)],[ca(3) ca(4)],'Color',colsplot(3,:),'Linewidth',lw)
plot([depths(4) depths(4)],[ca(4)/2 ca(4)],'Color',colsplot(4,:),'Linewidth',lw)
plot([depths(5) depths(5)],[ca(3) ca(4)],'Color',colsplot(5,:),'Linewidth',lw)
plot([depths(6) depths(6)],[ca(3) ca(4)],'Color',colsplot(6,:),'Linewidth',lw)
plot([depths(7) depths(7)],[ca(3) ca(4)],'Color',colsplot(7,:),'Linewidth',lw)
plot([depths(8) depths(8)],[ca(4)/2 ca(4)],'Color',colsplot(8,:),'Linewidth',lw)
plot([depths(9) depths(9)],[ca(3) ca(4)],'Color',colsplot(9,:),'Linewidth',lw)
plot([depths(10) depths(10)],[ca(3) ca(4)],'Color',colsplot(10,:),'Linewidth',lw)
plot([depths(11) depths(11)],[ca(3) ca(4)],'Color',colsplot(11,:),'Linewidth',lw)
xlabel('Centroid depth')
% 6.6120    6.5650    6.5860    6.5960    6.6770    6.5650    6.4740    6.6370    6.6800
axis([3.5 12.5 ca(3) ca(4)])


lw = 3; 
subplot(2,6,[5 6])
[nx,ny]= hist(mommagsam,100);
plot(ny,nx,'k','Linewidth',4)
set(gca,'YTickLabel',[]);
ca = [min([mommags ny]) max([mommags ny]) min(nx) max(nx)];
hold on
plot([mommags(1) mommags(1)],[ca(3) ca(4)],'Color',colsplot(1,:),'Linewidth',lw)
plot([mommags(2) mommags(2)],[ca(3) ca(4)/2],'Color',colsplot(2,:),'Linewidth',lw)
plot([mommags(3) mommags(3)],[ca(3) ca(4)/2],'Color',colsplot(3,:),'Linewidth',lw)
plot([mommags(4) mommags(4)],[ca(3) ca(4)],'Color',colsplot(4,:),'Linewidth',lw)
plot([mommags(5) mommags(5)],[ca(3) ca(4)],'Color',colsplot(5,:),'Linewidth',lw)
plot([mommags(6) mommags(6)],[ca(3) ca(4)/2],'Color',colsplot(6,:),'Linewidth',lw)
plot([mommags(7) mommags(7)],[ca(4)/2 ca(4)],'Color',colsplot(7,:),'Linewidth',lw)
plot([mommags(8) mommags(8)],[ca(3) ca(4)],'Color',colsplot(8,:),'Linewidth',lw)
plot([mommags(9) mommags(9)],[ca(4)/2 ca(4)],'Color',colsplot(9,:),'Linewidth',lw)
plot([mommags(10) mommags(10)],[ca(4)/2 ca(4)],'Color',colsplot(10,:),'Linewidth',lw)
plot([mommags(11) mommags(11)],[ca(3) ca(4)],'Color',colsplot(11,:),'Linewidth',lw)
xlabel('moment magnitude')
% 6.6120    6.5650    6.5860    6.5960    6.6770    6.5650    6.4740    6.6370    6.6800
axis(ca)
subplot(2,6,[7 8])
[nx,ny]= hist(samples(:,5),100);
plot(ny,nx,'k','Linewidth',4)
set(gca,'YTickLabel',[]);
ca = [min([strikes ny]) max([strikes ny]) min(nx) max(nx)];
hold on
plot([strikes(1) strikes(1)],[ca(3) ca(4)/5],'Color',colsplot(1,:),'Linewidth',lw)
plot([strikes(2) strikes(2)],[ca(3) ca(4)],'Color',colsplot(2,:),'Linewidth',lw)
plot([strikes(3) strikes(3)],[ca(4)/5 2*ca(4)/5],'Color',colsplot(3,:),'Linewidth',lw)
plot([strikes(4) strikes(4)],[ca(3) ca(4)/2],'Color',colsplot(4,:),'Linewidth',lw)
plot([strikes(5) strikes(5)],[ca(3) ca(4)],'Color',colsplot(5,:),'Linewidth',lw)
plot([strikes(6) strikes(6)],[2*ca(4)/5 3*ca(4)/5],'Color',colsplot(6,:),'Linewidth',lw)
plot([strikes(7) strikes(7)],[ca(3) ca(4)],'Color',colsplot(7,:),'Linewidth',lw)
plot([strikes(8) strikes(8)],[3*ca(4)/5 4*ca(4)/5],'Color',colsplot(8,:),'Linewidth',lw)
plot([strikes(9) strikes(9)],[ca(4)/2 ca(4)],'Color',colsplot(9,:),'Linewidth',lw)
plot([strikes(10) strikes(10)],[ca(3) ca(4)],'Color',colsplot(10,:),'Linewidth',lw)
plot([strikes(11) strikes(11)],[4*ca(4)/5 ca(4)],'Color',colsplot(11,:),'Linewidth',lw)

xlabel('strike [\circ]')
% 122 122 123 118 122 124 122 123 126
axis(ca)
subplot(2,6,[9 10])
[nx,ny]= hist(-samples(:,4),100);
plot(ny,nx,'k','Linewidth',4)
set(gca,'YTickLabel',[]);
ca = [min([dips ny]) max([dips ny]) min(nx) max(nx)];
hold on
plot([dips(1) dips(1)],[ca(3) ca(4)],'Color',colsplot(1,:),'Linewidth',lw)
plot([dips(2) dips(2)],[ca(3) ca(4)],'Color',colsplot(2,:),'Linewidth',lw)
plot([dips(3) dips(3)],[ca(3) ca(4)/3],'Color',colsplot(3,:),'Linewidth',lw)
plot([dips(4) dips(4)],[ca(3) ca(4)],'Color',colsplot(4,:),'Linewidth',lw)
plot([dips(5) dips(5)],[ca(3) ca(4)],'Color',colsplot(5,:),'Linewidth',lw)
plot([dips(6) dips(6)],[ca(4)/3 2*ca(4)/3],'Color',colsplot(6,:),'Linewidth',lw)
plot([dips(7) dips(7)],[ca(3) ca(4)],'Color',colsplot(7,:),'Linewidth',lw)
plot([dips(8) dips(8)],[ca(3) ca(4)],'Color',colsplot(8,:),'Linewidth',lw)
plot([dips(9) dips(9)],[ca(3) ca(4)],'Color',colsplot(9,:),'Linewidth',lw)
plot([dips(10) dips(10)],[2*ca(4)/3 ca(4)],'Color',colsplot(10,:),'Linewidth',lw)
plot([dips(11) dips(11)],[ca(3) ca(4)],'Color',colsplot(11,:),'Linewidth',lw)

xlabel('dip [\circ]')
% 88 93 94 79 93 93 91 92.3 93
axis(ca)
legend('our solution','GCMT',...
'NEIC',...
'NIED F-net',...
'Ito et al. (2006)',...
'Nishimura et al. (2006)',...
'Asano and Iwata (2006)',...
'Takenaka et al. (2006)',...
'Horikawa et al. (2006)',...
'Kobayashi et al. (2006)',...
'Sekiguchi et al. (2006)',...
'Ozawa et al. (2006)')