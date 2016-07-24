%% Data misfit with MAP 

addpath fig_extras
model = 1.0e+03 * [0.6048    3.7359    0.6139    3.7304    0.0174    0.0000   -0.0951   -0.0017         0         0]; 
model1 = change_okada(model);
magmodel = momentmag(model1); 
load 22837.dat
coastutm = FO_CrdTrans(X22837',1);
coastutm(2,:) = coastutm(2,:) -.3; 
coastutm(1,:) = coastutm(1,:) -.37; 

addpath ../../Noor/GPS_D17pt_D246/extras/
addpath ../../Noor/GPS_D17pt_D246/bin_util/

figure ; 

subplot(4,6,[1 2 7 8])
load D17_init_GPS_coseis_stack 
excl = [ 231   232   236   243   244  358  213  218 187];
D17.def(excl) = nan;
D17.pos.E(excl) = nan;
D17.pos.N(excl) = nan;
D17.pos.lon(excl) = nan;
D17.pos.lat(excl) = nan;

load FD
bim =  PlotAnyRGB(zeros(size(FD.dem,1),size(FD.dem,2)),FD.dem,jet,1.5,[-350 350]);
imagesc(FD.lonkm,FD.latkm,bim); hold on
scatter(D17.pos.E/1e3,D17.pos.N/1e3,100,D17.def,'.')
axis xy
caxis([min(D17.def) max(D17.def)])
plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',.5)
title('Data')
ca = axis ; 
axis equal; axis(ca)
%         set(gca,'FontSize',28);
%         xl=get(gca,'xlabel');
%         set(xl,'FontSize',28);
%         yl=get(gca,'ylabel');
%         set(yl,'FontSize',28);
%         tl=get(gca,'title');
%         set(tl,'FontSize',28);
ylabel('Northing [km]')

        
subplot(4,6,[3 4 9 10])
% calculate predicted data
%ind1 = find(isnan(FD.unw(:))==0); 
%ind2 = find(isnan(FD.unw(:))==1); 
%[x,y] = meshgrid(FD.lonkm(1:end),FD.latkm(1:end));
x = D17.pos.E/1e3;
y =  D17.pos.N/1e3;
coord = [x';y'];
dispmod = disloc(model1',coord,.25);
dispmod = FD.los'*dispmod;
%disall = zeros(size(FD.unw)); 
%disall(ind1) = dispmod; disall(ind2) = nan; 
bim =  PlotAnyRGB(zeros(size(FD.dem,1),size(FD.dem,2)),FD.dem,jet,1.5,[-350 350]);
imagesc(FD.lonkm,FD.latkm,bim); hold on
scatter(D17.pos.E/1e3,D17.pos.N/1e3,100,dispmod,'.')
axis xy
plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',.5)
caxis([min(D17.def) max(D17.def)])

title('predicted data from MAP')
set(gca,'YTickLabel',[]);
%         set(gca,'FontSize',28);
%         xl=get(gca,'xlabel');
%         set(xl,'FontSize',28);
%         yl=get(gca,'ylabel');
%         set(yl,'FontSize',28);
%         tl=get(gca,'title');
%         set(tl,'FontSize',28);
hold on; SurfProj(model1',1,0,'m');
axis xy
ca = axis ; 
axis equal; axis(ca)
        
subplot(4,6,[5 6 11 12])
err = D17.def - dispmod';
bim =  PlotAnyRGB(zeros(size(FD.dem,1),size(FD.dem,2)),FD.dem,jet,1.5,[-350 350]);
imagesc(FD.lonkm,FD.latkm,bim); hold on
scatter(D17.pos.E/1e3,D17.pos.N/1e3,100,err,'.')
axis xy
caxis([min(D17.def) max(D17.def)])
plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',.5)

title('residual')
set(gca,'YTickLabel',[]);
%         set(gca,'FontSize',28);
%         xl=get(gca,'xlabel');
%         set(xl,'FontSize',28);
%         yl=get(gca,'ylabel');
%         set(yl,'FontSize',28);
%         tl=get(gca,'title');
%         set(tl,'FontSize',28);
ca = axis ; 
axis equal; axis(ca)
        
subplot(4,6,[13 14 19 20])
load SD
bim =  PlotAnyRGB(SD.unw,SD.dem,jet,1.5,[-350 350],[min(D17.def) max(D17.def)]);
imagesc(SD.lonkm,SD.latkm,bim); hold on
axis xy
%scatter(SD,D17.pos.N/1e3,300,err,'.')

plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',.5)
%bim =  PlotAnyRGB(SD.unw,SD.dem);
%imagesc(SD.lonkm,SD.latkm,bim); 
%         set(gca,'FontSize',28);
%         xl=get(gca,'xlabel');
%         set(xl,'FontSize',28);
%         yl=get(gca,'ylabel');
%         set(yl,'FontSize',28);
%         tl=get(gca,'title');
%         set(tl,'FontSize',28);
axis xy
ca = axis ; 
axis equal; axis(ca)  
xlabel('Easting [km]')
ylabel('Northing [km]')

subplot(4,6,[15 16 21 22])
% calculate predicted data
ind1 = find(isnan(SD.unw(:))==0); 
ind2 = find(isnan(SD.unw(:))==1); 
[x,y] = meshgrid(SD.lonkm(1:end),SD.latkm(1:end));
coord = [x(ind1)';y(ind1)'];
dispmod = disloc(model1',coord,.25);
dispmod = SD.los'*dispmod;
disall = zeros(size(SD.unw)); 
disall(ind1) = dispmod; disall(ind2) = nan; 
bim =  PlotAnyRGB(disall,SD.dem,jet,1.5,[-350 350],[min(D17.def) max(D17.def)]);
imagesc(SD.lonkm,SD.latkm,bim); 
set(gca,'YTickLabel',[]);
%         set(gca,'FontSize',28);
%         xl=get(gca,'xlabel');
%         set(xl,'FontSize',28);
%         yl=get(gca,'ylabel');
%         set(yl,'FontSize',28);
%         tl=get(gca,'title');
%         set(tl,'FontSize',28);
hold on; SurfProj(model1',1,0,'m');
axis xy
ca = axis ; 
axis equal; axis(ca)
plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',.5)
xlabel('Easting [km]')

subplot(4,6,[17 18 23 24]); 
err = SD.unw - disall ; 
bim =  PlotAnyRGB(err,SD.dem,jet,1.5,[-350 350],[min(D17.def) max(D17.def)]);
imagesc(SD.lonkm,SD.latkm,bim); hold on
set(gca,'YTickLabel',[]);

%         set(gca,'FontSize',28);
%         xl=get(gca,'xlabel');
%         set(xl,'FontSize',28);
%         yl=get(gca,'ylabel');
%         set(yl,'FontSize',28);
%         tl=get(gca,'title');
%         set(tl,'FontSize',28);
axis xy
ca = axis ; 

plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',.5)
axis equal; axis(ca)
xlabel('Easting [km]')
