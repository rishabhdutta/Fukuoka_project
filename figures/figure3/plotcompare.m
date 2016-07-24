%% Data misfit with MAP 

addpath ../../..
model = 1.0e+03 * [0.6048    3.7359    0.6139    3.7304    0.0174    0.0000   -0.0951   -0.0017         0         0]; 
model1 = change_okada(model);
magmodel = momentmag(model1); 
load 22837.dat
coastutm = FO_CrdTrans(X22837',1);
coastutm(2,:) = coastutm(2,:) -.3; 
coastutm(1,:) = coastutm(1,:) -.37; 



addpath ../../../../data/
addpath ../../../../bin_util/

figure ; 

subplot(4,6,[1 2 7 8])
load D17_init_GPS_coseis_stack
excl = [ 231   232   236   243   244];  
D17.pos.E(excl) = []; D17.pos.N(excl) = []; D17.def(excl) = []; 
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
%%

% plot the samples and compare 
clear; close all; clc
cd /home/duttar/Desktop/insar/japan/Fukuoka/modeling/new_work/GPS_D17pt_D246/no_prior/
addpath ../..

load samples_nopriors_again
no_prior = unburnedsamples; 
okada_no_prior = change_okada(no_prior); 

cd ../prior_mag/
load samples_med_again2
mag_prior = unburnedsamples; 
okada_mag_prior = change_okada(mag_prior); 

cd ../prior_AS/
load samples_med_again
AS_prior = unburnedsamples; 
okada_AS_prior = change_okada(AS_prior); 

cd ../prior_mag_AS/
load samples_med_again.mat
with_prior = unburnedsamples; 
okada_with_prior = change_okada(with_prior); 

nos = 50;
figure; % figure for the paper 

subplot(4,6,[1 2])
[nx,ny] = hist(with_prior(:,1),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(no_prior(:,1),nos);
plot(ny1,nx1,'b','Linewidth',2)

[nx2,ny2] = hist(mag_prior(:,1),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(AS_prior(:,1),nos);
plot(ny3,nx3,'m','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
%ylabel('Longitude [km]')
%set(get(gca,'YLabel'),'Rotation',0)
xlabel('West X [km]','FontSize', 12)

subplot(4,6,[3 4])
[nx,ny] = hist(okada_with_prior(:,6),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(okada_no_prior(:,6),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(okada_mag_prior(:,6),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(okada_AS_prior(:,6),nos);
plot(ny3,nx3,'m','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Center X [km]','FontSize', 12)

subplot(4,6,[5 6])
[nx,ny] = hist(with_prior(:,3),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(no_prior(:,3),nos);
plot(ny1,nx1,'b','Linewidth',2)

[nx2,ny2] = hist(mag_prior(:,3),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(AS_prior(:,3),nos);
plot(ny3,nx3,'m','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('East X [km]','FontSize', 12)


subplot(4,6,[7 8])
[nx,ny] = hist(with_prior(:,2),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(no_prior(:,2),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(mag_prior(:,2),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(AS_prior(:,2),nos);
plot(ny3,nx3,'m','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])

xlabel('West Y [km]','FontSize', 12)


subplot(4,6,[9 10])
[nx,ny] = hist(okada_with_prior(:,7),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(okada_no_prior(:,7),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(okada_mag_prior(:,7),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(okada_AS_prior(:,7),nos);
plot(ny3,nx3,'m','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Center Y [km]','FontSize', 12)


subplot(4,6,[11 12])
[nx,ny] = hist(with_prior(:,4),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(no_prior(:,4),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(mag_prior(:,4),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(AS_prior(:,4),nos);
plot(ny3,nx3,'m','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('East Y [km]','FontSize', 12)

subplot(4,6,[13 14])
[nx,ny] = hist(okada_with_prior(:,1),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(okada_no_prior(:,1),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(okada_mag_prior(:,1),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(okada_AS_prior(:,1),nos);
plot(ny3,nx3,'m','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Length [km]','FontSize', 12)


subplot(4,6,[15 16])
[nx,ny] = hist(okada_with_prior(:,2),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(okada_no_prior(:,2),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(okada_mag_prior(:,2),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(okada_AS_prior(:,2),nos);
plot(ny3,nx3,'m','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Width [km]','FontSize', 12)

subplot(4,6,[17 18])
[nx,ny] = hist(okada_with_prior(:,3),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(okada_no_prior(:,3),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(okada_mag_prior(:,3),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(okada_AS_prior(:,3),nos);
plot(ny3,nx3,'m','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Depth [km]','FontSize', 12)

subplot(4,6,[19 20])
[nx,ny] = hist(okada_with_prior(:,5),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(okada_no_prior(:,5),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(okada_mag_prior(:,5),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(okada_AS_prior(:,5),nos);
plot(ny3,nx3,'m','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Strike [°]','FontSize', 12)


subplot(4,6,[21 22])
[nx,ny] = hist(okada_with_prior(:,4),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(okada_no_prior(:,4),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(okada_mag_prior(:,4),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(okada_AS_prior(:,4),nos);
plot(ny3,nx3,'m','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Dip [°]','FontSize', 12)



subplot(4,6,[23 24])
[nx,ny] = hist(okada_with_prior(:,8),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(okada_no_prior(:,8),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(okada_mag_prior(:,8),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(okada_AS_prior(:,8),nos);
plot(ny3,nx3,'m','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Strike-slip [m]','FontSize', 12)

legend('with priors','no priors','Fontsize',12)


%%

% plot CFS
load reliable_stress

addpath ../../../data/
load 22837.dat
coastutm = FO_CrdTrans(X22837',1);
coastutm(2,:) = coastutm(2,:) -.3; 
coastutm(1,:) = coastutm(1,:) -.37;

figure;


subplot(222)
coul = zeros(2251,81,81);
for i = 1:2251
    coul(i,:,:) = coulomb(:,:,i);
end
mean_coul = mean(coul);
mean_coul = mean_coul(:);
mean_coul= reshape(mean_coul,81,81);
imagesc(a1(1:81:end)/1e3,b1(1:81)/1e3,mean_coul)
axis xy
caxis([-5e5 5e5])
hold on
plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',2)


subplot(223)
coul = zeros(2251,81,81);
for i = 1:2251
    coul(i,:,:) = coulomb(:,:,i);
end
std_coul = std(coul);
std_coul = std_coul(:);
std_coul= reshape(std_coul,81,81);
imagesc(a1(1:81:end)/1e3,b1(1:81)/1e3,std_coul)
axis xy
caxis([0 5e5])
hold on
plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',2)


subplot(224)

covar = std_coul./mean_coul;
imagesc(a1(1:81:end)/1e3,b1(1:81)/1e3,abs(covar))
axis xy
caxis([0 1])
hold on
plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',2)
