 %% Figure 2 for the paper
 
 clear all; close all; clc
 addpath ../../Noor/GPS_D17pt_D246/extras/
 addpath ../../Noor/GPS_D17pt_D246/bin_util/ 
load D17_init_GPS_coseis_stack.mat
load SD.mat
load FD.mat
load 22837.dat

coastutm = FO_CrdTrans(X22837',1);
close all; 

excl = [ 231   232   236   243   244  358  213  218 187]; 
D17.def(excl) = nan;
D17.pos.E(excl) = nan; 
D17.pos.N(excl) = nan;
D17.pos.lon(excl) = nan;
D17.pos.lat(excl) = nan;

figure; 
subplot(4,4,[1 2 5 6])

szFD = size(FD.dem);
[rgbim] = PlotAnyRGB(nan(szFD(1),szFD(2)),FD.dem,jet,1.5,[-300 300],[],[1 2]);
imagesc(FD.lonkm,FD.latkm,rgbim)
axis xy; hold on; 
scatter(D17.pos.E/1e3,D17.pos.N/1e3,200,D17.def,'.')
coastutm(2,:) = coastutm(2,:) -.3; 
coastutm(1,:) = coastutm(1,:) -.37; 
clbar = colorbar; 
ylabel(clbar,'LOS displacement [m]')

rect1 =  [0.5992    0.6570    3.6930    3.7543]*1e3;
rectangle('Position',[rect1(1) rect1(3) rect1(2)-rect1(1) rect1(4)-rect1(3)],'Edgecolor',[48 102 117]/256,'Linewidth',2)

plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',.25)
caxis([min(D17.def) max(D17.def)]);
c = caxis;
xlabel('Easting [km]')
ylabel('Northing [km]')
axss1 = axis; 
axis equal
axis(axss1)


subplot(4,4,[3 4 7 8])
[rgbim] = PlotAnyRGB(SD.unw,SD.dem,jet,1.5,[-200 200],c,[1 2]);
 imagesc(SD.lonkm,SD.latkm,rgbim)
 axis xy 
  hold on

plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',.25)
xlabel('Easting [km]')
set(gca,'yaxislocation','right');
%ylabel('Northing [km]')
axss2 = axis; 
axis equal
axis(axss2)

subplot(4,4,[11 12 15  16])
 addpath error_cov
load error_cov
   %imagesc(SD.lonkm,SD.latkm,undef2);h = colorbar; 
   addpath ../
   [rgbim] = PlotAnyRGB(undef2,SD.dem,jet,1.5,[-200 200],c,[1 2]);
   
 imagesc(SD.lonkm,SD.latkm,rgbim)
  %ylabel(h,'LOS error [m]','Fontsize',12)
 hold on
 axis xy
 plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',.25)
 axss3 = axis; 
axis equal
axis(axss3)
 
 %colormap(slipcolor)
 
 axes('Position',[.5 .5 .2 .1])
 box on
 plot(lag2,cov22,'k.')
 hold on
 z6= [-5:.04:100];
 z5 = limiter(newparameters2,z6);
 plot(z6(:),z5(:),'r','Linewidth',2);
 axis([-5 50 -5 30])
 var2 = nanvar(undef2(:))*1e6;
 plot(0,var2,'b+','Linewidth',6)
 xlabel('lag distance [km]')
 ylabel('covariances InSAR [mm^2]')
 legend('covariogram','covariance function','variance')

  
%%  
%figure; % GPS
   
subplot(4,4,[9 10 13 14])
  
load FO_GPS_geonet
load FO_GPS_TN
addpath ../../Noor/GPS_D17pt_D246/

exclude = [];
%ii=find(FO.CGPSxy(1,:)>670);       % Here I find far-eastern stat
%mn = (mean(FO.CGPSenu(:,ii)'))';   % Calculate mean displ
% and remove the mean:
mn = [0 -0.9 -5.1]'*1e-3; % Shift vector
FO.CGPSenu = FO.CGPSenu-repmat(mn,1,size(FO.CGPSenu,2));
FN= ExcludeGPS(FO,exclude);
FN.CGPSll  = [FN.CGPSll FOtn.GPSll];
FN.CGPSenu = [FN.CGPSenu FOtn.GPSenu];
FN.CGPSerr = [FN.CGPSerr FOtn.GPSerr];
FN.CGPScov = BlockDiag(FN.CGPScov,FOtn.GPScov);
FN.CGPSxy  = [FN.CGPSxy FOtn.GPSxy];

   
denu    = FN.CGPSenu;
derr    = FN.CGPSerr;
dgps    = denu(:);
 
crdenu  = FN.CGPSxy;
tmp3    = repmat(crdenu,3,1);
crdgps  = zeros(2,size(crdenu,2)*3);
crdgps(:)= tmp3(:);
Ngps    = max(size(dgps));

ip=45; % number of station
ii=[3*ip-2:3*ip];

SARSTD     = 0.01;
DataCovGPS = FN.CGPScov;
DataCovGPS(ii,ii)=2^2.*DataCovGPS(ii,ii);


% disp('---> Plotting data...')
%  figure;
% % patch(FOqt.cx,FOqt.cy,dsar'); axis image; hold on  
%  DM_Quiver(crdenu,denu(:),DataCovGPS,0.2,[0 0],'k');
% 
% addpath ~/Desktop/softwares/arrows/
% arrows(crdgps(1,1:3:end),crdgps(2,1:3:end),denu(1,:)*100, denu(2,:)*100,[0.2,0.2,0.15,0.05],'Cartesian')
% hold on
% %quiver(crdgps(1,),crdgps(2,1:82),denu(1,:)*100, denu(2,:)*100)
% scatter(crdgps(1,1:3:end),crdgps(2,1:3:end),100,denu(3,:),'.')
%   
load Fukuokadem.mat
[xydem_x,xydem_y] = meshgrid(xdem,ydem);
utmdem = FO_CrdTrans([xydem_x(:)' ; xydem_y(:)'],1);

%imagesc(utmdem(1,1:1201:end),utmdem(2,1:1201),dem); axis xy
ca = [530         700        3660        3760];
%axis equal;axis(ca); 
[rgbim] = PlotAnyRGB(nan(size(dem)),dem,jet,1.5,[-200 200],[],[1 2]);
 imagesc(utmdem(1,1:1201:end),utmdem(2,1:1201),rgbim); axis xy; axis equal ; 
 %axis(ca)
 axis(axss1)
hold on
  plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',.25)  
   
addpath arrows/
% arrows([crdgps(1,1:3:end) 550],[crdgps(2,1:3:end) 3750],[denu(1,:) .1]*100,...
%    [denu(2,:) 0]*100,[.2,.2,.15,.01],'Cartesian','FaceColor','r','Ref',15,'EdgeColor','k')
[abc1] = arrows(FO.CGPSxy(1,:),FO.CGPSxy(2,:),FO.CGPSenu(1,:)*100,FO.CGPSenu(2,:)*100,[.2,.2,.15,.01],'Cartesian','FaceColor','r','Ref',15,'EdgeColor','k');
hold on
[abc2] = arrows(FOtn.GPSxy(1,:),FOtn.GPSxy(2,:),FOtn.GPSenu(1,:)*100,FOtn.GPSenu(2,:)*100,[.2,.2,.15,.01],'Cartesian','FaceColor','m','Ref',15,'EdgeColor','k');
arrows(550,3750,10,0,[.2,.2,.15,.01],'Cartesian','FaceColor','k','Ref',15,'EdgeColor','k')
%quiver(crdgps(1,),crdgps(2,1:82),denu(1,:)*100, denu(2,:)*100)
%scatter(crdgps(1,1:3:end),crdgps(2,1:3:end),300,denu(3,:),'.')
  
legend([abc1,abc2],'Continuous GPS stations','Campaign GPS stations')
%caxis([ min(denu(3,:))  max(denu(3,:))])
caxis([min(D17.def) max(D17.def)]);

rect1 =  [0.5992    0.6570    3.6930    3.7543]*1e3;
rectangle('Position',[rect1(1) rect1(3) rect1(2)-rect1(1) rect1(4)-rect1(3)],'Edgecolor',[48 102 117]/256,'Linewidth',2)

%rect2 = [0.5381    0.6643    3.6699    3.76]*1e3; 
%rectangle('Position',[rect2(1) rect2(3) rect2(2)-rect2(1) rect2(4)-rect2(3)],'Edgecolor',[39 119 53]/256,'Linewidth',3)

xlabel('Easting [km]')
 ylabel('Northing [km]')

 c= colorbar
 ylabel(c,'GPS vertical displacement [m]')
  hold on
  text(550,3750,'GPS [5cm]')
  
showaxes('on')


