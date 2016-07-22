% program to calcute error covariance matrix for FD and SD interferograms

addpath /home/duttar/Japan/Tohoku/util/
load ../data/FD 
load ../data/SD
load ../data/22837.dat

imagesc(FD.unw); 
undef1 = removearea(FD.unw); 
undef1 = RemovePlane(undef1,0);
close; 
imagesc(SD.unw); 
undef2 = removearea(SD.unw); 
undef2 = RemovePlane(undef2,0);
close; 

x_inc1 = FD.lonkm(2)- FD.lonkm(1); 
y_inc1 = FD.latkm(2)- FD.latkm(1);
[~,cov1,lag1] = variogram(undef1,50,x_inc1,y_inc1); 
cov12 = 1e6*cov1;
 figure('Color',[0.8 0.8 0.8]); subplot(411)
 imagesc(FD.londeg,FD.latdeg,undef1);h = colorbar; 
  ylabel(h,'LOS error [m]','Fontsize',16)
 hold on
 plot(X22837(:,1),X22837(:,2),'k','Linewidth',2)
 xlabel('Longitude [\circ E]')
 ylabel('Latitude [\circ N]')
 title('Masked Interferogram (Track 17)')
 set(gca,'FontSize',16);
xl=get(gca,'xlabel');
set(xl,'FontSize',16);
yl=get(gca,'ylabel');
set(yl,'FontSize',16);
tl=get(gca,'title');
set(tl,'FontSize',16);
 axis xy

 hold off
 
 subplot(412)
 plot(lag1,cov12,'k.')
 hold on 
 init = [20 10 400 600]; 
  [newparameters,error]=lsqcurvefit(@limiter,init,lag1(:)',cov12(:)');
  z4= [-5:.01:50];
  z5 = limiter(newparameters,z4);
 plot(z4(:),z5(:),'r','Linewidth',2);
 axis([-5 50 -5 45])
 var1 = nanvar(undef1(:))*1e6;
 plot(0,var1,'b+','Linewidth',6)
 xlabel('lag distance [km]')
 ylabel('covariances InSAR [mm^2]')
 legend('covariogram','covariance function','variance')
 set(gca,'FontSize',16);
xl=get(gca,'xlabel');
set(xl,'FontSize',16);
yl=get(gca,'ylabel');
set(yl,'FontSize',16);
tl=get(gca,'title');
set(tl,'FontSize',16);
 hold off



x_inc2 = SD.lonkm(2)- SD.lonkm(1); 
y_inc2 = SD.latkm(2)- SD.latkm(1);
[~,cov2,lag2] = variogram(undef2,50,x_inc2,y_inc2); 
cov22 = 1e6*cov2;
subplot(413);
imagesc(FD.londeg,FD.latdeg,undef2);h = colorbar; 
  ylabel(h,'LOS error [m]','Fontsize',16)
 hold on
 plot(X22837(:,1),X22837(:,2),'k','Linewidth',2)
 xlabel('Longitude [\circ E]')
 ylabel('Latitude [\circ N]')
 title('Masked Interferogram (Track 246)')
 set(gca,'FontSize',16);
xl=get(gca,'xlabel');
set(xl,'FontSize',16);
yl=get(gca,'ylabel');
set(yl,'FontSize',16);
tl=get(gca,'title');
set(tl,'FontSize',16);
 axis xy

 hold off
 
 subplot(414)

plot(lag2,cov22,'k.')
 hold on 
 init = [20 5 10 30]; 
  [newparameters,error]=lsqcurvefit(@limiter,init,cov22',lag2);
  z4= [-5:.04:100];
  z5 = limiter(newparameters,z4);
 plot(z4(:),z5(:),'r','Linewidth',2);
  axis([-5 50 -5 30])
 var2 = nanvar(undef2(:))*1e6;
 plot(0,var2,'b+','Linewidth',6)
 xlabel('lag distance [km]')
 ylabel('covariances InSAR [mm^2]')
 legend('covariogram','covariance function','variance')
 set(gca,'FontSize',16);
xl=get(gca,'xlabel');
set(xl,'FontSize',16);
yl=get(gca,'ylabel');
set(yl,'FontSize',16);
tl=get(gca,'title');
set(tl,'FontSize',16);
 
 hold off


np = [2.5e-5 20 1 2];
z4= [0:4:100];
  z5 = limiter(np,z4);
plot(z4,z5)




























