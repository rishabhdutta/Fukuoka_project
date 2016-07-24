
%% Data misfit with MAP 

addpath fig_extras
model = 1.0e+03 * [0.6048    3.7359    0.6139    3.7304    0.0174    0.0000   -0.0951   -0.0017         0         0]; 
model1 = change_okada(model);
magmodel = momentmag(model1); 

addpath ../../Noor/GPS_D17pt_D246/extras/
addpath ../../Noor/GPS_D17pt_D246/bin_util/

% figure ; 
% 
% subplot(4,6,[1 2 7 8])
 load D17_init_GPS_coseis_stack
 load FD
% excl = [ 231   232   236   243   244];
% D17.def(excl) = nan;
% bim =  PlotAnyRGB(zeros(size(FD.dem,1),size(FD.dem,2)),FD.dem,jet,1.5,[-50 50],[],[1 2]);
% imagesc(FD.lonkm,FD.latkm,bim); hold on
% scatter(D17.pos.E/1e3,D17.pos.N/1e3,300,D17.def,'.')
% axis xy
% caxis([min(D17.def) max(D17.def)])
% 
% title('Data')
%         set(gca,'FontSize',12);
%         xl=get(gca,'xlabel');
%         set(xl,'FontSize',12);
%         yl=get(gca,'ylabel');
%         set(yl,'FontSize',12);
%         tl=get(gca,'title');
%         set(tl,'FontSize',12);
% 
%         
% subplot(4,6,[3 4 9 10])
% % calculate predicted data
% %ind1 = find(isnan(FD.unw(:))==0); 
% %ind2 = find(isnan(FD.unw(:))==1); 
% %[x,y] = meshgrid(FD.lonkm(1:end),FD.latkm(1:end));
% x = D17.pos.E/1e3;
% y =  D17.pos.N/1e3;
% coord = [x';y'];
% dispmod = disloc(model1',coord,.25);
% dispmod = FD.los'*dispmod;
% %disall = zeros(size(FD.unw)); 
% %disall(ind1) = dispmod; disall(ind2) = nan; 
% bim =  PlotAnyRGB(zeros(size(FD.dem,1),size(FD.dem,2)),FD.dem,jet,1.5,[-50 50],[],[1 2]);
% imagesc(FD.lonkm,FD.latkm,bim); hold on
% scatter(D17.pos.E/1e3,D17.pos.N/1e3,300,dispmod,'.')
% axis xy
% caxis([min(D17.def) max(D17.def)])
% 
% title('predicted data from MAP')
% set(gca,'YTickLabel',[]);
%         set(gca,'FontSize',12);
%         xl=get(gca,'xlabel');
%         set(xl,'FontSize',12);
%         yl=get(gca,'ylabel');
%         set(yl,'FontSize',12);
%         tl=get(gca,'title');
%         set(tl,'FontSize',12);
% hold on; SurfProj(model1',1,0,'m');
% axis xy
%         
% subplot(4,6,[5 6 11 12])
% err = D17.def - dispmod';
% bim =  PlotAnyRGB(zeros(size(FD.dem,1),size(FD.dem,2)),FD.dem,jet,1.5,[-50 50],[],[1 2]);
% imagesc(FD.lonkm,FD.latkm,bim); hold on
% scatter(D17.pos.E/1e3,D17.pos.N/1e3,300,err,'.')
% axis xy
% caxis([min(D17.def) max(D17.def)])
% c=colorbar;
% ylabel(c,'[m]')
% set(gca,'FontSize',12);
% xl=get(gca,'xlabel');
% set(xl,'FontSize',12);
% yl=get(gca,'ylabel');
% set(yl,'FontSize',12);
% tl=get(gca,'title');
% set(tl,'FontSize',12);
% 
% title('residual')
% set(gca,'YTickLabel',[]);
% %         set(gca,'FontSize',28);
% %         xl=get(gca,'xlabel');
% %         set(xl,'FontSize',28);
% %         yl=get(gca,'ylabel');
% %         set(yl,'FontSize',28);
% %         tl=get(gca,'title');
% %         set(tl,'FontSize',28);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; 
subplot(4,6,[1 2 7 8])
load FD
bim =  PlotAnyRGB(FD.unw,FD.dem,jet,1.5,[-100 100],[min(D17.def) max(D17.def)],[1 2]);
imagesc(FD.lonkm,FD.latkm,bim); hold on
axis xy
%scatter(SD,D17.pos.N/1e3,300,err,'.')

%bim =  PlotAnyRGB(SD.unw,SD.dem);
%imagesc(SD.lonkm,SD.latkm,bim); 
        set(gca,'FontSize',12);
        xl=get(gca,'xlabel');
        set(xl,'FontSize',12);
        yl=get(gca,'ylabel');
        set(yl,'FontSize',12);
        tl=get(gca,'title');
        set(tl,'FontSize',12);
axis xy
ca = axis; 
axis equal
axis(ca)
        
subplot(4,6,[3 4 9 10])
% calculate predicted data
ind1 = find(isnan(FD.unw(:))==0); 
ind2 = find(isnan(FD.unw(:))==1); 
[x,y] = meshgrid(FD.lonkm(1:end),FD.latkm(1:end));
coord = [x(ind1)';y(ind1)'];
dispmod = disloc(model1',coord,.25);
dispmod = FD.los'*dispmod;
disall = zeros(size(FD.unw)); 
disall(ind1) = dispmod; disall(ind2) = nan; 
bim =  PlotAnyRGB(disall,FD.dem,jet,1.5,[-100 100],[min(D17.def) max(D17.def)],[1 2]);
imagesc(FD.lonkm,FD.latkm,bim); 
set(gca,'YTickLabel',[]);
        set(gca,'FontSize',12);
        xl=get(gca,'xlabel');
        set(xl,'FontSize',12);
        yl=get(gca,'ylabel');
        set(yl,'FontSize',12);
        tl=get(gca,'title');
        set(tl,'FontSize',12);
hold on; SurfProj(model1',1,0,'m');
axis xy
ca = axis; 
axis equal
axis(ca)
        
subplot(4,6,[5 6 11 12])
err = FD.unw - disall ; 
bim =  PlotAnyRGB(err,FD.dem,jet,1.5,[-100 100],[min(D17.def) max(D17.def)],[1 2]);
imagesc(FD.lonkm,FD.latkm,bim);
set(gca,'YTickLabel',[]);

        set(gca,'FontSize',12);
        xl=get(gca,'xlabel');
        set(xl,'FontSize',12);
        yl=get(gca,'ylabel');
        set(yl,'FontSize',12);
        tl=get(gca,'title');
        set(tl,'FontSize',12);
axis xy
ca = axis; 
axis equal
axis(ca)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        
subplot(4,6,[13 14 19 20])
load SD
bim =  PlotAnyRGB(SD.unw,SD.dem,jet,1.5,[-100 100],[min(D17.def) max(D17.def)],[1 2]);
imagesc(SD.lonkm,SD.latkm,bim); hold on
axis xy
%scatter(SD,D17.pos.N/1e3,300,err,'.')

%bim =  PlotAnyRGB(SD.unw,SD.dem);
%imagesc(SD.lonkm,SD.latkm,bim); 
        set(gca,'FontSize',12);
        xl=get(gca,'xlabel');
        set(xl,'FontSize',12);
        yl=get(gca,'ylabel');
        set(yl,'FontSize',12);
        tl=get(gca,'title');
        set(tl,'FontSize',12);
axis xy
ca = axis; 
axis equal
axis(ca)

        
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
bim =  PlotAnyRGB(disall,SD.dem,jet,1.5,[-100 100],[min(D17.def) max(D17.def)],[1 2]);
imagesc(SD.lonkm,SD.latkm,bim); 
set(gca,'YTickLabel',[]);
        set(gca,'FontSize',12);
        xl=get(gca,'xlabel');
        set(xl,'FontSize',12);
        yl=get(gca,'ylabel');
        set(yl,'FontSize',12);
        tl=get(gca,'title');
        set(tl,'FontSize',12);
hold on; SurfProj(model1',1,0,'m');
axis xy
ca = axis; 
axis equal
axis(ca)

        
subplot(4,6,[17 18 23 24])
err = SD.unw - disall ; 
bim =  PlotAnyRGB(err,SD.dem,jet,1.5,[-100 100],[min(D17.def) max(D17.def)],[1 2]);
imagesc(SD.lonkm,SD.latkm,bim);
set(gca,'YTickLabel',[]);

        set(gca,'FontSize',12);
        xl=get(gca,'xlabel');
        set(xl,'FontSize',12);
        yl=get(gca,'ylabel');
        set(yl,'FontSize',12);
        tl=get(gca,'title');
        set(tl,'FontSize',12);
axis xy
ca = axis; 
axis equal
axis(ca)

%%

% plot the samples and compare 
clear; close all; clc
cd /home/duttar/Desktop/insar/japan/Fukuoka/modeling/new_work/GPS_D17pt_D246/no_prior/
addpath ../..

load samples_nopriors_again
no_prior = unburnedsamples; 
okada_no_prior = change_okada(no_prior); 
moment_no_prior = momentmag(okada_no_prior);

cd ../prior_mag/
load samples_med_again2
mag_prior = unburnedsamples; 
okada_mag_prior = change_okada(mag_prior); 
moment_mag_prior = momentmag(okada_mag_prior);

cd ../prior_AS/
load samples_med_again
AS_prior = unburnedsamples; 
okada_AS_prior = change_okada(AS_prior); 
moment_AS_prior = momentmag(okada_AS_prior);

cd ../prior_mag_AS/
load samples_med_again.mat
with_prior = unburnedsamples; 
okada_with_prior = change_okada(with_prior); 
moment_with_prior = momentmag(okada_with_prior);

nos = 50;
figure; % figure for the paper 

subplot(5,6,[1 2])
[nx,ny] = hist(moment_with_prior,nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(moment_no_prior,nos);
plot(ny1,nx1,'b','Linewidth',2)

[nx2,ny2] = hist(moment_mag_prior+.02,nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(moment_AS_prior,nos);
plot(ny3,nx3,'y','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
%ylabel('Longitude [km]')
%set(get(gca,'YLabel'),'Rotation',0)
xlabel('Moment magnitude','FontSize', 12)
ylabel('a)','FontSize', 12)


subplot(5,6,[1 2]+6)
[nx,ny] = hist(with_prior(:,1),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(no_prior(:,1),nos);
plot(ny1,nx1,'b','Linewidth',2)

[nx2,ny2] = hist(mag_prior(:,1),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(AS_prior(:,1),nos);
plot(ny3,nx3,'y','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
%ylabel('Longitude [km]')
%set(get(gca,'YLabel'),'Rotation',0)
xlabel('West X [km]','FontSize', 12)
ylabel('b)','FontSize', 12)


subplot(5,6,[3 4]+6)
[nx,ny] = hist(okada_with_prior(:,6),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(okada_no_prior(:,6),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(okada_mag_prior(:,6),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(okada_AS_prior(:,6),nos);
plot(ny3,nx3,'y','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Center X [km]','FontSize', 12)
ylabel('c)','FontSize', 12)

subplot(5,6,[5 6]+6)
[nx,ny] = hist(with_prior(:,3),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(no_prior(:,3),nos);
plot(ny1,nx1,'b','Linewidth',2)

[nx2,ny2] = hist(mag_prior(:,3),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(AS_prior(:,3),nos);
plot(ny3,nx3,'y','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('East X [km]','FontSize', 12)
ylabel('d)','FontSize', 12)


subplot(5,6,[7 8]+6)
[nx,ny] = hist(with_prior(:,2),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(no_prior(:,2),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(mag_prior(:,2),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(AS_prior(:,2),nos);
plot(ny3,nx3,'y','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])

xlabel('West Y [km]','FontSize', 12)
ylabel('e)','FontSize', 12)

subplot(5,6,[9 10]+6)
[nx,ny] = hist(okada_with_prior(:,7),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(okada_no_prior(:,7),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(okada_mag_prior(:,7),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(okada_AS_prior(:,7),nos);
plot(ny3,nx3,'y','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Center Y [km]','FontSize', 12)
ylabel('f)','FontSize', 12)


subplot(5,6,[11 12]+6)
[nx,ny] = hist(with_prior(:,4),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(no_prior(:,4),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(mag_prior(:,4),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(AS_prior(:,4),nos);
plot(ny3,nx3,'y','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('East Y [km]','FontSize', 12)
ylabel('g)','FontSize', 12)


subplot(5,6,[13 14]+6)
[nx,ny] = hist(okada_with_prior(:,1),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(okada_no_prior(:,1),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(okada_mag_prior(:,1),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(okada_AS_prior(:,1),nos);
plot(ny3,nx3,'y','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Length [km]','FontSize', 12)
ylabel('h)','FontSize', 12)

subplot(5,6,[15 16]+6)
[nx,ny] = hist(okada_with_prior(:,2),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(okada_no_prior(:,2),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(okada_mag_prior(:,2),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(okada_AS_prior(:,2),nos);
plot(ny3,nx3,'y','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Width [km]','FontSize', 12)
ylabel('i)','FontSize', 12)


subplot(5,6,[17 18]+6)
[nx,ny] = hist(okada_with_prior(:,3),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(okada_no_prior(:,3),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(okada_mag_prior(:,3),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(okada_AS_prior(:,3),nos);
plot(ny3,nx3,'y','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Depth [km]','FontSize', 12)
ylabel('j)','FontSize', 12)


subplot(5,6,[19 20]+6)
[nx,ny] = hist(okada_with_prior(:,5),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(okada_no_prior(:,5),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(okada_mag_prior(:,5),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(okada_AS_prior(:,5),nos);
plot(ny3,nx3,'y','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Strike [\circ]','FontSize', 12)
ylabel('k)','FontSize', 12)


subplot(5,6,[21 22]+6)
[nx,ny] = hist(okada_with_prior(:,4),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(okada_no_prior(:,4),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(okada_mag_prior(:,4),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(okada_AS_prior(:,4),nos);
plot(ny3,nx3,'y','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Dip [\circ]','FontSize', 12)
ylabel('l)','FontSize', 12)


subplot(5,6,[23 24]+6)
[nx,ny] = hist(okada_with_prior(:,8),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(okada_no_prior(:,8),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(okada_mag_prior(:,8),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(okada_AS_prior(:,8),nos);
plot(ny3,nx3,'y','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Strike-slip [m]','FontSize', 12)
ylabel('m)','FontSize', 12)

legend('priors: Moment magnitude + aftershock locations','no priors','prior: moment magnitude','prior: Aftershock locations','Fontsize',12)


%%

% plot CFS
load reliable_stress
addpath ../../../bin_util/
load slipcolor

addpath ../../../data/
load 22837.dat
coastutm = FO_CrdTrans(X22837',1);
coastutm(2,:) = coastutm(2,:) -.3; 
coastutm(1,:) = coastutm(1,:) -.37;

addpath ../../
figure;
subplot(4,4,[1 2 5 6])
model = 1.0e+03 * [0.6048    3.7359    0.6139    3.7304    0.0174    0.0000   -0.0951   -0.0017         0         0]; 
model1 = change_okada(model);
magmodel = momentmag(model1);
 m = model1;
 m(:,1:3) = m(:,1:3).*1e3;
    m(:,6:7) = m(:,6:7)*1e3;
    
    figax = [0.5992    0.6570    3.6930    3.7543]*1e6; 
    inc = 200;
    [a1,b1] = meshgrid(figax(1):inc:figax(2),figax(4):-inc:figax(3));
    a2 = a1(:); b2 = b1(:); 
    
    c=-7000*ones(1,numel(b1));
    %make 100 km by 100km area with 1km increment around the epicenter
    x= [a2'; b2';c];
    xlen = size(a1,1); ylen = size(a1,2); 
    
    % calculate the stress changes
    [~,~,S]=disloc3d(m',x,30e9,.28);
    S1=S(1,:); S2=S(2,:); S3=S(3,:); 
    S4=S(4,:); S5=S(5,:); S6=S(6,:);

    nu=.4;  % effective friction coeff
    n=size(S,2);
    
    % the normal to the recipient fault is (unit vector)
    normalfault = [1/sqrt(2); 1/sqrt(2); 0];  
    alongfault = [1/sqrt(2);-1/sqrt(2) ; 0];
    parfor k = 1:n
        % stress change matrix for the location
        stressmat = [S1(k) S2(k) S3(k);S2(k) S4(k) S5(k); ...
            S3(k) S5(k) S6(k)];
        
        thet = -43.6; 
        %transformation vector
        at=[cosd(thet) sind(thet) 0; -sind(thet) cosd(thet) 0 ...
            ; 0 0 1];
        
        % stress mat aligned along recipient fault
        stressmatrec = at*stressmat*at';
       
        normalstr(k) = stressmatrec(5);
        shearstr(k) = -stressmatrec(2);  % minus sign means right-lat is positive
        
        coulombstr(k) = normalstr(k)*nu + shearstr(k) ; 
        
    end
coulombmap= reshape(coulombstr,xlen,ylen);
    shearmap= reshape(shearstr,xlen,ylen);
    normalmap= reshape(normalstr,xlen,ylen);
    
imagesc(a1(1:307:end)/1e3,b1(1:307)/1e3,coulombmap/1e6)
caxis([-.5 .5])
axis xy 
colorbar
hold on
plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',2)
c=colorbar;
ylabel(c,'\DeltaCFS [MPa]')

subplot(4,4,[3 4 7 8])
coul = zeros(2251,307,290);
for i = 1:2251
    coul(i,:,:) = coulomb(:,:,i);
end
mean_coul = mean(coul);
mean_coul = mean_coul(:);
mean_coul= reshape(mean_coul,307,290);

% addpath ../../../../data/
% load SD.mat 
% newdem = griddata(SD.lonkm,SD.latkm,SD.dem,a1/1e3,b1/1e3);
% load colormapcoulomb.mat
% [rgbim] = PlotAnyRGB(mean_coul,newdem,colormapcoulomb,1.5,[-100 100],[-.5e6 .5e6],[1 2]);
%  imagesc(SD.lonkm,SD.latkm,rgbim)
%  axis xy 
%   hold on
  
imagesc(a1(1:307:end)/1e3,b1(1:307)/1e3,mean_coul/1e6)
axis xy
caxis([-.5 .5])
hold on
plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',2)
colorbar
c=colorbar;
ylabel(c,'mean \DeltaCFS [MPa]')
set(gca,'YTickLabel',[]);


subplot(4,4,[9 10 13 14])
coul = zeros(2251,307,290);
for i = 1:2251
    coul(i,:,:) = coulomb(:,:,i);
end
std_coul = std(coul);
std_coul = std_coul(:);
std_coul= reshape(std_coul,307,290);
imagesc(a1(1:307:end)/1e3,b1(1:307)/1e3,log10(std_coul))
axis xy
caxis([3 7])
hold on
plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',2)
colorbar
c=colorbar;
ylabel(c,'log_{10} std \DeltaCFS [MPa]')




subplot(4,4,[11 12 15 16])

covar = std_coul./mean_coul;
imagesc(a1(1:307:end)/1e3,b1(1:307)/1e3,abs(covar))
axis xy
caxis([-0.3 .8])
hold on
plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',2)
%colormap(slipcolor)
colorbar
c=colorbar;
ylabel(c,'abs COV \DeltaCFS [MPa]')
set(gca,'YTickLabel',[]);

%% Figure 8: variability in some locations
 
locations = [130.391118 130.292938 129.78434 130.21081 130.23418; ...
    33.606178 33.667665 33.749971 33.649049 33.683618;];

addpath ../../../../data

loc_utm = FO_CrdTrans(locations,1);

plot(loc_utm(1,:),loc_utm(2,:),'.')

genkai = [614.3;3727]; 
fuku1 = [629.1;3719];
nishinoura = [612.3;3724];
shikano = [619.7;3726];

subplot(5,4,[1:16])
covar = std_coul./mean_coul;
imagesc(a1(1:307:end)/1e3,b1(1:307)/1e3,abs(covar))
ca = axis; 
axis equal
axis(ca)
axis xy
caxis([-0.3 .8])
hold on
plot(coastutm(1,:),coastutm(2,:),'k','Linewidth',2)
set(gca,'XAxislocation','top');
%colormap(slipcolor)
colorbar
c=colorbar;
ylabel(c,'abs COV \DeltaCFS [MPa]')
%set(gca,'YTickLabel',[]);

colk(1) = dsearchn(a1(1,:)',[genkai(1)*1e3]);
rowk(1) =  dsearchn(b1(1:307,1),[genkai(2)*1e3]);

scatter(a1(rowk(1),colk(1))/1e3,b1(rowk(1),colk(1))/1e3,2500,'b.')
text(a1(rowk(1),colk(1))/1e3,b1(rowk(1),colk(1))/1e3,'1')

colk(2) = dsearchn(a1(1,:)',[fuku1(1)*1e3]);
rowk(2) =  dsearchn(b1(1:307,1),[fuku1(2)*1e3]);

scatter(a1(rowk(2),colk(2))/1e3,b1(rowk(2),colk(2))/1e3,2500,'b.')
text(a1(rowk(2),colk(2))/1e3,b1(rowk(2),colk(2))/1e3,'2')

colk(3) = dsearchn(a1(1,:)',[nishinoura(1)*1e3]);
rowk(3) =  dsearchn(b1(1:307,1),[nishinoura(2)*1e3]);

scatter(a1(rowk(3),colk(3))/1e3,b1(rowk(3),colk(3))/1e3,2500,'b.')
text(a1(rowk(3),colk(3))/1e3,b1(rowk(3),colk(3))/1e3,'3')

colk(4) = dsearchn(a1(1,:)',[shikano(1)*1e3]);
rowk(4) =  dsearchn(b1(1:307,1),[shikano(2)*1e3]);

scatter(a1(rowk(4),colk(4))/1e3,b1(rowk(4),colk(4))/1e3,2500,'b.')
text(a1(rowk(4),colk(4))/1e3,b1(rowk(4),colk(4))/1e3,'4')

subplot(5,4,[17])
coul_loc1 = coulomb(rowk(1),colk(1),:);
coul_loc1 =reshape(coul_loc1,2251,1);
hist(coul_loc1/1e6,100)
axis([min(coul_loc1/1e6) max(coul_loc1/1e6) 0 70])
set(gca,'YTickLabel',[]);
title('Location 1')
xlabel('\DeltaCFS [MPa]')
h = findobj(gca,'Type','patch');
set(h,'FaceColor','k','EdgeColor','k')

subplot(5,4,[18])
coul_loc2 = coulomb(rowk(2),colk(2),:);
coul_loc2 =reshape(coul_loc2,2251,1);
hist(coul_loc2/1e6,100)
axis([min(coul_loc2/1e6) max(coul_loc2/1e6) 0 90])
set(gca,'YTickLabel',[]);
title('Location 2')
xlabel('\DeltaCFS [MPa]')
h = findobj(gca,'Type','patch');
set(h,'FaceColor','k','EdgeColor','k')

subplot(5,4,[19])
coul_loc3 = coulomb(rowk(3),colk(3),:);
coul_loc3 =reshape(coul_loc3,2251,1);
hist(coul_loc3/1e6,100)
axis([min(coul_loc3/1e6) max(coul_loc3/1e6) 0 110])
set(gca,'YTickLabel',[]);
title('Location 3')
xlabel('\DeltaCFS [MPa]')
h = findobj(gca,'Type','patch');
set(h,'FaceColor','k','EdgeColor','k')

subplot(5,4,[20])
coul_loc4 = coulomb(rowk(4),colk(4),:);
coul_loc4 =reshape(coul_loc4,2251,1);
hist(coul_loc4/1e6,100)
axis([min(coul_loc4/1e6) max(coul_loc4/1e6) 0 120])
set(gca,'YTickLabel',[]);
title('Location 4')
xlabel('\DeltaCFS [MPa]')
h = findobj(gca,'Type','patch');
set(h,'FaceColor','k','EdgeColor','k')




