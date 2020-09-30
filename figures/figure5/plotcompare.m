%% Figure5

% plot the samples and compare 
clear; close all; clc
addpath ../figure3/fig_extras/

load samples_nopriors_again
no_prior = unburnedsamples; 
okada_no_prior = change_okada(no_prior); 

load samples_priormag_again
mag_prior = unburnedsamples; 
okada_mag_prior = change_okada(mag_prior); 

load samples_priorAS_again
AS_prior = unburnedsamples; 
okada_AS_prior = change_okada(AS_prior); 

load ../figure4/samples_med_again2.mat
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


