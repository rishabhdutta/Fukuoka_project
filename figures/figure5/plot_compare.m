
% plot the samples and compare 
clear; close all; clc
addpath ../..

load samples_nopriors_again
no_prior = unburnedsamples; 
okada_no_prior = change_okada(no_prior); 
moment_no_prior = momentmag(okada_no_prior);

load samples_magprior0037_paper.mat
% nsamples = 100000;
% wsize = 10;
% X0= 1.0e+03 *  [  0.6026    3.7355    0.6141    3.7291   ...
%     0.0138    0.0010   -0.0910   -0.0016    0   0];
% samples = zeros(nsamples,numel(X0),wsize);
% samples = unburnedsamples(round(0.10*size(samples,1)):end,:,:); %burn in removed
% unburnedsamples = reshape(permute(samples, [1 3 2]),size(samples,1)*size(samples,3),size(samples,2));

mag_prior = unburnedsamples; 
okada_mag_prior = change_okada(mag_prior); 
moment_mag_prior = momentmag(okada_mag_prior);

load samples_ASprior_paper004.mat
nsamples = 100000;
wsize = 10;
X0= 1.0e+03 *  [  0.6026    3.7355    0.6141    3.7291   ...
    0.0138    0.0010   -0.0910   -0.0016    0   0];
samples = zeros(nsamples,numel(X0),wsize);
samples = unburnedsamples(round(0.10*size(samples,1)):end,:,:); %burn in removed
unburnedsamples = reshape(permute(samples, [1 3 2]),size(samples,1)*size(samples,3),size(samples,2));

AS_prior = unburnedsamples; 
okada_AS_prior = change_okada(AS_prior); 
moment_AS_prior = momentmag(okada_AS_prior);

load samples_allpriors_again_paper_mag0037_AS004.mat
% nsamples = 100000;
% wsize = 10;
% X0= 1.0e+03 *  [  0.6026    3.7355    0.6141    3.7291   ...
%     0.0138    0.0010   -0.0910   -0.0016    0   0];
% samples = zeros(nsamples,numel(X0),wsize);
% samples = unburnedsamples(round(0.10*size(samples,1)):end,:,:); %burn in removed
% unburnedsamples = reshape(permute(samples, [1 3 2]),size(samples,1)*size(samples,3),size(samples,2));

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

[nx2,ny2] = hist(moment_mag_prior,nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(moment_AS_prior,nos);
plot(ny3,nx3,'c','Linewidth',2)

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
plot(ny3,nx3,'c','Linewidth',2)

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
plot(ny3,nx3,'c','Linewidth',2)

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
plot(ny3,nx3,'c','Linewidth',2)

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
plot(ny3,nx3,'c','Linewidth',2)

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
plot(ny3,nx3,'c','Linewidth',2)

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
plot(ny3,nx3,'c','Linewidth',2)

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
plot(ny3,nx3,'c','Linewidth',2)

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
plot(ny3,nx3,'c','Linewidth',2)

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
plot(ny3,nx3,'c','Linewidth',2)

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
plot(ny3,nx3,'c','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Strike [N^\circE]','FontSize', 12)
ylabel('k)','FontSize', 12)


subplot(5,6,[21 22]+6)
[nx,ny] = hist(-okada_with_prior(:,4),nos);
plot(ny,nx,'r','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx1,ny1] = hist(-okada_no_prior(:,4),nos);
plot(ny1,nx1,'b','Linewidth',2)
[nx2,ny2] = hist(-okada_mag_prior(:,4),nos);
plot(ny2,nx2,'g','Linewidth',2)

[nx3,ny3] = hist(-okada_AS_prior(:,4),nos);
plot(ny3,nx3,'c','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Dip [\circ]','FontSize', 12)
ylabel('l)','FontSize', 12)


subplot(5,6,[23 24]+6)
[nx1,ny1] = hist(-okada_no_prior(:,8),nos);
plot(ny1,nx1,'b','Linewidth',2)
set(gca,'YTickLabel',[]);
hold on ;
[nx3,ny3] = hist(-okada_AS_prior(:,8),nos);
plot(ny3,nx3,'c','Linewidth',2)
[nx2,ny2] = hist(-okada_mag_prior(:,8),nos);
plot(ny2,nx2,'g','Linewidth',2)


[nx,ny] = hist(-okada_with_prior(:,8),nos);
plot(ny,nx,'r','Linewidth',2)

asx = [ny ny1 ny2 ny3];
asy = [nx nx1 nx2 nx3];
axis([min(asx) max(asx) min(asy) max(asy)])
xlabel('Strike-slip [m]','FontSize', 12)
ylabel('m)','FontSize', 12)

legend('no priors','prior: Aftershock locations','prior: Moment magnitude','priors: Moment magnitude + Aftershock locations','Fontsize',12)

