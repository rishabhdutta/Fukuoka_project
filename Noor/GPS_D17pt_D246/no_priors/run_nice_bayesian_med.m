clear all; close all; clc

addpath bin_util
addpath extras

NonlinFO_1
mbest= 1.0e+03 *  [  0.6026    3.7355    0.6141    3.7291   ...
    0.0138    0.0010   -0.0910   -0.0016    0   0];
varnum = numel(mbest);
nsamples = 100000;
wsize = 10;
X0 = mbest;
sqrflag=1;

SE =.01*[2 2 2.5 2 3 1 10 .5 0 0];    
 
proprnd=@(x) mvnrnd(x,SE);  % proposal random sampler

priorfrac = 1000;     % determines how strong is the prior

posterior=@(x) postfunc_med(x,priorfrac,d,coord,UD1,UD2,nu,W,fixind,fixpar,sqrflag,datind,parind);

%unburnedsamples = slicesample(mbest,100,'logpdf',posterior);
% 

tic

samples = zeros(nsamples,numel(X0),wsize);
acceptances = zeros(wsize,1);

parfor w =1:wsize
    [unburnedsamples(:,:,w),acceptances(w)] = mhsample(X0,nsamples,'pdf',...
        posterior,'proprnd',proprnd,'symmetric',1);
end
%warning('on');
%matlabpool close;

b = toc;
save('samples_no_priors_paper.mat','unburnedsamples','b','-v7.3')


disp('Total time sampling');
disp('Collected rows per second');
(wsize*nsamples) / b


disp('Average acceptance ratio (preferably 5-40% is usually good)');
b=mean(acceptances)

disp('Merging and burning chains');

samples = unburnedsamples(round(0.10*size(samples,1)):end,:,:); %burn in removed
unburnedsamples = reshape(permute(samples, [1 3 2]),size(samples,1)*size(samples,3),size(samples,2));

binss = 20; 
for i=1:8
    if i==1
        subplot(4,2,i)
        hist(unburnedsamples(:,i),binss);
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','r','EdgeColor','r')

        xlabel('X-loc left [km]')
    elseif i==2
        subplot(4,2,i)
        hist(unburnedsamples(:,i),binss);
        xlabel('Y-loc left [km]')
    elseif i==3
        subplot(4,2,i)
        hist(unburnedsamples(:,i),binss);
        xlabel('X-loc right [km]')
    elseif i==4
        subplot(4,2,i)
        hist(unburnedsamples(:,i),binss);
        xlabel('Y-loc right [degree]')
    elseif i==5
        subplot(4,2,i)
        hist(unburnedsamples(:,i),binss);
        xlabel('width [degree]')
    elseif i==6
        subplot(4,2,i)
        hist(unburnedsamples(:,i),binss);
        xlabel('depth [0.1km]')
    elseif i==7
        subplot(4,2,i)
        hist(unburnedsamples(:,i),binss);
        xlabel('dip [0.1km]')
    elseif i==8
        
        subplot(4,2,i)
        hist(unburnedsamples(:,i),binss);
        xlabel('strikeslip')
    end
    
end

figure;
for i=1:8
    subplot(4,2,i)
    plot(unburnedsamples(:,i))
end
% 
% disp('Estimating Geyer time tau, this is a somehwat conservative guess with chain of many million')
% disp('Probably you can well do with less thinning. But with a chain of tens of thousands, it might be too little.')
% 
% TF = max(geyer_imse(unburnedsamples(:,1:8)))



% map value 818813
% map =  1.0e+03 * [0.6048    3.7359    0.6139    3.7304    0.0174    0.0000   -0.0951   -0.0017         0         0]; 

