% plot the 2D kernel densities with the marginals 
close all 
addpath ~/Desktop/softwares/gkde2/
addpath ~/Desktop/softwares/mcmcdiag/
load ../figure5/samples_allpriors_again_paper_mag0037_AS004.mat

unburnedsamples2 = unburnedsamples;
unburnedsamples2(:,7:8) = -unburnedsamples2(:,7:8);

% shift x and y of endpoints 
unburnedsamples2(:,9) = unburnedsamples2(:,1);
unburnedsamples2(:,1) = unburnedsamples2(:,2);
unburnedsamples2(:,2) = unburnedsamples2(:,9);

unburnedsamples2(:,9) = unburnedsamples2(:,3);
unburnedsamples2(:,3) = unburnedsamples2(:,4);
unburnedsamples2(:,4) = unburnedsamples2(:,9);


mbest = [3.7361e+03 604.4676 3.7303e+03 613.945 ...
    16.25 0.2597 94.939 1.7175];

Tf = geyer_imse(unburnedsamples2);
Tf = sort(Tf);
downburned = unburnedsamples2(1:Tf(2)/10:end,:);

ind = [1 2;1 3;1 4;1 5; 1 6;1 7;1 8;2 3;2 4; 2 5;2 6;2 7;2 8;3 4;...
    3 5;3 6;3 7;3 8;4 5;4 6;4 7;4 8;5 6;5 7;5 8;6 7;6 8;7 8];


subind = [3 4; 5 6; 7 8; 9 10; 11 12;13  14; 15 16;...
   21 22;23 24 ;25 26;27 28;29 30;31 32;...
   39 40;41 42;43 44;45 46;47 48;...
   57 58;59 60;61 62;63 64;...
   75 76;77 78;79 80;...
   93 94;95 96;...
   111 112]+16;
%subind = [2 3 4 5 6 7 8 11 12 13 14 15 16 20 21 22 23 24 29 30 ...
    %31 32 38 39 40 47 48 56];
addpath ../../Noor/GPS_D17pt_D246/bin_util/
load slipcolor

figure; 
for i = 1:size(ind,1)
    subplot(8,16,subind(i,:))
    p = gkde2([downburned(:,ind(i,2)) downburned(:,ind(i,1))]);
    contour(p.x,p.y,p.pdf,1000)
    ca = axis; 
    if ind(i,2) == 8 
        axis([ca(1) 2.5 ca(3) ca(4)])
    end
    if ind(i,2) == 5
        axis([8 ca(2) ca(3) ca(4)])
    end
    hold on ; 
    scatter(mbest(ind(i,2)),mbest(ind(i,1)),'m*','Linewidth',5)
end
colormap(slipcolor)

indhist = [1 2;3 4; 5 6; 7 8; 9 10 ; 11 12; 13 14; 15 16];
bins = 3000;
lab = ['Y west [km]','X west [km]','Y east [km]','X east [km]',...
    'width [km]','depth [km]','dip [\circ]','strike-slip [m]'];

labind = [1 11;12 22;23 33; 34 44;45 54;55 64;65 75;76 90];
sz=12;
for i=1:length(indhist)
    subplot(8,16,indhist(i,:))
    
    [nx,ny]= hist(unburnedsamples2(:,i),bins,'color','k');
    plot(ny,nx)
    set(gca,'YTickLabel',[]);
    axis([min(ny) max(ny) min(nx) max(nx)])
    title(lab(labind(i,1):labind(i,2)))
    
    line([mbest(i) mbest(i)],[min(nx) max(nx)],'color','r','Linewidth',2)
    prc_samp = prctile(downburned(:,i),[2.5 97.5]);
    line([prc_samp(1) prc_samp(1)],[min(nx) max(nx)],'color','b','Linewidth',2)
    line([prc_samp(2) prc_samp(2)],[min(nx) max(nx)],'color','b','Linewidth',2)
    
    if i ==5
        ca = axis;
        axis([7 ca(2) ca(3) ca(4)])
    end
    
    if i ==8
        ca = axis;
        axis([ca(1) 3.6 ca(3) ca(4)])
    end
    
    set(gca,'FontSize',sz);
    xl=get(gca,'xlabel');
    set(xl,'FontSize',sz);
    yl=get(gca,'ylabel');
    set(yl,'FontSize',sz);
    tl=get(gca,'title');
    set(tl,'FontSize',sz);
end

addpath ../../Noor/GPS_D17pt_D246/bin_util/
%subplot(8,8,42)
okada_allprior = change_okada(unburnedsamples); 
Mw = momentmag(okada_allprior);

Mwnls = Mw(249627);

subplot(8,16,[68 69]); hold on; 
[nx,ny]=hist(Mw,bins);
plot(ny,nx)
set(gca,'YTickLabel',[]);
xlabel('Moment magnitude')
ca = axis;

prc_mag = prctile(Mw,[2.5 97.5]); 
line([Mwnls Mwnls],[min(nx) max(nx)],'Color','r','Linewidth',2)
line([prc_mag(1) prc_mag(1)],[min(nx) max(nx)],'Color','b','Linewidth',2)
line([prc_mag(2) prc_mag(2)],[min(nx) max(nx)],'Color','b','Linewidth',2)

axis([6.52 6.63 min(nx) max(nx)])

sz =12;
set(gca,'FontSize',sz);
xl=get(gca,'xlabel');
set(xl,'FontSize',sz);
yl=get(gca,'ylabel');
set(yl,'FontSize',sz);
tl=get(gca,'title');
set(tl,'FontSize',sz);

subplot(8,16,[99:102 115:118]); hold on 
addpath ../figure8/
sampconf1 = [downburned(:,2) downburned(:,1)];
sampconf11(:,1) = sampconf1(:,1)- mean(sampconf1(:,1));
sampconf11(:,2) = sampconf1(:,2)- mean(sampconf1(:,2));
nls_conf1 = [mbest(1) - mean(sampconf1(:,2));mbest(2) - mean(sampconf1(:,1)) ];
plot(nls_conf1(2),nls_conf1(1),'b.','MarkerSize',30)

sampconf2 = [downburned(:,4) downburned(:,3)];
sampconf22(:,1) = sampconf2(:,1)- mean(sampconf2(:,1));
sampconf22(:,2) = sampconf2(:,2)- mean(sampconf2(:,2));
nls_conf2 = [mbest(3) - mean(sampconf2(:,2));mbest(4) - mean(sampconf2(:,1)) ];
plot(nls_conf2(2),nls_conf2(1),'r.','MarkerSize',30)

[fp1,exy1] = confellipse2(sampconf11,.95);
[fp2,exy2] = confellipse2(sampconf22,.95);

fp1 = fill(exy1(:,1),exy1(:,2),'b');
fp2 = fill(exy2(:,1),exy2(:,2),'r');
legend([fp1 fp2],'West endpoint','East endpoint')



