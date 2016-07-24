% plot the 2D kernel densities with the marginals 

addpath ~/Desktop/softwares/gkde2/
addpath ~/Desktop/softwares/mcmcdiag/
load samples_med_again2

Tf = geyer_imse(unburnedsamples);
Tf = sort(Tf);
downburned = unburnedsamples(1:Tf(2)/10:end,:);

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
addpath ../../../bin_util/
load slipcolor

figure; 
for i = 1:size(ind,1)
    subplot(8,16,subind(i,:))
    p = gkde2([downburned(:,ind(i,2)) downburned(:,ind(i,1))]);
    contour(p.x,p.y,p.pdf,1000)
    
end
colormap(slipcolor)

indhist = [1 2;3 4; 5 6; 7 8; 9 10 ; 11 12; 13 14; 15 16];
bins = 500;
lab = ['X west [km]','Y west [km]','X east [km]','Y east [km]',...
    'width [km]','depth [km]','Dip [\circ]','strike-slip [m]'];

labind = [1 11;12 22;23 33; 34 44;45 54;55 64;65 75;76 90];
sz=10;
for i=1:length(indhist)
    subplot(8,16,indhist(i,:))
    
    [nx,ny]= hist(downburned(:,i),bins);
    plot(ny,nx)
    set(gca,'YTickLabel',[]);
    axis([min(ny) max(ny) min(nx) max(nx)])
    title(lab(labind(i,1):labind(i,2)))
    
    set(gca,'FontSize',sz);
xl=get(gca,'xlabel');
set(xl,'FontSize',sz);
yl=get(gca,'ylabel');
set(yl,'FontSize',sz);
tl=get(gca,'title');
set(tl,'FontSize',sz);
end

%subplot(8,8,42)

   nls(:,1) =  sqrt((unburnedsamples(:,3)-unburnedsamples(:,1)).^2 + (unburnedsamples(:,4)-unburnedsamples(:,2)).^2);
   nls(:,2) = unburnedsamples(:,5);
   nls(:,3) = unburnedsamples(:,6);
   nls(:,4) = unburnedsamples(:,7);
  
   
   nls(:,5) = 90 + atand(abs(unburnedsamples(:,4)-unburnedsamples(:,2))./(unburnedsamples(:,3)-unburnedsamples(:,1)));
       
   %    nls(5) = atand((unburnedsamples(3)-unburnedsamples(1))/(unburnedsamples(4)-unburnedsamples(2)));
   
   
   nls(:,6) = unburnedsamples(:,1)+ (unburnedsamples(:,3)-unburnedsamples(:,1))/2;
   nls(:,7) = unburnedsamples(:,2)+ (unburnedsamples(:,4)-unburnedsamples(:,2))/2;
   nls(:,8) = unburnedsamples(:,8);
   nls(:,9) = unburnedsamples(:,9);
   nls(:,10)= unburnedsamples(:,10);

unburnedsamples = nls;
moment=2.5*unburnedsamples(:,1).*unburnedsamples(:,2).*abs(unburnedsamples(:,8)).*1e23;
Mw=2/3*log10(moment)-10.7;
m =1e3*[0.0116 0.0119 0.0016 -0.0879 0.1196 0.6082 3.7326 -0.0021 0 0];
momnls= 2.5 *m(1)*m(2)*abs(m(8))*1e23;
Mwnls = 2/3* log10(momnls) - 10.7;

subplot(8,16,[83 84])
[nx,ny]=hist(Mw,bins);
plot(ny,nx)
set(gca,'YTickLabel',[]);
xlabel('Moment magnitude')

sz =14;
set(gca,'FontSize',sz);
xl=get(gca,'xlabel');
set(xl,'FontSize',sz);
yl=get(gca,'ylabel');
set(yl,'FontSize',sz);
tl=get(gca,'title');
set(tl,'FontSize',sz);



