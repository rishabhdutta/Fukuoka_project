function [Cpoints] = SurfProj(dis_geom,plotopt,im,col)
%   SurfProj       - Projects dislocation outlines onto surface
%
% Input:
% dis_geom = [len,width,depth,dip,strik,locE,locN,0,0,1];
% plotopt  - 1 plot, 0 no plotting (default is 1)
% im       - Image plot?  1 yes, 0 no (default).
% col      - (string), for color (default black)
%
% Output:
% k        - (2x5) matrix to plot a rectangle with " plot(k(1,:),k(2,:)) "

if nargin < 2; plotopt=1; im=0; col='k'; end
if nargin ==2; im=0; col='k'; end
if nargin ==3; col='k'; end
if nargin ==0 | nargin > 4; help SurfProj; end

Cpoints=[];
for kk=1:size(dis_geom,2)
  
len = dis_geom(1,kk); 
wid = dis_geom(2,kk); 
dep = dis_geom(3,kk);
dip = (dis_geom(4,kk))/180*pi;
str = dis_geom(5,kk); 
dx  = dis_geom(6,kk); 
dy  = dis_geom(7,kk);

p = [-len/2  len/2        	len/2       -len/2;
	  0      0 	 wid*cos(dip)  wid*cos(dip)];

a = (360-str+90)/180*pi;

R = [cos(a) -sin(a); sin(a) cos(a)];
k = R*p;
k = k + [dx;dy]*ones(1,4);
k = [k k(:,1)];

Cpoints = [Cpoints k ones(2,1)*NaN];

% If image plot we reverse the y-axis for plotting
if im==1; k(2,:)=-k(2,:); end

if plotopt==1
   hold on
   lw=get(gca,'DefaultLineLineWidth');
   plot(k(1,:),k(2,:),'Color',col)
       set(gca,'DefaultLineLineWidth',lw*5);

   if dip<0
    plot(k(1,1:2),k(2,1:2),'Color',col)
else
    plot(k(1,3:4),k(2,3:4),'Color',col)
end
   set(gca,'DefaultLineLineWidth',lw);
end

end
