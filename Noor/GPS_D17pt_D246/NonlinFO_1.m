% NON-LINEAR INVERSION OF the 2005 Fukuoka earthquake:
%clear all; close all; clc
% The procedure includes 4 scripts:
% NonlinSN_1:  Preparation of data and modeling
% NonlinSN_2:  Optimization of model parameters, storing
% NonlinSN_3:  Calculations of various things
% NonlinSN_4:  Plotting of various figures


% In this script we prepare the data and set all the parameters to run the
% Non-linear inversion to find the best parameters of the fault.  The
% inversion is performed in the next script 'NonlinSN_2' but here we do
% all the preperatory work:

% Script to find two optimal dislocation(s) for the 2005 SN event
% SJ. 17 Oct, 2005

amb = 0.004;  % Manual 5mm ambiguity:

% Exclude station vector, station 45 is on Genkaijima island
exclude = [45];
exclude = [];

%--------------------------------------------------------------------------
% LOAD IN DATA, and data covariance

% Load quadtree data:
load A
%load FD
load D17_init_GPS_coseis_stack
load SD
load intqtA
%load intqtFD
load intqtSD

% 
% % Load in SN image data
% load FOim
% FOim.unw = FOim.unw - amb;


% Load previous solution
load FO_Shinichi

% Load GPS data
load FO_GPS_geonet
load FO_GPS_TN
%ii=find(FO.CGPSxy(1,:)>670);       % Here I find far-eastern stat
%mn = (mean(FO.CGPSenu(:,ii)'))';   % Calculate mean displ
% and remove the mean:
mn = [0 -0.9 -5.1]'*1e-3; % Shift vector
FO.CGPSenu = FO.CGPSenu-repmat(mn,1,size(FO.CGPSenu,2));
FN=ExcludeGPS(FO,exclude);
FN.CGPSll  = [FN.CGPSll FOtn.GPSll];
FN.CGPSenu = [FN.CGPSenu FOtn.GPSenu];
FN.CGPSerr = [FN.CGPSerr FOtn.GPSerr];
FN.CGPScov = BlockDiag(FN.CGPScov,FOtn.GPScov);
FN.CGPSxy  = [FN.CGPSxy FOtn.GPSxy];
clear FO FOtn
    
%---------------------------------------------
% Unit vector from Mike Poland
UD1 = intqtA.los;
UD2 = intqtSD.los;


%------------------------------------------------------------------
% PREPARE DATA AND COORDINATES

% Prepare datavector.  The first entries are GPSenu, then Insar, then Insar
% south of the faults, these are three sets:
denu    = FN.CGPSenu;
derr    = FN.CGPSerr;
dgps    = denu(:);

dsar1    = intqtA.sv;% - 0.0792;
dsar2 =  D17.def';                      %intqtFD.sv;
dsar3 = intqtSD.sv;

d       = [dgps; dsar1'; dsar2'; dsar3'];
%d = dsar;
%disp('---> introduce shift in INSAR data...');

% prepare coordinates
%tmp     = SN.CGPSll-repmat([SNdif.corner_lon;SNdif.corner_lat],1,5);
crdenu  = FN.CGPSxy;
tmp3    = repmat(crdenu,3,1);
crdgps  = zeros(2,size(crdenu,2)*3);
crdgps(:)= tmp3(:);

crdsar1  = intqtA.cnt;
crdsar2 = [D17.pos.E'; D17.pos.N']/1e3;                      %intqtFD.cnt;
crdsar3 = intqtSD.cnt;

coord   = [crdgps crdsar1 crdsar2 crdsar3];
%coord   = crdsar;

% prepare data numbers:
N       = max(size(d));         % Total number of data
Ngps    = max(size(dgps));      % Number of GPS data (3xstations)
%Nres    = SNqtp.dnum(1);     % Number of Insar data north of faults
%Npat    = SNqtp.dnum(2);     % Number of Insar data on patch
%Nsar    = Nres + Npat;       % Number of all InSAR data
Nsar1    = max(size(dsar1));
Nsar2  = max(size(dsar2));
Nsar3 = max(size(dsar3));

% set up data-index vector
%dnum = Nsar; datind=dnum(1);
dnum    = [Ngps; Nsar1; Nsar2; Nsar3];
datind  = [dnum(1) sum(dnum(1:2)) sum(dnum(1:3)) sum(dnum(1:4))]'; 

%-------------------------------------------------------------
% DATA COVARIANCE MATRIX

disp('**** Warning: using simplified covariance matrix');

% Generate Data covariance matrix, InSAR variance 1e-4
% DataCov=eye(length(d))*(0.02).^2;
% DataCov([1:Ngps],[1:Ngps])=eye(Ngps).*repmat(SN.CGPSerr(:).^2,1,Ngps);

% GPS station on island is problematic, increase its sigma:
ip=45; % number of station
ii=[3*ip-2:3*ip];

SARSTD     = 0.01;
DataCovGPS = FN.CGPScov;
DataCovGPS(ii,ii)=2^2.*DataCovGPS(ii,ii);

DataCovSAR1 = eye(Nsar1)*10000000; 
DataCovSAR2 = zeros(Nsar2);
DataCovSAR3 = zeros(Nsar3); 

paramSAR2 = [67.1171e-6 2.4468 486.9097 647.0492];
varSAR2 = 3.9336e-5;
paramSAR3 = [6.821e-6 30.75 8.341 .721];
varSAR3 = 2.9401e-5;
% 
% for i = 1: numel(DataCovSAR2) 
%     if mod(i,Nsar2) == 0
%         row = Nsar2; 
%     else
%         row  = mod(i,Nsar2); 
%     end
%     
%     col = ceil(i/Nsar2); 
%     if row == col 
%         covSAR2(i) = varSAR2; 
%     else
%         co1  = crdsar2(:,row); co2 = crdsar2(:,col);
%         dist = sqrt((co1(1) - co2(1))^2 + (co1(2) - co2(2))^2);
%         covSAR2(i) = limiter(paramSAR2,dist);
%     end
% end
    DataCovSAR2 = sqrt(eye(Nsar2).*(D17.var*D17.var'));%reshape(covSAR2,Nsar2,Nsar2);
    
 
for i = 1: numel(DataCovSAR3) 
    if mod(i,Nsar3) == 0
        row = Nsar3; 
    else
        row  = mod(i,Nsar3); 
    end
    
    col = ceil(i/Nsar3); 
    if row == col 
        covSAR3(i) = varSAR3; 
    else
        co1  = crdsar3(:,row); co2 = crdsar3(:,col);
        dist = sqrt((co1(1) - co2(1))^2 + (co1(2) - co2(2))^2);
        covSAR3(i) = limiter(paramSAR3,dist);
    end
end
    DataCovSAR3 = reshape(covSAR3,Nsar3,Nsar3);   
%-------------------------------------------------------------
%edited by Rishabh (values taken from the variogram)
% varcov=zeros(62,62);           
% mn=mean(FOqt.sv(:));
% 
% x=zeros(1,62);
% y=zeros(1,62);
% 
% qt=FOqt.sv(:);
% varqt=var(qt);
% 
% x(1,:)=(FOqt.cx(1,:)+FOqt.cx(2,:))/2;
% y(1,:)=(FOqt.cy(3,:)+FOqt.cy(2,:))/2;
% 
% for i=1:62
%     for j=1:62
%         if j==i
%             varcov(i,j)=10e-5;                        %2*(FOqt.sv(i,1)-mn).^2;
%         else
%             %use of the obtained covariance function liable to change 
%             % C=1.4e-5*exp(-h/320)*cos((h+50)/120) where h is in 0.1km
%             well=sqrt((x(1,i)-x(1,j)).^2+(y(1,i)-y(1,j)).^2)*10;  
%             varcov(i,j)=1.4e-5*exp(-well./320)*cos((well+50)/120);
%         end
%     end
%     end

% edit by Rishabh ends 
%---------------------------------------------------------------



%DataCovSAR1 = eye(Nsar1) * (SARSTD).^2 ; %(earlier) DataCovSAR= eye(Nsar) * (SARSTD).^2;
%DataCovSAR2 = eye(Nsar2) * (SARSTD).^2 ;
%DataCovSAR3 = eye(Nsar3) * (SARSTD).^2 ;
%DataCov    = DataCovSAR;
DataCov    = BlockDiag(DataCovGPS, DataCovSAR1, DataCovSAR2, DataCovSAR3);

for i = 1:N 
    DataCov(i,i) = 8*DataCov(i,i);
end

disp(['---> Using InSAR std of : ',num2str(1e2*SARSTD),' cm']);
% Calculate inverso and Cholesky factorization:
COVinv     = inv(DataCov);
W          = chol(COVinv);
%exclude points 
excl = [ 231   232   236   243   244];    %removed1
 %excl = [243   244   272   277   281   282 ...
  %   290   291   293   298   299   300   301]; % removed2 (more than 4cm )
for i = 1:length(excl)
    W(1381+excl(i),1381+excl(i))=1e-3;
end

W2         = W'*W;
% 
% % Display the covariance matrix:
% load slipcolor
% figure;imagesc(DataCov);axis image
% colorbar; %colormap(slipcolor);


% ---------------------------------------------------------------------------------
% SET VARIABLES

% Poisson ratio, nu=0.31 is the average of all the undrained values for various
% materials in the book by Wang.  nu=0.20 is the average for drained state.
% 
% We use 0.28 as it is close to drained nu, because co-seismic data includes
% postseismic signal as well, or most of it (thus not use 0.27 (fully dr.), 
% nor 0.31 (undrained)).  For Sierra Negra we just use the standard 0.25

nu = 0.25;   

disp(['---> Using Poisson`s ratio of ',num2str(nu)]);



% %----------------------------------------------------------------------------------
% % PLOT DATA
% 
% disp('---> Plotting data...')
% figure;
% patch(FOqt.cx,FOqt.cy,dsar'); axis image; hold on  
% DM_Quiver(crdenu,denu(:),DataCovGPS,0.2,[0 0],'k');
% plot(FOim.clx(1,:),FOim.clx(2,:),'-k');
% title('Descending data'); xlabel('Easting [km]');
% ylabel('Northing [km]'); colorbar
% 
% % Put stars in the patch data, to be certain it is correct:
% %ii = sum(dnum(1:2))+1;
% plot(coord(1,1:dnum(1)),coord(2,1:dnum(1)),'ks');
% plot(coord(1,dnum(1)+1:end),coord(2,dnum(1)+1:end),'m*');
% 
% clx=FOim.clx;
% ii=find(clx(1,:)<670 & clx(1,:)>600 & clx(2,:)>3680 & clx(2,:)<3755 | isnan(clx(1,:))==1);
% clx=clx(:,ii);
% plot(clx(1,:),clx(2,:),'-k');



%--------------------------------------------------------
% PREPARE BOUNDS ON Sierra Negra fault model
%
% Make lower and upper bounds for possible dislocation parameters
%
% Length - Width - Depth - Dip - Strike - East - North - Sslip - Dslip - Op
% (km)     (km)    (km)    (deg) (deg)    (km)   (km)    (m)     (m)     (m)
%                  upper   neg-  
%                  edge    ative
%
% Tight bounds on strike-slip and opending = 0
%       Length - Width - Depth - Dip - Strike - East - North - Sslip - Dslip - Op


lb1 = [   590   3725    605     3715    5   0   -100    -3.0     0     0   ]';
ub1 = [   610   3745    630     3735    15  5   -80     -0.5     0     0   ]';

% Ambiguity for InSAR data north and south of the SN intra-caldera fault:
%lb4 = [  0.0   0.3 ]';
%ub4 = [  0.0   0.9 ]';
%lb4 = [  0.3 ]';
%ub4 = [  0.9 ]';

% Create bounds matrix
%bounds = [lb1 ub1; lb4 ub4];
bounds = [lb1 ub1];

% Parameter index, only used if want to keep two numbers the same 
% (e.g. center location of two faults the same), not used here.
parind = [1:10];

% Display bounding values on screen:
disp('---> Here are the bounding values in the non-linear optimization:'); disp(' ')
disp('X-loc-left  Y-loc-left   X-loc-right   Y-loc-right    width    Depth     Dip');
disp(num2str(bounds(1:7,:)'));disp(' ')
disp('Str-slip   Dip-slip   Opening   Ambiguity   Patch-Ambiguity');
disp(num2str(bounds(8:end,:)'));


%---------------------------------------------------------------
% OPTIMIZATION OPTIONS

% Options for Simulated annealing, see "help anneal"
annealopt = [4 3 4 2.5 0 0 1];

% Objective function that calculates simulated data from
% given model parameters, and calculates the residual
objfun='OneFault_obj_end';

% Find fully constrained parameters from specified bounds
parflag = bounds(:,1)-bounds(:,2);
fixind  = find(parflag==0);
fixpar  = bounds(fixind,1);
Npar    = size(max(parind))-size(max(fixind));

%----------------------------------------------------------------

disp(' ');
