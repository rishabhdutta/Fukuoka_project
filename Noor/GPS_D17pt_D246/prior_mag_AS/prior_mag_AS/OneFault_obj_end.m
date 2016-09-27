function [resiW,dpred,d_new,resi] = OneFault_obj_end(par,data,coord,unit1,unit2,nu,W,fixind,fixpar,sqrflag,datind,parind)
%   TwoFaults1721_obj     - Two faulst objective func. for june 17 and 21 eqs. for LEASTSQ and ANNEAL
%
% function [resiW,dpred,d_new,resi] = TwoFaults1721_obj(par,data,coord,unit1,unit2,nu,W2,fixind,fixpar,sqrflag,datind,parind)
%
% INPUT:
%		par	    - (24x1) vector of modelparameters where left and right
%		fault end parameters are input
%		data	- (nx1) vector of observations
%		coord	- (2xn) matrix of coordinates
%       unit1   - look vector for june 17
%       unit2   - look unit vector for june 21
%		nu	    - poisson's ratio
%		Sig     - Normalization vector, vector with variance values (like a diagonal of a COV matrix) 
%		fixind	- indices of fixed model parameters
%		fixpar	- the value of the fixed parameters
%		sqrflag - 1 - the norm of the residual vector is returned
%			      0 - the residual vector is returned (for LEASTSQ).
%       datind  - data index vector to separate different datasets
%       parind  - parameter index vector, if we want to keep fault parameters the same
%                 for the two faults
%
% OUTPUT:	
%       resiW	- either (nx1) residual vector between the input
%			      model and the data, or the norm of the vector
%       dpred   - predicted data vector
%       d_new   - corrected data (plane or amb. removed)
%       resi    - residual between the data and the model prediction


% Find total number of data points, plus number of 
% data points for individual data set
    N     = length(data);
	Ngps  = datind(1);
	N1    = datind(2)-datind(1);
	N2    = datind(3)-datind(2);
    N3    = datind(4)-datind(3);
	
% Put in the fixed parameters:
par(fixind)=fixpar;

% Ambiguities
amb1=0;
amb2=0;

%	amb1 = par(parind(11));
%	amb2 = par(parind(12));
    
%	desdx = par(parind(13));
%	desdy = par(parind(14));
	
% Plane for june 21 interferogram
%	abc21 = par(parind(22:24));
%	crd21 = coord(:,N17+1:N17+N21);	
%	plane = [ones(N21,1) crd21']*abc21;
	  
% Pick out June 17 igram coordinates, and GPS coordinates
   crdgps = coord(:,1:Ngps);
   tmp    = zeros(6,size(crdgps,2)/3);
   tmp(:) = crdgps(:);
   
   crdgps3=tmp(1:2,:);
   
   crd1   = coord(:,Ngps+1:Ngps+N1);
   crd2   = coord(:,Ngps+N1+1:Ngps+N1+N2);
   crd3   = coord(:,Ngps+N1+N2+1:end);
   
% Reconstruct the third fault, using parameters of the second fault
   %%%%define the par for disloc input
   
   par1(1) =  sqrt((par(3)-par(1)).^2 + (par(4)-par(2)).^2);
   par1(2) = par(5);
   par1(3) = par(6);
   par1(4) = par(7);
  
   if par(4)-par(2) <0
       par1(5) = 90 + atand(abs(par(4)-par(2))/(par(3)-par(1)));
   elseif par(4)-par(2) >= 0
       par1(5) = atand((par(3)-par(1))/(par(4)-par(2)));
   end
   
   par1(6) = par(1)+ (par(3)-par(1))/2;
   par1(7) = par(2)+ (par(4)-par(2))/2;
   par1(8) = par(8);
   par1(9) = par(9);
   par1(10)= par(10);
   
   par1 = par1';
   
   %par1(6)=par1(6)-0.5*par1(1)*sind(par1(5));  % parameter location shifted from end to centre
   %par1(7)=par1(7)-0.5*par1(1)*cosd(par1(5));
   
%   par2  = par(parind(1:10));
%   par2(6) = par2(6) - desdx;
%   par2(7) = par2(7) - desdy;
   

% Calculate model:   
   [ugps] = disloc(par1, crdgps3, nu);
   [u1]=disloc(par1, crd1, nu); 
   [u2]=disloc(par1, crd2, nu);
   [u3]=disloc(par1, crd3, nu);
   
% Add zeros for the June 21 igram, that didn't see u1
 %  u1b = [u1(:,1:N17) zeros(3,N21) u1(:,N17+1:N17+Ngps) ];
   
% Add up displacements from the two faults
  % u = u1;

% Project to radar range change, first for june 17 interferogram,
% and then for the june 21 igram
ugpsl   = ugps(:);
rngcng1 = unit1'*u1;
rngcng2 = unit2'*u2;
rngcng3 = unit2'*u3;
%rngcng1 = sum(u1.*unit1);
%rngcng2 = sum(u2.*unit2);

%	rngcng1 = sum(unit1'*u);
%	rngcng2 = sum(unit2.*u(:,datind(1)+1:datind(2)));
%ugps    =        u(:,datind(2)+1:datind(3));
%	ugps    = ugps(:);
	
    dpred   = [ugpsl;rngcng1';rngcng2';rngcng3'];%;rngcng2'];
	
% Take out individual data vectors
    dgps   = data(          1:datind(1));
    d1     = data(datind(1)+1:datind(2));
	d2     = data(datind(2)+1:datind(3));
    d3     = data(datind(3)+1:datind(4));
	
% Correct data for offset and/or plane
	d1_new = d1-amb1;
	d2_new = d2-amb2;
    d3_new = d3-amb2;
	
% Create new full and corrected data vector
	%d_new   = [dgps; d1_new; d2_new];
    d_new   = [dgps; d1_new; d2_new; d3_new];

% Calculate "weighted" residual between data and predicted data
	resi    = ( d_new - dpred );
 %    resiW   = resi ./ sqrt(Sig);
    resiW  = W*resi;
	
% Residual vector or sum or squares?
	if sqrflag==1
		resiW=resiW'*resiW;
	end
