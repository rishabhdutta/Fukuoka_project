function crdout = FO_CrdTrans(crdin,fwd);
% function crdout = FO_CrdTrans(crdin,fwd);
%
% Crd transformation function for Fukuoka, to 
% change between UTM and lat-lon.  Note, we have
% the UTM coords in km.
%
% Input:       crdin  - (2,N) matrix of N coordinates.
%              fwd    - 1 for Lat-Lon to UTM-km (default)
%                       0 for UTM-km to lat-lon
% Output       crdout - (2,N) of new coords.
%

if nargin == 1;
    fwd = 1;
end

if nargin > 2
    help FO_CrdTrans;
    return;
end

if size(crdin,1)~=2;
    disp(' ');disp('---> ERROR: Store coordinates columnwise');disp(' ');
    return;
end

% Define UTM
axesm utm

% Fukuoka UTM zone:
zone = '52S';

axesm('mapprojection','utm','zone',zone);

if fwd == 1
  % Transform coordinates to UTM
  % [x,y] = mfwdtran(lat,lon)
  [xx,yy] = mfwdtran(crdin(2,:),crdin(1,:));

  % Project to km
  xx = xx/1e3; 
  yy = yy/1e3;
  crdout = [xx;yy];
else
  % Transform from UTM to LatLon:
  % first change from km to meters:
  x2 = crdin(1,:) * 1e3;
  y2 = crdin(2,:) * 1e3;
  
  % [lat,lon] = minvtran(x,y):
  [lat,lon]=minvtran(x2,y2);
  
  crdout = [lon; lat];
end
