function [OutData,Plane,m,X,Y] = RemovePlane(InData,opt,Xv,Yv);
%   RemovePlane    - Removes a plane from an intensity image 
%
%   Input:   InData  - (n,m) matrix
%            opt     - (1,1) 0 - a+by+cy (default)
%                            1 - a+by+cy+dxy+ex^2+fy^2
%            Xv      - (m,1) vector with x values
%            Yv      - (n,1) vector with y values           
%
%   Output:  OutData - (n,m) matrix with the best plane removed
%            Plane   - (n,m) matrix with the plane
%            m       - (3,1) or (6,1) vector of model parameters

[lin,col] = size(InData);


if nargin>4,help RemovePlane;return;end;
if nargin<3; 
 [X,Y] = meshgrid([1:col],[1:lin]);
else
 [X,Y] = meshgrid(Xv,Yv);
end
if nargin<2;opt=0;end

% number of pts:
N = lin*col;

OutData = zeros(lin,col);

% Check for NaNs:
d0 = InData(:);
ind = find(isnan(d0)~=1);
nnv = length(ind);
if nnv<3;disp('Too few non-nan points');return;end;

% Create kernel and data vector
G = [ones(nnv,1) X(ind) Y(ind)];
G2= [ones(N,1)   X(:)   Y(:)];
d = InData(ind);

if opt==1
  G = [G  X(ind).*Y(ind) X(ind).^2 Y(ind).^2];
  G2= [G2 X(:).*Y(:)     X(:).^2   Y(:).^2];
end

% find best bi-linear plane
m = inv(G'*G)*G'*d;

% Calc. plane:
Plane=zeros(size(InData));
Plane(:) = G2*m;

% remove plane from data
drem = d0 - Plane(:);

% fill output matrix
OutData(:) = drem;

