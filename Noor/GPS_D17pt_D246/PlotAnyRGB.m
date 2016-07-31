function [rgbim,shad,clim] = PlotAnyRGB(colim,dem,colm,ampsc,shadsat,colsat,shadopt,file);
%    PlotTopoRGB    - Creates a color-topo-map with shading.\
%
% usage:  [rgbim,shad,clim] = PlotAnyRGB(colim,dem,colm,ampsc,shadsat,colsat,file);
%
% Input:
% colim      - (n,m) color image
% dem        - (n,m) dem matrix, for shading, can be higher res
% colm       - (63,3) colormap matrix, default is 'jet'
% ampsc      - (1,1) to scale the brightness (default=1), try 1.5!
% shadsat    - (2,1) vector with shading saturation values,
%                    default: max,min in shading (e.g. try [-50 50])
% colsat     - (2,1) vector with color saturation values,
%                    default: max,min of colorimage 
% shadopt    - (2,1) if (1,1)=1 then shad-calc, else file used directly, default=1.
%                    if (2,1)=1 the top-right shading, 2, top left (default=1);
% file       - (string) output jpeg file, default: no file
%
% Output:
% rgbim      - (n,m,3) rgb image of the topography
% shad       - (n,m) shading image of the topo, range [0 1]
% clim       - (2,1) vector with minimum and maximum colorvalues
%
% SJ, 7 July 2004

% Set default values, depending on number of input arguments
if nargin<3; colm=jet;  end
if nargin<4; ampsc=1;   end
if nargin<5; shadsat=1; end
if nargin<6; colsat=1;  end
if nargin<7; shadopt=[1,1]'; end

if length(shadopt)==1; shadopt=[shadopt(1); 1];end

%------------------------------------------------------------
% Create shading matrix:

% find number of lines and column in DEM-matrix
[ll,cc] = size(dem);
[ll2,cc2] = size(colim);

if ll~=ll2
  colim=imresize(colim,[ll cc],'nearest');
end

% create shading convolution matrix, sun shining from top-right:
if shadopt(2)==1
  con = [0 -1; 1 0];
elseif shadopt(2)==2
  con = [-1 0; 0 1];
elseif shadopt(2)==3
  con = [0 1; -1 0];
else
  con = [1 0; 0 -1];
end

if shadopt(1) == 1
  % Cross-correlate the two matrices:
  out = xcorr2(dem,con);

  % fix boundaries after convolution (this is pretty lame...)
  % convolution output is (ll+1,cc+1) in size, good output is
  % (ll-1,cc-1) because of wrapping:
  out = out(2:ll,2:cc);       % good part picked out
  out = [out; out(ll-1,:)];   % add a replicaof last line at bottom
  out = [out(:,1) out];       % add a replica of first column in front

  % Shading determines the amplitude, and is reported as output
  shad = out;
else
  shad = dem;
end

% Saturate shading, avoid extreme values
if shadsat ~=1
  shad(shad<shadsat(1))=shadsat(1);
  shad(shad>shadsat(2))=shadsat(2);

  % adjust shading file, so it only has values from 0-1
  amp = shad - shadsat(1);
  amp = amp/(shadsat(2)-shadsat(1));
else
  % adjust shading file, so it only has values from 0-1
  amp = shad - min(shad(:));
  amp = amp/max(amp(:));
end

% scale the amplitudes, for brighter image
if ampsc ~= 1
  amp = amp*ampsc;
  amp(amp>1) = 1;
end

%----------------------------------------------------------------------
% Prepare color image data

% Saturate color values, or find max,min of colorimage
if length(colsat)>1
  minc = colsat(1);
  maxc = colsat(2);
  colim(colim<minc)=minc;
  colim(colim>maxc)=maxc;
else
  minc = min(colim(:));
  maxc = max(colim(:));
end

% Report min,max elevation values
clim = [minc maxc]';

% adjust dem values, values between 0 and 63, as
% we'll keep top-bin empty for coloring the NaNs gray
d = colim-minc;
%d = ( d/max(d(:)) )*63;
d  = ( d/(maxc-minc) )*63;
d = round(d);            % has to be integers

% Here we modify the colortable, as the NaNs are colored
% by the top value, put gray for that
kjet = [colm(1:63,:); 0.9 0.9 0.9];

% Transform DEM image to rgb image using 'colm' colormap
rgbim = ind2rgb(d,kjet);

% Scale with amplitudes of RGB image using the scaled shading values
for k=1:3
  out  = rgbim(:,:,k);   % re-use 'out'
  out  = out.*amp;
  rgbim(:,:,k) = out;
end


% display image

%imagesc(rgbim);


% save image as jpeg:
if nargin == 8
  imwrite(rgbim,file,'JPEG','Quality',90);
end
