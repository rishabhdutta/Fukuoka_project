
function out = momentmag(in); 
% function out = momentmag(in);  
% in - okada parameters

nls = in;  

rig = 2.5e10; % shear modulus
slipmag = sqrt( nls(:,8).^2 + nls(:,9).^2 ); % slip magnitude
  
gmom = nls(:,1).* nls(:,2).* 1e6 .* slipmag;
  
moment = rig*gmom;
  
Mw = 2*log10(moment)/3 - 6.03;
out = Mw;
  
  