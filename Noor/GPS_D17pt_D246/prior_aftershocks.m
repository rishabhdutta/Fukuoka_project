function final = prior_aftershocks(trian,Lat,Lon,depth)

% function final = prior_aftershocks(trian,Lat,Lon,depth)
% output: 
% input : 
%   trian = 3 points on the plane [x1 y1 z1; x2 y2 z2; x3 y3 z3]
%   Lat = latitude aftershocks
%   Lon = longitude aftershocks
%   depth = depth aftershocks

% define plane from the 3 points
P1 = trian(1,:); P2 = trian(2,:); P3 = trian(3,:);

normal = cross(P1-P2, P1-P3);
syms x y z
P = [x,y,z];
planefunction = dot(normal, P-P1);

realdot = @(u, v) u*transpose(v);
planefunction = realdot(P-P1,normal);
planefunc = matlabFunction(planefunction);

D = planefunc(0,0,0); 
A = planefunc(1,0,0) - D ; 
B = planefunc(0,1,0) - D ;
C = planefunc(0,0,1) - D ; 

noEQ = length(Lat); 

for i = 1:noEQ 
    xpo = Lon(i); ypo= Lat(i); zpo = -depth(i);
    dist(i) = abs(A*xpo+B*ypo+C*zpo+D)/sqrt(A^2+B^2+C^2);
end

final = sqrt(sum(dist.^2));


    
    
