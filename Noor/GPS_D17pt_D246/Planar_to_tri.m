function [tri1,tri2,X,Y,Z] = Planar_to_tri(m,plot)
% [tri1,tri2,X,Y,Z] = Planar_to_tri(m)
% convert the okada parameters to triangular dislocation parameters 
% m is the okada parameter
% written by Rishabh Dutta

m = reshape(m,1,length(m));
sd = sin(m(4)*pi/180);
if sd == 0
    wtop = 0;
    wbot = m(2);
else
    wtop = (m(3)/sd) - m(2);
    wbot = m(3)/sd;
end
verts = [[wtop; m(1)/2; 0], ...
    [wbot; m(1)/2; 0], ...
    [wbot; -m(1)/2; 0], ...
    [wtop; -m(1)/2; 0]];
botmpt = [wbot; 0; 0];

thta = m(4)*pi/180;
phi = m(5)*pi/180;
R1 = [cos(thta) 0 sin(thta); 0 1 0; -sin(thta) 0 cos(thta)];
R2 = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
verts = R2 * R1 * verts;   
botmpt = R2 * R1 * botmpt;

verts = verts - repmat([botmpt(1:2); 0], 1, 4);    
verts = verts + repmat([m(6:7)'; 0], 1, 4);      

X = verts(1,:)';
Y = verts(2,:)';
Z = verts(3,:)'; 

if plot ==1
    plot3([X; X(1)],[Y; Y(1)],[Z; Z(1)],'r','Linewidth',6)
end
% triangle 1 ; 
X1 = X(1:3) ; Y1 = Y(1:3) ; Z1 = Z(1:3) ; 
%triangle 2 ; 
X2 = [X(3:4); X(1)] ; Y2 = [Y(3:4); Y(1)] ; Z2 = [Z(3:4); Z(1)] ; 
tri1 = [X1 Y1 Z1]; 
tri2 = [X2 Y2 Z2]; 







