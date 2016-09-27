function FN = ExcludeGPS(FO,exclude);
% A simple function to exclude certain GPS stations from
% the GPS data structure
%
% Input:     FO (struct)   - GPS data structure
%            exclude (nx1) - vector of stations to exclude
%
% Output:    FN (struct)   - new reduced data structure

% SJ 2 Nov 2005

% Get number of data
No = size(FO.CGPSenu,2);

if length(exclude) >= No
    disp(' ');disp('---> You want to exclude all stations?');disp(' ');
    return;
end

if length(exclude) > 0
 % Get index vector, put NaN for stations to exclude
 io = [1:No];
 it = io;
 it(exclude) = NaN;

 % New index vector:
 in = it(find(isnan(it)==0));

 % Long new vector:
 inl = [in*3-2; in*3-1; in*3];
 inl = inl(:);

 % Extract new information:
 FN.CGPSnam = FO.CGPSnam(in,:);
 FN.CGPSll  = FO.CGPSll(:,in);
 FN.CGPSenu = FO.CGPSenu(:,in);
 FN.CGPSerr = FO.CGPSerr(:,in);
 FN.CGPSxy  = FO.CGPSxy(:,in);
 FN.CGPScor = FO.CGPScor(inl);
 FN.CGPScov = FO.CGPScov(inl,inl);
else
 FN = FO;
end