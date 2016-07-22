function [vario,covar,lagi] = variogram(out,maxlag,x_inc,y_inc)
% function [variogram,covariogram] = variogram(out,x_inc,y_inc)

if nargin < 3
    x_inc = 1; 
    y_inc = x_inc; 
elseif nargin < 4
    x_inc = abs(x_inc);
    y_inc = x_inc;
end
x_inc = abs(x_inc); y_inc = abs(x_inc);


rows = size(out,1);    % y direction 
cols = size(out,2);    % x direction 

lagi = linspace(0,maxlag,200);
covar = zeros(length(lagi),1); vario = zeros(length(lagi),1);
for i = 1:length(lagi)
    lag = lagi(i); var = zeros(5000,1); cov = zeros(5000,1); 
    for j=1:5000
        k =1;
        while k == 1 
            row1 = ceil(rand*rows); col1 = ceil(rand*cols); 
            if isnan(out(row1,col1)) == 1
                continue
            else
                theta = rand * 360; 
                row2 = round(row1 + lag*sind(theta)/y_inc) ; 
                col2 = round(col1 + lag*cosd(theta)/x_inc) ; 
                if row2>rows || row2<1 || col2>cols || col2<1 || isnan(out(row2,col2)) ==1
                    continue
                else
                    k=k+1; 
                end
            end
        cov(j) = out(row2,col2)*out(row1,col1);
        var(j) = (out(row2,col2) - out(row1,col1))^2;
        end
    end
    covar(i) = sum(cov(:))/10000;
    vario(i) = sum(var(:))/10000;
end

figure; plot(lagi,covar,'r.') 
figure; plot(lagi,vario,'r.') 

    
                
                
        
        
        
end

