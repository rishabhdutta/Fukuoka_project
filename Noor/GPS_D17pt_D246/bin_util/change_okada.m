function [final] = change_okada(initial)
% function [final] = change_okada(initial)

nls1 = initial; 


   nls(:,1) =  sqrt((nls1(:,3)-nls1(:,1)).^2 + (nls1(:,4)-nls1(:,2)).^2);
   nls(:,2) = nls1(:,5);
   nls(:,3) = nls1(:,6);
   nls(:,4) = nls1(:,7);
  
   %if nls1(:,4)-nls1(:,2) <0
       nls(:,5) = 90 + atand(abs(nls1(:,4)-nls1(:,2))./(nls1(:,3)-nls1(:,1)));
   %elseif nls1(:,4)-nls1(:,2) >= 0
       %nls(:,5) = atand((nls1(:,3)-nls1(:,1))/(nls1(:,4)-nls1(:,2)));
   %end
   
   nls(:,6) = nls1(:,1)+ (nls1(:,3)-nls1(:,1))/2;
   nls(:,7) = nls1(:,2)+ (nls1(:,4)-nls1(:,2))/2;
   nls(:,8) = nls1(:,8);
   nls(:,9) = nls1(:,9);
   nls(:,10)= nls1(:,10);
   
    final = nls; 