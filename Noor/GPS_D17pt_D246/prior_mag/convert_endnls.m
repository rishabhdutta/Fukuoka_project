function [model_nls] = convert_endnls(model_end)

nls1 = model_end;

   nls(:,1) =  sqrt((nls1(:,3)-nls1(:,1)).^2 + (nls1(:,4)-nls1(:,2)).^2);
   nls(:,2) = nls1(:,5);
   nls(:,3) = nls1(:,6);
   nls(:,4) = nls1(:,7);
  
   
       nls(:,5) = 90 + atand(abs(nls1(:,4)-nls1(:,2))./(nls1(:,3)-nls1(:,1)));
   
       %nls(:,5) = atand((nls1(:,3)-nls1(:,1))/(nls1(:,4)-nls1(:,2)));
  
   
   nls(:,6) = nls1(:,1)+ (nls1(:,3)-nls1(:,1))/2;
   nls(:,7) = nls1(:,2)+ (nls1(:,4)-nls1(:,2))/2;
   nls(:,8) = nls1(:,8);
   nls(:,9) = nls1(:,9);
   nls(:,10)= nls1(:,10);
   
   model_nls = nls; 