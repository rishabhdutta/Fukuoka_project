function output = limiter(param,input)
a=param(1);
b=param(2);
c=param(3);
e=param(4);
i=1;
l=max(size(input));
while (i>=1 & i<=l)
    d(i)=cos((input(i)+e)/c);
    i=i+1;
end
output=a*d.*exp(-input./b);
