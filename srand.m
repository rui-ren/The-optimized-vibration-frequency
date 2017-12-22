function[a,s]=srand(low,high)
a=low+rand(3,4)*(high-low);
s=sumallelement(a);

function suma=sumallelement(M)
v=M(:);
suma=sum(v);

