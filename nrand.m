function [a,v]=nrand(low,high)    %????????
a=low+rand(3,4)*(high-low);
s=a(:);
v=sum(s);
end
