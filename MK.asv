syms x l
a=(x)/l;
N1=1-3*a^2+2*a^3;
N2=(a-2*a^2+a^3)*l;
N3=3*a^2-2*a^3;
N4=(-a^2+a^3)*l;
N=[N1 N2 N3 N4];
b=N'*N;
c=int(b,0,l);
d=c/420;
e=diff(N,x,2);
f=e'*e;
g=int(f,0,l);

%M
syms x p
a=(x)/p;
N1=1-3*a^2+2*a^3;
N2=(a-2*a^2+a^3)*p;
N3=3*a^2-2*a^3;
N4=(-a^2+a^3)*p;
N=[N1 N2 N3 N4];
b=diff(N',x,1)*N;
c=int(b,0,p)

m=N'*N;
f=int(m,0,p);
z=420*f

% ????  ??
fprintf('??Euler-Bernoulli ?????')
  












