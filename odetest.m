function xdot=odetest(t,x)
xdot=zero(2,1);
xdot(1)=x(2);
xdot(2)=-0.5*x(2)-2*x(1)-x(1)^2;

tspan=[0 20];
x0=[0 1];
[t,x]=ode23(@odetest,tspan,x0);
plot(t,x)
end
