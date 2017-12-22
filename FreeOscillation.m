function xdot=FreeOscillation(t,x,zeta,Alpha)
xdot=[x(2); -2.0*zeta*x(2)-x(1)-Alpha*x(1)^3];
end
