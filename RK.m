zeta=[0.1 1.0 5.0];
Alpha=[0.0, 0.0,0.0];
tspan=linspace(0,40,400);
lintype=char('-k','--k','-.k');
for i = 1: 3
    [t,x]=ode45(@FreeOscillation, tspan,[1 1],[],zeta(i),Alpha(i));
    figure(1);
    xlabel('Time(\tau)');
    ylabel('Displacement x(\tau)');
    axis([0 40 -1.5 1.5]);
    plot([0,40],[0,0],'k-');
    legend('\zeta=0.1','\zeta=1.0','\zeta=5.0');
    plot(t,x(:,1),lintype(i,:));
    hold on
    figure(2);
    xlabel('Displacement x(\tau)');
    ylabel('Velocity');
    title('Phase portrait');
    axis([-2.0 2.0 -2.0 2.0]);
    legend('\zeta=0.1','\zeta=1.0','\zeta=5.0');
    plot(x(:,1),x(:,2),lintype(i,:));
    hold on
end


