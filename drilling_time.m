 plot(1:166,A(1:166,1),'LineWidth',2)
 hold on
  plot(1:202,A(1:202,2), '-r', 'LineWidth',2)
  hold on
  plot(1:120, A(1:120,3), '.-k','LineWidth',2)
  hold on
  plot(1:162,A(1:162,4), '*g', 'LineWidth', 0.5)    % 有数据异常点
  hold on
  plot(1:170,A(1:170,5), '.m')
  hold on
  plot(1:204,A(1:204,6), 'ob')
  hold on
  plot(1:122,A(1:122,7), '*-r')
  hold on
  plot(1:144,A(1:144,8), '-.black')
  hold on
  plot(1:159,A(1:159,9), '*-b')
  hold on
  plot(1:140,A(1:140,10), '-.b')
  hold on
  plot(1:214,A(1:214,11), '*y')
  hold on
  plot(1:208,A(1:208,12), 'black', 'LineWidth', 1.5)
   set(gca,'xaxislocation','top')
   set(gca,'ydir','reverse')
    title('钻井周期表')
    gtext('磨溪009-4-X2' )   % 红线
    gtext('磨溪009-3-X2' )   % 黑线
     gtext('214天')
    xlabel('时间(d)')
    ylabel('井深(m)')
    
    
 plot(1:166,A(1:166,1),'LineWidth',2)
 hold on
  plot(1:202,A(1:202,2), '-r', 'LineWidth',1.5)
  legend('磨溪009-8-X1','磨溪009-2-H2','磨溪008-X16','磨溪008-20-H2','磨溪008-18X1','磨溪008-15-H1','磨溪008-11-X2','磨溪008-11-X1','磨溪008-7-X2')
  
  
 plot(1:246,D(1:246,1),'ob','LineWidth',2)  
 hold on
  plot(1:145,D(1:145,2),'r','LineWidth',2)
  hold on
   plot(1:177,D(1:177,1),'.-k','LineWidth',2)
   hold on
   set(gca,'xaxislocation','top')
   set(gca,'ydir','reverse')
   legend('高石001-H2井','高石001-X2井','高石001-X2井')
 gtext('246天')
  gtext('145天')
    gtext('高石001-2井')
   grid on
   
       xlabel('时间(d)')
    ylabel('井深(m)')
    
    
      