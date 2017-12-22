
    % 根据形函数、将移动载荷分解到节点等效力。
    % 对载荷位置进行初步判定
    timespace=gTimeEnd/gDeltaT;
    v=input('输入载荷移动的速度： ');
    Force_moving=input('输入移动载荷力的大小：');
    gDeltaT=input('输入时间步长：     ');
    f=zeros(gNode*2,gTimeEnd/gDeltaT);                                          % 定义力矩阵
   gDisp=zeros(gNode*2,gTimeEnd/gDeltaT);                                  % 定义位移矩阵                  
   gVelo=zeros(gNode*2,gTimeEnd/gDeltaT);                                  % 定义速度矩阵
   
    for i=1:gTimeEnd/gDeltaT
        cosin=mod(i*v,p)/p;
        N1=1-3*cosin^2+2*cosin^3;
        N2=(cosin-2*cosin^2+cosin^3)*l;
        N3=3*cosin^2-2*cosin^3;
        N4=(-cosin^2+cosin^3)*l;
        N=[N1, N2, N3, N4]';
        f(:,i)=Force_moving*N;
    end                                                                  % 施加移动力进行计算
    

   %采用振动Newmark方法进行振动分析
    gama = 0.5 ;
    beta = 0.25 ;                                                   % 采用平均加速度方法 Newmark- beta 方法，平均加速度法
    [~,N]=size(gK);
    alpha0 = 1/beta/gDeltaT^2;
    alpha1 = gama/beta/gDeltaT;
    alpha2 = 1/beta/gDeltaT;
    alpha3 = 1/2/beta - 1;
    alpha4 = gama/beta - 1;
    alpha5 = gDeltaT/2*(gama/beta-2);
    alpha6 = gDeltaT*(1-gama);
    alpha7 = gama*gDeltaT;
    K1 = gK + alpha0*gM + alpha1*gC;            %计算有效刚度矩阵
     
%-------------------------------------------------------------------------
      
  
% 对K1进行边界条件处理
   [bc1_number,dummy] = size(gBC1) ;
   K1im = zeros(Nnode, bc1_number);
    for ibc=1:1:bc1_number
        n=gBC1(ibc,1);
        d=gBC1(ibc,2);
        m=(n-1)*2+d;
        K1im(:,ibc)=K1(:,m);                                 %这是将原始边界条件储存到Klim中去，方便后面对力进行施加边界条件
        K1(:,m) = zeros( Nnode*2, 1 );                  %将有效刚度矩阵进行赋值
        K1(m,:) = zeros( 1, Nnode*2);                   % 化行、化列法对边界条件进行施加
        K1(m,m) = 1.0;                                          %施加边界条件
    end
    [KL,KU] = lu(K1);                                         % 进行LU分解，节省后面的求解时间
   
   %计算初始加速度
    gAcce(:,1) =gM\(-gK*gDisp(:,1)-gC*gVelo(:,1));         % 初始加速度一般都设为零的      
   
   %对每一个时间步计算
    for i=2:1:timestep
        if mod(i,100) == 0
            fprintf( '当前时间步：%d\n', i );          % 显示整数步长，模为零的情况
        end
        f1 =f(:,i)+gM*(alpha0*gDisp(:,i-1)+alpha2*gVelo(:,i-1)+alpha3*gAcce(:,i-1)) ...
                  + gC*(alpha1*gDisp(:,i-1)+alpha4*gVelo(:,i-1)+alpha5*gAcce(:,i-1)) ;            
       
        % 对f1进行边界条件处理, 施加力的边界条件
        [bc1_number,dummy] = size( gBC1 ) ;
        for ibc=1:1:bc1_number
            n = gBC1(ibc, 1 ) ;
            d = gBC1(ibc, 2 ) ;
            m = (n-1)*2 + d ;
        %如果是力  那么需要对力进行施加边界条件
            f1 = f1 - gBC1(ibc,3) * K1im(:,ibc) ;      % 这个是施加边界条件    采用化行化列法施加边界条件                                          
        end
        y = KL\f1 ;
        gDisp(:,i) = KU\y ;
        gAcce(:,i) = alpha0*(gDisp(:,i)-gDisp(:,i-1)) - alpha2*gVelo(:,i-1) - alpha3*gAcce(:,i-1) ;
        gVelo(:,i) = gVelo(:,i-1) + alpha6*gAcce(:,i-1) + alpha7*gAcce(:,i) ;
    end
    
   % 绘制时程曲线
    t = 0:gDeltaT:(gTimeEnd-gDeltaT);
    d = gDisp((floor(Nnode/2)*2)+1,:);
    subplot(2,1,1);
    plot( t, d );
    title( 'L/4处挠度时程曲线');
    xlabel( '时间(s)');
    ylabel( '挠度(cm)' );
    grid on
   
    % 对时程曲线进行FFT变换，获取频谱特性
    fd = fft( d ) ;
    df = 1/gTimeEnd ;
    f = (0:length(d)-1)*df ;
    subplot(2,1,2);
    plot(f,abs(fd)) ;
    set(gca,'xlim',[0,10]) ;
    title( 'L/4处挠度的频谱图' ) ;
    xlabel( '频率(Hz)') ;
    ylabel( '幅值' ) ;
    % 标注频率峰值
    fifi1 = diff(abs(fd));
    n = length(fifi1) ;
    d1 = fifi1(1:n-1);
    d2 = fifi1(2:n) ;
    indmax = find( d1.*d2<0 & d1>0 )+1;
    for i=1:length(indmax)
        if f(indmax(i)) > 10
            break;
        end
        text( f(indmax(i)+2), abs(fd(indmax(i)))*0.9, sprintf('f=%.3f',f(indmax(i))));
    end
   
