Matlab ???????
%???????
     %1. 这是么比 
     
     %2. ?????????
     
% ??Euler-Bernoulli Beam ????
     % ????????????????????????????????
syms gMaterial gNode L gElement gEigValue gEigVector 
%?????
fprintf('??Euler-Bernoulli Beam ???????????2\n') 
fprintf('??Timoshenko?????????3\n???????????????\n')

%?????????
O=input('??????(<6?: ');
Nel=input('?????? ');      
L=input('?????m?: '); 
E=input('??????????10^11pa??');
E=E*10^11;
S=input('?????????1?????2?????3??????');
%?????
if S==1
    D=input('???????m?:');
    I=D^4/64;
    A=pi*D^24;
end
if S==2
    a=input('???????m?: ');
    b=input('???????m?: ');
    I=a^3*b/12;
    A=a*b;
end
if S==3
    D=input('???????m?:');
    d=input('???????m?:');
    I=(D^4-d^4)/64;
    A=pi*(D^2-d^2)/4;
end

%???????????????
Nnode=Nel+1;        %???????
node=(1:Nnode)';    %?????
x=0:L/Nel:L;         
xx=x';              %??x????
yy=zeros(Nnode,1);
gNode=[ node,       xx,          yy];         

%????????
gMaterial=[  E          I           A          7800 ];    
gElement=[(1:Nel)',    (1:Nel)',     (2:Nnode)'];
gBco=[1,        1,          0
      1,        2,          0
      Nnode,    1,          0
      Nnode,    2,          0];
      %????  ???     ?
  gNF=[5,         2,      -8e3];
  ie=1;   
    E = gMaterial( 1, 1 ) ;
    I = gMaterial( 1, 2 ) ;
    A = gMaterial( 1, 3 ) ;
    xi = gNode( gElement( ie, 2),  2 ) ;
    yi = gNode( gElement( ie, 2),  3 ) ;
    xj = gNode( gElement( ie, 3),  2 ) ;
    yj = gNode( gElement( ie, 3),  3 ) ;
    l =( (xj-xi)^2 + (yj-yi)^2 )^(1/2) ;
   c=E*I/(l^3);
    k=[12        6*l     -12        6*l;
        6*l      4*l^2   -6*l       2*l^2;
       -12       -6*l      12       -6*l;
        6*l      2*l^2   -6*l       4*l^2];
    ke=c*k;
   mass=gMaterial(1,4)*A*l;
   m=[ 156        22*l      54        -13*l;
       22*l       4*l^2     13*l      -3*l^2;
       54         13*l      156       -22*l;
      -13*l      -3*l^2    -22*l      4*l^2];
   me=mass/420*m;
   gK=zeros(Nnode*2);
   gM=zeros(Nnode*2);
   for ie=1:Nel
     for i=1:2
       for j=1:2
           for p=1:2
               for q=1:2
                   m=(i-1)*2+p;
                   n=(j-1)*2+q;
                   M=(gElement(ie,i)-1)*2+m;
                   N=(gElement(ie,j)-1)*2+n;
                   gK(M,N)=gK(M,N)+ke(m,n);
                   gM(M,N)=gM(M,N)+me(m,n);
               end
           end
        end
     end
   end
    [bc1_number,~]=size(gBco);
    w2max = max( diag(gK)./diag(gM) ); 
    for ibc=1:1:bc1_number
        n = gBco(ibc, 1 ) ;   
        
d = gBco(ibc, 2 );
        m = (n-1)*2 + d ;       
        gK(:,m) = zeros( Nnode*2, 1 ) ;
        gK(m,:) = zeros( 1, Nnode*2 ) ;
        gK(m,m) = 1;  
    end
    for ibc=1:1:bc1_number
        n = gBco(ibc, 1 ) ;
        d = gBco(ibc, 2 ) ;
        m = (n-1)*2 + d ;       
        gM(:,m) = zeros( Nnode*2, 1 ) ;
        gM(m,:) = zeros( 1, Nnode*2 ) ;
        gM(m,m) = gK(m,m)/w2max/1e10 ;      
    end
        for i=1:Nnode*2
        for j=i:Nnode*2
            gK(j,i) = gK(i,j) ;
            gM(j,i) = gM(i,j) ;
        end
    end
    [gEigVector, gEigValue] = eigs(gK, gM, O, 'SM' );  
    for ibc=1:1:bc1_number
        n = gBco(ibc, 1 ) ;
        d = gBco(ibc, 2 ) ;
        m = (n-1)*2 + d ;
        gEigVector(m,:) = gBco(ibc,3) ;
    end
    
      fre_number=length(diag(gEigValue));
    fprintf( '\n\n ??   ????(??)  \n' ) ;
    for i=1:fre_number
        fprintf( '----------------') ;
    end
    fprintf( '\n' ) ;
    for i=1:fre_number
        fprintf( '  %6d        ', i ) ;
    end
    fprintf( '\n' ) ;
    for i=1:fre_number
        fprintf( '----------------') ;
    end
    fprintf( '\n' ) ;
    [dof,~]=size(gEigVector) ;
    for i=1:dof
        for j=fre_number:-1:1
            fprintf( '%15.7e ', gEigVector(i,j) ) ;
        end
        fprintf( '\n' ) ;
    end
    for i=1:fre_number
        fprintf( '----------------') ;
    end
    fprintf( '\n' ) ;

    % ?????
    fprintf( '\n\n\n\n ??   ???(??)??  \n' ) ;
    fprintf( '----------------------------------------------------------\n') ;
    fprintf( '   ??       ???       ??????(Hz)     ???(rad/s)\n') ;
    fprintf( '----------------------------------------------------------\n') ;
    for i=fre_number:-1:1
        fprintf( '%6d   %15.7e   %15.7e   %15.7e\n', fre_number-i+1, ...
            gEigValue(i,i), sqrt(gEigValue(i,i))/2/pi, sqrt(gEigValue(i,i)) ) ;
    end
    fprintf( '----------------------------------------------------------\n') ;
    for j=fre_number:-1:1
        figure
        x = gNode(:,2) ;
        y = gNode(:,3) ;
        dx = gEigVector(1:2:length(x)*2, j ) ;
        dy = gEigVector(2:2:length(x)*2, j );
        factor = max( [max(abs(x))/max(abs(dx)), max(abs(y))/max(abs(dy))] )* 0.05; 
        plot(x,y,'-', x+factor*dx, y+factor*dy,'r:');
        hold on
        title( sprintf( '?%d???: %.3f Hz', fre_number-j+1, sqrt(gEigValue(j,j))/2/pi ) ) ;
              grid on
    end
   
%---------------------------------------------------------------------
%?????????

% ????????
% ??????????????????
% gNode -------- ????
% gElement ----- ????
% gMaterial ---- ?????????????????????????
% gBC1 --------- ????
% gDeltaT ------ ????
% gTimeEnd ----- ??????
% gDisp -------- ??????
% gVelo -------- ??????
% gAcce -------- ???????

%??????

global gDeltaT gTimeEnd gDisp gVelo gAcce 
    %??????(Newmark?)
    gDeltaT=0.01;                 
    gTimeEnd = 4890*gDeltaT  ;    
    timestep = floor(gTimeEnd/gDeltaT) ;
   
 % ???????????
    gDisp = zeros( Nnode*2, timestep ) ;
    gVelo = zeros( Nnode*2, timestep ) ;
    gAcce = zeros( Nnode*2, timestep ) ;
    gForc = zeros( Nnode*2, timestep ) ;
    
    % ????
    gDisp(:,1) = zeros(Nnode*2, 1 ) ;    
    gVelo(:,1) = ones (Nnode*2, 1) ;     
    %gForce(:,1)                          
    gama = 0.5 ;
    beta = 0.25 ; 
    C=zeros(size(gK));
    [~,N]=size(gK);
    alpha0 = 1/beta/gDeltaT^2;
    alpha1 = gama/beta/gDeltaT;
    alpha2 = 1/beta/gDeltaT;
    alpha3 = 1/2/beta - 1;
    alpha4 = gama/beta - 1;
    alpha5 = gDeltaT/2*(gama/beta-2);
    alpha6 = gDeltaT*(1-gama);
    alpha7 = gama*gDeltaT;
    K1 = gK + alpha0*gM + alpha1*C;     
   [bc1_number,dummy] = size(gBco) ;
   K1im = zeros(N,bc1_number) ;
    for ibc=1:1:bc1_number
        n=gBco(ibc,1);
        d=gBco(ibc,2);
        m=(n-1)*2+d;
        K1im(:,ibc)=K1(:,m); 
        K1(:,m) = zeros( Nnode*2, 1 );
        K1(m,:) = zeros( 1, Nnode*2);
        K1(m,m) = 1.0;
    end
    [KL,KU] = lu(K1);   % ????????????????
    
   %???????
    gAcce(:,1) =gM\(-gK*gDisp(:,1)-C*gVelo(:,1));
    
   %?????????
    for i=2:1:timestep
        if mod(i,100) == 0
            fprintf( '??????%d\n', i );    % ??????
        end
        f1 =gM*(alpha0*gDisp(:,i-1)+alpha2*gVelo(:,i-1)+alpha3*gAcce(:,i-1)) ...
                  + C*(alpha1*gDisp(:,i-1)+alpha4*gVelo(:,i-1)+alpha5*gAcce(:,i-1)) ;
        
        % ?f1????????
        [bc1_number,dummy] = size( gBco ) ;
        for ibc=1:1:bc1_number
            n = gBco(ibc, 1 ) ;
            d = gBco(ibc, 2 ) ;
            m = (n-1)*2 + d ;
            f1 = f1 - gBco(ibc,3) * K1im(:,ibc) ;
            f1(m) = gBco(ibc,3) ;
        end
        y = KL\f1 ;
        gDisp(:,i) = KU\y ;
        gAcce(:,i) = alpha0*(gDisp(:,i)-gDisp(:,i-1)) - alpha2*gVelo(:,i-1) - alpha3*gAcce(:,i-1) ;
        gVelo(:,i) = gVelo(:,i-1) + alpha6*gAcce(:,i-1) + alpha7*gAcce(:,i) ;
    end
return

global  gElement gMaterial gBC1 gDisp gVelo gAcce gDeltaT gTimeEnd

    t = 0:gDeltaT:gTimeEnd-gDeltaT;
    d = gDisp((floor(Nnode/2)*2)+1,:);
    subplot(2,1,1);
    plot( t, d );
    title( 'L/4???????');
    xlabel( '??(s)');
    ylabel( '??(cm)' );
    grid on
    
    fd = fft( d ) ;
    df = 1/gTimeEnd ;
    f = (0:length(d)-1)*df ;
    subplot(2,1,2);
    plot(f,abs(fd)) ;
    set(gca,'xlim',[0,100]) ;
    title( 'L/4???????' ) ;
    xlabel( '??(Hz)') ;
    ylabel( '??' ) ;
    
    fifi1 = diff(abs(fd));
    n = length(fifi1) ;
    d1 = fifi1(1:n-1);
    d2 = fifi1(2:n) ;
    indmax = find( d1.*d2<0 & d1>0 )+1;
    for i=1:length(indmax)
        if f(indmax(i)) > 10 
            break ;
        end
        text( f(indmax(i)+2), abs(fd(indmax(i)))*0.9, sprintf('f=%.3f',f(indmax(i))));
    end
 
