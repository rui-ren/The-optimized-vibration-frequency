
% 模态分析理论    % 单位换算问题！！ unite is correct
% No.1 
%基本数据的输入
fprintf( '采用Euler-Bernoulli Beam 进行求解、节点自由度为2\n' );
L=input( '系统长度(m):  ' );
Nel=input( '微元体数量：' );
E=input( '输入材料的弹性模量(10^11pa):   ');
B=input( '输入钻头的外径(m)：   ');
D=input( '输入套管的外径(m)：   ');
d=input( '输入套管的内径(m):    ');
ef=input( '输入水泥浆的密度（kg/m^3):    ');
ep=input( '输入套管的密度 (kg/m^3)：     ');
                                                                  % mpa.s 和 pa.s 的单位换算
I=pi*(D^4-d^4)/64;                                                % 转动惯量的求解
A1=pi*(B^2)/4;                                                    % 井眼横截面积
A2=pi*(D^2)/4;                                                    % 套管外面积
A3=pi*(d^2)/4;                                                    % 套管内面积
E=E*10^11; 
cm=(B^2+D^2)/(B^2-D^2);
ma=cm*ef*A2;                                                      
mf=ef*A3;
mp=ep*(A2-A3);
m=mf+mp+ma;

%微元体节点进行编号
Nnode=Nel+1;                                                  % 节点总数
node=(1:Nnode);                                               % 生成节点向量
x=0:(L/Nel):L;                                                % 对节点进行坐标编号
xx=x';                                                        % 节点x坐标向量
yy=zeros(Nnode,1);                                            % 节点的node
                %节点编号            节点x坐标             节点y坐标
gNode=[   transpose(node)              xx                    yy];   %节点进行编号

 %No.3微元体和节点的关系矩阵

              %微元体编号             左端节点             右端节点
gElement=[    (1:Nel)',               (1:Nel)',            (2:Nnode)'];  


  %No.4 第一边界条件条件

	    %节点号       自由度号               边界值
 gBco=[    1,              1,                 0
	       1,              2,                 0
	       Nnode,          1,                 0
		   Nnode,          2,                 0];
         
   % No.5 微元体的长度
   xi=gNode(gElement(1,2),2);
   xj=gNode(gElement(1,3),2);
   yi=gNode(gElement(1,2),3);
   yj=gNode(gElement(1,3),3);
   p=sqrt((xi-xj)^2+(yi-yj)^2);    
   
   % 计算微元体的质量矩阵
     me=m/420*...            
	[156*p  22*p^2  54*p     -13*p^2;...
    22*p^2 4*p^3  13*p^2      -3*p^3;...
    54*p  13*p^2   156*p     -22*p^2;...
    -13*p^2  -3*p^3  -22*p^2  4*p^3];

	%套管的微元刚度矩阵
     Kea=E*I/(p^3)*...
	[12     6*p    -12     6*p;...
    6*p   4*p^2   -6*p    2*p^2;...
    -12    -6*p     12    -6*p;...
    6*p    2*p^2  -6*p    4*p^2];

 %微元的总体刚度矩阵
	   ke=Kea; 
   % No7. 质量矩阵和刚度矩阵进行组装 
       gK=zeros(Nnode*2);
       gM=zeros(Nnode*2);
     %按照微元体进行装配       
   for ie=1:Nel        % Nel 表示有多少个微元
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
   
   % No.8 采用第一类边界条件进行施加边界条件 
  [bc1_number, ~]=size(gBco);
  w2max = max( diag(gK)./diag(gM) ); 
    
   for ibc=1:1:bc1_number
        n = gBco(ibc, 1 );                                       %这里查找的是节点
        d = gBco(ibc, 2 );                                       %查找约束施加的自由度
        m = (n-1)*2 + d;                                         %计算约束自由度在总刚矩阵中占用的自由度
        gK(:,m) = zeros( Nnode*2, 1 );                           %列化成0
        gK(m,:) = zeros( 1, Nnode*2 );                           %行化成0
        gK(m,m) = 1;  
   end
    
   for ibc=1:1:bc1_number
        n = gBco(ibc, 1 );
        d = gBco(ibc, 2 );
        m = (n-1)*2 + d;      
        gM(:,m) = zeros( Nnode*2, 1 );
        gM(m,:) = zeros( 1, Nnode*2 ) ;
        gM(m,m) = gK(m,m)/w2max/1e10 ;         
   end
    
    for i=1:Nnode*2
           for j=i:Nnode*2
               gK(j,i) = gK(i,j);
               gM(j,i) = gM(i,j);                           % 进行对称化矩阵
           end 
    end
                     
   % 计算特征值和特征想
    [gEigVector, gEigValue] = eigs(gK, gM, 6, 'SM' );       %提取三阶特征值 
    fre_number=length(diag(gEigValue));  
 
	% 打印特征值
    fprintf( '\n\n\n\n 表二   特征值(频率)列表  \n' ) ;
    fprintf( '----------------------------------------------------------\n') ;
    fprintf( '   阶数       特征值       固有振动频率(Hz)     圆频率(rad/s)\n') ;
    fprintf( '----------------------------------------------------------\n') ;
    for i=fre_number:-1:1
        fprintf( '%6d   %15.7e   %15.7e   %15.7e\n', fre_number-i+1, ...
            gEigValue(i,i), sqrt(gEigValue(i,i))/2/pi, sqrt(gEigValue(i,i)) ) ;
    end
    fprintf( '----------------------------------------------------------\n') ;

