function gC=DampingMatrix_deformation(ie)
global  mf Ui gK  gM w1 w2 ma p  n_element  gElement Uo

%科氏阻尼矩阵
	cb=-(2*mf*Ui + ma*Uo)*...
    [0          -p/10       -1/2          p/10;...
    p/10       0           -p/10         p^2/60;...
    1/2        p/10           0            -p/10;...
    -p/10     -p^2/60     p/10         0];      

     
Cb=zeros(Nnode*2);
    for ie=1:n_element                                                    % Nel 表示有多少个微元
     for i=1:2
       for j=1:2
           for p=1:2
               for q=1:2
                   m=(i-1)*2+p;
                   n=(j-1)*2+q;
                   M=(gElement(ie,i)-1)*2+m;      
                   N=(gElement(ie,j)-1)*2+n;
                  Cb(M,N)=Cb(M,N)+cb(m,n);
               end
           end
        end
     end
    end
    
    %No.9 水泥浆引起的微元的粘性阻尼矩阵 
    dRatio=0.008;                                               % 结构阻尼比，钢材水泥选取0.008
    % Rayleigh Damping                                     % 粘性阻尼，采用比例阻尼方式
    alpha=2*(w1*w2)*dRatio/(w1+w2);            % w1、w2是管材的固有振动频率
    beta= 2*dRatio/(w1+w2);
    Ca=alpha*gM+beta*gK;                               % rayleigh 方法确定的结构阻尼矩阵
     % 阻尼矩阵由rayleigh 方法得出：
     % 总共的阻尼矩阵是：  
     gC=Cb+Ca;
     T = TransformMatrix( ie ) ;
    gC = T*gC*transpose(T) ;
     return
     
     % 已检查，没有问题
     