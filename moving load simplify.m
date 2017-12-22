
    %%  以下是移动载荷下的激振问题， 模拟移动载荷下的激振问题
       
    % 根据形函数、将移动载荷分解到节点等效力。
    % 对载荷位置进行初步判定
    t= input('输入载荷运动时间:    ');
    v=input('输入载荷移动的速度：');
    Force_moving=input('输入移动载荷力的大小：');
    delta_t=input('输入时间步长：');
    if (L/N) < (delta_t*v)
        printf('计算步长太大，请缩小步长');
    end
    delta_t=input('输入时间步长： ');
    tau=mod(t,delta_t);                                        %  这是tau的值，在A、B之间运动的时间
    S_moving=v*t;                                               %  t时间内，运动的位移
    j=floor(S_moving/p);                                    %  运动到的节点编号
    Element_x=mod(S_moving,p);                      %  移动载荷在微元体内的位置，即距离节点i 的距离。
    S_time=floor(t/delta_t)*v;                             %  整数步长内移动的距离
    S_element=gNode(j,2);                                  %   j 个节点所在的坐标
    % ---------------------------------------------------------------------------------------------------------------------------
    % 微元体内位置的判定, 以i 节点为分界点
  if S_time > S_element                                     %  判定A点是否在微元体内
       sigma_1=S_time -S_element;                     %  A 点到节点i 的距离
    % 微元体内位置的判定，以i+1节点为分界点
         if ceil(t/delta_t)*v>(j+1)*p
        sigma_2=p;                                                % 如果下一个时间步长后，移动力超出微元体i, 那么sigma_2 就等于微元体长度
         else
        sigma_2=sigma_1+delta_t*v;                    % 如果下一个时间步长后，移动内仍然在微元体j内，那么sigma_2级等于sigma_1 + delta_t *v
         end
    else 
        sigma_1=0;
        sigma_2 = ceil(t/delta_t)*v-i*p;
  end
    b1=1-3*sigma_1^2+2*sigma_1^3;
    b2=(6*sigma_1^2 - 6*sigma_1)*(sigma_2 - sigma_1);
    b3=(6*sigma_1 - 3)*(sigma_2 - sigma_1)^2;
    b4=2*(sigma_2 - sigma_1)^3;
    b5=(sigma_1 - 2*sigma_1^2 + sigma_1^3)*p;
    b6=(1-4*sigma_1+3*sigma_1^2)*(sigma_2 - sigma_1)*p;
    b7=(3*sigma_1 - 2 )*(sigma_2 - sigma_1)^2*p;
    b8=(sigma_2 - sigma_1)^3*p;
    b9=1-b1;
    b10=-b2;
    b11=-b3;
    b12=-b4;
    b13=(sigma_1^3 - sigma_1^2)*p;
    b14=(sigma_2 - sigma_1)*(3*sigma_1^2 - 2* sigma_1)* p;
    b15=(sigma_2-sigma_1)^2*(3*sigma_1 - 1)* p;
    b16=b8;
    
  % 根据Hermite 差值函数，将Force_moving进行离散。
   a0=[b1;b2;b3;b4];
   a1=[b5;b6;b7;b8];
   a2=[b9;b10;b11;b12];
   a3=[b13;b14;b15;b16];
   F=(a0+a1*(tau/delta_t)+a2*(tau/delta_t)^2 +a3*(tau/delta_t)^3)*Force_moving;       % 移动载荷离散后看成是初始的力，初始施加的力