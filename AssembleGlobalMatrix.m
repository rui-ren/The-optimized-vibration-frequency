function AssembleGlobalMatrix( ie, ke, me )
%  �ѵ�Ԫ�նȺ��������󼯳ɵ�����նȾ���
%      ie  --- ��Ԫ��
%      ke  --- ��Ԫ�նȾ���
%      me  --- ��Ԫ��������
%  ����ֵ:

  global gElement gK gM
    for i=1:1:2
        for j=1:1:2
            for p=1:1:2
                for q =1:1:2
                    m = (i-1)*2+p;
                    n = (j-1)*2+q ;
                    M = (gElement(ie,i)-1)*2+p ;
                    N = (gElement(ie,j)-1)*2+q ;
                    gK(M,N) = gK(M,N) + ke(m,n) ;
                    gM(M,N) = gM(M,N) + me(m,n) ;
                end
            end
        end
    end
return

% ���û������