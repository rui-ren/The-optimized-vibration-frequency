UTF-8
function AssembleGlobalMatrix( ie, ke, me )
%  參等啋試僅睿窐講撻淝摩傖善淕极試僅撻淝
%      ie  --- 等啋瘍
%      ke  --- 等啋試僅撻淝
%      me  --- 等啋窐講撻淝
%  殿隙硉:

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

% 潰脤羶衄恀枙
