function   me=MassMatrix_deformation(ie)
    me=m/420*...            
	[156*p   22*p^2  54*p    -13*p^2;...
    22*p^2  4*p^3  13*p^2   -3*p^3;...
    54*p     13*p^2    156*p    -22*p^2;...
    -13*p^2  -3*p^3  -22*p^2  4*p^3];
    T = TransformMatrix( ie ) ;
    me = T*me*transpose(T) ;
    return
    