To develop a distributed-parameter model for this horizontal casing, we use the extended Hamilton’s principle which can be expressed as follows(Ren R, et al, 2017):
∫_(t_1)^(t_2)▒〖δ(T_c-V_c)〗 dt+∫_(t_1)^(t_2)▒〖δW_c dt〗=0                       (1)
where, V_C is the total potential energy of the coupled system,W_C, is the virtual work done by Non-conservative forces,
The total kinetic energy of the coupled system
T_C is the total kinetic energy of the coupled system that can be calculated by: 
T_c=T_p+T_f                            (2)
     T_p=∫_0^L▒〖1/2 m_p 〖(w ̇)〗^2 dx〗                         (3)
T_f=∫_0^L▒〖1/2 m_f [(w ̇+U_fi  ∂w/∂x)^2+U_fi^2]dx〗                (4)
where, T_P and T_f are the kinetic energy of casing and the internal slurry respectively,  m_p and m_f is the unit mass of the casing and slurry, w is the lateral displacement of the casing, and the dot (•) donates the differentiation with respect to time t.
The total potential energy of the coupled system
The V_C term is the potential energy resulting from casing bending and slurry pressure:
V_c=V_P+V_N                                                (5)
where V_p the flexural potential energy due to the casing’s deformation. V_N is the axis average potential energy of casing caused by inner pressure. According to Euler-Bernoulli theory (Gere, J.M. 1997), curvature k is equal to w’’, thus, the strain of casing can be obtained:
        ε=-w^'' y                              (6)
where y is the distance from any point of cross-section to the center point of casing. Hence, V_P can be expressed as follows:
        V_P=∫_Ω▒〖1/2 E〖(-w^'' y)〗^2 dΩ=∫_0^L▒〖1/2 EI_y 〖(w^'')〗^2 dx〗〗            (7)
where Ω is the volume of the casing and EI_yis the flexural rigidity. L is the length of casing.
Taking into account Poisson coupling, the inner pressure of casing will change the axis potential energy of casing, which can be expressed as:
      V_N=∫_0^L▒〖1/2 NA〖(w^')〗^2 dx〗                        (8)
where N is the inner pressure, A is the inner cross-sectional area of casing. The prime ('') denotes differentiation with respect to axis coordinate x.
The virtual work from Non-conservative force
We assume that the velocity of slurry in the annulus is U_fo, hence, taking Fluid-Structure Interaction (FSI) into account, the coupling velocity of slurry in annulus is 
     V_((x,t))=w ̇+U_fo w^'                         (9)
In this case, according to Paidoussis. M.P. (1998); Paidoussis and Prabhakar, S. (2007), the lateral flow is identical to the 2-D potential flow that would result from the motion of the cylinder with velocity V_((x,t))through the outside fluid. The flow is supposed to have momentum m_a×V_((x,t)) and m_a is the unit mass added by the slurry. Under the consideration of the annulus (Wu T. et al, 1995), m_a can be expressed as:
  m_a=1/4 e_f πD^2  (B^2+D^2)/(B^2-D^2 )                         (10)
Where B is the size of drill bit, D is the diameter of casing.                             
The variation of this momentum gives rise to normal and tangential lateral forces on the cylinder, of which is the inviscid unit force:
              F_a=m_a (w ̈+U_fo w ̇^')                         (11)
   δW_c=∫_0^L▒〖(f(x)-cw ̇-F_a)δwdx〗                   (12)
where c is the damping coefficiency, δw is the virtual displacement, and f(x) is the external force.
Boundary Conditions
For the casing and slurry coupling system, the casing is emerged by slurry and centralizers are assumed as clamped support points, thus the boundary conditions for casing can be expressed as(Wang L, 2018 and Zalluhoglu U, 2019):
{█(w^' (0,t)=w^' (L,t)=0@w^'' (0,t)=w^'' (L,t)=0)┤                       (13)
Using the standard variation techniques (Rockafellar, et al. 2013) in accordance with the boundary conditions for a simply-supported casing, the following corresponding variation results are obtained:
{█(∫_(t_1)^(t_2)▒〖δT_P dt=-∫_(t_1)^(t_2)▒〖∫_0^L▒m_p  w ̈δwdxdt〗〗                         @∫_(t_1)^(t_2)▒〖δT_fi 〗 dt=-∫_(t_1)^(t_2)▒∫_0^L▒〖m_f (w ̈+U_fi^2 w^''+2w ̇^' U_fi )δwdxdt〗       @∫_(t_1)^(t_2)▒〖δV_P dt=∫_(t_1)^(t_2)▒∫_0^L▒〖EI_y w^'''' δwdxdt〗〗                         @∫_(t_1)^(t_2)▒〖δV_N dt=-∫_(t_1)^(t_2)▒∫_0^L▒〖NAw^'' δwdxdt〗〗                        @∫_(t_1)^(t_2)▒〖δW_c dt=∫_(t_1)^(t_2)▒∫_0^L▒〖(f(x)-cw ̇-m_a w ̈-m_a U_fo w ̇^')δwdxdt〗〗    )┤   (14)
Then, taking Eq. (14) into Eq. (1), according to variational analysis, the following reduced equation is obtained:
EI_y w^''''+(m_f U_fi^2-NA) w^''+(2m_f U_fi+m_a U_fo ) (w^' ) ̇+mw ̈+cw ̇=f(x)   (15)
  where, EI_y is bending stiffness of casing, w is the transverse displacement of casing, m_f is the mass of slurry per unit length, m_p is the mass of casing per unit length, m_a is added mass per unit length, U_fi is the velocity of inner slurry flow, U_fo is the velocity of slurry in annulus, N is the pumping pressure in casing, A is the flow area of casing, c is the viscous damping coefficiency, f(x) is the external force . Where m=m_f+m_p+m_a

