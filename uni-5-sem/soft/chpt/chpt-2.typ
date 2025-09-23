//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Elasticity
Elastic soft matter could be rubber or gel. The main property of these materials is that we can deform them, and after they'll recover they shape---these typically consist of polymers, in solution (cross-links) or not in solution.

== Stress and Strain
To deform a body we can use uniaxial tension and shear deformations---plus any combination of these. To capture this we define the Cauchy stress tensor $sigma_(alpha beta)$
$
  sigma_(alpha beta)^((d)) = F_alpha/A_beta
$
with $alpha, beta in {x,y,z}$ and $d<->"deviation"$---with the definition being fairly obvious, tension could be $sigma_(x x)$ while shear could be $sigma_(x y)$. Including ambient pressure the stress tensor becomes
$
  sigma_(alpha beta) equiv - P delta_(alpha beta) + sigma_(alpha beta)^((d))
$
what is the response to the applied stress---how does $r_alpha -> r'_alpha equiv r'_alpha (r_beta)$ due to the deformation? If we do a simple stretching of a body then
$
  r'_alpha = Lambda_(alpha alpha) r_alpha
$
where $Lambda_(alpha alpha)$ is some scaling---in the simplest case $Lambda_(alpha alpha) = L_alpha'\/L_alpha = lambda$ with $L_alpha$ being the length of our body before and after the deformation. For an incompressible material $V = V'$, if we're pulling along $x$ then
$
  L_x L_y L_z = L'_x L'_y L'_z
$
but $L'_y = L'_z$, so
$
  1 = V'/V = L'_x/L_x (L'_y/L_y)^2 = lambda_(parallel) (lambda_perp)^2 => lambda_perp = lambda_parallel^(-1\/2)
$
so in this case (volume preserving uniaxial deformation) $Lambda_(x x) = lambda => Lambda_(y y) = Lambda_(z z) = lambda^(-1\/2)$. For a simple shear deformation we can define $gamma = tan theta$---with $theta$ characterizing the deformation. If we shear along $x$ then $r'_z = r_z$ and $r'_y = r_y$, but $r'_x = r_x + gamma r_y$.

To generalize this we use the deformation gradient tensor $E_(alpha beta)$ defined by
$
  E_(alpha beta) equiv pdv(r'_alpha, r_beta)
$
so these are basically just the transformation coefficients. By definition we have
$
  dd(r'_alpha) = E_(alpha beta) dd(r_beta)
$
and
$
  dd(s^2) = dd(r_alpha, r_alpha)
$
so
$
  (dd(s)')^2 = dd(r'_alpha, r'_alpha) &= E_(alpha gamma) dd(r_gamma) E_(alpha delta) dd(r_delta) \
  &= E_(alpha gamma) E_(alpha delta) dd(r_gamma, r_delta) \
  &equiv C_(gamma delta) dd(r_gamma, r_delta)
$
with $C_(alpha beta)$ being the right Cauchy-Green tensor---note that $C_(alpha beta)=(E^T E)_(alpha beta) =E^T_(alpha gamma) E_(gamma beta)$. The previous examples give
$
  C_(alpha beta)^"scale" & = mat(Lambda_(x x)^2, 0, 0; 0, Lambda_(y y)^2, 0; 0, 0, Lambda_(z z)^2) \
  C_(alpha beta)^"shear" & = mat(1, gamma, 0; gamma, 1+ gamma^2, 0; 0, 0, 1)
$

We can now consider the extension
$
  (dd(s'))^2-(dd(s))^2 &= (C_(alpha beta)-delta_(alpha beta)) dd(r_alpha, r_beta) \
  &= 2 cal(E)_(alpha beta) dd(r_alpha, r_beta)
$
with $cal(E)_(alpha beta)$ being the Lagrangian strain tensor
$
  cal(E)_(alpha beta) equiv 1/2 (E_(gamma alpha) E_(gamma beta) - delta_(alpha beta))
$
this is nice because doing nothing gives $cal(E)_(alpha beta) = 0$ whereas $C_(alpha beta) = delta_(alpha beta)$.

Now we wan't to connect the Cauchy stress tensor (the force we apply) and the Lagrangian strain tensor (the response)---we assume a linear relationship, giving Hooke's law
$
  sigma_(alpha beta)^((d)) = K_(alpha beta gamma delta) cal(E)_(gamma delta)
$
so in the extreme case $K_(alpha beta gamma delta)$ is a four-index tensor---in what we're doing our material is typically homogeneous and isotropic such that
$
  K_(alpha beta gamma delta) = K delta_(alpha beta) delta_(gamma delta) + G (delta_(alpha beta) delta_(beta gamma) + delta_(alpha delta) delta_(beta gamma) - 2/3 delta_(alpha beta) delta_(gamma delta))
$
with $K$ being the bulk modulus and $G$ being the shear modulus---now we just have two material parameters instead of $81$. $K$ represents compressibily and is typically very large, meaning we assume our material is incompressible---any deformation is represented by $G$. For a simple shear
$
  sigma_(x y)^((d)) = sigma_(x y) = G gamma
$
and for simple uniaxial tension
$
  sigma_(x x) = - P + F_x/A_x => sigma_(x x) - sigma_(y y) = F_x/A_x = sigma_N = G(lambda^2 - lambda^(-1))
$

== Free energy
We can also look at the free energy density $f(E_(alpha beta))$ which should be invariant under coordinate transformations---given $B_(alpha beta) = (E E^T)_(alpha beta)$ the left Cauchy-Green tensor
$
  det (B_(alpha beta) - lambda delta_(alpha beta)) &= - lambda^3 + I_1 lambda^2 + I_2 lambda + I_3
$
with
$
  I_1 &= tr B = lambda_x^2 + lambda_y^2 + lambda_z^2 \
  I_2 &= 1/2 (tr (B)^2 - tr(B^2)) = lambda_x^2 lambda_y^2 + lambda_y^2 lambda_z^2 + lambda_x^2 lambda_z^2 \
  I_3 &= det B = lambda_x^2 lambda_y^2 lambda_z^2
$
these are all invariants by definition. So all dependence must lie in these
$
  f(E_(alpha beta)) = f(I_1,I_2,I_3)
$
for an incompressible body $I_3 = 1$ and drops out. We can Taylor expand this guy and find to lowest order
$
  f(I_1, I_2) = C_1 (I_1 - 3) + C_2 (I_2 - 3)
$
this is what characterizes a Mooney-Rivlin solid.

== Kuhn theory
We want to derive $G$ for rubber---polymer $+$ crosslinks. We assume loose strands follow $P(arrow(R))$ and $R'_alpha = E_(alpha beta) R_beta$. The free energy difference caused by the deformation is
$
  dd(F, d: Delta) &= rho_s [F(E_(alpha beta))- F(II)] \
  &= rho_s [1/2 (3 k_B T)/(b^2 N) (R'_alpha)^2 - 1/2 (3 k_B T)/(b^2 N) (R_alpha)^2]
$
with $rho_s$ being the density of strands.
