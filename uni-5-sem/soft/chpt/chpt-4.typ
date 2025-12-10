//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Colloids
A simple qualitative definition of a colloid are particles in solution that are relatively small with sizes from $tilde 10 "nm"$ to $tilde 1 mu"m"$. At this small scale thermal motion becomes highly relevant. As an example consider a $1 mu"m"$ spherical particle of some typical density. We can compute that $E_g tilde 5 E_"th"$ for such a particle. Given these are comparable we would call this particle a colloid.

As physicists we are interested in interactions between colloids. We now describe the typical repulsive and attractive interactions between colloids. We will find that combining these lead to a complex potential landscape that can be probed by experiment.

== Dispersion forces
Dispersion forces are forces between electrically symmetric particles caused by fluctuations in their electron distributions. Typically dispersion forces refer to the induced dipole-dipole interaction. We refer to this as the van der Waals interaction
$
  V_"Waals" = - C (a_0/r)^6
$
with $C$ and $a_0$ being constants. As a pair-potential this is quite weak due to the $r^(-6)$ scaling. We now show this scales as $h^(-1)$ for colloids.

We consider two spherical colloids with radii $R$ separated by a distance $h << R$. This assumption allows us to approximate the colloids as being planes. We now consider the interaction between all volumes $dd(V)$ in each colloid and find
$
  W = - integral underbracket(dd(bold(r)_1) n, "#molecules" #linebreak() "at" bold(r)_1) integral dd(bold(r)_2) n underbracket((C a_0^6)/abs(bold(r)_1-bold(r)_2)^6, "van der Waals" #linebreak() "interaction")
$
with $n$ being the number density within a colloid. We compute the integral using cylindrical coordinates
$
  abs(bold(r)_1-bold(r)_2)^2 = rho^2 + (z_1-z_2)^2
$
also defining $tilde(rho) equiv rho^2$ we find
$
  W/A &= -n^2 C a_0^6 integral_(-oo)^0 dd(z_1) integral_h^oo dd(z_2) integral_0^oo (pi dd(tilde(rho)))/[tilde(rho) + (z_1-z_2)^2]^3 \
  &= (n^2 C a_0^6 pi)/6 integral_(-oo)^0 dd(z_1) 1/(z_1-h)^3 \
  &= - (n^2 C a_0^6 pi)/(12 h^2) equiv w(h)
$
We define the Hamaker constant $A_H equiv n^2 C a_0^6 pi^2$ then
$
  w(h) = - A_H/(12 pi h^2)
$
so in the plane-plane approximation the potential scales as $h^(-2)$.

We can generalize the above to spherical colloids. Let $h$ be the minimal distance between two spherical colloids. The distance for any $rho$ is then
$
  tilde(h) (rho) & = h + 2(R-sqrt(R^2-rho^2)) \
                 & tilde.equiv h + rho^2/R
$
Each value of $rho$ defines an annulus with area $dd(A) = 2 pi rho dd(rho)$. Then the total potential is given by sweeping the surface
$
  U(tilde(h)) & = integral w(tilde(h)) dd(A) = pi R integral_h^(h+R) w (tilde(h)) dd(tilde(h)) tilde.eq^"Derjaguin" - R/12 A_H/h
$
so the potential scales as $h^(-1)$! The last step uses the Derjaguin approximation $h + R -> oo$.

== Depletion forces
Depletion forces are forces caused by smaller things in solution. A typical example would be polymers. Generally $n_"small" >> n_"big"$ for these types of interactions.

Consider small particles in solution with big colloids. Let  the size of the colloids be $tilde R$ and the size of the small particles $tilde r$. Clearly if two colloids are within a distance $d$ less than $r$ then small particles are depleted from the region between them. This depletion creates a density imbalance leading to an osmotic pressure since
$
  p prop N/V
$
so $dd(p, d: Delta) = n_"small" k_B T$. This pressure forces the colloids together.

Entropically this is also obvious. The small particles want to maximize their entropy and this is done by maximizing their available volume. For spherical colloids we define the excluded volume
$
  V_"excluded" = (4 pi)/3 (R + r)^3
$
Two separated colloids would have excluded volume $2 V_"excluded"$. But if the colloids come together and touch, then their exluded volumes will overlap leading to an effective excluded volume less than $2 V_"excluded"$. This makes it favorable for colloids to clump.

== Hard sphere interaction
The hard sphere interaction aries since colloids can not overlap.

The naive hard sphere potential is
$
  V = cases(oo"," #h(10pt) & r <= 2R, 0"," & "otherwise")
$
so an infinite barrier.

A more realistic potential would be
$
  V_"Pauli" prop 1/R^12
$
or
$
  V_"exp" prop e^(-r rho^(-1))
$
these allow for overlap between colloids. Though making it very unfavorable.

== Grafting
We can coat colloids with polymers. This process is called grafting.

Consider attaching polymers to a colloid. At sufficient polymer density the individual polymers get squished and we obtain a _brush_. We would like to determine the height of this brush.

We define the _grafting density_ $Gamma_p$ as the number of polymers per unit area. This defines a length scale $cal(l) = Gamma_p^(1\/2)$. We imagine each polymer occupies a cylinder of radius $cal(l)$ being a random walk inside it. We want to determine the height $h_p$ of these cylinders. Using these we write the volume fraction $phi$ as
$
  phi equiv V_"polymer"/(V_"cylinder") = (N v_c)/(cal(l)^2 h_p) = (N v_c Gamma_p)/h_p
$
with $N$ being the number of monomers in a polymer and $v_c$ being the volume of a monomer.

The free energy due to the polymers is
$
  f_E (h_p) = underbracket((3 k_B T h_p^2)/(2 N b^2), "single polymer") Gamma_p
$
with $b$ being the monomer length.

The free energy of mixing is
$
  F_"mix"/A &= underbracket(f_"mix", "by Flory-Huggins") h_p \
  &= (k_B T)/v_c [(1-phi) ln(1-phi) + underbracket(chi (1-phi) phi, "polymer solvent" #linebreak() "interaction")] h_p
$
assuming $phi$ is small we find
$
  f_"mix" (h_p) & = (k_B T h_p)/v_c [phi^2/2 -phi+ chi(phi-phi^2)] \
                & = (k_B T h_p phi)/v_c [(chi-1) + (1/2-chi) phi]
$
We find the total free energy
$
  f(h_p) = (3 k_B T Gamma_p)/(2 N b^2) h_p^2 + N k_B T Gamma_p [(chi-1) + (1/2 - chi) (N v_c Gamma_p)/h_p]
$
giving the brush height in equilibrium
$
  h_("eq",p)^3 = (N^3 v_c b^2 Gamma_p)/3 (1/2-chi)
$
which is nice! We see it depends on the length of the polymers, the grafting density, and their interaction with the solvent.

Consider two grafted colloids. Let the distance between them be $d$ and define $h = d\/2$. Then for $d$ less than $2 h_p^"eq"$ their polymer brushes will begin compressing. We compute the interaction potential as
$
  w(h) &= 2 [f(h) - f(h_p^"eq")] \
  & tilde.eq^"Taylor" 2 [f(h_p^"eq") + overbracket(evaluated(dv(f, h))_(h=h_p^"eq"), 0) (h-h_p^"eq") + 1/2 evaluated(dv(f, h, 2))_(h=h_p^"eq") (h-h_p^"eq")^2 + dots - f(h_p^"eq")] \
  &tilde.eq evaluated(dv(f, h, 2))_(h=h_p^"eq") (h-h_p^"eq")^2
$
we compute
$
  dv(f, h, 2) = (3 k_B T Gamma_p)/(N b^2) + (2 N k_B T v_c Gamma_p^2)/(h^3) (1/2-chi)
$
so
$
  evaluated(dv(f, h, 2))_(h=h_p^"eq") = (9 k_B T Gamma_p)/(N b^2)
$
The interaction potential is then
$
  w(h) tilde.eq (9 k_B T Gamma_p)/(N b^2) (h-h_p^"eq")^2
$
which is nice! We find a spring-like interaction depending on the grafting density.

== Charges
We can charge colloids and typically most colloids are negatively charged when made.

Recall the Coulomb force
$
  bold(f) = (q_1 q_2)/(4 pi epsilon) bold(hat(r))_(1 2)/abs(bold(r)_(1 2))^2
$
here $epsilon$ is important since $epsilon_"water"$ is relatively large meaning the interaction is weak.

We might expect ions in solution stick to colloids. However, this is not favorable due to entropy. Instead an ionic cloud develops around all colloids. At the surface we have
$
  E = sigma/epsilon
$
where we use the plane approximation. And in the bulk we have
$
  nabla dot bold(E) = rho/epsilon;"  "bold(E) = - grad V
$
This gets complicated since charges generate some potential and move under the mean field due to all other charges. This movement changes the charge density which changes the potential. We seek a steady state solution.

In one dimension we have
$
  dv(E, x) = - dv(V, x, 2) = rho/epsilon
$
where $rho = c_"ions" q_"ions"$. The concentration is determined by Boltzmann statistics
$
  c_"ions" = c_0 exp(- (q_"ions" V)/(k_B T))
$
as usual. A short proof is given below.

#proof[

  Consider a system with some potential difference. This gives rise to a flux
  $
    j_E = (q E)/zeta c
  $
  where $zeta$ is some friction coefficient. By Fick's law we also have a diffusive flux
  $
    j_D = - D dv(c, x)
  $
  By Einstein we have $D= k_B T zeta^(-1)$ so
  $
    j_"tot" = j_E + j_D = D (- dv(c, x) + (q E c)/(k_B T)) =^! 0
  $
  giving
  $
    integral 1/c dv(c, x) dd(x) & = integral (q E)/(k_B T) dd(x) \
                       ln c/c_0 & = - (q V)/(k_B T) \
                              c & = c_0 exp(- (q V)/(k_B T))
  $
  and we are done.

]
Let $q_"ions" = e$ and substitute to find
$
         dv(V, x, 2) & = - (e c_0)/epsilon exp(- (e V)/(k_B T)) \
  dv(tilde(V), x, 2) & = - 4 pi cal(l)_B c_0 e^(-tilde(V))
$
this is the Boltzmann-Poisson equation. We have defined $tilde(V)$ and the Bjerrum length $cal(l)_B$ by
$
  tilde(V) equiv (e V)/(k_B T)";  " cal(l)_B equiv e^2/(4 pi epsilon k_B T)
$
$cal(l)_B$ is the length scale at which entropic and electrical contributions are similar.

Consider a colloid in a neutral solution with many ions. Away from the colloid we have $n_(+0) = n_(-0) = n_0$ and generally we have $rho = e (n_+ - n_-)$. These densities are given by the Boltzmann distribution
$
  n_plus.minus = n_0 exp(minus.plus (e V)/(k_B T))
$
substituting we find
$
  rho & = - 2 e n_0 sinh (e V)/(k_B T)
$
Then
$
  dv(V, x, 2) & = (2 e n_0)/epsilon sinh (e V)/(k_B T) \
              & tilde.eq^(e V\/k_B T << 1) (2 e n_0)/epsilon (e V)/(k_B T)
$
where we use the Debye-Hueckel approximation $e V << k_B T$. We see $V$ is an exponential and physically we require the minus solution
$
  V = V_0 e^(- kappa x)
$
where we define
$
  kappa eq sqrt((2 e^2 n_0)/(epsilon k_B T)) equiv cal(l)_D^(-1)
$
with $cal(l)_D$ being the Debye length. This is a measure of the distance before the colloid appears neutral due to screening. By changing $n_0$ we can tune $cal(l)_D$!

Consider two colloids. Assuming they are far apart each will have an ionic cloud. We want to describe what happens when we bring them together and their ionic clouds combine. We seek the concentration of ions between the colloids when this happens. Taking the midpoint to be $x = 0$ we look for symmetric solutions to the Boltzmann-Poisson equation. Someone has already done this and a family of solutions of the form
$
  tilde(V) = B ln cos beta x
$
exists. We require
$
  B = 2";  " beta = sqrt(2 pi cal(l)_B c_0)
$
with $c_0 = c_+ (x=0)$. This is easily shown by substitution.

We also have the boundary condition
$
  evaluated(dv(V, x))_(x=D) = - sigma/epsilon
$
which gives the relation
$
  2 beta tan D beta = 4 pi cal(l)_B sigma/epsilon
$
this is nice! Since $beta prop sqrt(c_0)$ this is a function of the concentration which we can then solve for. We also have
$
  c = c_0 e^(- tilde(V)) = c_0/(cos beta x)^2
$
so the concentration is minimal at $x=0$. But at $x = 0$ we  have $dd(tilde(V))\/dd(x) = 0$ so they feel no electric field! We define $dd(c, d: Delta) = c_0 - c_oo$. This leads to an osmotic pressure keeping the colloids apart.
