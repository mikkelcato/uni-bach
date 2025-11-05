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
  sigma_(x x) = - P + F_x/A_x => sigma_(x x) - 1/2 (sigma_(y y)+sigma_(z z)) = F_x/A_x = sigma_N &= G(lambda^2 - lambda^(-1)) \
  &tilde.eq^(lambda = 1 + epsilon) 3 G epsilon = Y epsilon
$
note $sigma_N = lambda dd(f)\/dd(lambda)$---with $Y$ being Young's modulus. For an isotropic compression $V -> V'$ we have
$
  sigma_(x x) = sigma_(y y) = sigma_(z z) = -P +F/A = -P + Delta P
$
and we define
$
  Delta P = - K dd(V, d: Delta)/V = - 3 K epsilon_V
$
we'll treat $K -> oo$, but $G$ is important since soft matter can be deformed not compressed.

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
We want to derive $G$ for rubber---polymer (random walk) $+$ crosslinks gives a network (liquid $->$ elastic solid), and acts like a material we can deform. After the deformation strands have $R'_alpha = E_(alpha beta) R_beta$. To see how the material responds we find the free energy difference
$
  f(E_(alpha beta)) = tilde(f) (E_(alpha beta)) - tilde(f)_0
$
we denote the probability of a strand having $bold(R)$ with $N$ by $psi(bold(R), N)$. Then
$
  f(E_(alpha beta)) &= rho_s integral dd(N) integral dd(bold(R)) psi(bold(R), N) {F'_"strand"-F_"strand"} \ &= rho_s integral dd(N) integral dd(bold(R)) psi(bold(R), N) {3/2 (k_B T)/(b^2 N) bold(R') dot bold(R') -3/2 (k_B T)/(b^2 N) bold(R) dot bold(R)}
$
with $rho_s$ being the density of strands---we can then say $psi(bold(R), N) = P(bold(R),N) Phi_0 (N)$, so the probability of generating a given polymer multiplied by the probability of then making a given strand after cross-linking. We obtain
$
  f(E_(alpha beta)) &= rho_s 3/2 (k_B T)/b^2 integral dd(N) (Phi_0(N))/N integral dd(bold(R)) P(bold(R),N) (bold(R)' dot bold(R)'- bold(R) dot bold(R))
$
the $bold(R)$ integrand can be rewritten
$
  bold(R)' dot bold(R)' - bold(R) dot bold(R) &= E_(alpha beta) R_beta E_(alpha gamma) R_gamma - delta_(beta gamma) R_beta R_gamma \
  &= (E_(alpha beta) E_(alpha gamma) - delta_(beta gamma)) R_beta R_gamma \
  &= (C_(beta gamma) - delta_(beta gamma)) R_beta R_gamma \
  &= 2 cal(E)_(beta gamma) R_beta R_gamma
$
so the integral becomes
$
  integral dd(bold(R)) (dots) &= integral dd(x, y, z) (3/(2 pi b^2 N))^(3\/2) exp(- (3(x^2+y^2+z^2))/(2b^2 N)) (E_(alpha beta) E_(alpha gamma) - delta_(beta gamma)) R_beta R_gamma
$
with an implicit sum over $alpha, beta, gamma$, let $beta = x = gamma$ then
$
  I_(beta, gamma) &= sum_(beta, gamma) integral dd(x, y, z) (3/(2 pi b^2 N))^(3\/2) exp(- (3(x^2+y^2+z^2))/(2 b^2 N)) sum_alpha (E_(alpha beta) E_(alpha gamma) - delta_(beta gamma)) R_beta R_gamma \
  I_(x,x) &= [integral dd(y) (3/(2 pi b^2 N))^(1\/2) exp(- (3 y^2)/(2 b^2 N))]^2 \ & times integral dd(x) (3/(2 pi b^2 N))^(1\/2) exp(- (3 x^2)/(2b^2 N)) (E_(alpha x) E_(alpha x) - delta_(x x)) x^2 \
  &= integral dd(x) (3/(2 pi b^2 N))^(1\/2) exp(- (3 x^2)/(2b^2 N)) (E_(alpha x) E_(alpha x) - delta_(x x)) x^2 \
  &= (E_(alpha x) E_(alpha x) - delta_(x x)) integral dd(x) x^2 (3/(2 pi b^2 N))^(1\/2) exp(- (3 x^2)/(2 b^2 N)) \
  &= (E_(alpha x) E_(alpha x) - delta_(x x)) (b^2 N)/3
$
with $beta = y = gamma$ and $beta = z = gamma$ we get symmetric results. Take $beta = x, gamma = y$ then
$
  I_(x,y) &= integral dd(z) (dots) integral dd(y) (3/(2 pi b^2 N))^(1\/2) exp(- (3 y^2)/(2 b^2 N)) y integral dd(x) dots \
  &=^(integral dd(y)) 0
$
similarly for all other off-diagonal terms. So
$
  f(E_(alpha beta)) &= rho_s 3/2 (k_B T)/b^2 integral dd(N) (Phi_0 (N))/N {(E_(alpha beta) E_(alpha gamma) - delta_(beta gamma)) delta_(beta gamma) (b^2 N)/3} \
  &= (rho_s k_B T)/2 integral dd(N) Phi_0 (N) {E_(alpha beta) E_(alpha gamma)-delta_(beta gamma)) delta_(beta gamma)} \
  &= (rho_s k_B T)/2 (E_(alpha beta) E_(alpha gamma) - delta_(beta gamma)) delta_(beta gamma) integral dd(N) Phi_0 (N) \
  &= (rho_s k_B T)/2 (E_(alpha beta) E_(alpha gamma) - delta_(beta gamma) delta_(beta gamma)) \
  &= (k_B T rho_s)/2 (E_(alpha beta) E_(alpha beta) - 3)
$
for shear
$
  E_(alpha beta) = mat(1, gamma, 0; 0, 1, 0; 0, 0, 1)
$
we get
$
  f(E_(alpha beta)) = (k_B T rho_s gamma^2)/2
$
We apply a work when we shear
$
  dd(W) = sigma_( x y ) dd(gamma) = dd(f) => sigma_(x y) = k_B T rho_s gamma => G = k_B T rho_s
$
which is what we wanted to find---which is the obvious energy density one would guess.

\* percolation transition (liquid $->$ elastic solid $->$ solid), Kuhn theory assumes every strand is an entropic strand but in real life polymers are more constrained since they can't cross$-> G = G_K + G_E$.

#pagebreak()
= Surfaces and Interfaces
Surfactants---head (hydrophillic) and tail (hydrophobic) group, Colloids---small molecules in solution.

== Surface tension
Consider a water-air interface, or any liquid-air interface: water molecules in the bulk will be more content since it is energetically favorable---due to molecular interactions. Water molecules at the interface are less content since this is less energetically favorable, this leads to a surface tension $gamma = dd(E)\/dd(A)$. We can also consider a thin film with surface tension $gamma$, now we try to pull on this film increasing its length by $dd(x)$ using some force $F$. We obtain
$
  dd(W) = F dd(x) = 2 gamma a dd(x)
$
with $a$ being the width of our film and the factor two appearing since our film has two sides---so $gamma = F\/2 a$, and the surface tension counteracts our applied force. Similarly pressing liquid out of a syringe into air using some $P_0 + Delta P$, with $P_0$ being ambient pressure, will lead to an increased surface area, so $dd(W) = Delta P dd(V) = gamma dd(A)$. If the formed droplet is assumed spherical then
$
  dd(V) = 4 pi r^2 dd(r)",  " dd(A) = 8 pi r dd(r)
$
giving the Laplace pressure
$
  Delta P = (2 gamma)/r
$
for real droplets this is given by the more general Young-Laplace equation---but we won't cover this.

=== Wetting
Consider some surface and a droplet on the surface. In this case we have three different interfaces and three corresponding surface tensions: $gamma$ (liquid-vapor), $gamma_"SV"$ (solid-vapor) and $gamma_"SL"$ (solid-liquid). We want to know whether or not this droplet spreads. If the droplet has initial area $A$ then the energy before placing the droplet is $E_"before" = gamma_"SV" A$. Assuming the curvature ($A$ with $gamma$ same as $A$ with $gamma_"SL"$) of the drop is negligible we likewise have $E_"after" = (gamma + gamma_"SL") A$. Then we define the spreading coefficient
$
  gamma_"S" = (E_"before"-E_"after")/A = gamma_"SV" - gamma - gamma_"SL"
$
which is the energy difference per unit area---if $gamma_"S">0$ then the energy after is smaller and the droplet will spread as much as possible, if $gamma_"S" < 0$ the energy increases and the droplet stays.

Another measure that describes what happens is the contact angle $theta$. We consider what happens when our drop gets extended by $dd(x)$, this lengthens the liquid-vapor interface by $cos theta dd(x)$ and changes some solid-vapor to solid-liquid, and we obtain
$
  dd(E) = cos theta dd(x) gamma + (gamma_"SL"-gamma_"SV") dd(x)
$
in equilibrium
$
  dv(E, x) = 0 => gamma cos theta + gamma_"SL" - gamma_"SV" = 0
$
giving Young's equation
$
  cos theta = (gamma_"SV"-gamma_"SL")/gamma = (gamma_"S"+gamma)/gamma
$
what happens when we introduce gravity? The energy of the surface area is given by
$
  G_A = 4 pi gamma r^2
$
and the gravitational potential energy is
$
  E_"pot" = m g h = (4 pi rho g)/3 r^4
$
obviously $r^4 > r^2$ so for larger drops gravity wins and it will spread, but for small drops surface tension wins and it will keep its shape. They are equal at the scale
$
  G_A = E_"pot" => r^* tilde sqrt((gamma)/(rho g))
$
giving something like $r^*_"water" tilde 1.4 "mm"$.

Now consider placing a drop on a surface and letting it spread until it reaches equilibrium with height $h$. We assume that the volume $V$ of the drop is known, as well as all the surface tensions. Then ignoring curvature and using $h\/2$ as the height of the center of mass we obtain
$
  E_"pot" = V rho g h/2
$
and using $A = V\/h$
$
  G_A = (gamma + gamma_"SL") V/h - gamma_"SV" V/h = - gamma_"S" V/h
$
at equilibrium
$
  dv(E_"tot", h) = 0 => (V rho g)/2 + (gamma_S V)/(h^*)^2 = 0
$
giving
$
  h^* = sqrt((-2gamma_S)/(rho g))
$
relating this to $r^*$ can be done, and it simply gives $h^* tilde.eq r^* theta$.

=== Capillary effects
What happens when inserting a tube into some liquid? We denote the diameter of our tube by $2 a$, and the liquid rises within the tube to a height $h$ above the liquid outside. The volume within the tube is then $V_"tube"=pi a^2 h$, with center of mass $h\/2$, so
$
  E_"pot" = (pi a^2 rho g h^2)/2
$
and
$
  G_A = 2 pi a h (gamma_"SL"-gamma_"SV") = - 2 pi a h gamma cos theta
$
since it happens outside and inside, this then gives
$
  h^* = (2 gamma cos theta)/(a rho g) = (2 (r^*)^2 cos theta)/a
$

== Thermodynamics
=== The grand potential
The grand potential is a thermodynamic quantity depending on volume, temperature and the chemical potential. Importantly volume is extensive while temperature and the chemical potential are intensive so we can write
$
  G(V,T,mu) = V g(T,mu)
$
with $g(T,mu)$ being the grand potential volume density. Similarly
$
  G_A (A,T,mu) = A g_A (T,mu)
$
with $g_A (T,mu)$ being the grand potential surface density. One can show
$
  g = G/V = - P (T,mu)
$
#proof[
  We do a Legendre transform to get the free energy (Gibbs $->$ Helmholz)
  $
    G(V,T,mu_i) = F(V,T,N_i) - sum_i N_i mu_i
  $
  note
  $
    dd(F) = - P dd(V) - S dd(T) + sum mu_i dd(N_i)
  $
  then
  $
    dd(G) & = dd(F) - sum dd(N_i) mu_i - sum_i N_i dd(mu_i) \
          & = - P dd(V) - S dd(T) - sum_i N_i dd(mu_i)
  $
  but we also have
  $
    dd(G) = g dd(V) + V pdv(g, T) dd(T) + V sum_i pdv(g, mu_i) dd(mu_i)
  $
  with the obvious things being held constant. It immediately follows by comparison that
  $
        g & = - P \
        S & = V pdv(P, T) \
    N_i/V & = pdv(P, mu_i)
  $
]

Similarly we can write
$
  g_A = gamma(T, mu)
$
with a positive sign since surface tension is contractile.
#proof[
  For an area we have
  $
    G_A (A, T, mu_i) = A g_A (T, mu_i)
  $
  with
  $
    dd(G_A) & = gamma dd(A) -S_A dd(T) - sum_i^A N_i^A dd(mu_i) \
    dd(G_A) & = g_A dd(A) + A pdv(g_A, T) dd(T) + A sum_i pdv(g_A, mu_i) dd(mu_i)
  $
  everything follows as before giving
  $
        g_A & = gamma \
        S_A & = - A pdv(gamma, T) \
    N_i^A/A & = - pdv(gamma, mu_i)
  $
]
The sign difference stem from the surface wanting to contract $+gamma$, while both volume want to expand $-P$.

Using $G = - P V(T,mu_i)$ we can quickly find the Gibbs-Duhem relation
$
  dd(G) = - V dd(P) - P dd(V)
$
recall we had
$
  dd(G) = - S dd(T) - P dd(V) - sum_i N_i dd(mu_i)
$
subtracting these give
$
  sum_i N_i dd(mu_i) & = - S dd(T) + V dd(P)
$
we can do the same for $G_A = A gamma(T, mu_i)$ giving
$
  sum_i N_i^A dd(mu_i) = -S_A dd(T) - A dd(gamma)
$

=== Surfactants
We now consider a system with some air(I)-liquid(II) interface---in this case molecules can move freely (in principle) between the three (open) systems. We assume we know $V_I$, $V_(I I)$ and the interface area $A$. Then assuming thermal equilibrium we have
$
  G = G^I (V^I, T, mu^I) + G^(I I) (V^(I I), T, mu^(I I)) + G_A (A, T, mu^A)
$
We'd like to know how to determine $N_i^A$, to this end we define the surface excess
$
  Gamma_i = N_i^A/A
$
In the air phase the density of liquid molecules is obviously low, then around the interface we have some continuous gradient, but we want to define some zero---this is done by counting liquid in the air phase, and air in the liquid phase. This is an equivalent definition of the surface excess
$
  Gamma_0 = integral_(-oo)^0 dd(z) (n_0 (z)-n_(I I)) + integral_0^oo dd(z) (n_0 (z) - n_I)
$
with $i = 0$, and $n(z)$ being the number density. We define the zero by requiring this vanishes: $Gamma_0 = 0$, this is where the surface is.

If we add some surfactants (any molecule with hydrophillic head-group and hydrophobic tail-group) then $n_s (z)$ would be highly concentrated around the interface, since it is energetically favorable, then
$
  Gamma_s = integral_(-oo)^0 dd(z) (n_s (z) - n_(I I)) + integral_0^oo dd(z) (n_s (z) - n_I)
$
and $Gamma_s > 0$, an anti-surfactant would have $Gamma_(! s) < 0$, e.g. ions since it is energetically favorable to not be close to air. So we can write $N_i^A = A Gamma_i$.

Due to entropy there'll still be surfactants in the bulk, since they want to maximize their translational entropy. And at sufficient saturation, when the interface becomes _crowded_ then additional surfactants will start forming micelles, which begin to act like a single molecule---so entropy hates this, but it is energetically (lower enthalpy) favorable since the tails aren't exposed to liquid (water). Also consider some surface which is divided in two parts using some movable membrane. We add surfactant to one part of the surface leading to some surface tension $gamma$, while the part without surfactants has surface tension $gamma_0$, with $gamma < gamma_0$ since the surfactants disrupt the surface. This means that the surface with $gamma$ wants to expand and gives rise to some pressure.

We want to see what happens when we add surfactants to some liquid. So consider some surface with surfactant density $n$, this can be treated as a lattice since the surfactant molecules take up physical space, denote the fraction of occupied states by $theta$. We now imagine a surfactant molecule from the bulk hitting the surface, if the site is empty then it is absorbed, and if the site is filled then it bounces back. We denote the rate of _hitting_ the surface by $v_a$ (absorbance), explicitly $v_d = K_a (1- theta)$ with $K_a$ being the number of attempts per unit time. Similarly we have the rate of leaving the surface (desorbing) $v_d = K_d theta$ with $K_d$ being the number of attempts per unit time. In equilibrium these are equal since the rate of absorbtion and desorbtion should equal, $v_d = v_a$. This gives
$
  theta_"eq" = K/(1+K)",   " K=K_a/K_d
$
but the fraction of occupied sites is also
$
  theta_"eq" = Gamma/Gamma_s
$
and we can write
$
  K = n/n_s
$
with $s$ denoting saturation, so $Gamma_s$ is when every site is occupied and $n_s$ is the corresponding bulk density. Then
$
  Gamma = n/(n+n_s) Gamma_s
$
this is also called the Langmuir isotherm. Then what happens to the surface tension? We'd like to find $gamma(n)$ at constant $T$, from the Gibbs-Duhem relation we then obtain
$
              N^A dd(mu) & = -A dd(gamma) \
        => dv(gamma, mu) & = - Gamma = - (n Gamma_s)/(n+n_s) \
  dv(gamma, n) dv(n, mu) & = - (n Gamma_s)/(n+n_s)
$
now we need some $mu(n)$, since we want to be rid of the chemical potential. We can write
$
  mu(n) = mu_0 + k_B T ln n
$

#proof[
  Consider a gas with $n = N\/V$, we assume each molecule has some zero-energy $mu_0$. The partition function is
  $
    Z = (Z_"trans" Z_"part")^N/N!
  $
  each molecule has
  $
    E = 1/2 m v^2 + mu_0
  $
  the first term enters in $Z_"trans"$ and the second enters in $Z_"part"$. We know
  $
    Z_"trans" = V/lambda_"th"^3
  $
  and we can easily obtain
  $
    Z_"part" = exp(- mu_0/(k_B T))
  $
  then the free energy is
  $
    F & = - k_B T ln Z \
      & = - k_B T (- N ln N + N - (mu_0 N)/(k_B T) + N ln (V)/lambda_"th"^3)
  $
  by definition
  $
    mu = pdv(F, N) & = dots = mu_0 + k_B T ln (N lambda_"th"^3)/V \
                   & = mu_0 + k_B T ln(n lambda_"th"^3) \
                   & =^(lambda_"th" = "const") mu_0 + k_B T ln n
  $
  with the second term just being translational entropy.
]

Then we easily obtain
$
  dv(mu, n) = (k_B T)/n => dv(n, mu) = n/(k_B T)
$
so
$
  dv(gamma, n) = - (k_B T Gamma_s)/(n+n_s)
$
this can then be integrated to get $gamma(n)$ giving
$
  gamma(n) = gamma_0 - Gamma_s k_B T ln(1+n/n_s)
$
at low $n$ this is heavily dominated by $gamma_0$ as we'd expect and as $n$ increases $gamma$ lowers.

=== Micelles
Now we'd like to see what happens if we keep adding surfactants, since we'll eventually saturate the surface completely. Consider $m$ surfactant molecules in a box, either they'll stay apart or they'll form a micelle---this is a equilibrium process, meaning the chemical potential for both is equivalent. For a single surfactant we have
$
  mu_1 = mu_1^0 + k_B T ln n_1
$
and for the micelle we have
$
  mu_m = mu_m^0 + k_B T ln n_m
$
in equilibrium
$
  m mu_1 = mu_m => k_B T ln n_m/n_1^m & = m mu_1^0-mu_m^0 \
                            n_m/n_1^m & = exp((m mu_1^0-mu_m^0)/(k_B T))
$
this is essentially just the Boltzmann factor. We define
$
  n_c^(1-m) = exp((m mu_1^0-mu_m^0)/(k_B T)) => n_c = exp((m mu_1^0-mu_m^0)/((1-m)k_B T))
$
we can then write
$
  n_m/n^m_1 = n_c^(1-m) => n_m & = n_1^(m) n_c^(1-m) = (n_1/n_c)^m n_c \
$
now we use $n = n_1 + m n_m$ so
$
      n & = n_1 + m (n_1/n_c)^m n_c \
  n/n_c & = n_1/n_c + m (n_1/n_c)^m
$
we consider the case $n_1 < n_c$, then since $m$ is somewhat large we get
$
  n tilde.eq n_1
$
so below this critical value $n_c$ we have essentially no micelles. Consider also $n_1 > n_c$, then for a similar reason the second term now dominates
$
  n/n_c & tilde.eq m(n_1/n_c)^m \
        & => n_1 tilde.eq n_c (n/(m n_c))^(1\/m) tilde.eq n_c
$
since the $m$'th root for large $m$ just gives $1$. So $n_1 tilde.eq n_c$ with this we can obtain $n_m tilde.eq n_c$ aswell. Therefore $n_c$ is essentially the maximal concentration, for this reason we can't really take the limit $n_1\/n_c > 1$, since we just find $n_1\/n_c = 1$. We can also find
$
  n_m = n_c = n/m => n/m = n_c
$
by using
$
  n/n_c tilde.eq m (n_1/n_c)^m = m
$
so at $n_c$ micelles emerge, and all surfactant contribute to making micelles---after $n_1 = n_c$ then $n_1$ can't increase further forcing all surfactants to become micelles.
