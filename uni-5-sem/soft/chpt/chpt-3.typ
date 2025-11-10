//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Colloids
The second part of this course covers colloids and liquid crystals, starting with colloids.

It is very hard to define what a colloid is since much of soft matter can act in very different ways depending on many things. A simple qualitative definition could be particles in a solution that are small $10 "nm" -> 1 " "mu"m"$. For particles of this size thermal motion becomes relevant, so colloids are heavily affected by stuff like Brownian motion. For a spherical particle of $1 " "mu"m"$ with some typical density the gravitation energy $E_g tilde 20 dot 10^(-21) "J"$ while $E_"th" tilde k_B T tilde 4.5 dot 10^(-21) "J"$ at room temperature. These sizes are comparable with $E_g tilde 5 E_"th"$, and we'd call particles of this size colloids.

As physicists we are interested in how these colloids interact. We have distinguish between attractive and repulsive forces. For attractive forces we focus on dispersion (through Van der Waals forces) and depletion forces. For repulsive forces we focus on grafting, charges and the hard sphere interaction. Colloids are always attractive, these forces can't be avoided, so to avoid particles clumping and falling out of solution we need repulsive forces.

== Interactions
=== Dispersion forces
These consist of Van der Waals forces due to induced transient dipoles. The derivation of the VdW potential is done in other courses so we just posit the result
$
  V = - C (a_0/r)^6
$
with $C$ and $a_0$ being constants. As should be clear this is very short range, so we'd expect it not to really matter but as we'll now show it does matter.

We consider the situation in which two particles seperated by some $h << R$, this lets us approximate the particles as being parallel _walls_. We place the first at $z=0$ and the second at $z=h$, so for $z < 0$ we are inside the first particle while for $z > h$ we are inside the second. We consider all small volume elements $dd(V)$ and consider the energy between every one such element. If we have $n$ molecules per $dd(V)$ then
$
  W = - integral dd(bold(r)_1) n integral dd(bold(r)_2) n (C a_0^6)/abs(bold(r)_1-bold(r)_2)^6
$
we work in cylindrical coordinates so
$
  abs(bold(r)_1-bold(r)_2)^2 = (rho_2^2 + (z_1-z_2)^2)^3
$
then
$
  W &= - n^2 C a_0^6 integral_V_1 dd(bold(r)_1) integral_h^oo dd(z_2) integral_0^oo dd(rho_2) (2 pi rho_2)/(rho_2^2 + (z_1-z_2)^2)^3 \
  &=^(tilde(rho)=rho_2^2) -n^2 C a_0^6 integral_V_1 dd(bold(r)_1) integral_h^oo dd(z_2) integral_0^oo (pi dd(tilde(rho)))/(tilde(rho) + (z_1-z_2)^2)^3 \
  &= - (n^2 C a_0^6 pi)/2 integral_V_1 dd(bold(r)_1) integral_h^oo dd(z_2) 1/(z_1 - z_2)^4 \
  &= (n^2 C a_0^6 pi)/6 integral_V_1 dd(bold(r)_1) 1/(z_1-h)^3 \
  &= (n^2 C a_0^6 pi)/6 A integral_(-oo)^0 dd(z_1) 1/(z_1-h)^3 \
  &= - (n^2 C a_0^6 pi)/(12 h^2) A
$
which is typically written as
$
  W/A = w(h) = - A_H/(12 pi h^2)
$
where we define the Hamaker constant
$
  A_H = n^2 C a_0^6 pi^2
$
importantly we find a decay proportional to $h^(-2)$ instead of $r^(-6)$. We want to generalize this to two spheres instead of treating them as planes. To do this we let $h$ be the minimal distance between two spheres and $tilde(h) (rho)$---remember our $z$-axis goes into the spheres along $h$. One can show easily that
$
  tilde(h) (rho) & = h + 2(R-sqrt(R^2-rho^2)) \
                 & = h + 2 R (1 - sqrt(1 - rho^2/R^2)) \
                 & tilde.equiv h + rho^2/R
$
Now all values of $rho$ define an annulus with $dd(A) = 2 pi rho dd(rho)$ so to get the total interaction energy across the entire surface we calculate
$
  U(tilde(h)) & = integral dd(A) w(tilde(h)) \
              & = integral_0^R dd(rho) w (tilde(h)) 2 pi rho
$
we essentially sweep the surface. From before $rho = 0 -> tilde(h) = h$, $rho=R -> tilde(h) = h + R$ and $ tilde(h) = h + rho^2/R -> dd(tilde(h)) = (2 rho dd(rho))/R $ so
$
  U(tilde(h)) & = pi R integral_h^(h+R) w (tilde(h)) dd(tilde(h)) \
              & = - R/12 A_H [- 1/(h+R) + 1/h] \
              & =^"Derjaguin" - R/12 A_H/h
$
now the decay is just proportional to $h^(-1)$ even though the potential decays as $r^(-6)$ so we can't ignore it. In the last step we use the Derjaguin approximation to say that $1\/(h + R) -> 0$ since $h << R$, to see it explicitly
$
  -1/(h+R) + 1/h & = (-h+h+R)/(h(h+R)) \
                 & = - R/(h(h+R)) \
                 & tilde.equiv^(h << R) - R/(h R) \
                 & = - 1/h
$

=== Depletion forces
Depletion forces are forces caused by other things in solution with the colloids, this could for example be polymers---we have something big (colloids) with something small where $n_"small" >> n_"big"$ usually.

Approximating our colloids as cubes of length $D$ and the small things is solution as spheres with size $d$, then we obviously have two length scales. If two cubes are within a distance of $r < d$ then the small things don't fit between our colloids. This creates a density imbalance, leading to osmotic pressure forcing the colloids together. This can be seen by the ideal gas law
$
  p = N/V k_B T => p prop N/V
$
so $dd(p, d: Delta) = n_"small" k_B T$. We can also consider this in terms of entropy. The small things want to maximize their entropy, if we have a situation with $r < d$ then the small things can't be within this volume, and their entropy will be lower---to maximize their available volume, and thereby entropy, the colloids are forced together. We could quantify this by a potential of the form
$
  U = cases(0"," #h(50pt) & r > d, - dd(p, d: Delta) A r"," & r < d)
$
in the case of cubes---by Derjaguin this can be extended to spherical colloids.

This can also be quantified in terms of excluded volumes. For spherical colloids of radius $R$ and small things with radius $rho$ then the excluded volume is
$
  V_"excluded" = (4 pi)/3 (R + rho)^3
$
for a single colloid. If we have two colloids and they are seperated then the total excluded volume is $2 V_"excluded"$, but if they are touching or close together, then the effective excluded volume is $< 2 V_"excluded"$. To maximize entropy we want to minimize this effective excluded volume as much as possible.

=== Hydrophobic forces
Water molecules can form hydrogen bonds with eachother, but these are weak and are able to break at energies $tilde k_B T$, so any individual water molecule can make many possible bonds, leading to higher entropy. If we were to place a hydrophobic particle in water then the water molecules form a _cage_ (clathrate) around it, and they get stuck essentially. This naturally lowers the entropy since the stuck water molecules can no longer bond freely. For this reason it is beneficial for these particles to clump, to decrease the size of the clathrate, thereby leading to an attractive force. This also happens for large particles, such as colloids, even if a full clathrate isn't formed.

=== Hard sphere interaction
The hard sphere interaction is very repulsive, and essentially ensures that colloids don't overlap, which is unphysical. The naive hard sphere potential is
$
  V = cases(oo"," #h(10pt) & r <= 2R, 0"," & "otherwise")
$
a more realistic potential wouldn't be $oo$ but would be rapidly increasing---since we could have soft colloids. This interaction is always present.

=== Specific interactions
- \*phase transitions

- \*patchy colloids

- \*DNA coating---fine tuned grafting

- \*colloids with facets---none spherical colloids $tilde.equiv$ patchy colloids---entropic crystals

- \*active colloids---e.g. bacteria

=== Grafting
Grafting is the process of coating colloids with polymers (bonded through a covalent bond). At sufficient polymer density the individual polymers get _squished_ together, and we get a polymer brush, we'd like to know how tall this brush is---this also depends on the solvent.

We denote the number of polymer chains per unit area (the grafting density) by $Gamma_p$, this defines a characteristic length $cal(l) = Gamma_p^(1\/2)$. We imagine each polymer occupies a cylinder of radius $cal(l)$, being a random walk within it, we want to know the height of these cylinders $h_p$. We again define the volume fraction $phi = N v_c\/cal(l)^2 h_p$ with $N$ being the number of monomers in a polymer and $v_c$ being the volume of one such monomer, so $N v_c$ is the minimal polymer volume while $cal(l)^2 h_p$ is the maximal polymer volume. We can write this guy as
$
  phi = (v_c Gamma_p N)/h_p
$
Now we'll try to write the free energy per area. The elastic free energy is
$
  u_E (h_p) = underbrace((3 k_B T h_p^2)/(2 N b^2), "single polymer") Gamma_p
$
with $b$ being the monomer length. We also need the free energy of mixing per area
$
  F_"mix"/A = (f_"mix" V)/A = (f_"mix" A h_p)/A = f_"mix" h_p
$
where we've previously found
$
  f_"mix" = (k_B T)/v_c ((1-phi) ln(1-phi) + underbrace(chi (1-phi) phi, "polymer"<-->"solvent"))
$
assuming $phi$ is small, we find after computation
$
  U(h_p) = (3 k_B T)/(2 N b^2) h_p^2 Gamma_p + h_p (k_B T)/v_c [(chi-1) (v_c Gamma_p N)/h_p + (1/2 - chi) ((v_c Gamma_p N)/h_p^2)^2]
$
to find $h_p^"eq"$ we take the derivative
$
  h_p^"eq" = N root(3, Gamma_p (v_c b^2)/3 (1/2-chi))
$
this depends on the length of our polymers, the density of polymers and their interaction with the solvent---as we'd expect.

Now we imagine two colloids (in the Derjaguin approximation) with their own brush being brought close together. We denote the distance by $d$ and define $h = d\/2$, if $d < 2 h_p^"eq"$ then we begin compressing the polymer brushes. Then
$
  w(h) &= 2 (U(h) - U(h_p^"eq")) \
  & tilde.eq^"Taylor" 2 [U(h_p^"eq") + overbrace(evaluated(dv(U, h))_(h=h_p^"eq"), 0) (h-h_p^"eq") + 1/2 evaluated(dv(U, h, 2))_(h=h_p^"eq") (h-h_p^"eq")^2 + dots - U(h_p^"eq")] \
  &tilde.eq evaluated(dv(U, h, 2))_(h=h_p^"eq") (h-h_p^"eq")^2
$
we can find
$
  dv(U, h, 2) = (3 k_B T Gamma_p)/(N b^2) + (2 k_B T v_c Gamma_p^2 N)/(h^3) (1/2-chi)
$
giving
$
  evaluated(dv(U, h, 2))_(h=h_p^"eq") = (9 k_B T)/(N b^2) Gamma_p
$
so we obtain
$
  w(h) tilde.eq (9 k_B T Gamma_p)/(N b^2) (h-h_p^"eq")^2
$
which is a nice elastic interaction depending on the grafting density---this is positive meaning it is repulsive.

=== Charges
We now consider our colloids have surface charge---typically these are negatively charged, and usually all colloids have similar charge. We recall the Coulomb force
$
  bold(f) = (q_1 q_2)/(4 pi epsilon) bold(hat(r))_(1 2)/abs(bold(r)_(1 2))^2
$
here $epsilon$ is important, and electrical interactions in solution are generally weak since $epsilon_"water" tilde 100 epsilon_0$. One might expect that ions in solution stick to the colloids, but this is not advantageous due to entropy, instead a ionic cloud will develop around any colloid. At the surface we have (considering a flat wall)
$
  E = sigma/epsilon
$
and in the bulk we have
$
  nabla dot bold(E) = rho/epsilon;"  "bold(E) = - grad V
$
We have charges which generate some $V$, and we imagine that charges move under the mean field (or potential) due to all other charges, this makes the charge density change meaning the potential will change---this is called DVLO---we'd like some steady state or consistent state which is generally done by iteration.

In one dimension we have
$
  dv(E, x) = - dv(V, x, 2) = rho/epsilon
$
where $rho = c_"ions" q_"ions"$ the concentration is given by Boltzmann statistics
$
  c = c_0 exp(- (q V)/(k_B T))
$
To see this take a system with some potential difference then there'll be a flux
$
  j_E = (q E)/zeta c
$
where $zeta$ is some friction coefficient, similarly we have diffusive flux given by Fick's law
$
  j_D = - D dv(c, x)
$
by Einstein $D= k_B T zeta^(-1)$ so
$
  j_"tot" = j_E + j_D = D (- dv(c, x) + (q E c)/(k_B T)) =^! 0
$
giving
$
  integral 1/c dv(c, x) dd(x) & = integral (q E)/(k_B T) dd(x) \
                     ln c/c_0 & = - (q V)/(k_B T) \
                            c & = c_0 exp(- (q V)/(k_B T))
$
so $V$ leads to a concentration difference or a concentration difference leads to $V$. Now we obtain (for $q = e$) the Boltzmann-Poisson equation
$
         dv(V, x, 2) & = - e/epsilon c_0 exp(- (e V)/(k_B T)) \
  dv(tilde(V), x, 2) & = - 4 pi cal(l)_B c_0 e^(-tilde(V))
$
where we've defined
$
  tilde(V) equiv (e V)/(k_B T)
$
and the Bjerrum length
$
  cal(l)_B equiv e^2/(4 pi epsilon k_B T)
$
which is the characteristic length scale at which the entropic and electrical contributions are similar.

\*absence

We now consider a colloid in a neutral solution with many ions. Sufficiently far away from the colloid we have $n_(0+) = n_(0-) = n_0$ and generally we have $rho = e (n_+ - n_-)$ with the density given by the Boltzmann distribution
$
  n_plus.minus = n_0 exp(minus.plus (e V)/(k_B T))
$
so
$
  rho & = e n_0 (exp(- (e V)/(k_B T)) - exp(+ (e V)/(k_B T))) \
      & = - 2 e n_0 sinh (e V)/(k_B T)
$
giving
$
  dv(V, x, 2) & = - rho/epsilon \
              & = (2 e n_0)/epsilon sinh (e V)/(k_B T) \
              & tilde.eq^(e V\/k_B T << 1) (2 e n_0)/epsilon (e V)/(k_B T)
$
where we've used the Debye-Hueckel approximation $e V\/k_B T << 1$. So $V$ is an exponential, and physically we require the minus solution
$
  V = V_0 exp(- kappa x)
$
substituting gives
$
  kappa = sqrt((2 e^2 n_0)/(epsilon k_B T)) = cal(l)_D^(-1)
$
with $cal(l)_D$ being the Debye length, it is a measure of the distance before the colloid looks neutral (due to screening)---importantly it can be tuned by changing $n_0$ since $cal(l)_D prop n_0^(-1\/2)$.

Now we consider what happens if we have two colloids. If they are sufficiently far from each other then both will have their own ionic clouds. We'd like to know what happens when they come close and their ionic clouds begin combining. This leads to osmotic pressure since the positive charges between the colloids loose entropy---this pressure is given by $Pi = dd(n, d: Delta) k_B T = dd(c_+, d: Delta) k_B T$. We'd like to know the concentration of ions between the two colloids. In this case we have a symmetry about $x=0$, and it turns out the Boltzmann-Poisson equation has a family of solutions of the form
$
  tilde(V) = B ln(cos beta x)
$
which is symmetric around $x = 0$, but we require
$
     B & = 2 \
  beta & = sqrt(2 pi cal(l)_B c_0)
$
with $c_0 = c_+ (x=0)$. We also have the boundary condition, that at $x = D$, the surface of the colloid we require
$
  evaluated(dv(V, x))_(x=D) = - sigma/epsilon
$
after computing this gives
$
  2 beta tan D beta = 4 pi cal(l)_B sigma/epsilon
$
since $beta prop sqrt(c_0)$ this is a function of the concentration, which we can then solve for.

#proof[(problem 3)][
  We have
  $
    tilde(V) = B ln cos beta x
  $
  so the LHS
  $
       dv(tilde(V), x) & = - beta B tan beta x \
    dv(tilde(V), x, 2) & = - beta^2 B sec^2 beta x
  $
  and the RHS is
  $
    -4 pi cal(l)_B c_0 exp(- B ln cos beta x) &= - (4 pi cal(l)_B c_0)/(cos beta x)^B
  $
  they are the same if
  $
    -beta^2 B sec^2 beta x = - (4 pi cal(l)_B c_0)/(cos beta x)^B
  $
  meaning
  $
    B = 2;"  " beta = sqrt(2 pi cal(l)_B c_0)
  $
  we have boundary conditions $V(0)=0 => tilde(V)(0)=0$ which is trivial and
  $
    evaluated(dv(V, x))_(x=D) = - sigma/epsilon => dv(tilde(V), x) = - (e sigma)/(k_B T epsilon)
  $
  with
  $
    cal(l)_B = e^2/(4 pi epsilon k_B T)
  $
  this becomes
  $
    evaluated(dv(tilde(V), x))_(x=D) & = - (4 pi cal(l)_B sigma)/epsilon \
                 - 2 beta tan beta D & = - (4 pi cal(l)_B sigma)/epsilon \
                   2 beta tan beta D & = 4 pi cal(l)_B sigma/epsilon
  $
  which is what we had before.
]

We also have by back-substitution
$
  c = c_0 e^(- tilde(V)) = c_0/(cos beta x)^2
$
so the distribution is minimal at $x=0$. But we also have $dd(tilde(V))\/dd(x) = 0$ at $x=0$ so they feel no electric field, but an excess of charges lead to an osmotic pressure, since there is no electric field to hold them in place. Similarly at $x -> oo$ we have $E = 0$. We define $c_0 - c_oo = dd(c, d: Delta)$ which leads to an osmotic pressure.

\*colloids of opposite charge $->$ odd solution

== Problems on Colloids
=== P4.6
We consider a gas of rigid spherical molecules of radius $a$ confined in a volume $V$. If two rigid spheres of radius $R$ are introduced to the system with a gap $h$, each molecule is excluded by some region $R + a$. Let $dd(V(h), d: Delta)$ be the voume of this region.

a. By considering the free energy of the particle show that the average force acting on the sphere is given by
$
  F & = - n V k_B T pdv(ln(V-dd(V(h), d: Delta)), h) \
    & = n k_B T pdv(dd(V(h), d: Delta), h)
$
with $n$ being the number density of molecules.

For a gas with $N$ molecules the entropy goes as $(V-dd(V(h), d: Delta))^N$, so the free energy is
$
  A = - N k_B T ln(V-dd(V(h), d: Delta)) + "const"
$
the force is then
$
  F = pdv(A, h) & = -N k_B T pdv(ln(V-dd(V(h), d: Delta)), h) \
                & = (n V k_B T)/(V-dd(V(h), d: Delta)) pdv(dd(V(h), d: Delta), h)
$
for $V >> dd(V, d: Delta)$ we get the result.

b. Calculate $F(h)$.

We have
$
  dd(V(h), d: Delta) = (8 pi)/3 (R+a)^3 - V_"overlap" (d)
$
where $d = 2 R +h$ (center-center distance), if $h >= 2 a$ or $d >= 2 (R+a)$ then theres no overlap and
$
  dd(V(h), d: Delta) = (8 pi)/3 (R+a)^3
$
if $d <= 2 (R+a)$ ($h <= 2 a$) then we have overlap. The overlap can be bisected, we want the volume of one such cap. The height of the cap is $(R+a)-d\/2 = a - h\/2 = h^*$. The volume of a spherical cap is
$
  V = (pi h^*^2 (3 (R+a)- h^*))/3 => V_"overlap" = V/2
$
giving
$
  V_"overlap" (h) = pi/12 (2 a- h)^2 (6 R + 4 a + h)
$
so we've found
$
  dd(V(h), d: Delta) = cases((8pi)/3 (R+a)^3 - pi/12 (2a-h)^2 (6 R + 4 a + h)";  " h <= 2 a, (8 pi)/3 (R+a)^3";  " h >= 2 a)
$
we take the derivative
$
  pdv(dd(V(h), d: Delta), h) &= cases(-pi/12 [(2 a-h)^2 -2(6 R + 4 a + h)(2a - h)]";  " h <= 2a, 0";  " h>= 2a) \
  &= cases(-pi/12 [(2 a -h) (-6 a -3h - 12 R)]";  " h <= 2a, 0";  " h>= 2a) \
  &= cases(pi/4 [4 ( a^2 + 2R a)-h^2 - 4 R h]";  " h <= 2a, 0";  " h>= 2a) \
  &= cases(pi/4 [4 (R+a)^2 - (2 R + h)^2]";  " h <= 2a, 0";  " h>= 2a)
$
giving the force
$
  F(h) = cases((pi n k_B T)/4 [4 (R+a)^2 - (2 R + h)^2]";  " h <= 2 a, 0";  " h>= 2 a)
$

c. Take the limit $R >> h$.

We expand the force
$
  F(h) & = (pi n k_B T)/4 (4 a^2 + 8 R a - h^2 - 4 R h) \
       & tilde.eq n k_B T pi R (2 a - h)
$
we know
$
  w_"dep" (h) = n k_B T (h - 2 a)
$
for two parallel surfaces, since $R >> h$ we use the Derjaguin approximation
$
  U(h) &= pi R integral_h^oo omega(x) dd(x) => F(h) = -pdv(U(h), h) = - pi R omega_"dep" (h) = n k_B T pi R (2 a - h)
$
which is the same as above.


=== Debye-Hueckel
We know
$
  sigma = integral_0^oo rho (x) dd(x)
$
a. Use Debye-Hueckel to find a relation between $sigma$ and $cal(l)_D$.

From above we use
$
  rho & = - 2 e n_0 sinh (e V)/(k_B T) \
      & tilde.eq^"DB" - (2 e n_0 e V)/(k_B T) \
      & = - (2 e^2 n_0 V_0)/(k_B T) e^(- kappa x)
$
integrating gives with $kappa = cal(l)_D^(-1)$
$
  sigma & = (2 e^2 n_0 V_0)/(k_B T) cal(l)_D
$
and
$
  cal(l)_D^2 = (epsilon k_B T)/(2 e^2 n_0) => epsilon/cal(l)_D^2 = (2 e^2 n_0)/(k_B T)
$
so
$
  sigma = (epsilon V_0)/cal(l)_D
$


