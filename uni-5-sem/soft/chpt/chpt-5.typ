//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

#pagebreak()
= Liquid crystals
Liquid crystals are treated as a liquid phase between solids and _usual_ liquids. A molecules in a solid have fixed orientation and position, in a liquid neither is fixed. A liquid crystal is defined as a liquid with fixed orientation. The primary purpose of this course is describing the phase transitions
$
  "solid" -> "liquid crystal" <-->^! "liquid"
$
We have two ways to classify liquid crystals. The first is based on how the liquid crystal $<-->$ liquid phase transition is induced. We care about thermotropic and lyotropic liquid crystals, which are based on temperature and concentration respectively. The second is based on which symmetries the liquid crystal has. We mostly care about nematic liquid crystals whose symmetry is purely orientational. Examples of other types could be smectic, cholestric, columnar etc.

== Q-tensor
We want a to quantify the _orderedness_ of a liquid crystal. To this end we assign a unit vector $hat(u)$ to each molecule.

We consider the average
$
  expval(hat(u)) = integral hat(u) psi(hat(u)) dd(Omega)
$
where $psi(hat(u))$ is an angular probability distribution
$
  integral psi(hat(u)) dd(Omega) = 1
$
However this is a bad measure of orderedness since we always have $expval(hat(u)) = 0$ due to symmetry.

We instead consider the variance
$
  expval(u_alpha u_beta) "with" alpha, beta in {x,y,z}
$
For an isotropic distribution e.g. a pure liquid we trivially have $psi(hat(u)) = (4 pi)^(-1)$. We also have
$
  expval(u_alpha u_beta) = 1/3 delta_(alpha beta)
$
since the directions are independent and $sum abs(u_i)^2 =^! 1$. Therefore we define the $Q$-tensor by
$
  Q_(alpha beta) equiv expval(u_alpha u_beta - 1/3 delta_(alpha beta))
$
so for an isotropic distribution $Q_(alpha beta) = 0$ which is nice! To make our lives easier we can write this in terms of the nematic director $hat(n)$. This is the average direction of molecules in the liquid crystal. We write
$
  Q_(alpha beta) = S [n_alpha n_beta - 1/3 delta_(alpha beta)]
$
where we define the _scalar order parameter_
$
  S equiv 3/2 expval(underbrace((hat(n) dot hat(u))^2, cos^2 theta)-1/3)
$

#proof[
  We have
  $
    n_alpha n_beta Q_(alpha beta) &= n_alpha n_beta expval(u_alpha u_beta - 1/3 delta_(alpha beta)) = expval((n_alpha u_alpha)^2 - 1/3 underbracket(n_alpha n_alpha, 1)) = expval((n_alpha u_alpha)^2 - 1/3)
  $
  and
  $
    n_alpha n_beta S [n_alpha n_beta - 1/3 delta_(alpha beta)] &= 3/2 expval((n_alpha u_alpha)^2 - 1/3) underbracket([(n_alpha n_alpha)^2 - 1/3 n_alpha n_alpha], 2\/3) = expval((n_alpha u_alpha)^2 - 1/3)
  $
  so we are done.
]

For an isotropic distribution $S = 0$ and for perfect alignment with $psi(hat(u)) = delta(0)$ we find $S = 1$. Somewhere in between we have the nematic phase with $S tilde 0.5$ and $phi prop "Gaussian"$. The minimal value is for $psi = delta(pi/2)$ where $S=-1/2$.

== Maier-Saupe theory
We want to find the free energy $F = U - T S$ in a liquid crystal.

We define a inter-molecular pair-potential by
$
  w (hat(u),hat(u)') equiv - underbracket(tilde(U), "sets scale") underbracket((hat(u) dot hat(u)')^2, "favors" #linebreak() "alignment")
$
this is the Maier-Saupe interaction. The total potential energy is then
$
  U = underbracket((N z)/2, "all pairs") underbracket(integral dd(hat(u), hat(u)') psi(hat(u)) psi(hat(u)') w(hat(u),hat(u)'), "weighted energy")
$
where $z$ is a measure of how many interactions we include.

We use $S = k_B ln Omega$ for the entropy. Where
$
  Omega = N!/(N_1 ! dots N_i ! dots N_M !)
$
we think of each $N_i$ being the number of molecules within some angular region. Then
$
  S &tilde.eq^"Stirling's" N k_B (- sum_(i=1)^M N_i/N ln N_i/N) \
  &tilde.eq^("in limit" N_i\/N -> psi(hat(u))) -N k_B integral dd(hat(u)) psi(hat(u)) ln psi(hat(u))
$
this is the _orientational entropy_!

Then we obtain
$
  F = - (tilde(U) N z)/2 integral dd(hat(u), hat(u)') psi(hat(u)) psi(hat(u)') (hat(u)dot hat(u)')^2 + N k_B T integral dd(hat(u)) psi(hat(u)) ln psi(hat(u))
$
We want to minimize this with respect to $psi(hat(u))$ under the constraint
$
  integral dd(hat(u)) psi(hat(u)) = 1
$
meaning we minimize
$
  cal(L) = F - lambda (integral dd(hat(u)) psi(hat(u))-1)
$
so we compute
$
  dv(, psi, d: delta) {integral dd(hat(u)) [(N tilde(U)z)/2 integral dd(hat(u)') (hat(u) dot hat(u)')^2 psi(hat(u)) psi(hat(u)') - N k_B T psi(hat(u)) ln psi(hat(u)) - lambda psi(hat(u))] - lambda} = 0
$
the only non-trivial part is
$
  dv(, psi, d: delta) [integral dd(hat(u)') (hat(u) dot hat(u)')^2 psi(hat(u)) psi(hat(u)')] &= 2 integral dd(hat(u)') (hat(u) dot hat(u)')^2 psi(hat(u)')
$
we obtain
$
  psi &= exp[(tilde(U) z)/(k_B T) underbracket(integral dd(hat(u)') (hat(u) dot hat(u)')^2 psi(hat(u)'), hat(u) "in mean field of" hat(u)')] underbracket(C, "not dependent on" psi)
$
this is nice! We would however like an expression in terms of the scalar order paramter. We define
$
  w_"mf" (hat(u)) equiv - tilde(U) z integral dd(hat(u)') (hat(u) dot hat(u)')^2 psi(hat(u)')
$
under the mean field assumption we can write
$
  integral dd(hat(u)') (hat(u) dot hat(u)')^2 psi(hat(u)') &= u_alpha u_beta underbracket(expval(u_alpha u_beta), Q_(alpha beta) + 1/3 delta_(alpha beta)) \
  &= S u_z^2 + 1 - S
$
so we find
$
  psi(u) = overbracket(C, "constants") exp[(tilde(U) z S)/(k_B T) u_z^2]
$
$C$ is determined by
$
  C^(-1) = integral exp[(tilde(U) z S)/(k_B T) u_z^2] dd(u)
$
By defintion we also have
$
  S & = integral 3/2 [u_z^2 - 1/3] psi(hat(u)) dd(hat(u)) \
  & =^(u_x u_y "cancel") C^(-1) integral_0^1 dd(u_z) 3/2 [u_z^2 -1/3] exp[(tilde(U) z S)/(k_B T) u_z^2] \
  & = 1/(integral_0^1 exp(tilde(S) u_z^2) dd(u_z)) integral_0^1 3/2 [u_z^2 - 1/3] exp(tilde(S) u_z^2) dd(u_z) equiv I\( tilde(S) \)
$
where we have defined $tilde(S)$ by
$
  S\(tilde(S)\) equiv (k_B T tilde(S))/(tilde(U) z)
$
and we are done! We have found a set of simultaneous equation which $S$ must satisfy simultaenously. We see a clear temperature dependence since $S prop T$. We can show that they always coincide at $tilde(S) = 0$ and that this is the only solution at high temperatures. At lower temperatures a second solution appears.

So at high temperatures $S tilde 0$ and as we lower the temperature we eventually reach metastability where a second solution appears. As we keep lowering the temperature we find $S$ smoothly increasing until $S tilde 1$. This can be shown by finding $F(S)$ and the minima. We will not do this. The phase transition is defined by when the two minima have equal energy $F(0) = F\(S_"transition" \)$. At low enough temperatures the minima at $S = 0$ eventually becomes a maxima.

== Landau-de Gennes theory
There is a simpler way to find the Maier-Saupe behavior.

By Landau we can write the free energy near any phase transition as a Taylor series in the order parameter. Applying this to liquid crystals we write followig de Gennes
$
  F = F_0 + A S + B/2 S^2 + C/3! S^3 + dots
$
we do not know what the $A dots$ or $S$ correspond to.

We now try to determine what the terms correspond to. Consider the linear term $A S$. This must vanish since
$
  evaluated(dv(F, S))_(S=0) = A =^! 0
$
due to $F$ having a minima at $S=0$. We keep the cubic term since $S=-S$ are very different physically. Assuming all temperature dependence is $prop S^2$ we end up with
$
  F = F_0 + a (T-T^*) S^2 - b S^3 + c S^4
$
where anything not mentioned is for convenience. This is the free energy in Landau-de Gennes theory.

Consider
$
  0 & = pdv(F, S) \
  0 & = [2 a (T-T^*) - 3 b S + 4 c S^2] S \
$
so it is minimal for $S_"I" = 0$. The other extrema are
$
  S_"N" = (3b)/(8 c) [1 plus.minus sqrt(1 - (32 a c)/(9 b^2) (T-T^*))]
$
when the determinant vanishes the second extrema begins to appear
$
  T^(**) = T^* + (9 b^2)/(32 a c)
$
For $T$ greater than $T^(**)$ the only physical solution is $S_"I" = 0$. For $T=T^(**)$ we have two solutions
$
  S_I = 0";  " S^(**)= (3 b)/(8 c)
$
Consider the curvature
$
  pdv(F, S, 2) = 2 a (T-T^*) - 6 b S + 12 c S^2
$
For $S_I = 0$
$
  pdv(F, S, 2) = 2 a (T-T^*) cases(> 0 "for " T > T^*, < 0 "for " T < T^*)
$
so it is a minima for $T$ greater than $T^*$ and becomes a maxima for $T$ less than $T^*$. For $S^(**)$ we find
$
  pdv(F, S, 2) & = 2 a (T^(**)-T^*) - (9 b^2)/( 16 c) \
               & = 0
$
so it is a saddle point.

At the nematic-isotropic phase transition the free energy of the two minima are equal $F_"I" = F_"N"$ meaning
$
  S^2 (a(T-T^*) - b S + c S^2) = 0
$
we require this only has two solutions. The trivial solution is $S_I = 0$. For the other consider
$
  S_plus.minus = (b plus.minus sqrt(b^2 - 4 a c(T-T^*)))/(2 c)
$
the determinant must vanish so we find
$
  T_"NI" = T^* + b^2/(4 a c)
$
which is simple. All of the above matches what we claimed for Maier-Saupe even though we barely did any physics!


== Onsager theory
We use $F = - T S$ with the entropy being translational through _depletion_ and rotational.

We consider spherical particles with radius $r$ and volume $v_p$ enclosed in a box of volume $V$. Then
$
  v_p = (4 pi r^3)/3
$
and the excluded volume due to any particle is
$
  v_e = (4 pi (2 r)^3)/3
$
Assuming we have one particle then the volume available to it is $V$. Adding a second particle the volume available to it is $V-v_e$. Similarly the volume available for the $n$th particle is $V - (n-1) v_e$. The number of states is then
$
  W & = 1/n! product_(i=1)^n underbracket(A, "some constant") V (1-((i-1) v_e)/V) = (A^n V^n)/n! product_(i=0)^(n-1) (1-(i v_e)/V)
$
taking the logarithm gives
$
  ln W tilde.eq^"Stirling" n ln A + n ln V - n ln n + n + sum_(i=1)^n ln(1-(i v_e)/V)
$
we compute the sum as
$
  sum_(i=1)^n ln(1- (i v_e)/V) tilde.eq sum_(i=1)^n (- (i v_e)/V) eq - (v_e n^2)/(2 V)
$
then
$
  S & = k_B ln W \
  & = k_B [n ln (A V) - n ln n+n - (v_e n^2)/(2 V)] equiv underbracket(S(n,V), "mixing entropy")
$
We compare this to the separated system
$
  S_"mix" (n,V) - S_"unmixed" = Delta S_"mix"
$
and then let $F = - T Delta S_"mix"$. We compute
$
  Delta S_"mix" & = S_"mix" (n,V) - [underbracket(S_"mix" (0,V), "no particles" #linebreak() = 0) + underbracket(S_"mix" (n,n v_p), "all clumped")] \
  &= k_B [n ln V/(n v_p) - v_e/2 (n^2/V - n/v_p)]
$
so we obtain
$
  F = k_B T [- n ln V/(n v_p) + v_e/2 (n^2/V-n/v_p)]
$
We let $c equiv n\/V$ and consider
$
  F/n & = k_B T [ln c v_p + (v_e c)/2 + "constant"]
$
so this free energy only depends on the concentration and excluded volume! With the second term being is positive.

For rod-shaped particles the above is complicated. This is because orientation now matters with respect to $v_e$. As an example for rectangles the difference is
$
  dd(a, d: Delta) = a_perp - a_parallel = (l-w)^2
$
For three-dimensional rods one can find
$
  ln c v_p & --> integral ln(4 pi v_p c) psi(hat(u)) dd(hat(u)) \
  (v_e c)/2 & --> integral c underbracket(beta(hat(u), hat(u)'), "alignment") psi(hat(u)) psi(hat(u)') \
  &underbracket(+ integral psi(hat(u)) ln psi(hat(u)) dd(hat(u)), "same as Maier-Saupe")
$
we essentially introduce integrals to weigh over orientations and introduce rotational entropy. Onsager used
$
  beta = 2 d l^2 sin theta
$
with $theta$ being the angle between two rods. This can be used for any liquid crystal shape with $beta$ depending on the shape. Using this free energy one finds an isotropic-nematic phase transition as some critical $c_"NI"$. Entropy alone leads to liquid crystals!

== Frank elastic energy
Liquid crystals have elasticity, if we try to distort some particles in the nematic phase then there will be a restoring force. But we are able to globally rotate a liquid crystal without using energy. This is quantified by the Frank elastic energy. This is defined to have certain properties: if $hat(n) -> - hat(n)$ then it should be invariant, it should depend on $nabla n$ and it should be invariant under global rotations. Under these condition the terms that survive are
$
  underbrace((div hat(n))^2, "3D splay") " and " (curl hat(n))^2 = underbrace((hat(n) dot (curl hat(n)))^2, "3D twist") + underbrace((hat(n) times (curl hat(n)))^2, "in-plane bend")
$
We can then write the Frank energy density as
$
  cal(f)_"elastic" equiv 1/2 [k_1 (div hat(n))^2 + k_2 (hat(n) dot (curl hat(n)))^2 + k_3 (hat(n) times (curl hat(n)))^2]
$
with the $k_i tilde^"usual" 10 "pN"$ being elastic constants setting the scale.

Liquid crystals can also interact with external fields. The energy of an electric field in a medium is
$
  "Energy" prop bold(E) dot bold(D) "with" bold(D) = epsilon_0 epsilon bold(E)
$
due to anisotopy liquid crystals have two different dielectric constants $epsilon_parallel$ and $epsilon_perp$. This anisotropy also changes the refractive index $n$, the magnetic susceptibility $chi$, the viscosity $eta$, etc.. So we can write
$
  bold(D) = epsilon_0 epsilon_parallel bold(E)_parallel + epsilon_0 epsilon_perp bold(E)_perp
$
we now define the dielectric anisotropy $Delta epsilon equiv epsilon_parallel - epsilon_perp$ then
$
  epsilon_0 epsilon_parallel bold(E)_parallel &= epsilon_0 (Delta epsilon + epsilon_perp) bold(E)_parallel \
  &= epsilon_0 Delta epsilon bold(E)_parallel + epsilon_0 epsilon_perp bold(E)_parallel
$
so we can write
$
  bold(D) &= epsilon_0 Delta epsilon bold(E)_parallel + epsilon_0 epsilon_perp (bold(E)_parallel + bold(E)_perp) \
  &= epsilon_0 Delta epsilon bold(E)_parallel + epsilon_0 epsilon_perp bold(E)
$
The energy density is then
$
  "Energy density" &= - 1/2 bold(E) dot bold(D) \
  &= -1/2 epsilon_0 Delta epsilon (bold(E) dot bold(E)_parallel) - 1/2 epsilon_0 epsilon_perp E^2
$
we can write $bold(E)_parallel = (bold(E) dot hat(n)) hat(n)$ giving
$
  cal(f)_"field" &= underbrace(- 1/2 epsilon_0 Delta epsilon (bold(E) dot hat(n))^2, "all" hat(n) "dependence") + g(E^2)
$
if $Delta epsilon > 0$ then this contribution is negative in the free energy if the electric field is parallel with the crystal and vanishing if they are perpendicular.

Now consider a liquid crystal trapped between two glass slides coated with a conductive material. Given the liquid crystal is aligned in parallel with the slides we now apply an electric field across the crystal. When this is done the liquid crystal would like to rotate and align with the field. This necessarily requires energy due to boundary effects. To combat this the slides are coated with a polymer layer. The liquid crystals preferentially want to be aligned with the polymers. With this the liquid crystal is essentially fixed near the slides, but in the bulk of the liquid crystal it will align. So you gain energy by alignment with the electric field, but require elastic energy to fix the boundaries. For this reason if the field is turned off then the liquid crystal will again be parallel with the slides, since energy is needed to be perpendicular to the slides. Now due to the changing refractive index $n$, applying an electric field makes the liquid crystal opaque---this is one way to make a display.

The light we use must be polarized. To do this a cross polarizor is used. Now let the director $hat(n)$ make an angle $beta$ with one of the polarizers---in-plane angle. In this case the intensity is given by
$
  I = I_0 sin^2 (2 beta) sin^2 underbrace(Delta n, "birefringent")
$
so it is maximal for $beta = pi\/4$ and minimal for $beta = 0, pi\/2$. The director $hat(n)$ also has an out-of-plane angle or tilt-angle $phi$ with respect to the cross polarizer. This is what determines the size of $Delta n$, changing the electric field changes $phi$. One can make display with $phi = 0$, but in this case the electric field needs to be applied differently.

Assuming $k_i = k$ (one constant approximation) we can write the total energy density as
$
  cal(f) = 1/2 k [(div hat(n))^2 + (curl hat(n))^2] - 1/2 epsilon_0 Delta epsilon (hat(n) dot bold(E))^2
$
we want to minimize this under the constraint $abs(hat(n))=1$. This is done by writing
$
  hat(n) = vec(sin theta cos phi, sin theta sin phi, cos theta)
$
We consider the case described above: a liquid crystal parallel to two glass slides over which we apply an electric field. This is a two-dimensional problem and the angle only depends on $z$. Then $(phi = pi\/2)$
$
  hat(n) = vec(0, sin theta(z), cos theta(z))
$
The $cal(f)_"elastic"$ then becomes
$
  cal(f)_"elastic" = 1/2 k (dv(theta, z))^2
$
and $cal(f)_"field"$ becomes
$
  cal(f)_"field" = -1/2 epsilon_0 Delta epsilon E^2 cos^2 theta
$
so we obtain
$
  cal(f) = 1/2 k (dv(theta, z))^2 - 1/2 epsilon_0 Delta epsilon E^2 cos^2 theta
$
We minimize this using the Euler-Lagrange equation
$
  dv(cal(f), theta, d: delta) - dv(, z) dv(cal(f), dot(theta), d: delta) = 0
$
giving
$
  k dv(theta, z, 2) - epsilon_0 Delta epsilon E^2 cos theta sin theta = 0
$
or introducing
$
  xi^2 equiv k/(epsilon_0 Delta epsilon E^2)
$
we can write
$
  xi^2 dv(theta, z, 2) - cos theta sin theta = 0
$
this is solved by
$
  tan(theta/2 - pi/4) = exp(plus.minus z/xi)
$
so $xi$ sets the length scale. The above derivation ignores the energy loss in differing from the boundary condition, and completely ignores ions in the liquid crystal. We also completely ignore the isotropic-nematic phase transition---this assumes we are far from the transition.


