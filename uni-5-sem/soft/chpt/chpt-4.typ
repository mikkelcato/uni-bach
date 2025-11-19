//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

#pagebreak()
= Liquid crystals
Liquid crystals can be though of a new liquid phase between solids (crystals) and usual liquids. In a solid a material has both fixed (or uniform) orientation and position while a liquid has neither. In a liquid crystal we just have fixed orientation. We'd like to be able to describe the transition from solid $->$ liquid crystal $->$ liquid.

There are two main types of classification for liquid crystals. The first is based on how the transition liquid crystal $<->$ liquid is induced---we care about thermotropic (temperature) and lyotropic (concentration) liquid crystals. The second is based on the symmetry of the liquid crystal phase. The one we care about the most are nematic liquid crystals whose symmetry is purely orientational---one can also have other types e.g. smectic, cholestric, columnar, etc.

== Q-tensor
We want a number to quantify the _ordered-ness_ of a liquid crystal. To do this we assign a unit vector $hat(u)$ to each molecule. Consider the average
$
  expval(hat(u)) = integral_Omega u psi(u) dd(Omega)
$
where $psi(u)$ an angular probability distribution
$
  integral_Omega psi(u) dd(Omega) = 1
$
however we always have $expval(u) = 0$ due to apolar symmetry. Instead consider the second moment
$
  expval(u_alpha u_beta) "with" alpha, beta in {x,y,z}
$
for an isotropic distribution---as in a pure liquid---we have $psi(u) = (4 pi)^(-1)$, and it is easy to see that
$
  expval(u_alpha u_beta) = 1/3 delta_(alpha beta)
$
since the directions are independent and $sum abs(u_i)^2 = 1$. This leads us to construct the $Q$-tensor:
$
  Q_(alpha beta) = expval(u_alpha u_beta - 1/3 delta_(alpha beta))
$
so for an isotropic distribution $Q_(alpha beta) = 0$ by definition. To make our lives easier we write this in terms of the nematic director $hat(n)$ which is the average direction of the molecules. We do
$
  n_alpha n_beta Q_(alpha beta) = n_alpha n_beta expval(u_alpha u_beta - 1/3 delta_(alpha beta)) => Q_(alpha beta) = S (n_alpha n_beta - 1/3 delta_(alpha beta))
$
where we define the order parameter
$
  S = 3/2 expval(underbrace((hat(n) dot hat(u))^2, cos^2 theta)-1/3)
$
for an isotropic distribution $S = 0$ and for perfect alignment with $psi(u) = delta(0)$ we find $(hat(n) dot hat(u))^2 = 1 => S = 1$. Somewhere in between we have the nematic phase with $S tilde 0.5$ and $phi prop "Gaussian"$. The least possible value is for $psi = delta(pi/2)$ where $S=-1/2$.

== Maier-Saupe theory
We want to find the free energy $F = U - T S$.

We define the inter-molecular potential
$
  w (hat(u),hat(u)') equiv - underbrace(tilde(U), "sets scale") underbrace((hat(u) dot hat(u)')^2, "favors" #linebreak() "alignment")
$
then the total potential energy is
$
  U = underbrace((N z)/2, "all pairs") underbrace(integral dd(u, u') psi(u) psi(u') w(hat(u),hat(u)'), "weighted energy")
$
the parameter $z$ is a measure of how many interactions we include. For the entropy we use $ S = k_B ln Omega $ with
$
  Omega = N!/(N_1 ! N_2 ! dots N_M !)
$
where we think of each $1, 2, dots, M$ being a region of orientations---this is analogous to the usual translational entropy. Then we obtain
$
  S &tilde.eq^"Stirling's" k_B [N ln N - N - sum_(i=1)^M N_i (ln N_i - 1)] \
  &tilde.eq N k_B [ln N - 1 - sum_(i=1)^M N_i/N (ln N_i - 1)] \
  &tilde.eq N k_B (- sum_(i=1)^M N_i/N ln N_i/N) \
  &tilde.eq^("in limit" N_i\/N -> psi(u)) -N k_B integral dd(u) psi(u) ln psi(u)
$
so the free energy is
$
  F = - (tilde(U) N z)/2 integral dd(u, u') psi(u) psi(u') (hat(u)dot hat(u)')^2 + N k_B T integral dd(u) psi(u) ln psi(u)
$
we seek to minimize this with respect to $psi(u)$ under the constraint
$
  integral dd(u) psi(u) = 1
$
so we minimize
$
  cal(L) = F - lambda (integral dd(u) psi(u)-1)
$
with respect to $psi$:
$
  dv(cal(L), psi, d: delta)= 0
$
so
$
  dv(, psi, d: delta) {integral dd(u) [(N tilde(U)z)/2 integral dd(u') (u dot u')^2 psi(u) psi(u') - N k_B T psi(u) ln psi(u) - lambda psi(u)] - lambda} = 0
$
the seperate parts evaluate to
$
        dv(, psi, d: delta) [lambda] & = 0 \
  dv(, psi, d: delta) [- lambda psi] & = - lambda \
    dv(, psi, d: delta) [psi ln psi] & = ln psi + 1
$
$
  dv(, psi, d: delta) [integral dd(u') (u dot u')^2 psi(u) psi(u')] &= integral dd(u') (u dot u')^2 2 psi(u')
$
requiring the integrand vanishes then gives
$
  0 &= N tilde(U) z integral dd(u') (u dot u')^2 psi(u') - N k_B T ln (psi+1) - lambda \
  ln (psi + 1) &= 1/(N k_B T) [N tilde(U) z integral dd(u') (u dot u')^2 psi(u') - lambda] \
  psi &= exp[(tilde(U) z)/(k_B T) underbrace(integral dd(u') (u dot u')^2 psi(u'), u "in mean field of" u')] underbrace(C, "not dependent on" psi)
$
we define
$
  w_"mf" (u) equiv - tilde(U) z integral dd(u') (u dot u')^2 psi(u')
$
since it is a mean-field we can write (assuming $expval(u'_alpha u'_beta) = expval(u_alpha u_beta)$)
$
  integral dd(u') (u dot u')^2 psi(u') = u_alpha u_beta underbrace(expval(u_alpha u_beta), Q_(alpha beta) + 1/3 delta_(alpha beta))
$
this can be rewritten using the scalar order parameter
$
  S = 3/2 expval((u dot n)^2 - 1/3)
$
consider
$
  u_alpha u_beta expval(u_alpha u_beta) &= u_x^2 (Q_(x x)+1/3) + u_y^2 (Q_(y y)+1/3) + u_z^2 (Q_(z z)+1/3) \
  &= u_x^2 (-1/3 S +1/3) + u_y^2 (-1/3 S + 1/3) + u_z^2 (2/3 S +1/3) \
  &=^(abs(u)=1) 1/3 [(1-u_z^2) (-S+1) + u_z^2 (2 S + 1)] \
  &= S u_z^2 + 1 - S
$
so we can write
$
  psi(u) = overbrace(C, "constants") exp[(tilde(U) z S)/(k_B T) u_z^2]
$
we can determine $C$ by the constraint
$
  integral psi(u) dd(u) = 1 => C^(-1) = integral exp[(tilde(U) z S)/(k_B T) u_z^2] dd(u)
$
now consider
$
  S & = integral 3/2 (u_z^2 - 1/3) psi(u) dd(u) \
    & = (integral dd(u_z) 3/2 (u_z^2 -1/3) exp[(tilde(U) z S)/(k_B T) u_z^2]) 1/C
$
here we integrate over $u_z in [0,1]$ since the $u_x, u_y$ directions would just cancel. Now define
$
  tilde(S) equiv (S tilde(U) z)/(k_B T) => I(tilde(S)) = (integral_0^1 3/2 (u_z^2-1/3) exp(tilde(S) u_z^2) dd(u_z))/(integral_0^1 exp(tilde(S) u_z^2) dd(u_z))
$
now $S$ has to fulfill both of these so we try to find an intersection between the functions $f_1$ and $f_2$ defined by
$
  f_1 equiv I (tilde(S)) "and" f_2 equiv S(tilde(S)) = (tilde(S) k_B T)/(tilde(U) z)
$
here there is clear temperature dependence since it determines the slope of $f_2$. They always coincide at $tilde(S)=0$ but at high temperatures this is the only solution, while at lower temperatures we can find two solutions.

We find that at high temperatures $S tilde 0$ as we lower the temperature we reach metastability (essentially a non-continuous region) before $S$ smoothly increases to $S tilde 1$. In principle we can also go back and find the free energy $F(S)$---at high temperature this exhibits one minima (at $S = 0$) and at lower temperatures it exhibits two minima (at $S = 0$ and at $S_"transition"$) (if they are equal minima) with a barrier between them, if we keep lowering the temperature the minima at $S=0$ will eventually vanish with the other minima becoming deeper.

== Landau-de Gennes theory
There is a smarter way to do the above.

By Landau we can always write the free energy near a phase transition as a Taylor expansion in the order parameter---de Gennes applied this to liquid crystals. So we consider
$
  F = F_0 + A S + B/2 S^2 + C/3! S^3 + dots
$
here we do not know what $A, B, dots$ or $S$ is. We will now try to determine what the terms correspond to.

Consider the linear term $A S$, this must vanish, since:
$
  evaluated(dv(F, S))_(S=0) = A =^! 0
$
due to $F$ having a minima at $S=0$. We also know the cubic term is necessary, else $F$ would be even meaning $S=-S$ would correspond to the same physics. We assume all temperature dependence is $prop S^2$, so we write
$
  F = F_0 + a (T-T^*) S^2 - b S^3 + c S^4
$
where the $-$ sign is for convenience (since we want maxima and minima)---this is the free energy in Landau-de Gennes theory.

Consider
$
  0 & = pdv(F, S) \
  0 & = (2 a (T-T^*) - 3 b S + 4 c S^2) S \
$
so it is minimal for $S_"I" = 0$. The other physical extrema is
$
  S_"N" = (3b)/(8 c) [1 plus sqrt(1 - (32 a c)/(9 b^2) (T-T^*))]
$
this extrema appears at
$
  1 - (32 a c)/(9 b^2) (T-T^*) = 0 => T_C = T^* + (9 b^2)/(32 a c)
$
Before this $T > T_C$ and the only physical solution is $S_"I" = 0$. For $T=T_C$ we have solutions
$
  S = {S_"I" = 0, S_C= (3 b)/(8 c)}
$
To see if it these are minima or maxima consider
$
  pdv(F, S, 2) = 2 a (T-T^*) - 6 b S + 12 c S^2
$
For $S_"I" = 0$
$
  pdv(F, S, 2) = 2 a (T-T^*) cases(> 0 "for " T > T^*, < 0 "for " T < T^*)
$
so it is a minima for $T > T^*$ and becomes a maxima for $T < T^*$---meaning at $T < T^*$ it is very unlikely that the system is in the isotropic phase. For $S_C = 3 b \/ 8 c$
$
  pdv(F, S, 2) & = 2 a (T_C-T^*) - (9 b^2)/( 16 c) \
               & = 2 a (T^* + (9 b^2)/(32 a c) - T^*) - (9 b^2)/(16 c) \
               & = 0
$
so it is a saddle point.

At the nematic-isotropic transition the free energy of the two minima are the same $F_"I" = F_"N"$ so
$
  0 & = a(T-T^*) S^2 - b S^3 + c S^4 \
    & = (a(T-T^*) - b S + c S^2 ) S^2
$
so it is minimal for $S_"I" = 0$. For the quadratic
$
  a(T-T^*) &= b S_"N" - c S_"N"^2 \
  &= (12 b^2)/(32 c) + (12 b^2)/(32 c) sqrt(1-(32 a c)/(9 b^2) (T-T^*))] - (9 b^2)/(64 c) [1 + sqrt(1-(32 a c)/(9 b^2) (T-T^*))]^2
$
the second term
$
  -[dots]^2 &= -(9 b^2)/(64 c) [2 - (32 a c)/(9 b^2) (T-T^*) + 2 sqrt(1-(32 a c)/(9 b^2) (T-T^*))] \
  &= -(9 b^2)/(32 c) + a/2 (T-T^*) - (9 b^2)/(32 c) sqrt(1 - (32 a c)/(9 b^2) (T-T^*))
$
so we obtain
$
  [(16 a c)/(3 b^2) (T - T^*) -1 ]^2 & = 1 - (32 a c)/(9 b^2) (T-T^*) \
                              T_"NI" & = T^* + b^2/(4 a c)
$
in this case
$
  S_"NI" & = (3 b)/(8 c) [1 + sqrt(1 - (32 a c)/(9 b^2) (b^2/(4 a c)))] \
         & = (3 b)/(8 c) (1 + 1/3) \
         & = (b)/(2 c)
$

We could also say
$
  S^2 (a(T-T^*) - b S + c S^2) = 0
$
and require we only have two solutions having $F-F_0 = 0$, with one being $S = 0$. Consider then
$
  S_plus.minus = (b plus.minus sqrt(b^2 - 4 a c(T-T^*)))/(2 c)
$
for this to have one solution (meaning $F-F_0$ has two zeroes) the determinant must vanish giving
$
  0 & = b^2 - 4 a c (T-T^*) => T = T^* + b^2/(4 a c)
$
which is the same as above.

At the beginning we could have also written
$
  F = F_0 + A Tr Q + B/2 Tr Q^2 + C/3! Tr Q^3 + D/4! Tr Q^4
$
since $Tr Q = 0$ this term vanishes as before. The other terms are $prop S^2$, $prop S^3$, etc. just with different constants.

We can also apply Landau-de Gennes to other phase transitions. Consider for example a smectic-nematic transition. We define some order parameter describing this and denote it by $Psi$, one can show that $Psi = - Psi$ meaning we must have
$
  F = F_0 + a (T-T^*) Psi^2 + b Psi^4
$
then one could proceed as before.

== Onsager theory
Here we consider $F = - T S$ and rod-shaped particles, with the entropy being given by translational entropy (through depletion) and rotational entropy.
