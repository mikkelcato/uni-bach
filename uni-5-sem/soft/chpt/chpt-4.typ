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
  &tilde.eq^("in limit" N_i\/N -> psi(u)) - N k_B integral dd(u) psi(u) ln psi(u)
$
so the free energy is
$
  F = - (tilde(U) N z)/2 integral dd(u, u') psi(u) psi(u') (hat(u)dot hat(u)')^2 - N k_B T integral dd(u) psi(u) ln psi(u)
$
we seek to minimize this with respect to $psi(u)$ under the constraint
$
  integral dd(u) psi(u) = 1
$
so we use
$
  cal(L) = F - lambda (integral dd(u) psi(u)-1)
$
