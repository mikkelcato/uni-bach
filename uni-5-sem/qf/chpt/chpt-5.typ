//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Approximation Methods
== Nondegenerate time-independent pertubation theory
The central problem in quantum mechanics is solving the Schrödinger equation to get the spectrum of a given Hamiltonian, this is hard to do. To make our lives easier we write
$
  H = H_0 + V
$
where $H_0$ is a Hamiltonian for which we know the eigenkets and energies
$
  H_0 ket(n^((0))) = E_n^((0)) ket(n^((0)))
$
and $V$ is interpreted as a pertubation.

We want to solve
$
  (H_0 + V) ket(n) = E_n ket(n)
$
to track the pertubation we write
$
  (H_0 + lambda V) ket(n) = E_n ket(n)
$
where $lambda in [0,1]$ is a continuous parameter we use to keep track---in the limit $lambda -> 0$ we simply get back the free theory. What we really want is the difference in energy defined by
$
  Delta_n equiv E_n - E_n^((0))
$
so we have
$ E_n ket(n) = (Delta_n + E_n^((0))) ket(n) $
then we can rewrite the Schrödinger equation as
$
  (E_n^((0))-H_0) ket(n) = (lambda V - Delta_n) ket(n)
$
to get some insight we try to act with $bra(n^((0)))$ on this guy giving
$
  braket(n^((0)), (lambda V - Delta_n), n) = 0
$
this is nice because it tells us that $(lambda V - Delta_n) ket(n)$ has no component along $ket(n^((0)))$. We consider the operator
$
  phi_n equiv 1 - ketbra(n^((0))) = sum_(k eq.not n) ketbra(k^((0)))
$
which is just the difference between $ket(n^((0)))$ and $ket(n)$ and $phi_n ket(n^((0))) = 0$ by construction. Now we can write
$
  ket(n) = sum_k ket(k^((0))) braket(k^((0)), n) = braket(n^((0)), n) ket(n^((0))) + phi_n ket(n)
$
then
$
  (E_n^((0))-H_0) ket(n) &= underbrace((lambda V - Delta_n) ket(n^((0))), (E_n^((0))-H_0) ket(n^((0)))=0) braket(n^((0)), n) + (lambda V - Delta_n) phi_n ket(n) \
  &= (lambda V - Delta_n) phi_n ket(n)
$
we can then write
$
  ket(n) = 1/(E_n^((0))-H_0) (lambda V-Delta_n) phi_n ket(n)
$
this inversion is possible since we first act with $phi_n$ so we avoid problems like
$
  1/(E_n^((0))-E_n^((0))) = 1/0 -> "ill-defined"
$
now in the limit $lambda -> 0$ we have $Delta_n -> 0$ and $ket(n) -> ket(n^((0)))$ but in our case $ket(n) -> 0$ for this reason write
$
  ket(n) = 1/(E_n^((0))-H_0) (lambda V - Delta_n) phi_n ket(n) + c_n (lambda) ket(n^((0)))
$
notice that if we do $(E_n^((0))-H_0) ket(n)$ then we recover the previous since $(E_n^((0))-H_0) ket(n^((0))) = 0$ meaning this term acts like an arbitrary constant. In the same limit we see $c_n (lambda) -> 1$ and generally $c_n (lambda) = braket(n^((0)), n)$ which we fix by convention
$
  braket(n^((0)), n) = c_n (lambda) = 1
$
so $ket(n)$ is not normalized $braket(n, n) eq.not 1$. We finally obtain
$
  ket(n) = ket(n^((0))) + 1/(E_n^((0))-H_0) phi_n (lambda V- Delta_n) ket(n)
$
where we've switched the order of $phi_n$ and $(lambda V - Delta_n)$, since they commute. We can also easily find
$
  Delta_n = lambda braket(n^((0)), V, n)
$
this is the basic setup. So far everything is exact.

=== Pertubative solution
To solve this we assume that we can expand $ket(n)$ and $Delta_n$ in power series of $lambda$
$
  ket(n) &= ket(n^((0))) + lambda ket(n^((1))) + lambda^2 ket(n^((2))) + dots \
  Delta_n &= underbrace(Delta_n^((0)), 0) + lambda Delta_n^((1)) + lambda^2 Delta_n^((2)) + dots
$
now we can get expressions for the various $Delta_n^((i))$ by using
$
  Delta_n = lambda braket(n^((0)), V, n)
$
and comparing powers of $lambda$. Doing this we find
$
  Delta_n^((1)) & = braket(n^((0)), V, n^((0))) \
  Delta_n^((2)) & = braket(n^((0)), V, n^((1))) \
                & dots \
  Delta_n^((N)) & = braket(n^((0)), V, n^((N-1)))
$
we do the same with
$
  ket(n) = ket(n^((0))) + 1/(E_n^((0))-H_0) phi_n (lambda V-Delta_n) ket(n)
$
by substituting the expansions for $Delta_n$ and $ket(n)$ and comparing powers (and using $phi_n ket(n^((0)))=0$) we obtain $ket(n^((1)))$
$
  ket(n^((1))) = 1/(E_n^((0))-H_0) phi_n V ket(n^((0)))
$
which gives us $Delta_n^((2))$
$
  Delta_n^((2)) & = braket(n^((0)), V phi_n/(E_n^((0))-H_0) V, n^((0)))
$
which gives us $ket(n^((2)))$
$
  ket(n^((2))) = phi_n/(E_n^((0))-H_0) V phi_n/(E_n^((0))-H_0) V ket(n^((0))) - phi_n/(E_n^((0))-H_0) braket(n^((0)), V, n^((0))) phi_n/(E_n^((0))-H_0)V ket(n^((0)))
$
and we can keep going like this.

In practice these can be rewritten,
$
  Delta_n &= lambda Delta_1 + lambda^2 Delta_2 + dots \
  &= lambda V_(n n) + lambda^2 braket(n^((0)), V 1/(E_n^((0))-H_0) phi_n V, n^((0))) + dots \
  &= lambda V_(n n) + lambda^2 sum_(k eq.not n) braket(n^((0)), V 1/(E_n^((0))-H_0), k^((0))) braket(k^((0)), V, n^((0))) + dots \
  &= lambda V_(n n) + lambda^2 sum_(k eq.not n) abs(V_(n k))^2/(E_n^((0))-E_k^((0))) + dots
$
where we've defined
$
  V_(n k) equiv braket(n^((0)), V, k^((0)))
$
in a similar manner
$
  ket(n) &= ket(n^((0))) + lambda ket(n^((1))) + dots \
  &= ket(n^((0))) + lambda sum_(k eq.not n) V_(k n)/(E_n^((0))-E_k^((0))) ket(k^((0))) + dots
$
which are easy to work with.

Now it is obvious why this fails in the degenerate case, since then $E_n^((0)) = E_k^((0))$ making both of these ill-defined. This is why we want degenerate pertubation theory.

#pagebreak()
== Degenerate time-independent pertubation theory
To solve the problem above we seek a diagonal basis for $V$ such that $V_(n k) = 0$ for $n eq.not k$, thereby circumventing degeneracy.

Let $H = H_0 + lambda V$ then the unperturbed eigenkets satisfy
$
  H_0 ket(m^((0))) & = E_m^((0)) ket(m^((0)))
$
and the eigenkets of the full Hamiltonian satisfy
$
  H ket(l) = E_l ket(l)
$
as $lambda -> 0$ we recover $ket(l) = ket(l^((0)))$ with
$
  H_0 ket(l^((0))) = E_m^((0)) ket(l^((0)))
$
due to degeneracy we may have $ket(l^((0))) eq.not ket(m^((0)))$, but we can write
$
  ket(l^((0))) = sum_(m in D) c_m ket(m^((0)))
$
where we take $D$ to be a $g$-fold degenerate subspace. We want to diagonalize $V$ in this subspace, i.e. by going from $ket(m^((0))) -> ket(l^((0)))$. Assume
$
  V_(n k) = braket(m_n^((0)), V, m_k^((0)))
$
then diagonalizing and finding the eigenkets and eigenvalues give $ket(l^((0)))$ and $Delta_l^((1))$ respectively.

To see this consider
$
  H ket(l) = E_l ket(l) = (H_0 + lambda V) ket(l) => (E_l - H_0 - lambda V) ket(l) = 0
$
solve iteratively with $ket(l) tilde.eq ket(l^((0)))$:
$
  (E_l - H_0 - lambda V) ket(l^((0))) = (E_l - E_l^((0))-lambda V) ket(l^((0))) = (Delta_l^((1))-lambda V) ket(l^((0))) = 0
$
so the first order correction is
$
  Delta_l^((1)) = braket(l^((0)), V, l^((0)))
$

=== Linear Stark effect
We consider the hydrogen atom with energies
$
  E_n = - (alpha e^2 hbar c)/(2a_0) 1/n^2
$
for $n > 1$ these are degenerate with degeneracy $n^2$. Consider the $n=2$ case, here we have four states with the same energy: $l=0, m=0$ or $l=1, m=plus.minus 1, 0$. Now we apply a uniform electric field in the $z$-direction:
$
  V = -e z abs(bold(E))
$
we must diagonalize this in the $(n l m)$ basis. To simplify things we know $braket(alpha, x, beta) = 0$ if $ket(alpha)$ and $ket(beta)$ share parity. The parity of $ket(n l m)$ is $(-1)^l$ so for non-zero matrix elements $braket(2 l' m', z, 2 l m)$ we must have $l'=0$ and $l=1$ or $l'=1$ and $l=0$. We can do further simplifications since $z$ acts like a spherical tensor:
$
  T_q^((k)) = Y_(l = k)^(m = q) (hat(n)) "with" hat(n) -> bold(V)
$
given $k = l = 1$ and $hat(n)_z = z/r -> V_z$:
$
  Y_l^0 = sqrt(3/(4 pi)) cos theta = sqrt(3/(4 pi)) z/r -> T_0^((1)) = sqrt(3/(4 pi)) V_z
$
now by Winger-Eckart
$
  braket(alpha' j' m', T_q^((k)), alpha j m) = 0 "unless" m' = q +m
$
so we need $m' = m = 0$. Giving
$
  V = mat(0, braket(2s, V, 2p\, m=0), 0, 0; braket(2 p\, m=0, V, 2 s), 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0)
$
explictly
$
  braket(2 s, V, 2 p\, m=0) = braket(2p\, m=0, V, 2 s) = 3 e a_0 abs(bold(E))
$
the upper left block looks like $sigma_x$ so we can immediately find
$
  Delta_(plus.minus)^((1)) = plus.minus 3 e a_0 abs(bold(E)) underbrace(prop abs(bold(E)), "linear" #linebreak() "stark effect")
$
and
$
  ket(plus.minus) = 1/sqrt(2) (ket(2 s\,m=0) plus.minus ket(2p\,m=0))
$
which is quite nice.

In a real hydrogen atom $2 s$ and $2 p$ are not degenerate meaning we can just use degenerate perturbation theory---given $V$ is much smaller than splitting. In this case
$
  Delta_n = - e abs(bold(E)) z_(n n) + e^2 abs(bold(E))^2 sum_(k eq.not n) abs(z_(n k))^2/(E_n^((0))-E_k^((0))) + dots
$
in non-degenerate case $z_(n n) = 0$ since $ket(n^((0)))$ is a parity eigenstate. So
$
  Delta_n = e^2 abs(bold(E))^2 sum_(k eq.not n) abs(z_(n k))^2/(E_n^((0))-E_k^((0))) underbrace(prop abs(bold(E))^2, "quadratic" #linebreak() "stark effect")
$


#pagebreak()
== Hydrogen-like atoms
We consider a nucleus with $Z$ protons orbited by a single electron $e^-$.

=== Relativistic correction
The non-relativistic Hamiltonian is
$
  H_0 = p^2/(2 m_e) - (Z e^2)/r
$
the kinetic energy changes due to relativity
$
  K = p^2/(2m_e) -> K_"rel" &= sqrt(p^2 c^2 + m_e^2 c^4) - underbrace(m_e c^2, "rest energy") \
  &tilde.eq^"Taylor" p^2/(2m_e) underbrace(- (p^2)^2/(8 m_e^3 c^2), "perturbation") + dots
$
so
$
  H = H_0 + V "with" V=- (p^2)^2/(8 m_e^3 c^2)
$
nicely there is no time-dependence, but we still need degenerate perturbation theory since $ket(n l m)$ is degenerate---so we need to find a basis where $V$ is diagonal. Consider
$
  [L_l, p^2] = [L_l,p_i]p_i + p_i [L_l,p_i]
$
and
$
  [L_i,p_l] = [epsilon_(i j k) x_j p_k, p_l] = epsilon_(i j k) delta_(j l) p_k
$
so $p_i [L_l,p_i] = 0$. This implies $[bold(L),V] = 0$ so it is already diagonal in the $ket(n l m)$ basis. Then the first order shift is
$
  Delta_(n l)^((1)) = braket(n l m, V, n l m)
$
we rewrite
$
  (p^2)^2/(8 m_e^3 c^2) = 1/(2 m_e c^2) (p^2/(2 m_e))^2 = 1/(2 m_e c^2) (H_0 + (Z e^2)/r)^2
$
giving
$
  Delta_(n l)^((1)) = - 1/(2 m_e c^2) [(E_n^((0)))^2 + 2 E_n^((0)) braket(n l m, (Z e^2)/r, n l m)+braket(n l m, (Z e^2)^2/r^2, n l m)]
$
by Weyl
$
  expval(K) = -1/2 expval(U) => underbrace(expval(E_"tot"), expval(H_0) = E_n^((0))) = 1/2 expval(U)
$
for a stationary state. For $H_0$
$
  E_n^((0)) = -1/2 expval((Z e^2)/r) => expval((Z e^2)/r) = -2 E_n^((0))
$
for a general radial Hamiltonian
$
  H = p^2/(2 m) + V(r)
$
we had
$
  [- hbar^2/(2 m^2) dv(, r) (r^2 dv(, r)) + underbrace((l(l+1) hbar^2)/(2 m r^2) + V(r), V_"eff" (r))] R_(E l) (r) = E R_(E l) (r)
$
consider
$
  (l(l+1)hbar^2)/(2 m r^2) ->^"relativistic correction" (l(l+1)hbar^2)/(2 m r^2) + underbrace((Z e^2)^2/r^2, "since" prop r^(-2)) = (tilde(l) (tilde(l)+1) hbar^2)/(2 m r^2)
$
by shifting $l(l+1) + Z e^4 (2 m) hbar^(-2)$ giving
$
  tilde(l) = l - 1/2 Z^2 e^4 (2 m)/(hbar^2 (l+1/2))
$
then
$
  E_n^((0)) = 1/2 m c^2 (Z^2 alpha^2)/n^2
$
with $n equiv N+l+1 -> tilde(n) equiv N + tilde(l)+1$ giving
$
  E_tilde(n)^((0)) &= - 1/2 m c^2 (Z^2 alpha^2)/(n-1/2 Z^2 e^4 (2m)/(hbar^2(l+1/2)))^2 \
  &tilde.eq - 1/2 m c^2 (Z^2 alpha^2)/n^2 (1 + underbrace(Z^2 e^4 (2 m)/(n hbar^2 (l+1/2)), "correction")) \
  &= E_n^((0)) - 4 E_n^((0)) n/(l+1/2)
$
so
$
  expval((Z e^2)^2/r^2) = - 4 E_n^((0)) n/(l+1/2)
$
collecting terms we obtain
$
  Delta_(n l)^((1)) &= -1/(2 m_e c^2) [(E_n^((0)))^2-4 (E_n^((0)))^2- 4 E_n^((0)) n/(l+1/2)] \
  &= E_n^((0)) [(Z^2 alpha^2)/n^2 (-3/4 + n/(l+1/2))] \
  &= -1/2 m_e c^2 Z^4 alpha^4 [- 3/(4 n^4) + 1/(n^3 (l+1/2))]
$

=== Spin-orbit coupling
We consider a nucleus of $Z$ protons with a single electron $e^-$ in the outermost shell.

This electron lives in a potential which is no longer a pure Coulomb potential:
$
  V_c (r) = e underbrace(phi(r), "cloud of inner" e^- #linebreak() "affect this guy") => bold(E) = - 1/e nabla V_c (r)
$
and because it is moving it feels a magnetic field
$
  bold(B) = - bold(v)/c times bold(E)
$
the electron has magnetic moment
$
  bold(mu) = (e bold(S))/(m_e c)
$
and this gives rise to the spin-orbit coupling
$
  V_"LS" = - 1/2 bold(mu) dot bold(B) = 1/(2 m_e^2 c^2) 1/r dv(V, r) (bold(L) dot bold(S))
$
we treat this as a pertubation
$
  H = H_0 + V_"LS" "with" H_0 = p^2/(2m) + V_c (r)
$
now we have two choices with respect to choice of basis: eigenkets of ${L^2, L_z, S^2, S_z}$ with $ket(l s m m')$ or eigenkets of ${L^2, S^2, J^2, J_z}$ with $ket(j=l+s\, m)$. We want the one where $V_"LS"$ is diagonal so the one where $bold(L) dot bold(S)$ commutes with everything---this is the second set with basis $ket(j=l+s\, m)$. Then is it simple to find the correction, we use
$
  psi_(n l m) = R_(n l) underbrace(cal(Y)_l^(j = l plus.minus 1/2\, m), "spin-angular function")
$
to obtain
$
  Delta_(n l j)^((1)) &= 1/(2 m_e^2 c^2) expval(1/r dv(V, r))_(n l) hbar^2/2 cases(l "for" j = l +1/2 "(spin-up)", -(l+1) "for" j = l -1/2 "(spin-down)") \
  expval(1/r dv(V_c, r))_(n l) &equiv integral_0^oo R_(n l) 1/r dv(V_c, r) R_(n l) r^2 dd(r)
$
this is called Lande's interval rule. To find the above we used
$
  integral cal(Y)^dagger bold(S) dot bold(L) cal(Y) dd(Omega) = 1/2 [j(j+1)-l(l+1)-3/4]hbar^2
$
in the case of a pure Coulomb potential
$
  expval(1/r dv(V_c, r))_(n l) & = expval((Z e^2)/r^3)_(n l)
$
to solve this consider
$
  expval([H_0,A]) = 0
$
for any operator $A$ in the case of $A = p_r$ we find
$
  0 & = expval([(l(l+1)hbar^2)/(2 m_e r^2) - (Z e^2)/r, p_r]) \
    & = expval(-(l(l+1)hbar^2)/(m_e r^3) + (Z e^2)/r^2)
$
then using the previous we can write
$
  expval((Z e^2)/r^3)_(n l) &= m_e/(l(l+1)hbar^2) expval((Z e^2)^2/r^2)_(n l) \
  &= - (2 m_e^2 c^2 Z^2 alpha^2)/(n l(l+1)(l+1/2)hbar^2) E_n^((0))
$


=== Summary of corrections
Due to relativistic effects we had
$
  Delta_(n l)^((1)) = E_n [(Z^2 alpha^2)/n^2 (-3/4 + n/(l+1/2))] prop cal(O)(alpha^4)
$
and due to the electron not being a point-like particle we need to smear the potential over $lambda_c = hbar\/m_e c$
$
  V(bold(r)) -> 1/V integral V(bold(r)+"smear")
$
using $V(bold(r) + "smear") = V(bold(r)) + "smear" times nabla V(bold(r)) + dots$---this gives the Darwin term
$
  Delta_(n 0)^((1)) = - E_n^((0)) (Z alpha)^2/n prop cal(O)(alpha^4)
$
this is equivalent to the above for $l = 0$. We also had spin-orbit coupling (fine-structure) as above.

Everything is of order $cal(O)(alpha^4)$ so terms naturally cancel and the total correction is
$
  E_n^((0)) -> E_n^((1)) = E_n^((0)) [1 + (Z alpha^2)/n (-3/4 + 1/(j+1/2))]
$
so $l$ drops out!

There are other higher order corrections, e.g. Lamb shift ($cal(O)(alpha^5)$) due to virtual photons and hyperfine splitting due to the nucleus not being spherically symmetric.


#pagebreak()
== Variational methods
Essentially just the principle
$
  expval(H) >= E_"gs"
$
for any trial state $ket(psi)$.

We make a guess as to what the ground state looks like $ket(tilde(0)) eq.not ket(0)$ and define
$
  tilde(H) = braket(tilde(0), H, tilde(0))/braket(tilde(0))
$
then clearly $tilde(H) >= E_"gs"$. We can see this by writing
$
  ket(tilde(0)) = sum_k ket(k) braket(k, tilde(0))
$
with $H ket(k) = E_k ket(k)$. Then
$
  tilde(H) &= (sum_k abs(braket(k, tilde(0)))^2 E_k)/(sum_k abs(braket(k, tilde(0)))^2) \
  &= underbrace((sum_k abs(braket(k, tilde(0)))^2 (E_k-E_0))/(sum_k abs(braket(k, tilde(0)))^2), >0 "since" E_k >= E_0) + E_0 \
  & >= E_0
$
so $tilde(H) >= E_0$.

If $ket(tilde(0)) = ket(0)$ then $braket(k, tilde(0)) = 0$ for all $k eq.not 0$. Assume $braket(k, tilde(0)) tilde cal(O)(epsilon)$ with $epsilon$ small, i.e. our guess is good. Then
$
  tilde(H) - E_0 tilde underbrace(cal(O)(epsilon^2), "from" Sigma)
$
which is relatively precise. We write our guess as
$
  ket(tilde(0)) = ket(0) + delta ket(tilde(0))
$
then $tilde(H) = 0 + cal(O) (delta^2)$ and then we can minimize with respect to $delta$.

=== Infinite square well
We consider a square well of width $2a$ centered at $x=0$. Then we know the solution
$
  braket(x, 0) = 1/sqrt(a) cos((pi x)/(2 a)) " and " E_n = hbar/(2 m) (pi^2/(4 a^2))
$
We try:
$
  braket(x, tilde(0)) = abs(a)^lambda - abs(x)^lambda
$
then we can find $tilde(H)$ and minimize to find $lambda$. We compute
$
  tilde(H) &= braket(tilde(0), H, tilde(0))/braket(tilde(0)) \
  &= (integral_(-a)^a dd(x) braket(tilde(0), x) braket(x, H, x) braket(x, tilde(0)))/(integral_(-a)^a dd(x) abs(braket(tilde(0), x))^2) \
  &= (- hbar^2/(2 m)) (integral_(-a)^a dd(x) (abs(a)^lambda-abs(x)^lambda) pdv(, x, 2) (abs(a)^lambda-abs(x)^lambda))/(integral_(-a)^a dd(x) (abs(a)^lambda - abs(x)^lambda)^2) \
  &= [((lambda+2)(2 lambda+1))/(2 lambda -1)] hbar^2/(4 m a^2)
$
minimizing we obtain
$
  pdv(tilde(H), lambda) = 0 => lambda = (1+sqrt(6))/2 tilde.eq 1.72
$
giving the energy
$
  tilde(H)_"min" = (5 + 2 sqrt(6))/pi^2 E_0 tilde.eq 1.063 E_0
$
which is pretty close. One can also show that
$
  abs(braket(0, tilde(0)))^2 >= 0.99963
$
meaning
$
  braket(0, tilde(0)) = cos theta => theta <= 1.1 degree
$
from
$
  H_"min" = sum_k abs(braket(k, tilde(0)))^2 E_k^2 >= abs(braket(0, tilde(0)))^2 E_0 + dots
$

#pagebreak()
== Time-dependent potentials
We will only consider time-dependence in the form
$
  H = H_0 + V(t)
$
where we take the solution to $H_0$ as being known $H_0 ket(n) = E_n ket(n)$. Then we can write any state at $t = 0$ as
$
  ket(alpha) = sum_n c_n (t=0) ket(n)
$
since we take $V(t=0) = 0$. We want to know what happens for $t > 0$. We write
$
  ket(alpha\, t) = sum_n overbrace(c_n (t), "from perturbation") underbrace(exp[(-i E_n t)/h], "time-evolution" #linebreak() "of unperturbed") ket(n)
$
so all time-dependence from $V(t)$ is carried by $c_n (t)$ with $c_n (0)$ corresponding to $V(t=0)=0$. Then the problem we want to solve is finding these expansion coefficients.

=== The interaction picture
From the Schrödinger picture we have
$
  ket(alpha\, t)_S = exp[(-i H t)/hbar] ket(alpha\, t=0)_S
$
meaning the states evolve in time. Then
$
  braket(alpha\, t, A_S, alpha\, t)_S &= braket(alpha\, t=0, exp[(i H t)/hbar] A_S exp[(- i H t)/hbar], alpha\, t=0) \
  &= braket(alpha, A_H (t), alpha)_H
$
meaning
$
  A_H = exp[(i H t)/hbar] A_S exp[(-i H t)/hbar]
$
and
$
  ket(alpha)_H = exp[(i H t)/hbar] ket(alpha\, t)_S
$
In the interaction picture we define
$
  ket(alpha\, t)_I = exp[(i H_0 t)/hbar] ket(alpha\, t)_S
$
if we have no interactions and $H = H_0$ then this state is the same as the Heisenberg state, but the interaction part lives in the Schrödinger picture. Then
$
  braket(alpha\, t, A_S, alpha\, t)_S &= braket(alpha\, t, exp[(- i H_0 t)/hbar] exp[(i H_0 t)/hbar] A_S exp[(-i H_0 t)/hbar] exp[(i H_0 t)/hbar], alpha\, t)_S \
  &= braket(alpha\, t, A_I (t), alpha\,t)_I
$
where
$
  A_I (t) equiv exp[(i H_0 t)/hbar] A_S exp[(-i H_0 t)/hbar]
$
so the operators evolve under the free theory, while the states evolve under the interaction. This means we get two equations. Consider
$
  i hbar pdv(, t) ket(alpha\, t)_I &= i hbar pdv(, t) (exp[(i H_0 t)/hbar] ket(alpha\, t)_S) \
  &= - H_0 exp[(i H_0 t)/hbar] ket(alpha\, t)_I + exp[(i H_0 t)/hbar] underbrace((i hbar pdv(, t) ket(alpha\,t)_S), "Schrödinger equation") \
  &= - H_0 exp[(i H_0 t)/hbar] ket(alpha\, t)_S + exp[(i H_0 t)/hbar] (H_0 + V) ket(alpha\, t)_S \
  &= exp[(i H_0 t)/hbar] V exp[(-i H_0 t)/hbar] exp[(i H_0 t)/hbar] ket(alpha\, t)_S \
  &= V_I ket(alpha\,t)_I
$
so we obtain
$
  i hbar pdv(, t) ket(alpha\, t)_I = V_I ket(alpha\, t)_I
$
and by definition
$
  dv(A_I, t) = 1/(i hbar) [A_I, H_0]
$

Before we had
$
  ket(alpha\, t)_S = sum_n c_n (t) exp[(-i E_n t)/hbar] ket(n)
$
giving
$
  ket(alpha\, t)_I & = exp[(i H_0 t)/hbar] ket(alpha\, t)_S \
                   & = sum_n c_n (t) ket(n)
$
which is nice. We act with $bra(n)$ on this guy and take the derivative
$
  i hbar pdv(, t) underbrace(braket(n, alpha\, t)_I, c_n (t)) &= sum_m braket(n, V_I, m) braket(m, alpha\, t)_I
$
consider
$
  braket(n, V_I, m) &= braket(n, exp[(i H_0 t)/hbar] V (t) exp[(-i H_0 t)/hbar], m) \
  &= underbrace(braket(n, V, m), V_(n m)) exp[(i (E_n - E_m) t)/hbar] \
  &= V_(n m) e^(i omega_(n m) t)
$
giving
$
  i hbar pdv(, t) c_n (t) & = sum_m V_(n m) e^(i omega_(n m) t) c_m (t)
$
=== Time-dependent two state problem
We have the Hamiltonian
$
  H_0 = E_1 ketbra(1) + E_2 ketbra(2) " with " E_2 > E_1
$
we add the perturbation $ V(t) = gamma e^(i omega t) ketbra(1, 2)+gamma e^(-i omega t) ketbra(2, 1) $
which acts to mix the states. We compute the four matrix elements:
$
  V_12 & = gamma e^(i omega t) => V_21 = gamma e^(- i omega t) \
  V_11 & = V_22 = 0
$
we assume $c_1 (0) = 1$ and $c_2 (0) = 0$. This gives
$
  dot(c)_1 & = 1/(i hbar) V_12 e^(i omega_12 t) c_2 \
  dot(c)_2 & = 1/(i hbar) V_21 e^(i omega_21 t) c_1
$
this can be solved to give
$
  abs(c_2 (t))^2 & = (gamma^2 \/hbar^2)/Omega^2 sin^2 (Omega t) \
  abs(c_1 (t))^2 & = 1 - abs(c_2 (t))^2
$
with
$
  Omega^2 = gamma^2/hbar^2 + (omega - omega_21)/4
$
as  $omega -> omega_21$ we find $abs(c_2 (t))^2 -> 1$ which is resonance.

=== Sudden approximation
We turn on a perturbation so fast that the state does not have time to adjust, so we are still in the unperturbed state (e.g. a step-function).

We define some characteristic time $T = Omega^(-1)$. Consider
$
  i hbar pdv(, t) U(t,t_0) = H U (t,t_0)
$
if we let $t = S Omega^(-1)$ then
$
  i pdv(, S) U(t,t_0) = H/(hbar Omega) U(t,t_0)
$
we take $T -> 0$ then $U(t,t_0) -> bb(1)$, so we are free to assume that the time-evolution becomes trivial. This is a fair approximation when $hbar Omega >> E_(n m) = hbar omega_(n m)$ or $T << 2 pi omega_(n m)^(-1)$.

=== Adiabatic approximation
This is the opposite of the sudden approximation, so a very slow time-variation meaning we are always in the perturbed state. The idea is to just use time-independent perturbation theory but make the coefficients time-dependent.

Consider the Schrödinger equation
$
  i hbar pdv(, t) ket(alpha\, t) = H (t) ket(alpha\, t)
$
we make the ansatz
$
  ket(alpha\, t) = sum_n c_n (t) e^(i theta_n) ket(n\, t)
$
with
$
  theta_n (t) equiv - 1/hbar integral_0^t E_n (t') dd(t')
$
which reduces to the time-independent case. Substituting this guy gives
$
  sum_n e^(i theta_n) [dot(c)_n (t) + c_n (t) pdv(, t)] ket(n\, t) = 0
$
acting with $bra(m\, t)$ we find
$
  dot(c)_m (t) &= - sum_n c_n (t) e^(i(theta_m - theta_n)) braket(m\,t, pdv(, t), n\, t)
$
Consider
$
  H(t) ket(n\, t) = E_n (t) ket(n\, t)
$
taking the derivative and acting with $bra(m\, t)$ we find
$
  braket(m\, t, dot(H), n\, t) = (E_n (t)-E_m (t)) braket(m\, t, pdv(, t), n\, t)
$
for $n eq.not m$ where $braket(m, dot(E), n)=0$. Combining the above we find
$
  dot(c)_m (t) = underbrace(c_m (t) braket(m\, t, pdv(, t), m\, t), "diagonal term") + sum_(m eq.not n) c_n (t) e^(i (theta_m-theta_n)) braket(m\, t, dot(H), n\, t)/(E_n - E_m)
$
we take $dot(H)$ to be small, or the first term should be larger then the sum:
$
  E_m >> expval(dot(H))/E_(n m) equiv 1/tau
$
then we can write
$
  c_n (t) = e^(i gamma_n (t)) c_n (0)
$
where
$
  gamma_n (t) equiv i integral_0^t braket(n\, t', pdv(, t'), n\, t') dd(t')
$
is the Berry phase.
