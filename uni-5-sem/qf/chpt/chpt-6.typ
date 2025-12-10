//**** init-ting
#import "@preview/physica:0.9.7": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Scattering theory
== The $S$-matrix
Let
$
  H = H_0 + V(bold(r)) " where " H_0 = p^2/(2 m)
$
then
$
  H_0 ket(bold(k)) = E_bold(k) ket(bold(k)) " with " E_bold(k) = (hbar^2 bold(k)^2)/(2 m)
$
so $ket(bold(k))$ are free plane-waves, we imagine these come in form negative infinity and scatter. Recall
$
  ket(alpha\, t)_I = U_I (t,t_0) ket(alpha\,t_0)_I
$
which satisfies
$
  i hbar pdv(, t) U_I (t,t_0) = V_I (t) U_I (t,t_0)
$
with
$
      V_I (t) & = exp((i H_0 t)/hbar) V exp(-(i H_0 t)/hbar) \
  U_I (t,t_0) & = 1 - i/hbar integral_(t_0)^t V_I (t') U_I (t',t_0) dd(t')
$
then
$
  braket(n, U_I (t,t_0), i) &= delta_(n i) - i/hbar sum_m braket(n, V, m) integral_(t_0)^t e^(i omega_(m n) t') braket(m, U_I (t',t_0), i) dd(t')
$
We regularize our problem by enclosing it in a box of length $L$. Then the free plane-wave solution becomes
$
  psi_k (x) = 1/(2 pi)^(3\/2) e^(i bold(k) dot bold(x)) -> psi_(k_n) (x) = 1/L^(3\/2) e^(i bold(k)_n dot bold(x))
$
with
$
  integral dd(x, 3) psi_(k_n)^* (x) psi_(k_m) (x) = delta_(n m) "etc."
$
we use these later.

For now to first order the Dyson series becomes ($U -> delta_(m i)$)
$
  braket(n, U_I (t,t_0), i) &= delta_(n i) - i/hbar braket(n, V, i) integral_(t_0)^t e^(i omega_(n i) t') dd(t')
$
we consider $t_0 -> -oo$ and $t -> oo$ (and as before $V(t) = e^(eta t) V)$. Define
$
  braket(n, V, i) equiv T_(n i) underbrace(e^(epsilon.alt t), "inside" integral)
$
then
$
  braket(n, U_I (t,t_0), i) &= delta_(n i) - i/hbar T_(n i) integral_(t_0)^t e^(i omega_(n i) t' + epsilon.alt t') dd(t')
$
By assumption $epsilon.alt > 0$, and for $epsilon.alt -> 0$ we have $t << epsilon.alt^(-1)$. Then $lim_(t-> oo) lim_(epsilon.alt -> 0) e^(epsilon.alt t) = 1$ and $t_0 = -oo$ gives $e^(epsilon.alt t_0) -> 0$. Taking the limits in the correct order we find
$
  S_(n i) &equiv lim_(t-> oo) [lim_(epsilon.alt -> 0) braket(n, U_I (t,-oo), i)] \
  &= delta_(n i) - i/hbar T_(n i) integral_(-oo)^oo e^(i omega_(n i) t') dd(t') \
  &= delta_(n i) - 2 pi i delta(E_n - E_i) T_(n i)
$
which is the $S$-matrix. These are related to the probability of scattering.

== The cross-section
The transition rate is
$
  omega_(i -> n) & = dv(, t) abs(c_n (t))^2 \
                 & = dv(, t) abs(braket(n, U_I (t,-oo), i))^2
$
for $i eq.not n$ we compute
$
  braket(n, U_I (t,-oo), i) &= -i/hbar T_(n i) integral_(-oo)^t e^(i omega_(n i) t' + epsilon t') dd(t') \
  &= -i /hbar T_(n i) e^(i omega_(n i) t+ epsilon t)/(i omega_(n i)+t)
$
so
$
  omega_(i -> n) &= dv(, t) [1/hbar^2 abs(T_(n i))^2 e^(2 epsilon.alt t)/(omega_(n i)^2 + epsilon.alt^2)] \
  &= 1/hbar^2 abs(T_(n i))^2 (2 epsilon.alt e^(2 epsilon.alt t))/(omega_(n i)^2 + epsilon.alt^2)
$
we take the limit
$
  lim_(t -> 0) (epsilon.alt #h(0.5em) overbrace(e^(2 epsilon.alt t), -> 1))/(omega_(n i)^2 + epsilon.alt^2) &= pi delta(omega_(n i)) \
  &= pi hbar delta (E_n-E_i)
$
so we recover Fermi's golden rule
$
  omega_(i -> n) = (2 pi )/hbar abs(T_(n i))^2 delta(E_n-E_i)
$

We want to relate this to the cross-section. We consider the probability of any scattering occuring. Recall the density of states
$
  rho (E_n) = dv(n, E_n, d: Delta) = (m k)/(hbar^2) (L/(2 pi))^3 dd(Omega)
$
with
$
  E_bold(n) = hbar^2/(2 m) ((2 pi)/L)^3 abs(bold(n))^2 " and " dd(n, d: Delta) = 4 pi abs(n)^2 Delta abs(n) dd(Omega)/(4 pi)
$
we then compute
$
  integral dd(n) omega_(i -> n) &= integral dd(E) rho(E) omega_(i -> n) \
  &= integral (m k)/hbar^2 (L/(2 pi))^3 omega_(i -> n) dd((Delta E)) dd(Omega) \
  &= (m k L^3)/((2 pi)^2 hbar^3) abs(T_( n i))^2 dd(Omega)
$
this is the probability of going from $i -> "any state in" dd(Omega)$. We now consider shooting some target having some cross-section with a flux of particles. The size of the cross-section is obviously related to how probable it is to hit the target---it can be treated as an effective interaction area. Then
$
  "flux" times "cross-section" = "transition rate"
$
In the plane-wave approximation
$
  v = (hbar k)/m";  " t = L/v";  " A = L^2
$
so
$
  "flux" = 1/(A t) = v/L^3 = (hbar k)/(m L^3)
$
then the differential cross-section is
$
  dv(sigma, Omega) & = ( (m L^3)/(2 pi hbar^2) )^2 abs(T_(n i))^2
$
integrating this guy will then give the cross-section $sigma$.

== Lippmann-Schwinger
By the above we have
$
  braket(n, U_I (t,-oo), i) &=^(cal(O)(V)) delta_(n i) + 1/hbar T_(n i) e^(i omega_(n i) t + epsilon.alt t)/(- omega_(n i) + i epsilon.alt) \
  braket(n, U_I (t,-oo), i) &= delta_(n i) - i /hbar sum_m V_(n m) integral_(-oo)^t e^(i omega_(n m) t') braket(m, U_I (t',-oo), i) dd(t')
$
using the first in the second gives three terms. The first is just $delta_(n i)$ the second is
$
  1/hbar V_(n i) e^(i omega_(n i) t + epsilon.alt t)/(-omega_(n i) + i epsilon.alt)
$
and the third is
$
  T_3 &= - i/hbar^2 sum_m V_(n m) T_(m i)/(-omega_(m i) + i epsilon.alt) integral_(-oo)^t e^(i omega_(n m) t' + i omega_(m i) t' + epsilon.alt t') dd(t') \
  &= - i/hbar^2 sum_m V_(n m) T_(m i)/(- omega_(m i) + i epsilon.alt) integral_(-oo)^t e^(i omega_(n i) t' + epsilon.alt t') dd(t') \
  &= - i/hbar^2 sum_m V_(n m) T_(m i)/(-omega_(m i) + i epsilon.alt) e^(i omega_(n i) t + epsilon.alt t)/(i omega_(n i) + epsilon.alt) \
  &= 1/hbar^2 e^(i omega_(n i) t + epsilon.alt t)/(-omega_(n i) + i epsilon.alt) sum_m V_(n m) T_(m i)/(-omega_(m i) + i epsilon.alt) \
  &= 1/hbar e^(i omega_(n i) t + epsilon.alt t)/(-omega_(n i) + i epsilon.alt) sum_m V_(n m) T_(m i)/(E_i- E_m + i hbar epsilon.alt)
$
we compare this with the $cal(O)(V)$ expression
$
  braket(n, U_I (t,-oo), i) &= delta_(n i) + 1/hbar V_(n i) e^(i omega_(n i) t + epsilon.alt t)/(-omega_(n i) + i epsilon.alt) + 1/hbar e^(i omega_(n i) t + epsilon.alt t)/(-omega_(n i) + i epsilon.alt) sum_m V_(n m) T_(m i)/(E_i- E_m + i hbar epsilon.alt) \
  &= delta_(n i) + 1/hbar e^(i omega_(n i) t + epsilon.alt t)/(-omega_(n i) + i epsilon.alt) [V_(n i) + sum_m V_(n m) T_(m i)/(E_i-E_m + i hbar epsilon.alt)] \
  &= delta_(n i) + 1/hbar e^(i omega_(n i) t +epsilon.alt t)/(- omega_(n i) + i epsilon.alt) T_(n i)
$
this gives a relation between $T_(n i)$ and $V_(n i)$
$
  T_(n i) = V_(n i) + sum_m V_(n m) T_(m i)/(E_i-E_m + i hbar epsilon.alt)
$
It is now convenient to define the Lippmann-Schwinger state $ket(psi^((+)))$ by
$
  T_(n i) = sum_j braket(n, V, j) braket(j, psi^((+))) = braket(n, V, psi^((+)))
$
by the above
$
  braket(n, V, psi^((+))) = braket(n, V, i) + sum_m braket(n, V, m) braket(m, V, psi^((+)))/(E_i-E_m + i hbar epsilon.alt)
$
this must be true for all $ket(n)$ giving
$
  ket(psi^((+))) &= ket(i) + sum_m ket(m) braket(m, V, psi^((+)))/(E_i-E_m+i hbar epsilon.alt) \
  &= ket(i) + sum_m 1/(E_i-H_0 + i hbar epsilon.alt) ket(m) braket(m, V, psi^((+))) \
  &= ket(i) + 1/(E_i-H_0 + i hbar epsilon.alt) V ket(psi^((+)))
$
this is the Lippmann-Schwinger equation. Using $ket(psi^((+)))$ we can write
$
  dv(sigma, Omega) = ((m L^3)/(2 pi hbar^2))^2 abs(braket(n, V, psi^((+))))^2
$
We can also define an operator $T$ with matrix elements $T_(n i) = braket(n, T, i)$ by $T ket(i) = V ket(psi^((+)))$. Acting with $V$ on $ket(psi^((+)))$ gives the equation
$
  T ket(i) = V ket(i) + V 1/(E_i-H_0 + i hbar epsilon.alt) T ket(i)
$
which can be written as an operator equation by
$
  T = V + V 1/(E_i-H_0+i hbar epsilon.alt) T
$
given $V$ is _weak_ an expansion obviously follows by
$
  T = V + V 1/(E_i-H_0+i hbar epsilon.alt) V + V 1/(E_i-H_0+i hbar epsilon.alt) V 1/(E_i-H_0+i hbar epsilon.alt) V + dots
$

== The scattering amplitude
Above we had
$
  dv(sigma, Omega) = ((m L^3)/(2 pi hbar^2))^2 abs(braket(n, V, psi^((+))))^2
$
we will use this later. We define $psi^((-))$ to be the inverse of $psi^((+))$ then
$
  braket(bold(x), psi^((plus.minus))) = braket(bold(x), i)+ integral dd(x, 3) braket(bold(x), 1/(E_i - H_0 plus.minus i epsilon.alt), bold(x)') braket(bold(x)', V, psi^((plus.minus)))
$
Consider the propagator
$
  G_plus.minus (bold(x),bold(x)') &= hbar^2/(2 m) braket(bold(x), 1/(E-H_0 plus.minus i epsilon.alt), bold(x)') \
  &= hbar^2/(2 m) sum_k' sum_(k'') underbracket(braket(bold(x), bold(k)'), prop e^(i bold(k)' dot bold(x))) underbracket(braket(bold(k)', 1/(E-H_0 plus.minus i epsilon.alt), bold(k)''), delta_(k' k'') [E-hbar^2 k''^2\/2 m plus.minus i epsilon.alt]^(-1)) braket(bold(k)'', bold(x)') \
  &= 1/L^3 sum_k' e^(i bold(k)' dot (bold(x)-bold(x)'))/(k^2 - k'^2 plus.minus i epsilon.alt) \
  &=^(L^(-3) sum ->^(L -> oo) integral dd(k, 3)\/(2 pi)^3) integral dd(k, 3)/(2 pi)^3 e^(i bold(k) dot (bold(x)-bold(x)'))/(k^2-k'^2 plus.minus i epsilon.alt) \
  &= 1/(2 pi)^3 integral_0^oo dd(k') k'^2 integral_1^(-1) dd(cos theta) integral_0^(2 pi) dd(phi) [e^(i k' abs(bold(x)-bold(x)') cos theta)/(k^2-k'^2 plus.minus i epsilon.alt)] \
  &= 1/(8 pi^2) 1/(i abs(bold(x)-bold(x)')) integral dd(k')k' [(e^(-i k' abs(bold(x)-bold(x)'))-e^(i k' abs(bold(x)-bold(x)')))/(k^2-k'^2 plus.minus i epsilon.alt)]
$
so we need to compute
$
  I_- &= integral_(-oo)^oo dd(k') k' e^(-i k' abs(bold(x)-bold(x)'))/(k^2-k'^2 plus.minus i epsilon.alt) \
$
we use a contour along $RR$ closing in lower half-plane. The closing part vanishes since $k' = Re (k') + i Im (k')$, meaning for $Im (k') -> oo$ the exponential dies. Then by the residue theorem
$
  I_- &= integral.cont_C_L dd(k') k' e^(-i k' abs(bold(x)-bold(x)'))/(k^2-k'^2 plus.minus i epsilon.alt) \
  &= 2 pi i Res(k' e^(-i k' abs(bold(x)-bold(x)'))/(k^2-k'^2 plus.minus i epsilon.alt))
$
since we pick up the pole at $k' tilde.eq -k$. We use
$
  Res(f) & = lim_(z->z_0) (z-z_0) f(z) \
         & =^(f=g h^(-1)) lim_(z -> z_0) (z g(z) - z_0 g(z))/h(z) \
         & = lim_(z->z_0) (z g'(z) + g(z) - z_0 g'(z))/(h'(z)) \
         & = g(z_0)/(h'(z_0))
$
to obtain
$
  I_- & = 2 pi i ((-k e^(i k abs(bold(x)-bold(x)')))/(2 k)) \
      & = - pi i k e^(i k abs(bold(x)-bold(x)'))
$
for the other integral we pick up the pole at $k' tilde.eq k$ giving
$
  I_+ & = - pi i k e^(-i k abs(bold(x)-bold(x)'))
$
instead taking $k -> -k$ gives $I_+ = I_-$ so
$
  G_plus.minus = -1/(4 pi) e^(plus.minus i k abs(bold(x)-bold(x)'))/(abs(bold(x)-bold(x)'))
$
This is the Green's function for
$
  (nabla^2 + k^2) G_plus.minus (bold(x),bold(x)') = delta^((3)) (bold(x)-bold(x)')
$
We can now compute
$
  braket(bold(x), psi^((plus.minus))) = braket(bold(x), i) - (2 m)/hbar^2 integral dd(x', 3) e^(plus.minus i k abs(bold(x)-bold(x)'))/(4 pi abs(bold(x)-bold(x)')) braket(bold(x)', V, psi^((plus.minus)))
$
we take the potential to be local
$
  braket(bold(x)', V, bold(x)'') = V (bold(x)) delta^((3)) (bold(x)'-bold(x)'')
$
giving
$
  braket(bold(x)', V, psi^((plus.minus))) &= integral dd(x'', 3) braket(bold(x)', V, bold(x)'') braket(bold(x)', psi^((plus.minus))) \
  &= V(bold(x)') braket(bold(x)', psi^((plus.minus)))
$
so
$
  braket(bold(x), psi^((plus.minus))) &= braket(bold(x), i) - (2 m)/hbar^2 integral dd(bold(x)', 3) e^(plus.minus i k abs(bold(x)-bold(x)'))/(4 pi abs(bold(x)-bold(x)')) V(bold(x)') braket(bold(x)', psi^((plus.minus)))
$
Now we use $abs(bold(x)) >> abs(bold(x)')$ and introduce $r = abs(bold(x))$, $r' = abs(bold(x)')$ then
$
  abs(bold(x)-bold(x)') tilde.eq r - hat(r) dot bold(r)'
$
with $bold(k)' equiv k hat(r)$ we can write
$
  e^(plus.minus i k abs(bold(x)-bold(x)')) &tilde.eq e^(plus.minus i k r) e^(minus.plus i bold(k)' dot bold(x)') \
  1/abs(bold(x)-bold(x)') &tilde.eq 1/r
$
Taking the incoming mode to be a plane-wave we find
$
  braket(bold(x), psi^((plus))) &=^"large r" braket(bold(x), bold(k)) - (2 m)/(4 pi hbar^2) e^(i k r)/r integral dd(x', 3) e^(-i bold(k)' dot bold(x)') V(bold(x)') braket(bold(x)', psi^((+))) \
  &= 1/L^(3\/2) [underbrace(e^(i bold(k) dot bold(x)), "plane-wave") + underbrace(e^(i k r)/r f(bold(k)',bold(k)), "spherical wave")]
$
with
$
  f(bold(k)',bold(k)) &= - 1/(4 pi) (2 m)/hbar^2 L^3 integral dd(x', 3) underbrace(e^(-i bold(k)' dot bold(x)')/(L^(3\/2)), "plane-wave") V(bold(x)') braket(bold(x)', psi^((+))) \
  &= - (2 m L^3)/(4 pi hbar^2) braket(bold(k)', V, psi^((+)))
$
for $psi^((-))$ we have
$
  braket(bold(x), psi^((-))) = 1/L^(3\/2) [e^(bold(k) dot bold(x)) + e^(-i k r)/r f(-bold(k)',bold(k))]
$
Then the cross-section can be written as
$
  dv(sigma, Omega) = abs(f(bold(k)',bold(k)))^2
$

== The optical theorem
For forward scattering $f(theta = 0) equiv f(bold(k),bold(k))$ then
$
  Im f(theta = 0) = (k sigma_"tot")/(4 pi)
$
with
$
  sigma_"tot" equiv integral dd(Omega) dv(sigma, Omega)
$

== The Born approximation
We defined
$
  T_(n i) = braket(n, V, psi^((+)))";  " braket(n, T, i) = T_(n i)
$
we can write this as
$
  braket(bold(k)', T, bold(k)) = braket(bold(k)', V, psi^((+)))
$
the first order Born approximation takes
$
  T tilde.eq V
$
i.e. only keeping the first term in the expansion for $T$.

We also derived
$
  dv(sigma, Omega) = abs(f(bold(k)',bold(k)))^2
$
with
$
  f(bold(k)',bold(k)) & = - (m L^3)/(2 pi hbar^2) braket(bold(k)', V, psi^((+))) \
  & =^"by above" - (m L^3)/(2 pi hbar^2) braket(bold(k)', T, bold(k))
$
to first order
$
  f^((1)) (bold(k)',bold(k)) &= - (m L^3)/(2 pi hbar^2) braket(bold(k)', V, bold(k)) \
  &= - (m L^3)/(2 pi hbar^2) integral dd(bold(x)', 3) integral dd(bold(x)'', 3) underbracket(braket(bold(k)', bold(x)'), "plane-wave") braket(bold(x)', V, bold(x)'') braket(bold(x)'', bold(k)) \
  & = -(m)/(2 pi hbar^2) integral dd(bold(x)', 3) integral dd(bold(x)'', 3) e^(-i bold(k)' dot bold(x)') V(bold(x)'') delta^((3)) (bold(x)'-bold(x)'') e^(i bold(k) dot bold(x)'') \
  &= - m/(2 pi hbar^2) integral dd(bold(x)', 3) e^(i (bold(k)-bold(k)') dot bold(x)') V(bold(x)')
$
so $f^((1)) (bold(k)',bold(k))$ is the Fourier transform of $V(bold(x))$ with respect to $bold(q) = bold(k)-bold(k)'$.

By momentum conservation $abs(bold(k)) = abs(bold(k)') equiv k$ meaning
$
  q = abs(bold(k)-bold(k)') = 2 k sin theta\/2
$
with $theta$ being the angle between $bold(k)$ and $bold(k)'$. Then for a spherically symmetric potential $V(bold(x)) = V(r)$ we find
$
  f^((1)) (theta) &= - m/(2 pi hbar^2) integral dd(r) r^2 integral dd(cos tilde(theta)) integral dd(phi) e^(i q r cos tilde(theta)) V(r) \
  &= - m/(i q hbar^2) integral dd(r) r V(r) [e^(i q r)-e^(-i q r)] \
  &= - (2 m)/(q hbar^2) integral dd(r) r V(r) sin(q r)
$
which is quite simple.

As an example consider a finite well with
$
  V(r) = cases(V_0 &" for" r <= a, 0 & " otherwise")
$
we find
$
  f^((1)) (theta) & = - (2 m)/(q hbar^2) integral_0^oo dd(r) r V(r) sin (q r) \
  & = - (2 m V_0)/(q hbar^2) integral_0^a dd(r) r sin(q r) \
  &= - (2 m)/hbar^2 (V_0 a^3)/(q a)^2 [(sin q a)/(q a) - cos q a]
$
the zeroes of this guy can be used to determine $a$.

Another example is the Yukawa potential
$
  V(r) = (V_0 e^(- mu r))/(mu r)
$
which gives
$
  f^((1)) (theta) = - ((2 m V_0)/(mu hbar^2)) 1/(q^2 + mu^2)
$
by the usual trick $sin q r = Im (e^(i q r))$. By definition we have
$
  q^2 = 4 k^2 sin^2 theta/2 = 2 k^2 (1- cos theta)
$
so we can write
$
  dv(sigma, Omega) tilde.eq ((2 m V_0)/(mu hbar^2))^2 1/[2 k^2 (1- cos theta)+mu^2]^2
$
The Yukawa potential reduces to the Coulomb potential in the limit $mu -> 0$ provided the ratio $V_0 mu^(-1)$ is held fixed. Let $V_0 mu^(-1) = Z Z' e^2$ then
$
  dv(sigma, Omega) tilde.eq ((2 m)^2 (Z Z' e^2)^2)/hbar^4 1/(16 k^4 sin^4 theta\/2)
$
taking $p = hbar k$ with $E = p^2\/2 m$ we obtain the classical Rutherford scattering
$
  dv(sigma, Omega) tilde.eq 1/16 ((Z Z' e^2)/E)^2 1/(sin^4 theta\/2)
$

When is the Born approximation valid? We had
$
  ket(psi^((+))) tilde.eq ket(bold(k))
$
and
$
  braket(bold(x), psi^((+))) = braket(bold(x), bold(k)) - underbracket((2 m)/hbar integral dd(bold(x)', 3) e^(i k abs(bold(x)-bold(x)'))/(4 pi abs(bold(x)-bold(x)')) V(bold(x)') braket(bold(x)', psi^((+))), "must be small!")
$
To get a proper criterion we take $V tilde.eq V_0$ within some range $a$ and let $r' = abs(bold(x)-bold(x)')$. Then
$
  abs((2 m)/hbar^2 ((4 pi)/3 a^3) e^(i k r')/(4 pi a) V_0 e^(i bold(k) dot bold(x)')/L^(3\/2)) << abs(e^(i bold(k) dot bold(x))/L^(3\/2))
$
for $k a << 1$ (low energy) we find
$
  (m abs(V_0) a^2)/hbar^2 << 1
$
as an example for the Yukawa potential we take $a = mu^(-1)$ giving
$
  (m abs(V_0))/(hbar^2 mu^2) << 1
$
For $k a >> 1$ (high energy) we find
$
  (2 m)/hbar^2 (abs(V_0) a)/k ln k a << 1
$
we see that the approximation gets better and better for higher and higher energies.

To second order we use
$
  T tilde.eq V + V 1/(E - H_0 + i epsilon.alt) V
$
giving
$
  f(bold(k),bold(k)') tilde.eq f^((1)) (bold(k),bold(k)') + f^((2)) (bold(k),bold(k)')
$
where
$
  f^((2)) (bold(k),bold(k)') &= - (m L^3)/(2 pi hbar^2) braket(bold(k)', V 1/(E-H_0 + i epsilon.alt) V, bold(k)) \
  &=^"continuous limit" - 1/(4 pi) (2 m)/hbar^2 (2 pi)^3 integral dd(bold(x)', 3) integral dd(bold(x)'', 3) braket(bold(k)', bold(x)') V (bold(x)') \
  & #h(6em) times bra(bold(x)') 1/(E-H_0 +i epsilon.alt) ket(bold(x)'') V(bold(x)'') braket(bold(x)'', bold(k)) \
  &= - 1/(4 pi) (2 m)/hbar^2 integral dd(bold(x)', 3) integral dd(bold(x)'', 3) e^(-i bold(k)' dot bold(x)') V(bold(x)') [(2 m)/hbar^2 G_+ (bold(x)',bold(x)'')] V(bold(x)'') e^(i bold(k) dot bold(x)'')
$
physically this can be viewed as a two-part process where we first scatter at $bold(x)''$ after which we propagate to and scatter at $bold(x)'$.

== Partial waves
By a Fourier transform we can expand a function in terms of plane-waves
$
  f(x) = integral dd(k, 3)/(2 pi)^3 f(bold(k)) e^(i bold(k) dot bold(x))
$
but using
$
  e^(i bold(k) dot bold(x)) = 4 pi sum_(l,m) i^l j_l (k x) Y_(l m) (hat(x)) Y_(l m)^* \(hat(k))
$
we can expand the function in terms of the spherical harmonics
$
  f(x) = integral_0^oo dd(k) sum_(l, m) f_(l m) (k) sqrt(2/pi) k j_l (k x) Y_(l m) (theta, phi)
$
with
$
  f_(l m) (k) equiv k i^l integral dd(Omega_k) f(k) Y_(l m)^* \(hat(k))
$
we treat this as the _spherical Fourier transform_. We can do something analogously in quantum mechanics.

=== Free states
Consider a free particle with
$
  H_0 = K
$
then $[H_0, L^2] = [H_0, L_z] = 0$. Then we can expand $ket(bold(k))$ in terms of $ket(E\,l\,m)$ which we call the spherical wave state. We then claim
$
  braket(bold(k), E\,l\,m) = g_(l E) Y_(l m) \(hat(k))
$
with $g_(l E)$ to be determined.

Let $bold(k) = k hat(z)$ then $L_z ket(k hat(z)) = 0$ trivially meaning $m = 0$ for this state. Then by Wigner-Eckart we have
$
  braket(E'\,l'\,m', k hat(z)) = 0
$
for $m' eq.not 0$. So we can write
$
  ket(k hat(z)) = sum_l' integral dd(E') ket(E'\,l'\, m' = 0) braket(E'\,l'\,m' = 0, k hat(z))
$
this is nice because now we can find any $ket(bold(k))$ by rotating this guy
$
  ket(bold(k)) = cal(D) (phi, theta) ket(k hat(z))
$
We find
$
  braket(E\,l\,m, bold(k)) &= sum_l' integral dd(E') braket(E\,l\,m, D(phi,theta), E'\,l'\,m'=0) braket(E'\,l'\,m'=0, k hat(z)) \
  &= sum_l' integral dd(E') cal(D)_(m 0)^((l')) (phi,theta) delta_(l l') delta(E-E') braket(E'\,l'\,m'=0, k hat(z)) \
  &= cal(D)_(m 0)^((l)) (phi,theta) braket(E\,l\, m = 0, k hat(z))
$
recall
$
  cal(D)_(m 0)^((l)) (theta, phi) = sqrt((4 pi)/(2 l+1)) Y_(l m)^* (theta, phi)
$
so
$
  braket(bold(k), E\,l\,m) = g_(l E) (k) Y_(l m) \(hat(k))
$
where we define
$
  g_(l E) equiv sqrt((4 pi)/(2 l+1)) braket(E\, l\, m = 0, k hat(z))
$
Consider
$
  (H_0 - E) ket(E\,l\,m) = 0
$
and
$
  bra(bold(k)) (H_0-E) = ((hbar^2 k^2)/(2 m) - E) bra(bold(k))
$
so
$
  braket(bold(k), H_0 -E, E\,l\,m) = ((hbar^2 k^2)/(2 m)- E) underbracket(braket(bold(k), E\,l\,m), "non-vanishing if" #linebreak() E = dots) = 0
$
meaning
$
  g_(l E) (k) = N delta((hbar^2 k^2)/(2 m)- E)
$
To determine $N$ we assume $braket(E'\,l'\,m', E\, l\, m) = delta(E-E') delta_(l l') delta_(m m')$ which eventually gives
$
  N = hbar/sqrt(m k)
$
and we finally obtain
$
  braket(bold(k), E\,l\,m) = hbar/sqrt(m k) delta((hbar^2 k^2)/(2 m)- E) Y_(l m) \(hat(k))
$
so
$
  ket(bold(k)) &= sum_(l, m) integral dd(E) ket(E\,l\, m) braket(E\, l\, m, bold(k)) \
  &= sum_(l=0)^oo sum_(m=-l)^l evaluated(ket(E\,l\,m))_(E = hbar^2 k^2\/2 m) [hbar/sqrt(m k) Y_(l m)^* \(hat(k)\) ]
$
this is nice!

We can also write this in position space. We write
$
  braket(bold(x), E\,l\,m) &= integral dd(bold(k), 3) braket(bold(x), bold(k)) braket(bold(k), E\,l\, m) \
  &= integral dd(k) k^2 4 sum_(l', m') i^l' j_l' (k x) Y_(l'm') (hat(x)') hbar/sqrt(m k) delta((hbar^2 k^2)/(2 m)- E) \
  &#h(1em) times integral dd(Omega_k) Y_(l' m')^* \(hat(k)) Y_(l m) \(hat(k)) \
  &= (i^l)/hbar sqrt((2 m k)/pi) j_l (k x) Y_(l m) \(hat(x))
$

=== The expansion
We assume the potential is spherically symmetric meaning $T$ commutes with $L^2$ and $bold(L)$. Then we can write
$
  braket(E'\,l'\,m', T, E\,l\,m) = T_l (E) delta_(l l') delta_(m m') delta(E-E')
$
and we can compute
$
  f(bold(k),bold(k)) &= - 1/(4 pi) (2 m L^3)/hbar^2 braket(bold(k)', T, bold(k)) \
  &= - 1/(4pi) (2 m)/hbar^2 (2pi)^3 sum_(l m) sum_(l' m') integral dd(E) integral dd(E') braket(bold(k)', E'\,l'\,m') \
  &#h(7em) times braket(E'\,l'\,m', T, E\,l\,m) braket(E\,l\,m, bold(k)) \
  &= - (4 pi^2)/k sum_(l m) T_l (E) evaluated(Y_(l m) \(hat(k)'\) Y_(l m)^* \(hat(k)\))_(E = hbar^2 k^2\/2m)
$
We pick $bold(k) = k hat(z)$ then
$
  Y_(l m) \(hat(k)\) = sqrt((2 l +1)/(4 pi)) delta_(m 0)
$
so only $m = 0$ contributes. Let $theta$ be the angle between $bold(k)$ and $bold(k)'$ then
$
  Y_(l m) \(hat(k)') = sqrt((2 l + 1)/(4 pi)) P_l (cos theta)
$
We define the _partial wave amplitude_
$
  f_l (k) equiv - (pi T_l (E))/k
$
then
$
  f(bold(k)',bold(k)) = sum_(l=0)^oo (2 l +1) f_l (k) P_l (cos theta)
$
this is the partial wave expansion---see Sakurai for further details. We can write
$
  braket(bold(x), psi^((+))) =^("large" r) 1/(2 pi)^(3\/2) sum_l (2 l + 1) P_l/(2 i k) [(1 + 2 i k f_l (k)) e^(i k r)/r - e^(-i (k r - l pi))/r]
$
we see that scattering only changes the coefficient of the ougoing wave by
$
  1 -> 1 + 2 i k f_l (k)
$
with the incoming wave unaffected.

== Phase shifts
Recall the Schrödinger equation
$
  i hbar pdv(, t) psi = - hbar^2/(2m) nabla^2 psi + V psi
$
define $rho = abs(psi)^2$ then
$
  i hbar pdv(, t) (psi^* psi) &= - hbar^2/(2 m) nabla (psi^* nabla psi - (nabla psi^*) psi)
$
by continuity
$
  pdv(rho, t) + nabla dot bold(j) = 0
$
so we define
$
  bold(j) & = (i hbar)/(2 m) (psi^* nabla psi - (nabla psi^*) psi) = hbar/m Im (psi^* nabla psi)
$
We make the ansatz
$
  psi = sqrt(rho) e^(i S\/hbar)
$
giving
$
  bold(j) = (rho nabla S)/m
$
in the classical limit one can show that $S$ is Hamilton's principle function. The above $bold(j)$ ensure probability is conserved. This is what we mean by unitarity.

Take $rho$ to be time-independent. Then
$
  div bold(j) = 0
$
or by Gauss' theorem
$
  integral_"spherical surface" bold(j) dot dd(S) = 0
$
so the flux is zero. By angular momentum conservation this must be true for each partial wave separately. This means the magnitude of the coefficiens for $e^(-i k r)\/r$ and $e^(i k r)\/r$ must be the same. We define
$
  S_l equiv 1 + 2 i k f_l (k)
$
then we require
$
  abs(S_l (k)) =^! 1
$
this is the unitarity relation. We write
$
  S_l = e^(2 i delta_l)
$
giving
$
  f_l (k) = (e^(2 i delta_l) - 1)/(2 i k) = (e^(i delta_l) sin delta_l)/k
$
The scattering amplitude is then
$
  f(theta) = 1/k sum_(l=0)^oo (2 l +1) e^(i delta_l) sin delta_l P_l (cos theta)
$
and we derived it by assuming rotational invariance and unitarity. With the latter being a sacred principle. We can also find the cross-section
$
  sigma_"tot" & = integral abs(f(theta))^2 dd(Omega) \
              & = (4 pi)/k^2 sum_l (2 l +1) sin^2 delta_l
$
this is nice! And is also in agreement with the optical theorem.

As an example we consider hard sphere scattering
$
  V = cases(oo &" for" r < R, 0 &" otherwise")
$
physically the expansion should be a plane-wave for $r > R$ but for $r <= R$ is should vanish. The current expansion we have
$
  e^(i bold(k) dot bold(x)) = sum_l (2l + 1) i^l j_l (k r) P_l (cos theta)
$
does not vanish for $r <= R$. So we need some other solution that vanishes and that satisfies the proper boundary condition.

We look for solutions using the other Bessel function $n_l$. We write
$
  braket(bold(x), psi^((+))) &= 1/(2 pi)^(3\/2) sum_l i^l ( 2l + 1) A_l (k r) P_l (cos theta)
$
we need $evaluated(A_l (k r))_(r = R) = 0$ and
$
  A_l equiv c_l^((1)) j_l (k r) + c_l^((2)) n_l (k r)
$
in the large $r$ limit we had
$
  braket(bold(x), psi^((+))) &= 1/(2 pi)^(3\/2) sum_l (2l +1) P_l/(2 i k) [(1 + 2 i k f_l (k)) e^(i k r)/r - e^(-i(k r- l pi))/r]
$
by comparing with the previous we see
$
  evaluated(A_l)_("large" r) = (e^(2 i delta_l) e^(i k r))/(2 i k r) - e^(-i (k r - l pi))/(2 i k r)
$
expanding $j_l$ and $n_l$ we find
$
  A_l = e^(i delta_l) underbracket([cos delta_l j_l (k r) - sin delta_l n_l (k r)], -> 0 "at" r=R)
$
by the boundary condition
$
  tan delta_l = (j_l (k R))/(n_l (k R))
$
and we are done! For $l = 0$ we find
$
  tan delta_0 = - tan k R
$
so $delta_0 = - k R$ giving
$
  A_(l=0) prop (sin k r)/(k r) cos delta_0 + (cos k r)/(k r) sin delta_0 = 1/(k r) sin(k r + delta_0)
$
For $k r << 1$ we find
$
  tan delta_l tilde.eq (- (k R)^(2l + 1))/((2 l + 1) [(2l-1)!!]^2)
$
and for $l eq.not 0$ this is tiny and we only care about $l = 0$. Then
$
  sigma_"tot" & = (4 pi)/k^2 sum_l (2l + 1) sin^2 delta_l \
              & = (4 pi)/k^2 sin^2 delta_0 \
              & tilde.eq 4 pi R^2
$
which is four times the geometric cross-section!

== Eikonal approximation

