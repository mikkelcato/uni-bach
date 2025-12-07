//**** init-ting
#import "@preview/physica:0.9.5": *
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
