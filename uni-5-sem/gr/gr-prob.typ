//**** init-ting
#import "@preview/physica:0.9.5": *
#import "temp.typ": *


#show: thmrules.with(qed-symbol: $square$)
#show: note.with(
  title: [
    *problems in general relativity and cosmology*
  ],
  authors: (
    (
      name: "mkh",
    ),
  ),
  abstract: [
    Solved problem sets.
  ],
)

= Lecture 1 & 2
== Relativistic rocket
A spaceship propels itself in a straight line by giving portions of its mass a constant (backward) velocity $u$ relative to its instantaneous rest frame. It continues to do so until it reaches a velocity $v$ relative to its initial rest-frame. Show that
$
  m_i/m_f = ((c+v)/(c-v))^(c/(2u))
$
Use the result to find the classical result and the first order relativistic correction.

We'll focus on the lab frame, wherein the rocket at some time has velocity $v'$. In the instantaneous rest frame of the rocket the exhaust has velocity $-u$, so in the lab frame it would have velocity
$
  overline(u) = (v-u)/(1-(u v)/c^2)
$
take the rocket to have $p^mu$ and mass $m + dd(m)$ before giving up some mass and $p^(mu) + dd(p)^mu$ after giving up the mass $dd(m)$---the exhaust has momentum $p^mu_"ex"$ after, thus by momentum conservation,
$
  p^mu = p^mu + dd(p)^mu + p_"ex"^mu => dd(p)^mu + p_"ex"^mu = 0
$
where
$
   dd(p)^mu & = vecrow(dd((gamma m c)), dd((gamma m v'))) \
  p_"ex"^mu & = vecrow(gamma_u dd(m) c, gamma_u dd(m) overline(u))
$
so we get two equations
$
   c dd((gamma m)) & = gamma_u c dd(m) \
  dd((gamma m v')) & = gamma_u dd(m) overline(u)
$
the first gives $gamma_u dd(m) = dd((gamma m))$, so
$
                   dd((gamma m v')) & = dd((gamma m)) overline(u) \
  v' dd((gamma m)) + gamma m dd(v') & = dd((gamma m)) overline(u) \
     dd((gamma m)) (overline(u)-v') & = gamma m dd(v') \
            dd((gamma m))/(gamma m) & = dd(v')/(overline(u)-v')
$
this can be integrated
$
  integral_(gamma_i m_i)^(gamma_f m_f) dd((gamma m))/(gamma m) &= integral_0^v dd(v')/(overline(u)-v') \
  ln((gamma_f m_f)/(gamma_i m_i)) &= integral_0^v dd(v')/((v'-u)/(1-u v'\/c^2)-v') \
  ln((gamma m_f)/m_i) &= 1/u integral_0^v ((c^2-u v') dd(v)')/( v'^2-c^2) \
  ln((gamma m_f)/m_i) &= 1/u [- c^2 integral_0^v dd(v)'/(c^2 - v'^2) +u integral_0^v (v' dd(v)')/(c^2 - v'^2)]
$
we solve each integral, starting with the second. We let $w = c^2 - v'^2 => -dd(w)\/2 = v' dd(v)'$ giving
$
  I_2 & = -1/(2) integral dd(w)/w \
      & = -1/2 ln(w_2/w_1) = -1/2 ln (1 - v^2/c^2) \
      & = ln (1-v^2\/c^2)^(-1\/2) = ln gamma
$
and the first
$
  I_1 & = integral_0^v dd(v)'/(c^2-v'^2) \
      & = 1/c^2 integral_0^v dd(v)'/(1-v'^2\/c^2)
$
letting $v'\/c = x => dd(x) = dd(v)'\/c$, giving
$
  I_1 & = 1/c integral dd(x)/(1-x^2) \
      & = 1/(2c) integral ( 1/(1+x) + 1/(1-x)) dd(x) \
      & = 1/(2c) ( ln((1+x_2)/(1+x_1)) + ln((1-x_1)/(1-x_2))) \
      & = 1/(2c) (ln ((1+v\/c)/(1))+ ln((1)/(1-v\/c))) \
      & = 1/(2 c) ln ((1+ v\/c)/(1-v\/c))
$
So we obtain,
$
  ln (m_f/m_i) & = - c/(2 u) ln ((1+v\/c)/(1-v\/c)) \
  ln (m_i/m_f) & = c/(2 u) ln ((1+v\/c)/(1-v\/c)) \
       m_i/m_f & = ((1+v\/c)/(1-v\/c))^(c\/2u)
$
which is what we were asked to find.

Now to get the classical limit we assume $v\/c, u\/c << 1$,
$
  ln (m_i/m_f) = c/(2 u) ln ((1 + v\/c)/(1 - v\/c))
$
and we use
$
  ln ((1 + x)/(1 - x)) = 2 x + (2 x^3)/3 + (2 x^5)/5 + dots
$
giving
$
  ln((1+ v\/c)/(1-v\/c)) & = 2 v/c + 2/3 v^3/c^3 + 2/5 v^5/c^3 + dots \
                         & tilde.eq 2 v/c
$
so
$
  ln(m_i/m_f) & tilde.eq c/(2u) 2 v/c \
              & = v/u => m_i/m_f = exp(v/u)
$
which is the classical rocket equation. The first correction would be
$
  ln((1+v\/c)/(1-v\/c)) tilde.eq 2 v/c + 2/3 v^3/c^3
$
giving
$
  ln(m_i/m_f) & tilde.eq v/u + 1/3 v/(u) v^2/c^2 \
              & = v/u (1+1/3 v^2/c^2)
$
so
$
  m_i/m_f & tilde.eq exp[v/u (1 + 1/3 v^2/c^2)] \
          & =^"or" e^(v\/u)exp(v^3/(3c^2u))
$
so the first relativistic correction is a factor
$
  "correction" = exp(v^3/(3 c^2 u))
$

== Rotating metric
Consider the metric defined by
$
  dd(tau)^2 = - c^2 dd(t)^2 + dd(x)^2 + dd(y)^2 + dd(z)^2 tilde eta_(mu nu) = eta^(mu nu)
$
find the metric in the rotating coordinate system defined by the transformation
$
  t' & = t \
  x' & = sqrt(x^2 + y^2) cos(phi - omega t) \
  y' & = sqrt(x^2+y^2) sin(phi - omega t) \
  z' & = z
$
where $phi = arctan(y\/x)$ and $omega$ is a constant.

We use
$
  g'^(alpha beta) = pdv(x'^alpha, x^gamma) pdv(x'^beta, x^delta) g^(gamma delta)
$
with $g^(gamma delta) = eta^(gamma delta)$. We need to find all the derivatives,
$
  pdv(x'^0, x^0) & = 1 \
  pdv(x'^0, x^1) & = pdv(x'^0, x^2) = pdv(x'^0, x^3) = 0
$
$
  pdv(x'^3, x^0) & = pdv(x'^3, x^1) = pdv(x'^3, x^2) = 0 \
  pdv(x'^3, x^3) & = 1
$
now the annoying ones
$
  pdv(x'^1, x^0) &= sqrt(x^2+y^2) sin(phi-omega t) omega \
  pdv(x'^1, x^1) &= x (x^2+y^2)^(-1\/2) cos(phi - omega t) - sqrt(x^2+y^2) sin(phi - omega t) 1/(1+y^2\/x^2) (-y/x^2) \
  &= x/sqrt(x^2+y^2) cos(phi - omega t) + y/sqrt(x^2+y^2) sin(phi-omega t) \
  pdv(x'^1, x^2) &= y/(sqrt(x^2+y^2)) cos(phi-omega t) -sqrt(x^2+y^2) sin(phi - omega t) x/(x^2+y^2) \
  &= y/sqrt(x^2+y^2) cos(phi - omega t) - x/sqrt(x^2+y^2) sin(phi-omega t) \
  pdv(x'^1, x^3) &= 0 \
  pdv(x'^2, x^0) &= -sqrt(x^2+y^2) cos(phi - omega t) omega \
  pdv(x'^2, x^1) &= x/sqrt(x^2+y^2) sin(phi - omega t) - sqrt(x^2+y^2) cos(phi-omega t) y/(x^2+y^2) \
  &= x/sqrt(x^2+y^2) sin(phi-omega t)- y/sqrt(x^2+y^2) cos(phi-omega t) \
  pdv(x'^2, x^2) &= y/sqrt(x^2+y^2) sin(phi-omega t) + sqrt(x^2+y^2) cos(phi-omega t) x/(x^2+y^2) \
  &= y/sqrt(x^2+y^2) sin(phi-omega t) + x/sqrt(x^2+y^2) cos(phi-omega t) \
  pdv(x'^2, x^3) &= 0
$
Now we can find the metric elements.
$
  g'^(alpha beta) = pdv(x'^alpha, x^gamma) pdv(x'^beta, x^delta) g^(gamma delta)
$
so
$
  g'^(00) & = pdv(x'^0, x^gamma) pdv(x'^0, x^delta) g^(gamma delta) \
  & = (pdv(x'^0, x^0))^2 g^(00) = -1 \
  g'^(01) & = pdv(x'^0, x^gamma) pdv(x'^1, x^delta) g^(gamma delta) \
  &= pdv(x'^0, x^0) (pdv(x'^1, x^0)g^(00)+pdv(x'^1, x^1)g^(01)+pdv(x'^1, x^2)g^(02)+pdv(x'^1, x^3)g^(03)) \
  &= - pdv(x'^1, x^0) = -sqrt(x^2+y^2) sin(phi-omega t) omega \
  g'^(02) &= pdv(x'^0, x^gamma) pdv(x'^2, x^delta) g^(gamma delta) \
  &= pdv(x'^0, x^0) (pdv(x'^2, x^0)g^(00)+pdv(x'^2, x^1)g^(01)+pdv(x'^2, x^2)g^(02)+pdv(x'^2, x^3)g^(03)) \
  &= - pdv(x'^2, x^0) = sqrt(x^2+y^2) cos(phi - omega t) omega \
  g'^(03) &= - pdv(x'^3, x^0) = 0 \
  g'^(11) &= pdv(x'^1, x^gamma) pdv(x'^1, x^delta) g^(gamma delta) \
  &= -(pdv(x'^1, x^0))^2 + (pdv(x'^1, x^1))^2 + (pdv(x'^1, x^2))^2 + (pdv(x'^1, x^3))^2 \
  &= -(x^2+y^2) sin^2(phi-omega t) omega^2 + (x/sqrt(x^2+y^2) cos(phi - omega t) + y/sqrt(x^2+y^2) sin(phi-omega t))^2 \
  &+ (y/sqrt(x^2+y^2) cos(phi - omega t) - x/sqrt(x^2+y^2) sin(phi-omega t))^2 \
  &= 1-omega^2(x^2+y^2)sin^2 (phi-omega t) \
  g'^(12) &= pdv(x'^1, x^gamma) pdv(x'^2, x^delta) g^(gamma delta) \
  &= - pdv(x'^1, x^0) pdv(x'^2, x^0) + pdv(x'^1, x^1) pdv(x'^2, x^1)+ pdv(x'^1, x^2) pdv(x'^2, x^2)+pdv(x'^1, x^3) pdv(x'^2, x^3) \
  &= omega^2(x^2+y^2) sin(phi-omega t)cos(phi - omega t)\
  &+ x^2/(x^2+y^2) cos sin - y^2/(x^2+y^2) cos sin - (x y)/(x^2+y^2) cos^2 + (x y)/(x^2+y^2)sin^2 \
  &+ y^2/(x^2+y^2) cos sin - x^2/(x^2+y^2) cos sin - (x y)/(x^2+y^2) sin^2 + (x y)/(x^2+y^2) cos^2 \
  &= omega^2(x^2+y^2) sin(phi-omega t)cos(phi - omega t) \
  g'^(13) &= 0 \
  g'^(22) &= -(pdv(x'^2, x^0))^2 + (pdv(x'^2, x^1))^2 + (pdv(x'^2, x^2))^2 \
  &= - omega^2(x^2+y^2) cos^2(phi - omega t) \
  &+ x^2/(x^2+y^2)sin^2 + y^2/(x^2+y^2)cos^2 - (2 x y)/(x^2+y^2) sin cos + y^2/(x^2+y^2) sin^2 + x^2/(x^2+y^2) cos^2 + (2 x y)/(x^2+y^2) sin cos \
  &= 1 - omega^2(x^2+y^2)cos^2(phi-omega t) \
  g'^(23) &= 0 \
  g'^(33) & = (pdv(x'^3, x^3))^2 g^(33) = 1
$
now we have every element, and:
$
  g'^(mu nu) &= mat(-1, - omega sqrt(x^2+y^2)sin(phi-omega t), omega sqrt(x^2+y^2)cos(phi-omega t), 0; - omega sqrt(x^2+y^2) sin(phi-omega t), 1-omega^2 (x^2 + y^2)sin^2 (phi-omega t), omega^2 (x^2+y^2) sin(phi-omega t) cos(phi-omega t), 0; omega sqrt(x^2+y^2)cos(phi-omega t), omega^2(x^2+y^2)sin(phi-omega t) cos(phi-omega t), 1-omega^2(x^2+y^2)cos^2(phi-omega t), 0; 0, 0, 0, 1) \
  &= mat(-1, - omega r sin(phi'), omega r cos(phi'), 0; - omega r sin(phi'), 1-omega^2 r^2 sin^2 (phi'), omega^2 r^2 sin(phi') cos(phi'), 0; omega r cos(phi'), omega^2 r^2 sin(phi') cos(phi'), 1-omega^2 r^2 cos^2(phi'), 0; 0, 0, 0, 1)
$
where we define $r^2 = x^2+y^2$ and $phi' = phi - omega t$.

You get the same using differentials:
$
  t' & = t \
  x' & = r cos(phi - omega t) \
  y' & = r sin(phi-omega t) \
  z' & = z
$
with $r = sqrt(x^2+y^2)$ and $phi = arctan(y\/x)$, we express $dd(x), dd(y)$ in terms of $dd(r),dd(phi),dd(t)$,
$
  dd(x) & = cos phi dd(r) - r sin phi dd(phi) \
  dd(y) & = sin phi dd(r) + r cos phi dd(phi)
$
and
$
  dd(x)^2+dd(y)^2 = dd(r)^2 + r^2 dd(phi)^2
$
we have in the rotating frame
$
  phi' = phi-omega t => phi = phi' + omega t
$
and $r' = r$, $t' = t$, $z' = z$,
$
  dd(phi) = dd(phi)' + omega dd(t)",  " dd(r)=dd(r)'",  " dd(z)=dd(z)'",  " dd(t)=dd(t)'
$
so
$
  dd(s)^2 &= -c^2 dd(t)^2 + dd(r)^2 + r^2 dd(phi)^2 + dd(z)^2 \
  &= -c^2 dd(t)^2 + dd(r)^2 + r^2 (dd(phi)' + omega dd(t))^2 + dd(z)^2 \
  &=^(c=1) - (1 - r'^2 omega^2) dd(t)'^2 + 2 r'^2 omega dd(phi)' dd(t)' + dd(r)'^2 + r'^2 dd(phi)'^2 + dd(z)'^2 \
  &= -(1 - omega^2 (x'^2 + y'^2)) dd(t)'^2 + 2 omega(x' dd(y)' - y' dd(x)') dd(t)' + dd(x)'^2 + dd(y)'^2 + dd(z)''
$
this gives
$
  g_(mu nu) &= mat(-(1-omega^2(x'^2+y'^2)), -omega y', omega x', 0; - omega y', 1, 0, 0; omega x', 0, 1, 0; 0, 0, 0, 1) \
  g^(mu nu) &= mat(-1, -omega y', omega x', 0; - omega y', 1-omega y'^2, omega^2 x' y', 0; omega x', omega^2 x' y', 1- omega^2 x'^2, 0; 0, 0, 0, 1)
$
which is equivalent to what was found before.

#pagebreak()
= Lecture 3 & 4
== short index problems
Evaluate
$
  g^(mu nu) g_(mu nu) = tensor(delta, mu, -mu) = 4
$

Evaluate
$
  g^(i mu) g_(mu i) = tensor(delta, i, -i) = 3
$

Prove
$
  partial_mu g^(nu rho) = - g^(nu kappa) g^(rho lambda) partial_mu g_(lambda kappa)
$
we start with
$
  g^(mu alpha) g_(alpha nu) = tensor(delta, mu, -nu)
$
take the derivative
$
  partial_rho (g^(mu alpha)) g_(alpha nu) + g^(mu alpha) partial_rho g_(alpha nu) &= 0 \
  partial_rho (g^(mu alpha)) g_(alpha nu) g^(nu beta) &= - g^(mu alpha) partial_rho (g_(alpha nu)) g^(nu beta) \
  partial_rho (g^(mu alpha)) tensor(delta, beta, -alpha) &= -g^(mu alpha) g^(beta nu) partial_rho g_(alpha nu) \
  partial_rho g^(mu beta) &= - g^(mu alpha) g^(beta nu) partial_rho g_(alpha nu) \
  partial_rho g^(mu nu) &= - g^(mu alpha) g^(nu beta) partial_rho g_(alpha beta)
$

Evaluate $T^((alpha beta)) omega_([alpha beta])$.

We have
$
  T^((alpha beta)) omega_([alpha beta]) &= 1/4 [T^(alpha beta) + T^(beta alpha)][omega_(alpha beta) - omega_(beta alpha)] \
  &= 1/4 (T^(alpha beta) omega_(alpha beta) - T^(alpha beta)omega_(beta alpha) + T^(beta alpha) omega_(alpha beta) - T^(beta alpha) omega_(beta alpha)) \
  &= 0
$
or
$
  T^((alpha beta)) omega_([alpha beta]) &= T^((beta alpha)) omega_([beta alpha]) = - T^((alpha beta)) omega_([alpha beta]) => T^((alpha beta)) omega_([alpha beta]) = 0
$

Evaluate $omega_([alpha beta]) U^alpha U^beta$.

$
  omega_([alpha beta]) U^alpha U^beta = omega_([beta alpha]) U^beta U^alpha = - omega_([alpha beta]) U^alpha U^beta => omega_([alpha beta]) U^a U^beta = 0
$
we always have
$
  omega_(alpha beta) U^alpha U^beta &= omega_((alpha beta)) U^alpha U^beta + omega_([alpha beta]) U^alpha U^beta \
  &= omega_((alpha beta)) U^alpha U^beta
$

We define the trace of $T_(alpha beta)$ as $g^(alpha beta) T_(alpha beta)$. The trace-free part is $T^"TF"_(alpha beta) = T_(alpha beta) - 1/D g_(alpha beta) T$ where $T$ is the trace of $T_(alpha beta)$. Verify that this is trace-free.

We have
$
  g^(alpha beta) T^"TF"_(alpha beta) &= g^(alpha beta) (T_(alpha beta) - 1/D g_(alpha beta) T) \
  &= T - 1/D g^(alpha beta) g_(alpha beta) T \
  &= T - 1/D tensor(delta, alpha, -alpha) T \
  &= T - 1/D D T \
  &= 0
$

#pagebreak()
== long index problem
Show that ($g equiv det g_(mu nu)$)
$
  g^(mu nu) pdv(g_(mu nu), x^lambda) = pdv(ln(-g), x^lambda)
$
use this to show that for a scalar $phi.alt$ we have
$
  square_g phi.alt = 1/sqrt(-g) partial_mu (sqrt(-g) g^(mu sigma) partial_sigma phi.alt)
$
where $square_g$ is the generally covariant generalization of the d'Alembertian $square equiv partial_mu partial^mu$.

We start with deriving the Jacobi formula
$
  pdv(det M, x) = det M Tr(M^(-1) pdv(M, x))
$
to do this we choose to start with the identity
$
  ln det M = tr ln M
$
taking the derivative of both sides
$
  1/(det M) pdv(det M, x) & = pdv(, x) tr ln M \
                          & = tr pdv(, x) (ln M) \
                          & = tr (M^(-1) pdv(M, x))
$
and we are done.

Now consider the RHS of
$
  g^(mu nu) pdv(g_(mu nu), x^lambda) = pdv(ln(-g), x^lambda)
$
rewriting
$
  pdv(ln (-g), x^lambda) &= 1/(-g) pdv((-g), x^lambda) \
  &= 1/g pdv(det g_(mu nu), x^lambda) \
  &= 1/g [det g_(mu nu) Tr (g^(mu nu) pdv(g_(mu nu), x^lambda)) ] \
  &= g^(mu nu) pdv(g_(mu nu), x^lambda)
$
so
$
  pdv(ln(-g), x^lambda) = g^(mu nu) pdv(g_(mu nu), x^lambda)
$

For the second part we need to calculate $nabla^mu nabla_mu phi.alt$. For scalars $nabla_mu phi.alt = partial_mu phi.alt$, so
$
  square_g phi.alt &= nabla^mu nabla_mu phi.alt \
  & = g^(mu nu) nabla_nu partial_mu phi.alt \
  & = g^(mu nu) partial_nu partial_mu phi.alt - g^(mu nu) tensor(Gamma, rho, -mu nu) partial_rho phi.alt
$
now we want to derive an expression for the Christoffel symbol. We have
$
  g^(mu nu) tensor(Gamma, rho, -mu nu) &= 1/2 g^(mu nu) g^(rho sigma) [ pdv(g_(nu sigma), x^mu) + pdv(g_(mu sigma), x^nu) - pdv(g_(mu nu), x^sigma)] \
  &= g^(rho sigma) g^(mu nu) pdv(g_(nu sigma), x^mu) - 1/2 g^(rho sigma) g^(mu nu) pdv(g_(mu nu), x^sigma)
$
for the first term we use
$
  g^(mu alpha) g_(alpha nu) &= tensor(delta, mu, -nu) \
  pdv(, x^rho) (g^(mu alpha) g_(alpha nu)) &= 0 \
  pdv(g^(mu alpha), x^rho) g_(alpha nu) g^(nu beta) + g^(mu alpha) pdv(g_(alpha nu), x^(rho)) g^(nu beta) &= 0 \
  pdv(g^(mu beta), x^rho) &= - g^(mu alpha) g^(nu beta) pdv(g_(alpha nu), x^rho) \
  pdv(g^(mu beta), x^rho) &= - g^(mu alpha) g^(beta nu) pdv(g_(alpha nu), x^rho) \
  pdv(g^(mu nu), x^rho) &= - g^(mu alpha) g^(nu beta) pdv(g_(alpha beta), x^rho)
$
so
$
  g^(mu nu) tensor(Gamma, rho, -mu nu) = - pdv(g^(rho mu), x^mu) - 1/2 g^(rho sigma) g^(mu nu) pdv(g_(mu nu), x^sigma)
$
for the second term we use the derived identity
$
  g^(mu nu) pdv(g_(mu nu), x^sigma) = pdv(ln (-g), x^sigma)
$
giving
$
  g^(mu nu) tensor(Gamma, rho, -mu nu) &= - pdv(g^(rho mu), x^mu) - 1/2 g^(rho sigma) pdv(ln(-g), x^sigma) \
  -g^(mu nu) tensor(Gamma, rho, -mu nu) &= pdv(g^(rho mu), x^mu) + g^(rho sigma) pdv(ln(sqrt(-g)), x^sigma) \
  &= pdv(g^(rho mu), x^mu) +g^(rho sigma)/sqrt(-g) pdv(sqrt(-g), x^sigma) \
  &= pdv(g^(mu rho), x^(mu)) + 1/sqrt(-g) g^(mu rho) pdv(sqrt(-g), x^mu) \
  &= 1/sqrt(-g) pdv(, x^mu) (g^(mu rho) sqrt(-g))
$
so
$
  square_g phi.alt &= g^(mu rho) partial_rho partial_mu phi.alt + 1/sqrt(-g) partial_mu (g^(mu rho) sqrt(-g)) partial_rho phi.alt \
  square_g phi.alt &= 1/sqrt(-g) partial_mu (sqrt(-g) g^(mu rho) partial_rho phi.alt)
$
