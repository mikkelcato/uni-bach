
#import "chpt-temp.typ": *
#show: chpt-note.with()

// 1st lecture
= Quantum Electrodynamics
The ABC theory described above is a perfectly legitimate QFT. However, it does not describe the real world. This is partly due to particles being described differently depending on their spin in real QFTs:
$
     "spin" 0 & tilde "Klein-Gordon equation" \
  "spin" 1\/2 & tilde "Dirac equation" \
     "spin" 1 & tilde "Proca equation"
$
Once we have established the Feynman rules, however, then the underlying field equation _becomes obsolete_. We would like to understand quantum electrodynamics meaning we need to understand the dynamics of $e^-$ and $gamma$.

== The Dirac equation
We can "derive" the Klein-Gordon equation by considering
$
  p^mu p_mu - m^2 c^2 = 0
$
and promoting $p_mu$ to an operator
$
  p_mu -> i hbar partial_mu
$
Then we obtain
$
  partial^mu partial_mu psi + (m^2 c^2)/hbar^2 psi = 0
$
or with $square equiv partial_mu partial^mu$ we have
$
  [square + (m^2 c^2)/hbar^2] psi = 0
$
which is the Klein-Gordon equation. A more correct derivation uses the Klein-Gordon Lagrangian
$
  cal(L)_"KG" = 1/2 partial_mu phi.alt partial^mu phi.alt - 1/2 m^2 phi.alt^2
$
which describes a free massive scalar field $phi.alt$.

When first conceived the Klein-Gordon equation appeared very problematic.#footnote[When trying to use it to describe $e^-$ people ran into problems. This is due to the Klein-Gordon equation describing particles with spin $0$. The Klein-Gordon equation is also second order in time leading to an apparent breaking of Born's statistical interpretation.] This lead Dirac to seek an equation consistent with
$
  p^mu p_mu - m^2 c^2 = 0
$
while also being first order in time. The idea is to "factor" the above
$
  p^mu p_mu - m^2 c^2 &= (beta^kappa p_kappa + m c)(gamma^lambda p_lambda - m c) \
  &= beta^kappa gamma^lambda p_kappa p_lambda - m c( beta^kappa - gamma^kappa) p_kappa - m^2 c^2
$
Then we need $beta^kappa = gamma^kappa$ and#footnote[Since $p_kappa p_lambda$ is symmetric.]
$
  p^mu p_mu &= gamma^kappa gamma^lambda p_kappa p_lambda \
  &= 1/2 (gamma^kappa gamma^lambda + gamma^lambda gamma^kappa) p_kappa p_lambda \
  &=^! eta^(kappa lambda) p_kappa p_lambda
$
which implies the $gamma^mu$ are matrices satisfying the Clifford algebra
$
  {gamma^mu,gamma^nu} = 2 eta^(mu nu)
$
We will use the $4 times 4$ matrices
$
  gamma^0 = mat(bb(1)_2, 0; 0, -bb(1)_2)";  " gamma^i = mat(0, sigma^i; -sigma^i, 0)
$
However, many other $gamma^mu$ also do the job. Then
$
  gamma^mu p_mu - m c = 0
$
and promoting $p_mu$ to an operator we obtain
$
  (i hbar gamma^mu partial_mu - m c) psi = 0
$
which is the Dirac equation. We will call $psi$ the Dirac spinor. Note, the components of $psi$ do not transform like a four-vector!

== Solutions to the Dirac equation
Consider $psi(bold(x), t) = psi(t)$.#footnote[Then the Dirac equation describes a state with $bold(p) = 0$ i.e. a particle at rest.] Then the Dirac equation becomes
$
  i hbar gamma^0 partial_t psi - m c^2 psi = 0
$
or
$
  mat(bb(1)_2, 0; 0, -bb(1)_2) vec(partial_t psi_A, partial_t psi_B) = -i (m c^2)/hbar vec(psi_A, psi_B)
$
implying
$
  partial_t psi_A & = -i ((m c^2)/hbar) psi_A \
  partial_t psi_B & = i ((m c^2)/hbar) psi_B
$
with solutions
$
  psi_A & = psi_A (0) exp[-i ((m c^2)/hbar) t] \
  psi_B & = psi_B (0) exp[+i ((m c^2)/hbar) t]
$
We recognize the familiar
$
  "time dependence" tilde exp(-(i E t)/hbar)
$
Then $psi_A$ corresponds to a particle with $E = m c^2$ as expected for a particle at rest. However, $psi_B$ corresponds to a particle with negative energy $E = - m c^2$! We take these negative energy solutions to represent antiparticles with positve energy. Then $psi_A$ describes $e^-$ while $psi_B$ describes $e^+$. Both being two-component spinors which can nicely describe a spin $1/2$ system. We have found four independent solutions with $bold(p)=0$
$
  psi^((1)) & = exp[-i ((m c^2)/hbar) t] vec(1, 0, 0, 0) tilde "spin-up" e^- \
  psi^((2)) & = exp[-i ((m c^2)/hbar) t] vec(0, 1, 0, 0) tilde "spin-down" e^- \
  psi^((3)) & = exp[+i ((m c^2)/hbar) t] vec(0, 0, 1, 0) tilde "spin-down" e^+ \
  psi^((4)) & = exp[+i ((m c^2)/hbar) t] vec(0, 0, 0, 1) tilde "spin-up" e^+
$
Consider plane-wave solutions
$
  psi(x) = a e^(-i k^mu x_mu) u(k)
$
implying
$
  partial_mu psi = - i k_mu psi
$
so
$
  (hbar gamma^mu k_mu - m c) u = 0
$
Using
$
  gamma^mu k_mu = mat(k^0, - bold(k) dot bold(sigma); bold(k) dot bold(sigma), -k^0)
$
we find
$
  u_A & = 1/(k^0-m c\/hbar) (bold(k) dot bold(sigma)) u_B \
  u_B & = 1/(k^0 + m c\/hbar) (bold(k) dot bold(sigma)) u_A
$
or
$
  u_A &= 1/((k^0)^2 - (m c \/hbar)^2) (bold(k) dot bold(sigma))^2 u_A \
  &=^((bold(k) dot bold(sigma))^2 = bold(k)^2) bold(k)^2/((k^0)^2 - (m c\/hbar)^2) u_A
$
implying#footnote[or the trivial solution $u_A = u_B = 0$]
$
  k^mu k_mu = ((m c)/hbar)^2
$
Then $hbar k^mu$ must be a four-vector whose square is $m^2 c^2$ implying
$
  hbar k^mu = plus.minus p^mu
$
with $+$ being associated with particle states and $-$ being associated with antiparticle states. We can now construct our four independent solutions:
$
  u_A &= vec(1, 0) => u_B = (bold(p) dot bold(sigma))/(p^0 + m c) vec(1, 0) = c/(E + m c^2) vec(p_z, p_z+i p_y) \
  u_A &= vec(0, 1) => u_B= c/(E+m c^2) vec(p_x-i p_y, -p_z) \
  u_B &= vec(1, 0) => u_A = c/(E+m c^2) vec(p_z, p_x+i p_y) \
  u_B &= vec(0, 1) => u_A = c/(E+m c^2) vec(p_x-i p_y, -p_z)
$
We need the $+$ sign in the first two, otherwise the $u_B$ would blow up as $bold(p) -> 0$. These are again our particle states. We will normalize these by
$
  u^dagger u = (2 E)/c
$
We find the normalisation
$
  N equiv sqrt((E+m c^2)/c)
$
and the canonical solutions become
$
  u^((1)) &= N vec(1, 0, (c p_z)/(E+m c^2), (c(p_x+i p_y))/(E+m c^2))";  " u^((2)) = N vec(0, 1, (c(p_x-i p_y))/(E+m c^2), (-c p_z)/(E+m c^2)) \
  v^((1)) &= N vec((c(p_x-i p_y))/(E+m c^2), (-c p_z)/(E+m c^2), 0, 1)";  " v^((2)) = - N vec((c p_z)/(E+m c^2), (c(p_x+i p_y))/(E+m c^2), 1, 0)
$
with
$
  psi & = a exp[-(i p_mu x^mu)/hbar] u tilde "particles" \
  psi & = a exp[+(i p_mu x^mu)/hbar] v tilde "antiparticles"
$
The particle states satisfy
$
  (gamma^mu p_mu - m c)u = 0
$
However, the antiparticle states satisfy
$
  (gamma^mu p_mu + m c) v = 0
$

== Bilinear covariants
We would like to construct a scalar from $psi$. This is done by introducing the adjoint $overline(psi)$ defined by
$
  overline(psi) equiv psi^dagger gamma^0
$
Then
$
  overline(psi) psi = psi^dagger gamma^0 psi tilde "scalar"
$
We can also construct a pseudoscalar by introducing $gamma^5$
$
  gamma^5 equiv i gamma^0 gamma^1 gamma^2 gamma^3
$
with
$
  {gamma^mu, gamma^5} = 0
$
Then
$
  overline(psi) gamma^5 psi tilde "pseudoscalar"
$
Likewise,
$
          overline(psi) gamma^mu psi & tilde "vector" \
  overline(psi) gamma^mu gamma^5 psi & tilde "pseudovector" \
     overline(psi) sigma^(mu nu) psi & tilde "antisymmetric tensor"
$
with
$
  sigma^(mu nu) equiv i/2 (gamma^mu gamma^nu - gamma^mu gamma^nu)
$
Any other _biliniear_ will always reduce to one of the above. One can further show that
$
  j^mu = overline(psi) gamma^mu psi
$
is conserved and leads to charge!

== $gamma$ and Maxwell's equations
We define the four-potential
$
  A^mu = vecrow(V, bold(A))
$
and the four-current
$
  J^mu = vecrow(c rho, bold(J))
$
Then Maxwell's equations become#footnote[The inhomogeneous ones.]
$
  partial_mu F^(mu nu) = (4 pi)/c J^nu
$
where we have introduced the antisymmetric Faraday tensor
$
  F^(mu nu) = partial^mu A^nu - partial^nu A^mu
$
Note, this implies
$
  partial_mu J^mu = 0
$
We see Maxwell's equations are unchanged under gauge transformations
$
  A_mu -> A_mu + partial_mu lambda
$
We refer to this as having gauge freedom. We can use this to impose
$
  partial_mu A^mu = 0
$
which is called the Lorentz gauge. With this gauge Maxwell's equations become
$
  square A^mu = (4 pi)/c J^mu
$
This still does not uniquely specify $A^mu$ since any gauge transformation with
$
  square lambda = 0
$
will still leave our equations unchanged. We will remedy this by imposing
$
  A^0 = 0
$
in empty space where $J^mu = 0$. Then the Lorentz gauge becomes
$
  div bold(A) = 0
$
which is called the Coulomb gauge. Whenever we single out a component like $A^0$ we break the otherwise manifest Lorentz invariance. Then whenever we perform a Lorentz transformation we are forced to also perform a gauge transformation in order to restore the Coulomb gauge.

$A^mu$ becomes the wave function of $gamma$ in quantum electrodynamics. The free $gamma$ satisfies the above with $J^mu = 0$
$
  square A^mu = 0
$
which we see is the Klein-Gordon equation for a particle with $m = 0$. Consider plane-wave solutions
$
  A^mu (x) = a exp[-(i p_mu x^mu)/hbar] epsilon.alt^mu (p)
$
where $p^mu = vecrow(E\/c, bold(p))$ and $epsilon.alt^mu$ is the _polarization vector_.#footnote[$epsilon.alt^mu$ characterises the spin of $gamma$.] We find the constraint
$
  p^mu p_mu = 0
$
as we would expect. We also have by the Lorentz gauge
$
  p^mu epsilon.alt_mu = 0
$
and in the Coulomb gauge
$
  epsilon.alt^0 = 0
$
implying
$
  bold(epsilon.alt) dot bold(p) = 0
$
Then $bold(epsilon.alt)$ is perpendicular to the direction of propagation. We say a free $gamma$ is _transversely polarized_. There are always two linearly independent $bold(epsilon.alt)$ perpendicular to $bold(p)$. As an example if $bold(p)$ is along $hat(z)$ then
$
  bold(epsilon.alt)^((1)) = vecrow(1, 0, 0)";  " bold(epsilon.alt)^((2)) = vecrow(0, 1, 0)
$
We would expect a $gamma$ has three spin states. However, only massive particles with spin $s$ admit $2 s + 1$ different spin orientations. All massless particles only have two.#footnote[expect for $s=0$ which only has one.] Then along the direction of motion the $gamma$ can have
$
  m_s = plus.minus s tilde h = plus.minus 1
$

// 2nd lecture
== Feynman rules for QED
We briefly summarize what we have found above.

We found that free $e^-$ and $e^+$ with $p = vecrow(E\/c, bold(p))$ are represented by
$
  psi(x) & = a exp[-(i p_mu x^mu)/hbar] u^((s)) (p) tilde e^- \
  psi(x) & = a exp[+(i p_mu x^mu)/hbar] v^((s)) (p) tilde e^+
$
with $s = 1,2$. The spinors $u^((s))$ and $v^((s))$ satisfy
$
  (gamma^mu p_mu - m c) u = 0";  " (gamma^mu p_mu + m c) v = 0
$
and their adjoints satisfy
$
  overline(u) (gamma^mu p_mu - m c)= 0";  " overline(v) ( gamma^mu p_mu + m c) = 0
$
These are orthogonal
$
  overline(u)^((1)) u^((2)) = 0";  " overline(v)^((1)) v^((2)) = 0
$
normalised
$
  overline(u) u = 2 m c";  " overline(v) v = - 2 m c
$
and complete
$
  sum u^((s)) overline(u)^((s)) = gamma^mu p_mu + m c";  " sum v^((s)) overline(v)^((s)) = gamma^mu p_mu - m c
$
We have already determined a convenient set ${u^((i)), v^((i))}$.

We also found that free $gamma$ with $p = vecrow(E\/c, bold(p))$ are represented by
$
  A_mu (x) = a exp[-(i p_mu x^mu)/hbar] epsilon.alt_mu^((s))
$
where $s = 1,2$. The $epsilon.alt_mu^((s))$ satisfy
$
  p^mu epsilon.alt_mu = 0
$
are orthogonal
$
  epsilon.alt_mu^((1) *) epsilon.alt^((2) mu) = 0
$
and normalised#footnote[Note, these are arbitrary choices.]
$
  epsilon.alt^(mu *) epsilon.alt_mu = 1
$
The Coulomb gauge has
$
  epsilon.alt^0 = 0";  " bold(epsilon.alt) dot bold(p) = 0
$
and the $bold(epsilon.alt)$ are complete
$
  sum epsilon.alt_i^((s)) epsilon.alt_j^((s)*) = delta_(i j) - hat(p)_i hat(p)_j
$
We have already determined a convenient pair ${epsilon.alt^((1)), epsilon.alt^((2))}$.


We now consider a general Feynman diagram in QED
#figure(
  feynman(orientation: "vertical", (
    vertex("i1"),
    vertex("i2"),
    vertex("i3"),
    vertex("node", shape: "blob", size: 1),
    vertex("f4"),
    vertex("f5"),
    vertex("f6"),
    edge("i1", "node", label: $p_1$),
    edge("i2", "node", type: "boson", label: $p_2$),
    edge("i3", "node", label: $p_3$),
    edge("node", "f4", label: $p_4$),
    edge("node", "f5", type: "boson", label: $p_5$),
    edge("node", "f6", label: $p_6$),
  )),
)
Then the Feynman rules are:

1. Label the incoming- and outgoing momenta $p_1, dots, p_n$ as above and label internal momenta $q_1, dots$.#footnote[Note, the external momenta should go forward in time.]

2. Each external line carries a factor: $   e^- & tilde cases("incoming:" u, "outgoing:" overline(u)) \
    e^+ & tilde cases("incoming:" overline(v), "outgoing:" v) \
  gamma & tilde cases("incoming:" epsilon.alt_mu, "outgoing:" epsilon.alt_mu^*) $

3. Each vertex carries a factor $i g_e gamma^mu$. Where $g_e$ is related to the fundamental charge $e$ by $ g_e = e sqrt((4 pi)/(hbar c)) = sqrt(4 pi alpha) $

4. Each internal line carries a propagator: $ e^(plus.minus) & tilde (i gamma^mu q_mu + m c)/(q^2 - m^2 c^2) \
           gamma & tilde (- i eta_(mu nu))/q^2 $

5. Each vertex carries a $delta$ function $ (2 pi)^4 delta^((4)) (k_1 + k_2 + k_3) $ with $k_i$ being negative for outward momenta.

6. Each internal line carries an integral $ integral dd(q, 4)/(2 pi)^4 $

7. The result will carry a factor $ (2 pi)^4 delta^((4)) (p_1+p_2+ dots - p_n) $ replace this by $i$.

8. We include a $-$ between diagrams differing only by the interchange of two incoming or outgoing $e^(plus.minus)$#footnote[or of an incoming $e^(minus.plus)$ and outgoing $e^(plus.minus)$]

== $e^- mu^-$ scattering
Consider $e^- mu^-$ scattering. We have one diagram contributing to second-order
#figure(
  feynman(
    (
      vertex("i1"),
      vertex("i2"),
      vertex("t"),
      vertex("b"),
      vertex("f1"),
      vertex("f2"),
      edge("i1", "t", label: $mu^-$),
      edge("t", "f1", label: $mu^-$),
      edge("i2", "b", label: $e^-$),
      edge("b", "f2", label: $e^-$),
      edge("t", "b", type: "boson", label: $q$),
    ),
    orientation: "vertical",
  ),
)
with $q$ upwards. We find
$
  scr(M) &tilde integral dd(q, 4)/(2 pi)^4 overbracket([overline(u)^((s_3)) (p_3) i g_e gamma^mu u^((s_1)) (p_1)], e^- "line") (-i eta_(mu nu))/q^2 underbracket([overline(u)^((s_4)) (p_4) i g_e gamma^nu u^((s_2))(p_2)], mu^- "line") \ & times (2 pi)^4 delta^((4)) (p_2 - p_4 + q) (2 pi)^4 delta^((4)) (p_1 - p_3 - q)
$
Above we have constructed the "$e^-$ line" and "$mu^-$ line" _backwards_ and _connected_ them with the $gamma$ propagator. We do this to ensure the matrix multiplication is consistent. Then
$
  scr(M) &= [overline(u)^((s_3)) (p_3) i g_e gamma^mu u^((s_1)) (p_1)] (eta_(mu nu))/(p_1 - p_3)^2 [overline(u)^((s_4)) (p_4) i g_e gamma^nu u^((s_2)) (p_2)] \
  &= - (g_e^2)/(p_1-p_3)^2 [overline(u)^((s_3)) (p_3) gamma^mu u^((s_1)) (p_1)] [overline(u)^((s_4)) (p_4) gamma_mu u^((s_2)) (p_2)]
$

== $e^- e^-$ scattering
We have the same diagram as above. However, now we also have a second diagram where the $e^-$ emerging with $(p_3,s_3)$ comes from the $e^-$ with $(p_2,s_2)$ instead of the $e^-$ with $(p_1,s_1)$. We essentially twist the diagram above. The total amplitude is then
$
  scr(M) &= -g_e^2/(p_1-p_3)^2 [overline(u)(3) gamma^mu u(1)][overline(u)(4) gamma_mu u(2)] \ & + underbracket(g_e^2/(p_1-p_4)^2 [overline(u) (4) gamma^mu u(1)] [overline(u)(3) gamma_mu u(2)], 3 <-> 4)
$
where $-$ comes from rule $8$.

== $e^- e^+$ scattering
We again have two diagrams. One is similar to $e^- mu^-$ scattering
#figure(
  feynman(
    (
      vertex("i1"),
      vertex("i2"),
      vertex("t"),
      vertex("b"),
      vertex("f1"),
      vertex("f2"),
      edge("i1", "t", type: "antifermion", label: $e^+$),
      edge("t", "f1", type: "antifermion", label: $e^+$),
      edge("i2", "b", label: $e^-$),
      edge("b", "f2", label: $e^-$),
      edge("t", "b", type: "boson", label: $q$),
    ),
    orientation: "vertical",
  ),
)
We find
$
  scr(M)_1 &tilde integral dd(q, 4)/(2 pi)^4 [overline(u)(3) i g_e gamma^mu u(1)] (- i eta_(mu nu))/q^2 underbracket([overline(v)(2) i g_e gamma^nu v(4)], e^+ "line") \
  &times (2 pi)^4 delta^((4)) (p_1-p_3-q) (2 pi)^4 delta^((4)) (p_2 + q - p_4)
$
Note, when going _backwards_ along an antiparticle we are going forwards in time.#footnote[The order is always $"adjoint" times "interaction" times "spinor"$ with the momenta going _backwards_ from left to right.] Then
$
  scr(M)_1 &= - g_e^2/(p_1-p_3)^2 [overline(u)(3) gamma^mu u(1)][overline(v)(2) gamma_mu v(4)]
$
The second diagram is
#figure(
  feynman(
    (
      vertex("i1"),
      vertex("i2"),
      vertex("t"),
      vertex("b"),
      vertex("f1"),
      vertex("f2"),
      edge("i1", "t", type: "fermion", label: $e^-$),
      edge("b", "f1", type: "fermion", label: $e^-$),
      edge("i2", "t", type: "antifermion", label: $e^+$),
      edge("b", "f2", type: "antifermion", label: $e^+$),
      edge("t", "b", type: "boson", label: $q$),
    ),
  ),
)
We find
$
  scr(M)_2 &tilde integral dd(q, 4)/(2 pi)^4 [overline(v)(2) i g_e gamma^mu u(1)] (-i eta_(mu nu))/q^2 [overline(u)(3) i g_e gamma^nu v(4)] \
  &times (2pi)^4 delta^((4)) (p_1 + p_2 - q) (2 pi)^4 delta^((4)) (q - p_3-p_4) \
  scr(M)_2 &= - g_e^2/(p_1+p_2)^2 [overline(u)(3) gamma^mu v(4)][overline(v)(2) gamma_mu u(1)]
$
Also, exchanging the incoming $e^+$ with the outgoing $e^-$ in this diagram we recover the first diagram. Then by rule $8$ we find
$
  scr(M) &= - g_e^2/(p_1-p_3)^2 [overline(u)(3) gamma^mu u(1)][overline(v)(2) gamma_mu v(4)] \
  &+ g_e^2/(p_1+p_2)^2 [overline(u)(3) gamma^mu v(4)][overline(v)(2) gamma_mu u(1)]
$

== Compton scattering
We again have two diagrams. The first diagram is
#figure(feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("t"),
    vertex("b"),
    vertex("f1"),
    vertex("f2"),
    edge("i1", "t", type: "boson"),
    edge("t", "f1"),
    edge("i2", "b"),
    edge("b", "f2", type: "boson"),
    edge("b", "t", type: "fermion", label: $q$),
  ),
  orientation: "vertical",
))
We find
#let feyn(body) = math.cancel(angle: 15deg, body)


$
  scr(M)_1 &tilde integral dd(q, 4)/(2 pi)^4 epsilon.alt_mu (2) [overline(u)(4) i g_e gamma^mu (i (feyn(q) + m c))/(q^2-m^2 c^2) i g_e gamma^nu u(1)] epsilon.alt_nu (3)^* \
  &times (2 pi)^4 delta^((4)) (p_1-p_3-q) (2 pi)^4 delta^((4)) (p_2 + q - p_4)
$
Note, we contract each $epsilon.alt_mu$ with the index of the vertex where it was _created_. Then
$
  scr(M)_1 &= g_e^2/((p_1-p_3)^2 - m^2 c^2) [overline(u)(4) feyn(epsilon.alt) (2) (feyn(p)_1 - feyn(p)_3+ m c) feyn(epsilon.alt)(3)^* u(1)]
$
where we have defined
$
  feyn(a) equiv gamma^mu a_mu
$
Note, we will write $feyn(epsilon.alt^*) = gamma^mu (epsilon.alt_mu^*)$. The second diagram is

#figure(feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("t"),
    vertex("b"),
    vertex("f1"),
    vertex("f2"),
    edge("i1", "t"),
    edge("b", "f1", type: "boson"),
    edge("i2", "t", type: "boson"),
    edge("b", "f2"),
    edge("t", "b", label: $q$),
  ),
))

and gives
$
  scr(M)_2 = g_e^2/((p_1+p_2)^2 - m^2 c^2) [overline(u)(4) feyn(epsilon.alt)(3)^* (feyn(p)_1+feyn(p)_2+m c) feyn(epsilon.alt)(2) u(1)]
$
Then
$
  scr(M) = scr(M)_1 + scr(M)_2
$

== Computing $abs(scr(M))^2$
We typically do not know $s_i$ and $s_f$ exactly. We would then need
$
  expval(abs(scr(M))^2) equiv "average over" s_i "and sum over" s_f "of" abs(scr(M) (s_i -> s_f))^2
$
We would prefer computing $expval(abs(scr(M))^2)$ directly. Consider $e^- mu^-$ scattering. Then
$
  abs(scr(M))^2 = g_e^4/(p_1-p_3)^4 [overline(u)(3) gamma^mu (1)][overline(u)(4) gamma_mu u(2)][overline(u)(3) gamma^nu u(1)]^* [overline(u)(4) gamma_nu u(2)]^*
$
We see terms of the form
$
  G equiv [overline(u)(a) Gamma_1 u(b)][overline(u)(a) Gamma_2 u(b)]^*
$
with $a$ and $b$ representing momenta and spin, while $Gamma_1$ and $Gamma_2$ are $4 times 4$ matrices. We have
$
  [overline(u)(a) Gamma_2 u(b)]^* & = [u(a)^dagger gamma^0 Gamma_2 u(b)]^dagger \
  & = u(b)^dagger Gamma_2^dagger gamma^0 u(a) \
  & = u(b)^dagger gamma^0 gamma^0 Gamma_2^dagger gamma^0 u(a) \
  & = overline(u)(b) overline(Gamma)_2 u(a)
$
with
$
  overline(Gamma)_2 equiv gamma^0 Gamma_2^dagger gamma^0
$
Then
$
  G = [overline(u)(a) Gamma_1 u(b)][overline(u)(b) overline(Gamma)_2 u(a)]
$
We sum over $b$ orientations
$
  sum_(s_b) G &= overline(u) (a) Gamma_1 [sum_(s_b) u^((s_b)) (p_b) overline(u)^((s_b)) (p_b)] overline(Gamma)_2 u(a) \
  &=^"completeness" overline(u)(a) Gamma_1 (feyn(p)_b + m_b c) overline(Gamma)_2 u(a) \
  &= overline(u)(a) Q u(a)
$
with
$
  Q equiv Gamma_1 (feyn(p)_b + m_b c) overline(Gamma)_2
$
We sum over $a$ orientations
$
  sum_(s_a) sum_(s_b) G & = sum_(s_a) overline(u)^((s_a)) (p_a) Q u^((s_a)) (p_a) \
  &= sum_(s_a) sum_(i, j =1)^4 overline(u)^((s_a)) (p_a)_i Q_(i j) u^((s_a)) (p_a)_j \
  &= sum_(i,j=1)^4 Q_(i j) [sum_(s_a) u^((s_a)) (p_a) overline(u)^(s_a)(p_a)]_(j i) \
  &= sum_(i,j = 1)^4 Q_(i j) (feyn(p)_a + m_a c)_(j i) \
  &= tr [Q (feyn(p)_a+m_a c)]
$
implying
$
  sum_"spins" [overline(u)(a) Gamma_1 u(b)][overline(u)(a) Gamma_2 u(b)]^* = tr [Gamma_1 (feyn(p)_b + m_b c) overline(Gamma)_2 (feyn(p)_a + m_a c)]
$
which is sometimes called _Casimir's trick_. When either $u$ is replaced by a $v$ the corresponding $m_i$ switches sign.

Again for $e^- mu^-$ scattering we have $Gamma_2 = gamma^nu$ so
$
  overline(Gamma)_2 = gamma^0 gamma^(nu dagger) gamma^0 = gamma^nu
$
so applying Casimir's trick twice we find
$
  expval(abs(scr(M))^2) = g_e^4/(4 (p_1-p_3)^4) tr [gamma^mu (feyn(p)_1 + m c) gamma^nu (feyn(p)_3 + m c)] tr [gamma_mu (feyn(p)_2 + M c) gamma_nu (feyn(p)_4 + M c)]
$
where $m = m_e$ and $M = m_mu$. We include $1/4$ since we average over all initial spins, i.e. from $2 "particles" times 2 "orientations"$.

Casimir's trick reduces the computation of $abs(scr(M))^2$ to computing traces. Note,
$
  tr (A + B) & = tr A + tr B \
  tr alpha A & = alpha tr A \
      tr A B & = tr B A
$
Recall,
$
  eta_(mu nu) eta^(mu nu) & = 4 \
     {gamma^mu, gamma^nu} & = 2 eta^(mu nu) => {feyn(a),feyn(b)} = 2 a dot b
$
These imply the _contraction theorems_,
$
  gamma_mu gamma^mu &= 4 \
  gamma_mu gamma^nu gamma^mu &= - 2 gamma^nu => gamma_mu feyn(a) gamma^mu = - 2 feyn(a) \
  gamma_mu gamma^nu gamma^lambda gamma^mu &= 4 eta^(nu lambda) => gamma_mu feyn(a) feyn(b) gamma^mu = 4 a dot b \
  gamma_mu gamma^nu gamma^lambda gamma^sigma gamma^mu &= -2 gamma^sigma gamma^lambda gamma^nu => gamma_mu feyn(a)feyn(b)feyn(c) gamma^mu = - 2 feyn(a) feyn(b) feyn(c)
$
and the _trace theorems_,
$
  tr "odd" \# gamma^mu & = 0 \
  tr 1 & = 4 \
  tr gamma^mu gamma^nu &= 4 eta^(mu nu) => tr feyn(a) feyn(b) = 4 a dot b \
  tr gamma^mu gamma^nu gamma^lambda gamma^sigma &= 4 (eta^(mu nu) eta^(lambda sigma) - eta^(mu lambda) eta^(nu sigma) + eta^(mu sigma) eta^(nu lambda)) \
  => tr feyn(a)feyn(b)feyn(c)feyn(d) &= 4 [(a dot b )(c dot d) - (a dot c) (b dot d) + (a dot d) (b dot c)]
$
With $gamma^5 = i gamma^0 gamma^1 gamma^2 gamma^3$ we find
$
  tr gamma^5 &= 0 \
  tr gamma^5 gamma^mu gamma^nu &= 0 => tr gamma^5 feyn(a) feyn(b) = 0 \
  tr gamma^5 gamma^mu gamma^nu gamma^lambda gamma^sigma &= 4 i epsilon.alt^(mu nu lambda sigma) => tr gamma^5 feyn(a)feyn(b)feyn(c)feyn(d) = 4 i epsilon.alt^(mu nu lambda sigma) a_mu b_nu c_lambda d_sigma
$
As an example we compute the traces in $e^- mu^-$ scattering
$
  tr [gamma^mu (feyn(p)_1 + m c) gamma^nu (feyn(p)_3 + m c)] & = tr gamma^mu feyn(p)_1 gamma^nu feyn(p)_3 + m c underbracket([tr gamma^mu feyn(p)_1 gamma^nu + tr gamma^mu gamma^nu feyn(p)_3], 0 "by" \# gamma^mu "being odd") + (m c)^2 underbracket(tr gamma^mu gamma^nu, 4 eta^(mu nu))
$
We compute the first term by
$
  tr gamma^mu feyn(p)_1 gamma^nu feyn(p)_3 &= overbracket((p_1)_lambda (p_3)_sigma, "numbers") tr overbracket(gamma^mu gamma^lambda gamma^nu gamma^sigma, "matrices") \
  &= (p_1)_lambda (p_3)_sigma 4 (eta^(mu lambda) eta^(nu sigma) - eta^(mu nu) eta^(lambda sigma) + eta^(mu sigma) eta^(lambda nu)) \
  &= 4 [p_1^mu p_3^nu - eta^(mu nu) (p_1 dot p_3) + p_3^mu p_1^nu]
$
implying
$
  tr[gamma^mu (feyn(p)_1 + m c) gamma^nu (feyn(p)_3 + m c)] &= 4 {p_1^mu p_3^nu + p_3^mu p_1^nu + eta^(mu nu) [(m c)^2 - (p_1 dot p_3)]}
$
The second trace is similar so
$
  expval(abs(scr(M))^2) &= (4 g_e^4)/(p_1-p_3)^4 {p_1^mu p_3^nu + p_3^mu p_1^nu + eta^(mu nu) [(m c)^2 - (p_(1 kappa) p_3^kappa)]} times {p_(2 mu) p_(4 nu) + p_(4 mu) p_(2 nu) + eta_(mu nu) [(M c)^2 - (p_2 dot p_4)]} \
  &= (8 g_e^4)/(p_1-p_3)^4 [(p_1 dot p_2)(p_3 dot p_4) + (p_1 dot p_4) (p_2 dot p_3)- (p_1 dot p_3) (M c)^2 - (p_2 dot p_4) (m c)^2 + 2 (m M c^2)^2]
$

== Computing $Gamma$ and $sigma$
Consider $e^- mu^-$ scattering with $M >> m$. We assume the recoil of $mu^-$ can neglected, and would like to determine the cross-section in the lab frame with $M$ at rest. We know
$
  dv(sigma, Omega) = (hbar/(8 pi M c))^2 expval(abs(scr(M))^2)
$
We have
$
  p_1 = vecrow(E/c, bold(p)_1)";  " p_2 = vecrow(M c, bold(0))";  " p_3 = vecrow(E/c, bold(p)_3)";  " p_4 = vecrow(M c, bold(0))
$
and $abs(bold(p)_1) = abs(bold(p)_3) = abs(bold(p))$ with $bold(p)_1 dot bold(p)_3 = bold(p)^2 cos theta$. Then
$
                 (p_1-p_3)^2 & =- (bold(p)_1-bold(p)_3)^2
                               = - 4 bold(p)^2 sin^2 theta/2 \
               (p_1 dot p_3) & = m^2 c^2 + 2 bold(p)^2 sin^2 theta/2 \
  (p_1 dot p_2)(p_3 dot p_4) & = (p_1 dot p_4)(p_2 dot p_3) = (M E)^2 \
               (p_2 dot p_4) & = (M c)^2
$
implying
$
  expval(abs(scr(M))^2) = ((g_e^2 M c)/(bold(p)^2 sin^2 theta\/2))^2 [(m c)^2 + bold(p)^2 cos^2 theta/2]
$
Then
$
  dv(sigma, Omega) = ((alpha hbar)/(2 bold(p)^2 sin^2 theta\/2))^2 [(m c)^2 + bold(p)^2 cos^2 theta/2]
$
which is called the _Mott formula_. Assuming $bold(p)^2 << (m c)^2$ we recover the _Rutherford formula_
$
  dv(sigma, Omega) = (e^2/(2 m v^2 sin^2 theta\/2))^2
$
which is quite nice.

However, in QED decays do not exist! We may consider processes like $pi^0 -> gamma + gamma$ as decays, but these are also scatterings $q + overline(q) -> gamma+ gamma$. We will consider _pair annihilation_ of the form $e^+ + e^- -> 2 gamma$. We find two diagrams with
$
  scr(M)_1 &= g_e^2/((p_1-p_3)^2 - m^2 c^2) overline(v)(2) feyn(epsilon.alt)_4 (feyn(p)_1-feyn(p)_3 + m c) feyn(epsilon.alt)_3 u(1) \
  scr(M)_2 &= g_e^2/((p_1-p_4)^2-m^2 c^2) overline(v)(2) feyn(epsilon.alt)_3 (feyn(p)_1-feyn(p)_4 + m c) feyn(epsilon.alt)_4 u(1)
$
and $scr(M) = scr(M)_1 + scr(M)_2$. We assume $e^(plus.minus)$ are $tilde$ at rest initially. Then
$
  p_1 = m c vecrow(1, bold(0))";  " p_2 = m c vecrow(1, bold(0))";  " p_3 = m c vecrow(1, 0, 0, 1)";  " p_4 = m c vecrow(1, 0, 0, -1)
$
and
$
  (p_1-p_3)^2 - m^2 c^2 = (p_1-p_4)^2 - m^2 c^2 = - 2 (m c)^2
$
Also,
$
  feyn(p)_1 feyn(epsilon.alt)_3 &= - feyn(epsilon.alt)_3 feyn(p)_1 + 2 underbracket(p_1 dot epsilon.alt_3, 0 "in Coulomb gauge") \
  feyn(p)_3 feyn(epsilon.alt)_3 &= -feyn(epsilon.alt)_3 feyn(p)_3 + 2 underbracket(p_3 dot epsilon.alt_3, 0 "by Lorenz gauge")
$
implying
$
  (feyn(p)_1 - feyn(p)_3 + m c) feyn(epsilon.alt)_3 u(1) &= feyn(epsilon.alt)_3 (-feyn(p)_1 + feyn(p)_3 + m c) u(1) \
  &=^"by Dirac" feyn(epsilon.alt)_3 feyn(p)_3 u(1)
$
similarly
$
  (feyn(p)_1 - feyn(p)_4 + m c) feyn(epsilon.alt)_4 u(1) = feyn(epsilon.alt)_4 feyn(p)_4 u(1)
$
Then
$
  scr(M) = - g_e^2/(2 (m c)^2) overline(v)(2) [feyn(epsilon.alt)_4 feyn(epsilon.alt)_3 feyn(p)_3 + feyn(epsilon.alt)_3 feyn(epsilon.alt)_4 feyn(p)_4] u(1)
$
Note,
$
  feyn(p)_3 = m c( gamma^0 - gamma^3)";  " feyn(p)_4 = m c (gamma^0 + gamma^3)
$
and
$
  feyn(epsilon.alt) = - bold(epsilon.alt) dot bold(gamma) = - mat(0, bold(sigma) dot bold(epsilon.alt); - bold(sigma) dot bold(epsilon.alt), 0)
$
so
$
  feyn(epsilon.alt)_3 feyn(epsilon.alt)_4 = - mat((bold(sigma) dot bold(epsilon.alt)_3)(bold(sigma) dot bold(epsilon.alt)_4), 0; 0, (bold(sigma) dot bold(epsilon.alt)_3)(bold(sigma) dot bold(epsilon.alt)_4))
$
with
$
  (bold(sigma) dot bold(a))(bold(sigma) dot bold(b)) = bold(a) dot bold(b) + i bold(sigma) dot (bold(a) times bold(b))
$
implying
$
  (feyn(epsilon.alt)_4 feyn(epsilon.alt)_3 + feyn(epsilon.alt)_3 feyn(epsilon.alt)_4) &= - 2 bold(epsilon.alt)_3 dot bold(epsilon.alt)_4 \
  (feyn(epsilon.alt)_4 feyn(epsilon.alt)_3 - feyn(epsilon.alt)_3 feyn(epsilon.alt)_4) &= 2 i (bold(epsilon.alt)_3 times bold(epsilon.alt)_4) dot bold(Sigma)
$
with
$
  bold(Sigma) = mat(bold(sigma), 0; 0, bold(sigma))
$
Then
$
  scr(M) = g_e^2/(m c) overline(v)(2) [(bold(epsilon.alt)_3 dot bold(epsilon.alt)_4)gamma^0 + i (bold(epsilon.alt)_3 times bold(epsilon.alt)_4) dot bold(Sigma) gamma^3] u(1)
$
We consider the _singlet state_ of $e^- e^+$. Then
$
  scr(M)_"singlet" = 1/sqrt(2) (scr(M)_(arrow.t arrow.b) - scr(M)_(arrow.b arrow.t))
$
where $scr(M)_(arrow.t arrow.b)$ is found by using $u^((1))$ and $v^((2))$ given by
$
  u(1) = sqrt(2 m c) vec(1, 0, 0, 0)";  " overline(v)(2) = sqrt(2 m c) vecrow(0, 0, 1, 0)
$
i.e. spin-up for $e^-$ and spin-down for $e^+$. Then
$
  overline(v)(2) gamma^0 u(1) = 0";  " overline(v)(2) bold(Sigma) gamma^3 u(1) = - 2 m c hat(z)
$
so
$
  scr(M)_(arrow.t arrow.b) = - 2 i g_e^2 (bold(epsilon.alt)_3 times bold(epsilon.alt)_4)_z
$
Likewise,
$
  scr(M)_(arrow.b arrow.t) = - scr(M)_(arrow.t arrow.b)
$
implying the amplitude for annihilation of a stationary $e^- e^+$ pair into two $gamma$ emerging in $plus.minus hat(z)$ is
$
  scr(M)_"singlet" = -2 sqrt(2) i g_e^2 (bold(epsilon.alt)_3 times bold(epsilon.alt)_4)_z
$
Note, for the _triplet state_ we would have
$
  scr(M)_"triplet" = 1/sqrt(2) (scr(M)_(arrow.t arrow.b) + scr(M)_(arrow.b arrow.t)) = 0
$
implying it is forbidden! We have
$
  bold(epsilon.alt)_+ = - 1/sqrt(2) vec(1, i, 0)";  " bold(epsilon.alt)_- = 1/sqrt(2) vec(1, -i, 0)
$
Also, the $gamma$ spins must be opposite so we can have
$
  bold(epsilon.alt)_3 = bold(epsilon.alt)_+";  " bold(epsilon.alt)_4 = bold(epsilon.alt)_-
$
so
$
  bold(epsilon.alt)_3 times bold(epsilon.alt)_4 = i hat(k)
$
or in the other case
$
  bold(epsilon.alt)_3 times bold(epsilon.alt)_4 = -i hat(k)
$
We need the antisymmetric combination#footnote[Where $arrow.t$ and $arrow.b$ refer to $gamma$ polarization.]
$
  scr(M)_"singlet" &= 1/sqrt(2) (scr(M)_(arrow.t arrow.b) - scr(M)_(arrow.b arrow.t)) \
  &= - 4 g_e^2
$
since again $scr(M)_"triplet" = 0$.

Then in the CM frame
$
  dv(sigma, Omega) = ((hbar c)/(8 pi (E_1+E_2)))^2 abs(bold(p)_f)/abs(bold(p)_i) abs(scr(M))^2
$
with
$
  E_1 = E_2 = m c^2";  " abs(bold(p)_f) = m c";  " abs(bold(p)_i) = m v
$
we find
$
  dv(sigma, Omega) = 1/(c v) ((hbar alpha)/m)^2 => sigma = (4 pi)/(c v) ((hbar alpha)/m)^2
$
We can find $tau$ by considering
$
  dv(sigma, Omega) = 1/L dv(N, Omega)
$
or $N = L sigma$ with
$
  L = rho v
$
The $e^-$ density for a _single atom_ is $abs(psi(0))^2$ so
$
  Gamma &equiv N = rho v sigma = v sigma abs(psi(0))^2 = (4 pi)/c ((hbar alpha)/m)^2 abs(psi(0))^2
$
The ground state has
$
  abs(psi(0))^2 = 1/pi ((alpha m c)/(2 hbar))^2
$
so
$
  tau = Gamma^(-1) = (2 hbar)/(alpha^5 m c^2) tilde 10^(-10) "s"
$
