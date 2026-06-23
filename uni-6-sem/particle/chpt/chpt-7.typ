#import "../../temp.typ": *
#show: chpt-note.with()

// 1st lecture
= Quantum Electrodynamics
== The Dirac equation
We can "derive" the Klein-Gordon equation#footnote[Which describes spin zero particles.] by considering
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
which is the Klein-Gordon equation. The Klein-Gordon equation is problematic since the equation is second order in time and only describes spin zero particles.

This lead Dirac to try and find an equation consistent with
$
  p^mu p_mu - m^2 c^2 = 0,
$
which is first order in time. We "factor" the above
$
  p^mu p_mu - m^2 c^2 &= (beta^kappa p_kappa + m c)(gamma^lambda p_lambda - m c) = beta^kappa gamma^lambda p_kappa p_lambda - m c( beta^kappa - gamma^kappa) p_kappa - m^2 c^2
$
We see $beta^kappa = gamma^kappa$ and
$
  p^mu p_mu &= gamma^kappa gamma^lambda p_kappa p_lambda = 1/2 (gamma^kappa gamma^lambda + gamma^lambda gamma^kappa) p_kappa p_lambda =^! eta^(kappa lambda) p_kappa p_lambda
$
implying the $gamma^mu$#footnote[We could write $gamma^mu_(alpha beta)$ with $alpha, beta$ being spinor indices.] are satisfy the Clifford algebra
$
  {gamma^mu,gamma^nu} = 2 eta^(mu nu)
$
We will use the following $gamma^mu$
$
  gamma^0 = mat(bb(1)_2, 0; 0, -bb(1)_2)",  " gamma^i = mat(0, sigma^i; -sigma^i, 0),
$
where the $sigma^i$ are the Pauli matrices. However, many other $gamma^mu$ also work.

We find
$
  gamma^mu p_mu - m c = 0,
$
and promoting $p_mu$ to an operator we obtain
$
  (i hbar gamma^mu partial_mu - m c) psi = 0,
$
which is the Dirac equation. We will call $psi$ the Dirac spinor. Note, the components of $psi$ do not transform like a four-vector!

== $bold(p)=0$ solutions
Now, consider $psi(bold(x), t) = psi(t)$ or $bold(p)=0$. We find
$
  i hbar gamma^0 partial_t psi - m c^2 psi = 0,
$
or
$
  mat(bb(1)_2, 0; 0, -bb(1)_2) vec(partial_t psi_A, partial_t psi_B) = -i (m c^2)/hbar vec(psi_A, psi_B),
$
implying
$
  partial_t psi_A & = -i ((m c^2)/hbar) psi_A",  " partial_t psi_B & = i ((m c^2)/hbar) psi_B,
$
which has solutions
$
  psi_A & = psi_A (0) exp[-i ((m c^2)/hbar) t]",  " psi_B & = psi_B (0) exp[+i ((m c^2)/hbar) t].
$
Now, consider the usual time-dependence
$
  "time dependence" tilde exp(-(i E t)/hbar),
$
implying $psi_A$ corresponds to a particle with $E = m c^2$. However, $psi_B$ corresponds to a particle with negative energy $E = - m c^2$! We say these correspond to antiparticles. We have found four independent solutions ($bold(p)=0$)
$
  psi^((1)) & = exp[-i ((m c^2)/hbar) t] vec(1, 0, 0, 0) tilde "spin-up" e^- \
  psi^((2)) & = exp[-i ((m c^2)/hbar) t] vec(0, 1, 0, 0) tilde "spin-down" e^- \
  psi^((3)) & = exp[+i ((m c^2)/hbar) t] vec(0, 0, 1, 0) tilde "spin-down" e^+ \
  psi^((4)) & = exp[+i ((m c^2)/hbar) t] vec(0, 0, 0, 1) tilde "spin-up" e^+
$

== Plane-wave solutions
Now, consider plane-wave solutions
$
  psi(x) = a e^(-i k^mu x_mu) u(k),
$
implying
$
  (hbar gamma^mu k_mu - m c) u = 0.
$
or
$
  u_A & = 1/(k^0-m c\/hbar) (bold(k) dot bold(sigma)) u_B",  " u_B & = 1/(k^0 + m c\/hbar) (bold(k) dot bold(sigma)) u_A.
$
We can find a condition on $k_mu$ by
$
  u_A &= (bold(k) dot bold(sigma))^2/((k^0)^2 - (m c \/hbar)^2) u_A &=^((bold(k) dot bold(sigma))^2 = bold(k)^2) bold(k)^2/((k^0)^2 - (m c\/hbar)^2) u_A
$
implying#footnote[or the trivial solution $u_A = u_B = 0$]
$
  k^mu k_mu = ((m c)/hbar)^2
$
We find $hbar k^mu$ is a four-vector whose square is $m^2 c^2$ implying
$
  hbar k^mu = plus.minus p^mu.
$
We say $hbar k^mu = +p^mu$ corresponds to particle states and $hbar k^mu = -p^mu$ corresponds to antiparticle states.#footnote[We have $tilde e^(minus.plus i k_mu x^mu)$ then acting with $p^mu = i hbar partial^mu$ makes this obvious.]

We find four solutions:
$
  u_A &= vec(1, 0) => u_B = (bold(p) dot bold(sigma))/(p^0 + m c) vec(1, 0) = c/(E + m c^2) vec(p_z, p_z+i p_y) \
  u_A &= vec(0, 1) => u_B= c/(E+m c^2) vec(p_x-i p_y, -p_z) \
  u_B &= vec(1, 0) => u_A = c/(E+m c^2) vec(p_z, p_x+i p_y) \
  u_B &= vec(0, 1) => u_A = c/(E+m c^2) vec(p_x-i p_y, -p_z)
$
Where $hbar k^mu = + p^mu$ in the first two. We normalise these by
$
  u^dagger u = (2 E)/c.
$
With
$
  N equiv sqrt((E+m c^2)/c),
$
we find the "canonical" solutions
$
  u^((1)) &= N vec(1, 0, (c p_z)/(E+m c^2), (c(p_x+i p_y))/(E+m c^2))",  " u^((2)) = N vec(0, 1, (c(p_x-i p_y))/(E+m c^2), (-c p_z)/(E+m c^2)) \
  v^((1)) &= N vec((c(p_x-i p_y))/(E+m c^2), (-c p_z)/(E+m c^2), 0, 1)",  " v^((2)) = - N vec((c p_z)/(E+m c^2), (c(p_x+i p_y))/(E+m c^2), 1, 0)
$
Now, we can write the plane-wave solutions#footnote[With motion along $z$ the $u^((1)), v^((1))$ are spin up, while the $u^((2)), v^((2))$ are spin down.]
$
  psi & = a exp[-(i p_mu x^mu)/hbar] u tilde "particles" \
  psi & = a exp[+(i p_mu x^mu)/hbar] v tilde "antiparticles"
$
Note, the particle states satisfy
$
  (gamma^mu p_mu - m c)u = 0.
$
However, the antiparticle states satisfy
$
  (gamma^mu p_mu + m c) v = 0.
$

== Bilinear covariants
We define the adjoint $overline(psi)$ by
$
  overline(psi) equiv psi^dagger gamma^0
$
With $overline(psi)$ we can make scalars
$
  overline(psi) psi = psi^dagger gamma^0 psi tilde "scalar"
$
We can also make pseudoscalars by defining $gamma^5$
$
  gamma^5 equiv i gamma^0 gamma^1 gamma^2 gamma^3,
$
which satisfies
$
  {gamma^mu, gamma^5} = 0
$
We can then make
$
  overline(psi) gamma^5 psi tilde "pseudoscalar"
$
Likewise, we can make more complicated stuff
$
          overline(psi) gamma^mu psi & tilde "vector" \
  overline(psi) gamma^mu gamma^5 psi & tilde "pseudovector" \
     overline(psi) sigma^(mu nu) psi & tilde "antisymmetric tensor"
$
with
$
  sigma^(mu nu) equiv i/2 (gamma^mu gamma^nu - gamma^mu gamma^nu)
$
All other biliniears will always reduce to one of these.

Note, one can show
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
We can then write Maxwell's equations as
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
which is called the Lorenz gauge. With this gauge Maxwell's equations become
$
  square A^mu = (4 pi)/c J^mu
$
This still does not uniquely specify $A^mu$ since any gauge transformation with
$
  square lambda = 0
$
will still leave our equations unchanged. We remedy this by imposing
$
  A^0 = 0
$
in empty space where $J^mu = 0$. Then the Lorenz gauge becomes
$
  div bold(A) = 0
$
which is called the Coulomb gauge. Note, now we have broken the manifest Lorentz invariance.

$A^mu$ becomes the wave function of the photon $gamma$ in quantum electrodynamics. The free photon $gamma$ satisfies the above with $J^mu = 0$
$
  square A^mu = 0,
$
which we see is the Klein-Gordon equation for a particle with $m = 0$. Now, consider plane-wave solutions
$
  A^mu (x) = a exp[-(i p_mu x^mu)/hbar] epsilon.alt^mu (p),
$
where $p^mu = vecrow(E\/c, bold(p))$ and $epsilon.alt^mu$ is the polarization vector.#footnote[$epsilon.alt^mu$ characterises the spin of the photon $gamma$.] We find $p^mu p_mu = 0$.

With the Lorenz gauge
$
  p^mu epsilon.alt_mu = 0,
$
and in the Coulomb gauge
$
  epsilon.alt^0 = 0,
$
implying
$
  bold(epsilon.alt) dot bold(p) = 0.
$
We find $bold(epsilon.alt)$ is perpendicular to $bold(p)$. We say a free photon $gamma$ is transversely polarised. We can always find two $bold(epsilon.alt)$ perpendicular to $bold(p)$. We take $bold(p) = p_z hat(z)$ then
$
  bold(epsilon.alt)^((1)) = vecrow(1, 0, 0)",  " bold(epsilon.alt)^((2)) = vecrow(0, 1, 0).
$

We find the photon $gamma$ has two "spin" states.#footnote[Which all massless particles do, except for $s=0$ which only has one.] So the photon $gamma$ can have
$
  m_s = plus.minus s tilde h = plus.minus 1
$
Note, the linear polarisation $bold(epsilon.alt)^((1,2))$ above are superpositions of states with definite $h$. We could have used
$
  bold(epsilon.alt)^((+)) = 1/sqrt(2) vec(1, i, 0)",  " bold(epsilon.alt)^((-)) = 1/sqrt(2) vec(1, -i, 0)
$

// 2nd lecture
== Summary
We have found particles $e^-$ and antiparticles $e^+$ with $p = vecrow(E\/c, bold(p))$ are represented by
$
  psi(x) & = a exp[-(i p_mu x^mu)/hbar] u^((s)) (p) tilde e^- \
  psi(x) & = a exp[+(i p_mu x^mu)/hbar] v^((s)) (p) tilde e^+
$
with $s = 1,2$. We know the spinors $u^((s))$ and $v^((s))$ satisfy
$
  (gamma^mu p_mu - m c) u &= 0",  " (gamma^mu p_mu + m c) v = 0 \
  overline(u) (gamma^mu p_mu - m c)&= 0",  " overline(v) ( gamma^mu p_mu + m c) = 0
$
and (orthogonal and normalised)
$
  overline(u)^((1)) u^((2)) = 0",  " & overline(v)^((1)) v^((2)) = 0 \
          overline(u) u = 2 m c",  " & overline(v) v = - 2 m c
$
We also have the completeness relations
$
  sum_s u^((s)) overline(u)^((s)) = gamma^mu p_mu + m c",  "& sum_s v^((s)) overline(v)^((s)) = gamma^mu p_mu - m c
$
with an example being the "canonical" ${u^((s)), v^((s))}$ above.

We have found photons $gamma$ with $p = vecrow(E\/c, bold(p))$ are represented by
$
  A_mu (x) = a exp[-(i p_mu x^mu)/hbar] epsilon.alt_mu^((s))
$
with $s = 1,2$.

The $epsilon.alt_mu^((s))$ satisfy (Lorenz gauge)
$
  p^mu epsilon.alt_mu = 0
$
and (orthogonal and normalised)
$
  epsilon.alt_mu^((1) *) epsilon.alt^((2) mu) & = 0 \
            epsilon.alt^(mu *) epsilon.alt_mu & = 1
$
Also, with the Coulomb gauge we have
$
  epsilon.alt^0 = 0";  " bold(epsilon.alt) dot bold(p) = 0
$
and the $bold(epsilon.alt)$ are complete
$
  sum epsilon.alt_i^((s)) epsilon.alt_j^((s)*) = delta_(i j) - hat(p)_i hat(p)_j
$
with an example being the ${bold(epsilon.alt)^((1)), bold(epsilon.alt)^((2))}$ above.

== Feynman rules for QED
Now, consider a general Feynman diagram in QED
#figure(
  feynman.feynman(orientation: "vertical", (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("i3"),
    feynman.vertex("node", shape: "blob", size: 1),
    feynman.vertex("f4"),
    feynman.vertex("f5"),
    feynman.vertex("f6"),
    feynman.edge("i1", "node", label: $p_1$),
    feynman.edge("i2", "node", type: "boson", label: $p_2$),
    feynman.edge("i3", "node", label: $p_3$),
    feynman.edge("node", "f4", label: $p_4$),
    feynman.edge("node", "f5", type: "boson", label: $p_5$),
    feynman.edge("node", "f6", label: $p_6$),
  )),
)
The Feynman rules are:

1. Label the incoming- and outgoing momenta $p_1, dots, p_n$ as above and label internal momenta $q_1, dots$.#footnote[Note, the external momenta should go forward in time.]

2. Each external line carries a factor: $   e^- & tilde cases("incoming:" u, "outgoing:" overline(u)) \
    e^+ & tilde cases("incoming:" overline(v), "outgoing:" v) \
  gamma & tilde cases("incoming:" epsilon.alt_mu, "outgoing:" epsilon.alt_mu^*) $

3. Each vertex carries a factor $i g_e gamma^mu$. Where $g_e$ is related to the fundamental charge $e$ by $ g_e = e sqrt((4 pi)/(hbar c)) = sqrt(4 pi alpha) $

4. Each internal line carries a propagator: $ e^(plus.minus) & tilde (i (gamma^mu q_mu + m c))/(q^2 - m^2 c^2) \
           gamma & tilde (- i eta_(mu nu))/q^2 $

5. Each vertex carries a $delta$ function $ (2 pi)^4 delta^((4)) (k_1 + k_2 + k_3) $ with $k_i$ being negative for outward momenta.

6. Each internal line carries an integral $ integral dd(q, 4)/(2 pi)^4 $

7. The result will carry a factor $ (2 pi)^4 delta^((4)) (p_1+p_2+ dots - p_n) $ replace this by $i$.

8. We include a $-$ between diagrams differing only by the interchange of two incoming or outgoing $e^(plus.minus)$#footnote[or of an incoming $e^(minus.plus)$ and outgoing $e^(plus.minus)$]

== $e^- mu^-$ scattering
Now, consider $e^- mu^-$ scattering. We have one diagram contributing to second-order
#figure(
  feynman.feynman(
    (
      feynman.vertex("i1"),
      feynman.vertex("i2"),
      feynman.vertex("t"),
      feynman.vertex("b"),
      feynman.vertex("f1"),
      feynman.vertex("f2"),
      feynman.edge("i1", "t", label: $mu^-$),
      feynman.edge("t", "f1", label: $mu^-$),
      feynman.edge("i2", "b", label: $e^-$),
      feynman.edge("b", "f2", label: $e^-$),
      feynman.edge("t", "b", type: "boson", label: $q$),
    ),
    orientation: "vertical",
  ),
)
with $q$ upwards. We find
$
  scr(M) &tilde integral dd(q, 4)/(2 pi)^4 overbracket([overline(u)^((s_3)) (p_3) i g_e gamma^mu u^((s_1)) (p_1)], e^- "line") (-i eta_(mu nu))/q^2 underbracket([overline(u)^((s_4)) (p_4) i g_e gamma^nu u^((s_2))(p_2)], mu^- "line") \ & times (2 pi)^4 delta^((4)) (p_2 - p_4 + q) (2 pi)^4 delta^((4)) (p_1 - p_3 - q)
$
Where the "$e^-$ line" and "$mu^-$ line" have been constructed backwards and "connected" with the photon $gamma$ propagator. We find
$
  scr(M) &= [overline(u)^((s_3)) (p_3) i g_e gamma^mu u^((s_1)) (p_1)] (eta_(mu nu))/(p_1 - p_3)^2 [overline(u)^((s_4)) (p_4) i g_e gamma^nu u^((s_2)) (p_2)] \
  &= - (g_e^2)/(p_1-p_3)^2 [overline(u)^((s_3)) (p_3) gamma^mu u^((s_1)) (p_1)] [overline(u)^((s_4)) (p_4) gamma_mu u^((s_2)) (p_2)]
$

== $e^- e^-$ scattering
Now, consider $e^- e^-$ scattering. We have the same diagram as above contributing. However, we also have a second diagram where the $e^-$ emerging with $(p_3,s_3)$ comes from the $e^-$ with $(p_2,s_2)$ instead of the $e^-$ with $(p_1,s_1)$. We essentially "twist" the diagram above. We find
$
  scr(M) &= underbracket(-, "by rule" 8)g_e^2/(p_1-p_3)^2 [overline(u)(3) gamma^mu u(1)][overline(u)(4) gamma_mu u(2)] \ & + underbracket(g_e^2/(p_1-p_4)^2 [overline(u) (4) gamma^mu u(1)] [overline(u)(3) gamma_mu u(2)], 3 <-> 4)
$

== $e^- e^+$ scattering
Now, consider $e^- e^+$ scattering. We again have two diagrams contributing. One is similar to $e^- mu^-$ scattering
#figure(
  feynman.feynman(
    (
      feynman.vertex("i1"),
      feynman.vertex("i2"),
      feynman.vertex("t"),
      feynman.vertex("b"),
      feynman.vertex("f1"),
      feynman.vertex("f2"),
      feynman.edge("i1", "t", type: "antifermion", label: $e^+$),
      feynman.edge("t", "f1", type: "antifermion", label: $e^+$),
      feynman.edge("i2", "b", label: $e^-$),
      feynman.edge("b", "f2", label: $e^-$),
      feynman.edge("t", "b", type: "boson", label: $q$),
    ),
    orientation: "vertical",
  ),
)
We find
$
  scr(M)_1 &tilde integral dd(q, 4)/(2 pi)^4 [overline(u)(3) i g_e gamma^mu u(1)] (- i eta_(mu nu))/q^2 underbracket([overline(v)(2) i g_e gamma^nu v(4)], e^+ "line") \
  &times (2 pi)^4 delta^((4)) (p_1-p_3-q) (2 pi)^4 delta^((4)) (p_2 + q - p_4)
$
Note, when going "backwards" along an antiparticle we are going forwards in time.#footnote[The order is always $"adjoint" times "interaction" times "spinor"$ when we follow "fermion-flow".] We find
$
  scr(M)_1 &= - g_e^2/(p_1-p_3)^2 [overline(u)(3) gamma^mu u(1)][overline(v)(2) gamma_mu v(4)]
$
The second contributing diagram is
#figure(
  feynman.feynman(
    (
      feynman.vertex("i1"),
      feynman.vertex("i2"),
      feynman.vertex("t"),
      feynman.vertex("b"),
      feynman.vertex("f1"),
      feynman.vertex("f2"),
      feynman.edge("i1", "t", type: "fermion", label: $e^-$),
      feynman.edge("b", "f1", type: "fermion", label: $e^-$),
      feynman.edge("i2", "t", type: "antifermion", label: $e^+$),
      feynman.edge("b", "f2", type: "antifermion", label: $e^+$),
      feynman.edge("t", "b", type: "boson", label: $q$),
    ),
  ),
)
We find
$
  scr(M)_2 &tilde integral dd(q, 4)/(2 pi)^4 [overline(v)(2) i g_e gamma^mu u(1)] (-i eta_(mu nu))/q^2 [overline(u)(3) i g_e gamma^nu v(4)] \
  &times (2pi)^4 delta^((4)) (p_1 + p_2 - q) (2 pi)^4 delta^((4)) (q - p_3-p_4) \
  scr(M)_2 &= - g_e^2/(p_1+p_2)^2 [overline(u)(3) gamma^mu v(4)][overline(v)(2) gamma_mu u(1)]
$
Here, exchanging the incoming positron $e^+$ with the outgoing electron $e^-$ we recover the first diagram. We find
$
  scr(M) &= - g_e^2/(p_1-p_3)^2 [overline(u)(3) gamma^mu u(1)][overline(v)(2) gamma_mu v(4)] \
  &underbracket(+, "by rule" 8) g_e^2/(p_1+p_2)^2 [overline(u)(3) gamma^mu v(4)][overline(v)(2) gamma_mu u(1)]
$

== Compton scattering
Now, consider $e^- gamma$ scattering. We again have two diagrams contributing. The first diagram is
#figure(feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("t"),
    feynman.vertex("b"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),
    feynman.edge("i1", "t", type: "boson"),
    feynman.edge("t", "f1"),
    feynman.edge("i2", "b"),
    feynman.edge("b", "f2", type: "boson"),
    feynman.edge("b", "t", type: "fermion", label: $q$),
  ),
  orientation: "vertical",
))
We find
#let feyn(body) = math.cancel(angle: 15deg, body)


$
  scr(M)_1 &tilde integral dd(q, 4)/(2 pi)^4 overbracket(epsilon.alt_mu (2) underbracket([overline(u)(4) i g_e gamma^mu (i (feyn(q) + m c))/(q^2-m^2 c^2) i g_e gamma^nu u(1)], e^- "line") epsilon.alt_nu (3)^*, "sandwich") \
  &times (2 pi)^4 delta^((4)) (p_1-p_3-q) (2 pi)^4 delta^((4)) (p_2 + q - p_4)
$
Here, we contract each $epsilon.alt_mu$ with the index of the vertex where it was "created". We find
$
  scr(M)_1 &= g_e^2/((p_1-p_3)^2 - m^2 c^2) [overline(u)(4) feyn(epsilon.alt) (2) (feyn(p)_1 - feyn(p)_3+ m c) feyn(epsilon.alt)(3)^* u(1)]
$
where we define
$
  feyn(a) equiv gamma^mu a_mu
$
Note, we will write $feyn(epsilon.alt^*) = gamma^mu (epsilon.alt_mu^*)$.

The second contributing diagram is

#figure(feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("t"),
    feynman.vertex("b"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),
    feynman.edge("i1", "t"),
    feynman.edge("b", "f1", type: "boson"),
    feynman.edge("i2", "t", type: "boson"),
    feynman.edge("b", "f2"),
    feynman.edge("t", "b", label: $q$),
  ),
))

We find
$
  scr(M)_2 = g_e^2/((p_1+p_2)^2 - m^2 c^2) [overline(u)(4) feyn(epsilon.alt)(3)^* (feyn(p)_1+feyn(p)_2+m c) feyn(epsilon.alt)(2) u(1)]
$

== Computing $abs(scr(M))^2$
Now, computing $abs(scr(M))^2$ requires knowing the spins. We usually do not know the incoming spins $s_i$ and outgoing spins $s_f$ exactly. We instead use the "spin-averaged" amplitude
$
  expval(abs(scr(M))^2) equiv 1/N_i sum_(s_i) sum_(s_f) abs(scr(M) (s_i -> s_f))^2,
$
where we only average over the initial spins $s_i$!#footnote[When our experimental setup is unaware of the initial spin state and our detector cannot measure outgoing spin.]

Now, consider $e^- mu^-$ scattering. We found
$
  abs(scr(M))^2 = g_e^4/(p_1-p_3)^4 [overline(u)(3) gamma^mu u (1)] [overline(u)(3) gamma^nu u(1)]^* [overline(u)(4) gamma_mu u(2)][overline(u)(4) gamma_nu u(2)]^*
$
Here, "all" terms are of the form
$
  G eq [overline(u)(a) Gamma_1 u(b)][overline(u)(a) Gamma_2 u(b)]^*
$
where $Gamma_1 = gamma^mu$ and $Gamma_2 = gamma^nu$ above. We have
$
  [overline(u)(a) Gamma_2 u(b)]^* & = [u(a)^dagger gamma^0 Gamma_2 u(b)]^dagger = u(b)^dagger Gamma_2^dagger gamma^0 u(a) = u(b)^dagger underbracket(gamma^0 gamma^0, bb(1)) Gamma_2^dagger gamma^0 u(a) = overline(u)(b) overline(Gamma)_2 u(a)
$
with $overline(Gamma)_2$ being
$
  overline(Gamma)_2 equiv gamma^0 Gamma_2^dagger gamma^0
$
implying
$
  G = [overline(u)(a) Gamma_1 u(b)][overline(u)(b) overline(Gamma)_2 u(a)]
$
Now, we sum over $s_b$
$
  sum_(s_b) G &= overline(u) (a) Gamma_1 [sum_(s_b) u^((s_b)) (p_b) overline(u)^((s_b)) (p_b)] overline(Gamma)_2 u(a) \
  &=^"completeness" overline(u)(a) Gamma_1 (feyn(p)_b + m_b c) overline(Gamma)_2 u(a) \
  &= overline(u)(a) Q u(a)
$
with $Q$ being
$
  Q equiv Gamma_1 (feyn(p)_b + m_b c) overline(Gamma)_2
$
Now, we sum over $s_a$
$
  sum_(s_a) sum_(s_b) G & = sum_(s_a) overline(u)^((s_a)) (p_a) Q u^((s_a)) (p_a) \
  &= sum_(s_a) sum_(alpha, beta =1)^4 overline(u)_alpha^((s_a)) (p_a) Q_(alpha beta) u_beta^((s_a)) (p_a) \
  &= sum_(alpha, beta=1)^4 Q_(alpha beta) [sum_(s_a) u_beta^((s_a)) (p_a) overline(u)_alpha^((s_a)) (p_a)] \
  &= sum_(alpha, beta = 1)^4 Q_(alpha beta) (feyn(p)_a + m_a c)_(beta alpha) \
  &= tr [Q (feyn(p)_a+m_a c)]
$
implying
$
  sum_(s_a, s_b) [overline(u)(a) Gamma_1 u(b)][overline(u)(a) Gamma_2 u(b)]^* = tr [Gamma_1 (feyn(p)_b + m_b c) overline(Gamma)_2 (feyn(p)_a + m_a c)]
$
which is called Casimir's trick.#footnote[Here, we had an expression of the form $dots overline(u)(a) [bullet] u (b) overline(u)(b) [bullet] u(a) dots$. We should think of this expression as being circular i.e. the $dots$ connect. We can immediately replace $u overline(u) -> (feyn(p) + m c)$ in this case or $v overline(v) -> (feyn(p) - m c)$ when we have "mixed" or "antiparticle" cases. ] When we have done the sums over $s_i$ and $s_f$ we need to remember the $N_i^(-1)$.

Now, with $e^- mu^-$ scattering we had $Gamma_2 = gamma^nu$ so
$
  overline(Gamma)_2 = gamma^0 gamma^(nu dagger) gamma^0 = gamma^nu
$
We apply Casimir's trick to find
$
  expval(abs(scr(M))^2) = underbracket(1/4, N_i^(-1)) g_e^4/((p_1-p_3)^4) underbracket(tr [gamma^mu (feyn(p)_1 + m_e c) gamma^nu (feyn(p)_3 + m_e c)], "from" sum_(s_1, s_3)) overbracket(tr [gamma_mu (feyn(p)_2 + m_mu c) gamma_nu (feyn(p)_4 + m_mu c)], "from" sum_(s_2,s_4))
$
We recall,
$
  tr (A + B) & = tr A + tr B",  " tr alpha A & = alpha tr A",  " tr A B & = tr B A
$
and,
$
  eta_(mu nu) eta^(mu nu) & = 4 \
     {gamma^mu, gamma^nu} & = 2 eta^(mu nu) => {feyn(a),feyn(b)} = 2 a dot b
$
With these one can derive the contraction theorems
$
  gamma_mu gamma^mu &= 4 \
  gamma_mu gamma^nu gamma^mu &= - 2 gamma^nu => gamma_mu feyn(a) gamma^mu = - 2 feyn(a) \
  gamma_mu gamma^nu gamma^lambda gamma^mu &= 4 eta^(nu lambda) => gamma_mu feyn(a) feyn(b) gamma^mu = 4 a dot b \
  gamma_mu gamma^nu gamma^lambda gamma^sigma gamma^mu &= -2 gamma^sigma gamma^lambda gamma^nu => gamma_mu feyn(a)feyn(b)feyn(c) gamma^mu = - 2 feyn(a) feyn(b) feyn(c)
$
and the trace theorems
$
  tr "odd" \# gamma^mu & = 0 \
  tr 1 & = 4 \
  tr gamma^mu gamma^nu &= 4 eta^(mu nu) => tr feyn(a) feyn(b) = 4 a dot b \
  tr gamma^mu gamma^nu gamma^lambda gamma^sigma &= 4 (eta^(mu nu) eta^(lambda sigma) - eta^(mu lambda) eta^(nu sigma) + eta^(mu sigma) eta^(nu lambda)) \
  => tr feyn(a)feyn(b)feyn(c)feyn(d) &= 4 [(a dot b )(c dot d) - (a dot c) (b dot d) + (a dot d) (b dot c)]
$
With $gamma^5 = i gamma^0 gamma^1 gamma^2 gamma^3$ we also have
$
  tr gamma^5 &= 0 \
  tr gamma^5 gamma^mu gamma^nu &= 0 => tr gamma^5 feyn(a) feyn(b) = 0 \
  tr gamma^5 gamma^mu gamma^nu gamma^lambda gamma^sigma &= 4 i epsilon.alt^(mu nu lambda sigma) => tr gamma^5 feyn(a)feyn(b)feyn(c)feyn(d) = 4 i epsilon.alt^(mu nu lambda sigma) a_mu b_nu c_lambda d_sigma
$
Now, consider one of the traces in $e^- mu$ scattering
$
  tr [gamma^mu (feyn(p)_1 + m_e c) gamma^nu (feyn(p)_3 + m_e c)] & = tr gamma^mu feyn(p)_1 gamma^nu feyn(p)_3 + m_e c underbracket([tr gamma^mu feyn(p)_1 gamma^nu + tr gamma^mu gamma^nu feyn(p)_3], 0 "by" \# gamma^mu "being odd") + (m_e c)^2 underbracket(tr gamma^mu gamma^nu, 4 eta^(mu nu))
$
We compute
$
  tr gamma^mu feyn(p)_1 gamma^nu feyn(p)_3 &= overbracket((p_1)_lambda (p_3)_sigma, "numbers") tr overbracket(gamma^mu gamma^lambda gamma^nu gamma^sigma, "matrices") \
  &= (p_1)_lambda (p_3)_sigma 4 (eta^(mu lambda) eta^(nu sigma) - eta^(mu nu) eta^(lambda sigma) + eta^(mu sigma) eta^(lambda nu)) \
  &= 4 [p_1^mu p_3^nu - eta^(mu nu) (p_1 dot p_3) + p_3^mu p_1^nu]
$
implying
$
  tr[gamma^mu (feyn(p)_1 + m_e c) gamma^nu (feyn(p)_3 + m_e c)] &= 4 {p_1^mu p_3^nu + p_3^mu p_1^nu + eta^(mu nu) [(m_e c)^2 - (p_1 dot p_3)]}
$
We find
$
  expval(abs(scr(M))^2) &= (4 g_e^4)/(p_1-p_3)^4 {p_1^mu p_3^nu + p_3^mu p_1^nu + eta^(mu nu) [(m_e c)^2 - (p_(1 kappa) p_3^kappa)]} {p_(2 mu) p_(4 nu) + p_(4 mu) p_(2 nu) + eta_(mu nu) [(m_mu c)^2 - (p_2 dot p_4)]} \
  &= (8 g_e^4)/(p_1-p_3)^4 [(p_1 dot p_2)(p_3 dot p_4) + (p_1 dot p_4) (p_2 dot p_3)- (p_1 dot p_3) (m_mu c)^2 - (p_2 dot p_4) (m_e c)^2 + 2 (m_e m_mu c^2)^2]
$

== $Gamma$ for $e^- mu$ scattering
Now, we can compute $Gamma$ and $sigma$ using $expval(abs(scr(M))^2)$.

We consider $e^- mu^-$ scattering assuming $m_mu >> m_e$ and would like $sigma$ in the lab frame with $m_mu$ at rest. We have
$
  dv(sigma, Omega) = (hbar/(8 pi M c))^2 expval(abs(scr(M))^2)
$
with
$
  p_1 = vecrow(E/c, bold(p)_1)";  " p_2 = vecrow(m_mu c, bold(0))";  " p_3 = vecrow(E/c, bold(p)_3)";  " p_4 = vecrow(m_mu c, bold(0))
$
and $abs(bold(p)_1) = abs(bold(p)_3) = abs(bold(p))$ with $bold(p)_1 dot bold(p)_3 = bold(p)^2 cos theta$. We have
$
                 (p_1-p_3)^2 & =- (bold(p)_1-bold(p)_3)^2
                               = - 4 bold(p)^2 sin^2 theta/2 \
               (p_1 dot p_3) & = m_e^2 c^2 + 2 bold(p)^2 sin^2 theta/2 \
  (p_1 dot p_2)(p_3 dot p_4) & = (p_1 dot p_4)(p_2 dot p_3) = (m_mu E)^2 \
               (p_2 dot p_4) & = (m_mu c)^2
$
implying
$
  expval(abs(scr(M))^2) = ((g_e^2 m_mu c)/(bold(p)^2 sin^2 theta\/2))^2 [(m_e c)^2 + bold(p)^2 cos^2 theta/2]
$
and
$
  dv(sigma, Omega) = ((alpha hbar)/(2 bold(p)^2 sin^2 theta\/2))^2 [(m_e c)^2 + bold(p)^2 cos^2 theta/2]
$
which is called the Mott formula. With $bold(p)^2 << (m_e c)^2$ we find the Rutherford formula
$
  dv(sigma, Omega) = (e^2/(2 m v^2 sin^2 theta\/2))^2
$
which is quite nice.

== Pair annihilation
Now, we consider pair annihilation of the form $e^+ + e^- -> 2 gamma$. We have two diagrams with#footnote[One is "twisted".]
$
  scr(M)_1 &= g_e^2/((p_1-p_3)^2 - m^2 c^2) overline(v)(2) feyn(epsilon.alt)_4 (feyn(p)_1-feyn(p)_3 + m c) feyn(epsilon.alt)_3 u(1) \
  scr(M)_2 &= g_e^2/((p_1-p_4)^2-m^2 c^2) overline(v)(2) feyn(epsilon.alt)_3 (feyn(p)_1-feyn(p)_4 + m c) feyn(epsilon.alt)_4 u(1)
$
and $scr(M) = scr(M)_1 + scr(M)_2$.

We will work in the CM-frame and assume the $e^(plus.minus)$ are initially $tilde$ at rest. Here, we will be unable to average over $s_i$ since the composite system $e^- e^+$ is in the singlet or triplet configuration. When $e^+$ and $e^-$ are at rest the photons $gamma$ come out back-to-back. We take the $z$-axis to coincide with the direction of the first photon. We have#footnote[We use $E_i = 2 m c^2 =^! 2 E_gamma = 2 abs(bold(p)) c$]
$
  p_1 = vec(m c, bold(0))",  " p_2 = vec(m c, bold(0))",  " p_3 = vec(m c, 0, 0, m c)",  " p_4 = vec(m c, 0, 0, -m c)
$
implying
$
  (p_1-p_3)^2 - m^2 c^2 = (p_1-p_4)^2 - m^2 c^2 = - 2 (m c)^2
$
Now, we use
$
  feyn(p)_1 feyn(epsilon.alt)_3 &= - feyn(epsilon.alt)_3 feyn(p)_1 + 2 underbracket(p_1 dot epsilon.alt_3, 0 "in Coulomb gauge") \
  feyn(p)_3 feyn(epsilon.alt)_3 &= -feyn(epsilon.alt)_3 feyn(p)_3 + 2 underbracket(p_3 dot epsilon.alt_3, 0 "by Lorenz gauge")
$
implying
$
  (feyn(p)_1 - feyn(p)_3 + m c) feyn(epsilon.alt)_3 u(1) &= feyn(epsilon.alt)_3 (-feyn(p)_1 + feyn(p)_3 + m c) u(1) \
  &=^"by Dirac" feyn(epsilon.alt)_3 feyn(p)_3 u(1) \
  (feyn(p)_1 - feyn(p)_4 + m c) feyn(epsilon.alt)_4 u(1) &= feyn(epsilon.alt)_4 feyn(p)_4 u(1)
$
We find
$
  scr(M) = - g_e^2/(2 (m c)^2) overline(v)(2) [feyn(epsilon.alt)_4 feyn(epsilon.alt)_3 feyn(p)_3 + feyn(epsilon.alt)_3 feyn(epsilon.alt)_4 feyn(p)_4] u(1)
$
with
$
  feyn(p)_3 = m c( gamma^0 - gamma^3)";  " feyn(p)_4 = m c (gamma^0 + gamma^3)
$
implying
$
  [feyn(epsilon.alt)_4 feyn(epsilon.alt)_3 feyn(p)_3 + feyn(epsilon.alt)_3 feyn(epsilon.alt)_4 feyn(p)_4] &= m c [(feyn(epsilon.alt)_4 feyn(epsilon.alt)_3 + feyn(epsilon.alt)_3 feyn(epsilon.alt)_4) gamma^0 - (feyn(epsilon.alt)_4 feyn(epsilon.alt)_3 - feyn(epsilon.alt)_3 feyn(epsilon.alt)_4) gamma^3]
$
Now, we use
$
  feyn(epsilon.alt) = - bold(epsilon.alt) dot bold(gamma) = - mat(0, bold(sigma) dot bold(epsilon.alt); - bold(sigma) dot bold(epsilon.alt), 0)
$
or
$
  feyn(epsilon.alt)_3 feyn(epsilon.alt)_4 = - mat((bold(sigma) dot bold(epsilon.alt)_3)(bold(sigma) dot bold(epsilon.alt)_4), 0; 0, (bold(sigma) dot bold(epsilon.alt)_3)(bold(sigma) dot bold(epsilon.alt)_4))
$
with
$
  (bold(sigma) dot bold(a))(bold(sigma) dot bold(b)) = bold(a) dot bold(b) + i bold(sigma) dot (bold(a) times bold(b))
$
implying
$
  (feyn(epsilon.alt)_4 feyn(epsilon.alt)_3 + feyn(epsilon.alt)_3 feyn(epsilon.alt)_4) &= - 2 bold(epsilon.alt)_3 dot bold(epsilon.alt)_4",  " (feyn(epsilon.alt)_4 feyn(epsilon.alt)_3 - feyn(epsilon.alt)_3 feyn(epsilon.alt)_4) &= 2 i (bold(epsilon.alt)_3 times bold(epsilon.alt)_4) mat(bold(sigma), 0; 0, bold(sigma))
$
We find
$
  scr(M) = g_e^2/(m c) overline(v)(2) [(bold(epsilon.alt)_3 dot bold(epsilon.alt)_4)gamma^0 + i (bold(epsilon.alt)_3 times bold(epsilon.alt)_4) mat(bold(sigma), 0; 0, bold(sigma)) gamma^3] u(1)
$
Now, we consider the singlet configuration of $e^- e^+$ or symbolically
$
  scr(M)_"singlet" = 1/sqrt(2) (scr(M)_(arrow.t arrow.b) - scr(M)_(arrow.b arrow.t))
$
where $scr(M)_(arrow.t arrow.b)$ is found by using $u^((1))$ (spin up $e^-$) and $v^((2))$ (spin down $e^+$) given by
$
  u(1) = sqrt(2 m c) vec(1, 0, 0, 0)",  " overline(v)(2) = sqrt(2 m c) vecrow(0, 0, 1, 0)
$
We find
$
  overline(v)(2) gamma^0 u(1) = 0",  " overline(v)(2) mat(bold(sigma), 0; 0, bold(sigma)) gamma^3 u(1) = - 2 m c hat(z)
$
implying
$
  scr(M)_(arrow.t arrow.b) = - 2 i g_e^2 (bold(epsilon.alt)_3 times bold(epsilon.alt)_4)_z
$
Likewise, we find
$
  scr(M)_(arrow.b arrow.t) = - scr(M)_(arrow.t arrow.b)
$
implying
$
  scr(M)_"singlet" = -2 sqrt(2) i g_e^2 (bold(epsilon.alt)_3 times bold(epsilon.alt)_4)_z
$
Note, the triplet configuration would have
$
  scr(M)_"triplet" = 1/sqrt(2) (scr(M)_(arrow.t arrow.b) + scr(M)_(arrow.b arrow.t)) = 0
$
implying the decay is forbidden in this case!

Now, we need the corresponding $bold(epsilon.alt)$. Recall the "spin-up" ($m_s = +1$) and "spin-down" ($m_s = -1$) polarisations
$
  bold(epsilon.alt)_+ = - 1/sqrt(2) vec(1, i, 0)";  " bold(epsilon.alt)_- = 1/sqrt(2) vec(1, -i, 0)
$
We know the value of $L_z$ should vanish, so the photon spins should be antiparallel $arrow.t arrow.b$ or $arrow.b arrow.t$. When $arrow.t arrow.b$ we have#footnote[Here, the photon 4 has $h = +1$ along $-z$ i.e. $epsilon.alt_+ (+z) = epsilon.alt_- (+z)$]
$
  bold(epsilon.alt)_3 = bold(epsilon.alt)_+",  " bold(epsilon.alt)_4 = bold(epsilon.alt)_-
$
implying
$
  bold(epsilon.alt)_3 times bold(epsilon.alt)_4 = i hat(k)
$
When $arrow.b arrow.t$ we have
$
  bold(epsilon.alt)_3 times bold(epsilon.alt)_4 = -i hat(k)
$
We need the antisymmetric combination#footnote[Where $arrow.t$ and $arrow.b$ refer to photon polarisation.]
$
  scr(M)_"singlet" &= 1/sqrt(2) (scr(M)_(arrow.t arrow.b) - scr(M)_(arrow.b arrow.t)) = - 4 g_e^2
$
since $scr(M)_"triplet" = 0$.

Now, we have in the CM-frame
$
  dv(sigma, Omega) =1/2 ((hbar c)/(8 pi (E_1+E_2)))^2 abs(bold(p)_f)/abs(bold(p)_i) abs(scr(M))^2
$
with
$
  E_1 = E_2 = m c^2",  " abs(bold(p)_f) = m c",  " abs(bold(p)_i) = (m v)/2
$
where $v$ is the relative velocity between $e^-$ and $e^+$. We find
$
  dv(sigma, Omega) = 1/(c v) ((hbar alpha)/m)^2 => sigma = (4 pi)/(c v) ((hbar alpha)/m)^2
$
Now, we can compute $tau_(e^- e^+)$. We can write
$
  dv(sigma, Omega) = 1/scr(L) dv(N, Omega)
$
or $N = scr(L) sigma$ with
$
  scr(L) = rho v
$
where $rho$ is the number density of incident particles. The $e^-$ density of a single "atom" is $abs(psi(0))^2$ and $N = Gamma$ implying
$
  Gamma = v sigma abs(psi(0))^2 = (4 pi)/c ((hbar alpha)/m)^2 abs(psi(0))^2
$
Now, the ground state has
$
  abs(psi(0))^2 = 1/pi ((alpha m c)/(2 hbar))^2
$
implying
$
  tau = Gamma^(-1) = (2 hbar)/(alpha^5 m c^2) tilde 10^(-10) "s"
$
