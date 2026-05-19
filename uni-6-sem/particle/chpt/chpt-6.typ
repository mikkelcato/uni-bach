#import "chpt-temp.typ": *
#show: chpt-note.with()

= Computing $Gamma$ and $sigma$
Within particle physics we can probe interactions through processes like decays and scattering. We would therefore like to be able to compute various physical quantities such as the lifetime of a particle.

== $Gamma$ and $sigma$
We define the _decay rate_ $Gamma$ by#footnote[This is simply the probability per unit time that a particle will decay.]
$
  dv(N, t) = - Gamma N
$
implying
$
  N = N_0 e^(-Gamma t)
$
The mean lifetime $tau$ is then
$
  tau = 1/Gamma
$
We typically have multiple decay modes for a given particle. With multiple decays we define the total decay rate $Gamma_"tot"$ by
$
  Gamma_"tot" = sum_i Gamma_i
$
and the lifetime is
$
  tau = 1/Gamma_"tot"
$
We can also define the _branching ratio_ by
$
  "branching ratio for" i^"th" "mode" = Gamma_i/Gamma_"tot"
$
Which is a measure of how many particles decay through each mode. Then to compute both $tau$ and the branching ratios we need to calculate $Gamma_i$.

When discussing scattering we are interested in the _cross section_ $sigma$. All processes have their own cross section $sigma_i$ and we define the total cross section $sigma_"tot"$ by
$
  sigma_"tot" = sum_i sigma_i
$
We now consider a particle encountering some potential and scattering off at an angle $theta$. This _scattering angle_ will be a function of the _impact parameter_ $b$. Classical _hard-sphere scattering_ gives
$
  b =^"hard-sphere" R cos theta/2
$
and _Rutherford scattering_ gives
$
  b=^"Rutherford" (q_1 q_2)/(2 E) cot theta/2
$
We imagine a particle coming in with an impact parameter between $b$ and $b + dd(b)$. This particle with then scatter with a scattering angle between $theta$ and $theta + dd(theta)$. We can more generally consider a particle passing through an infinitesimal area $dd(sigma)$ and scattering into a corresponding solid angle $dd(Omega)$. They are related as
$
  dd(sigma) = D(theta) dd(Omega)
$
with $D$ being the _differential cross section_. We have
$
  D & = dv(sigma, Omega) = abs(b/(sin theta) dv(b, theta)) \
    & =^"hard-sphere" R^2/4 \
    & =^"Rutherford" ((q_1 q_2)/(4 E sin^2 (theta\/2)))^2
$
The cross section can then be found as
$
  sigma & = integral dd(sigma) \
        & =^"hard-sphere" pi R^2 tilde "circle" \
        & =^"Rutherford" oo tilde "infinite range"
$
Any particles within this area will scatter while any outside will pass by unaffected.

Consider a beam of particles with uniform _luminosity_ $L$. Then $dd(N) = L dd(sigma)$ is the number of particles per unit time passing through the area $dd(sigma)$ implying
$
  dd(N) = L dv(sigma, Omega) dd(Omega)
$
or if $L$ and $dd(Omega)$ are fixed then
$
  dv(sigma, Omega) = dd(N)/(L dd(Omega))
$

== Fermi's _Golden Rule_
Consider a particle decaying as
$
  1 -> 2 + 3+ dots + n
$
Then the decay rate is given by
$
  Gamma &= S/(2 hbar m_1) integral abs(scr(M))^2 (2 pi)^4 delta^((4)) (p_1-p_2 - dots - p_n) \
  & times product_(j=2)^n 2 pi delta(p_j^2-m_j^2 c^2) Theta(p_j^0) dd(p_j, 4)/(2 pi)^4
$
where $S$ is the _symmetry factor_, for each group of $s$ identical particles $S$ gets a factor of $(1\/s!)$. The dynamics of the process are contained within the _amplitude_ $scr(M) (p_1, dots, p_n)$ everything else is referred to as _phase space_.#footnote[The units of $scr(M)$ are $(m c)^(4-n)$ where $n = "#external lines"$.] We are told to integrate over all outgoing four-momenta with the constraints:

1. All outgoing particles lie on their mass shell $p_j^2 = m_j^2 c^2$.

2. All outgoing energies are positive $p_j^0 > 0$.

3. Energy and momentum are conserved $p_1 = p_2 + dots + p_n$.

These are enforced by the $delta$ and $Theta$ functions.

We can write $dd(p_j, 4) = dd(p_j^0, bold(p)_j, [1,3])$ and
$
  delta(p_j^2 -m_j^2 c^2) &= delta((p_j^0)^2 - (bold(p)_j^2 + m_j^2 c^2)) \
  &= 1/(2 sqrt(bold(p)_j^2 + m_j^2 c^2)) [delta(p_j^0 - sqrt(bold(p)_j^2 + m_j^2 c^2)) + delta(p_j^0 + sqrt(bold(p)_j^2 + m_j^2 c^2))]
$
The $Theta$ function kills the spike at $p_j^0 = -sqrt(bold(p)_j^2 + m_j^2 c^2)$ implying
$
  Theta(p_j^0) delta(p_j^2 - m_j^2 c^2) = 1/(2 sqrt(bold(p)_j^2 + m_j^2 c^2)) delta(p_j^0 - sqrt(bold(p)_j^2 + m_j^2 c^2))
$
We can now perform the integral over $p_j^0$ using the $delta$ function and obtain#footnote[Since $Theta$ is unity at $p_j^0 = sqrt(bold(p)_j^2 +m_j^2 c^2)$]
$
  Gamma = S/(2 hbar m_1) integral abs(scr(M))^2 (2 pi)^4 delta^((4)) (p_1 - dots - p_n) product_(j=2)^n 1/(2 sqrt(bold(p)_j^2 + m_j^2 c^2)) dd(bold(p)_j, 3)/(2 pi)^3
$
with
$
  p_j^0 -> sqrt(bold(p)_j^2 + m_j^2 c^2)
$
everywhere.

As an example consider the decay
$
  1 -> 2 + 3
$
Then
$
  Gamma = S/(32 pi^2 hbar m_1) integral abs(scr(M))^2 (delta^((4)) (p_1 - p_2 - p_3))/(sqrt(bold(p)_2^2 + m_2^2 c^2) sqrt(bold(p)_3^2 + m_3^2 c^2)) dd(bold(p)_2, bold(p)_3, [3,3])
$
Also,
$
  delta^((4)) (p_1 - p_2 - p_3) = delta(p_1^0-p_2^0-p_3^0) delta^((3)) (bold(p)_1 - bold(p)_2 - bold(p)_3)
$
We assume $bold(p)_1 = 0$ and $p_1^0 = m_1 c$
$
  Gamma &= S/(32 pi^2 hbar m_1) integral abs(scr(M))^2 (delta(m_1 c - sqrt(bold(p)_2^2 + m_2^2 c^2) - sqrt(bold(p)_3^2 +m_3^2 c^2)))/(sqrt(bold(p)_2^2 + m_2^2 c^2) sqrt(bold(p)_3^2 + m_3^2 c^2)) delta^((3)) (bold(p)_2 + bold(p)_3) dd(bold(p)_2, bold(p)_3, [3,3]) \
  &= S/(32 pi^2 hbar m_1) integral abs(scr(M))^2 (delta(m_1 c - sqrt(bold(p)_2^2 + m_2^2 c^2)- sqrt(bold(p)_2^2 + m_3^2 c^2)))/(sqrt(bold(p)_2^2 + m_2^2 c^2) sqrt(bold(p)_2^2 + m_3^2 c^2)) dd(bold(p)_2, 3)
$
We adopt spherical coordinates with $dd(bold(p)_2, 3) = r^2 sin theta dd(r, theta, phi)$#footnote[This is momentum space $r = abs(bold(p)_2)$.]
$
  Gamma &= S/(32 pi^2 hbar m_1) integral abs(scr(M))^2 (delta(m_1 c - sqrt(r^2 + m_2^2 c^2) - sqrt(r^2 + m_3^2 c^2)))/(sqrt(r^2+m_2^2 c^2) sqrt(r^2 + m_3^2 c^2)) r^2 sin theta dd(r, theta, phi)
$
We had $scr(M) = scr(M)(p_1, p_2, p_3)$ initially and after everything we now have $scr(M) = scr(M)(bold(p)_2)$. However, $scr(M)$ is a scalar implying $scr(M) = scr(M) (bold(p)_2^2 = r^2)$. Then
$
  Gamma = S/(8 pi hbar m_1) integral_0^oo abs(scr(M)(r))^2 (delta(m_1 c - sqrt(r^2 +m_2^2 c^2)-sqrt(r^2+m_3^2 c^2)))/(sqrt(r^2+m_2^2 c^2) sqrt(r^2+m_3^2 c^2)) r^2 dd(r)
$
We define
$
  u = sqrt(r^2 + m_2^2 c^2) + sqrt(r^2 + m_3^2 c^2)
$
Then
$
  Gamma = S/(8 pi hbar m_1) integral_((m_2 +m_3) c)^oo abs(scr(M)(r))^2 delta(m_1 c - u) r/u dd(u)
$
Using the $delta$ function we send $u -> m_1 c$ and $r$ to
$
  r_0 = c/(2 m_1) sqrt(m_1^4+m_2^4+m_3^4 - 2 m_1^2 m_2^2 - 2 m_1^2 m_3^2 - 2 m_2^2 m_3^2)
$
which is the value of $abs(bold(p)_2)$ consistent with conservation of energy. Then#footnote[Being able to perform all the integrals as above without knowing the form of $scr(M)$ is typically not possible. Also, we assume $m_1 > m_2 + m_3$ otherwise $Gamma = 0$ as expected.]
$
  Gamma = (S abs(bold(p)))/(8 pi hbar m_1^2 c) abs(scr(M))^2
$
where $abs(bold(p))$ is the magnitude of either outgoing momentum.

== Scattering
Consider the process
$
  1+2 -> 3+4+dots+n
$
Then the cross section is given by
$
  sigma &= (S hbar^2)/(4 sqrt((p_1 dot p_2)^2 - (m_1 m_2 c^2)^2)) integral abs(scr(M))^2 (2 pi)^4 delta^((4)) (p_1+p_2-p_3 - dots - p_n) \
  & times product_(j=3)^n 2 pi delta(p_j^2 - m_j^2 c^2) Theta(p_j^0) dd(p_j, 4)/(2 pi)^4
$
We see everything is almost the same as before and we can again perform the $p_j^0$ integrals
$
  sigma &= (S hbar^2)/(4 sqrt((p_1 dot p_2)^2 - (m_1 m_2 c^2)^2)) integral abs(scr(M))^2 (2 pi)^4 delta^((4)) (p_1 + p_2 - p_3 - dots - p_n) \
  & times product_(j=3)^n 1/(2 sqrt(bold(p)_j^2 + m_j^2 c^2)) dd(bold(p)_j, 3)/(2 pi)^3
$
with
$
  p_j^0 = sqrt(bold(p)_j^2 + m_j^2 c^2)
$
everywhere.

As an example consider the process
$
  1 + 2 -> 3 + 4
$
in the CM frame $bold(p)_2 = - bold(p)_1$ and
$
  sqrt((p_1 dot p_2)^2 - (m_1 m_2 c^2)^2) = ((E_1 + E_2) abs(bold(p)_1))/c
$
Then
$
  sigma &= (S hbar^2 c)/(64 pi^2 (E_1 + E_2) abs(bold(p)_1)) integral abs(scr(M))^2 (delta^((4)) (p_1 +p_2 -p_3-p_4))/(sqrt(bold(p)_3^2 + m_3^2 c^2) sqrt(bold(p)_4^2 + m_4^2 c^2)) dd(bold(p)_3, bold(p)_4, [3,3])
$
We write
$
  delta^((4)) (p_1 + p_2 - p_3 - p_4) = delta((E_1+E_2)/c - p_3^0-p_4^0) delta^((3)) (bold(p)_3 + bold(p)_4)
$
and perform the $bold(p)_4$ integral using the $delta$ function
$
  sigma &= (hbar/(8 pi))^2 (S c)/((E_1 + E_2) abs(bold(p)_1)) integral abs(scr(M))^2 (delta((E_1 + E_2)\/c - sqrt(bold(p)_3^2 + m_3^2 c^2)- sqrt(bold(p)_3^2 + m_4^2 c^2)))/(sqrt(bold(p)_3^2 + m_3^2 c^2) sqrt(bold(p)_3^2 +m_4^2 c^2)) dd(bold(p)_3, 3)
$
This time we cannot do the angular integrals since $scr(M)$ depends on the direction of $bold(p)_3$ as well as its magnitude.#footnote[We have three invariants $bold(p)_1^2$, $bold(p)_3^2$ and $bold(p)_1 dot bold(p)_3 = abs(bold(p)_1)abs(bold(p)_3) cos theta$. However, $bold(p)_1$ is fixed so the only integration variables $scr(M)$ can depend on are $abs(bold(p)_3)$ and $theta$.] We still write
$
  dd(bold(p)_3, 3) = r^2 dd(r, Omega)
$
Then
$
  dv(sigma, Omega) &= (hbar/( 8 pi))^2 (S c)/((E_1 + E_2) abs(bold(p)_1)) underbracket(integral_0^oo abs(scr(M))^2 (delta((E_1+E_2)\/c - sqrt(r^2 + m_3^2 c^2) - sqrt(r^2+m_4^2 c^2)))/(sqrt(r^2+m_3^2 c^2) sqrt(r^2+m_4^2 c^2)) r^2 dd(r), "as before")
$
We do the integral as before to obtain
$
  dv(sigma, Omega) &= ((hbar c)/(8 pi))^2 (S abs(scr(M))^2)/(E_1+E_2)^2 abs(bold(p)_f)/abs(bold(p)_i)
$
where $abs(bold(p)_f)$ and $abs(bold(p)_i)$ are the magnitudes of either- outgoing and incoming momentum respectively.

== ABC theory
We now want to begin computing ampltitudes $scr(M)$. This is done by writing all the relevant Feynman diagrams for a given process and then using the _Feynman rules_. We start by considering a simple _toy_ theory which will serve to illustrate the method.

We imagine a world with three particles ${A,B,C}$.#footnote[We will assume $m_A > m_B + m_C$.] We assume they are spin zero and each is its own antiparticle. There is one primitive vertex
#figure(
  feynman(
    (
      vertex("A"),
      vertex("B"),
      vertex("C"),
      vertex("M"),
      edge("A", "M", label: $A$, type: "line"),
      edge("M", "B", label: $B$, type: "line"),
      edge("M", "C", label: $C$, type: "line"),
    ),
    orientation: "vertical",
  ),
)
This diagram contributes to $A -> B + C$. We will also consider other processes such as $A + A -> B + B$ and $A + B -> A + B$ both of which happen through the exchange of a $C$.

We now consider a general Feynman diagram
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
    edge("i2", "node", label: $p_2$),
    edge("i3", "node", label: $p_3$),
    edge("node", "f4", label: $p_4$),
    edge("node", "f5", label: $p_5$),
    edge("node", "f6", label: $p_6$),
  )),
)
Then the Feynman rules are:

1. Label the incoming- and outgoing momenta $p_1, dots, p_n$ as above and label internal momenta $q_1, dots$.

2. Each vertex carries a factor $-i g$. Where $g$ is referred to as the _coupling_.

3. Each internal line carries a _propagator_#footnote[We have $q_j^2 eq.not m_j^2 c^2$ since _virtual particles_ are off shell.] $ i/(q_j^2 - m_j^2 c^2) $

4. Each vertex carries a $delta$ function $ (2 pi)^4 delta^((4)) (k_1+k_2+k_3) $ with $k_i$ being negative for outward momenta.

5. Each internal line carries an integral $ integral dd(q_j, 4)/(2 pi)^4 $

6. The result will carry a factor $ (2 pi)^4 delta^((4)) (p_1+p_2+ dots - p_n) $ replace this by $i$.

As an example we consider the decay
$
  A -> B + C
$
This diagram has no internal lines and we obtain
$
  (-i g) (2 pi)^4 delta^((4)) (p_1-p_2-p_3)
$
which upon $(2pi)^4 delta^((4)) -> i$ becomes
$
  scr(M) = g + cal(O)(g^2)
$
Then
$
  Gamma = (g^2 abs(bold(p)))/(8 pi hbar m_A^2 c)
$
where
$
  abs(bold(p)) = c/(2 m_A) sqrt(m_A^4 + m_B^4 +m_C^4 - 2 m_A^2 m_B^2 - 2 m_A^2 m_C^2 - 2m_B^2 m_C^2)
$
We have found the lifetime of $A$
$
  tau_A = (8 pi hbar m_A^2 c)/(g^2 abs(bold(p)))
$

As a second example we consider the process
$
  A + A -> B + B
$
The lowest-order contribution is
#figure(
  feynman(
    orientation: "vertical",
    (
      vertex("A1", label: $A$),
      vertex("A2", label: $A$),
      vertex("M1"),
      vertex("M2"),
      vertex("B1", label: $B$),
      vertex("B2", label: $B$),
      edge("A1", "M1", label: $p_2$, type: "fermion"),
      edge("M1", "B1", label: $p_4$, type: "fermion"),
      edge("A2", "M2", label: $p_1$, type: "fermion"),
      edge("M2", "B2", label: $p_3$, type: "fermion"),
      edge("M1", "M2", label: $C$, type: "antifermion"),
      edge("M2", "M1", type: "line", label: $q$),
    ),
  ),
)
We obtain
$
  &(-i g)^2 integral dd(q, 4)/(2 pi)^4 i/(q^2 - m_C^2 c^2) (2 pi)^4 delta^((4))(p_1 - q- p_3) (2 pi)^4 delta^((4))(p_2 + q- p_4) \
  &-i g^2 (2 pi)^4 (delta^((4)) (p_1 + p_2 - p_3 - p_4))/((p_4 - p_2)^2 - m_C^2 c^2)
$
which upon $(2 pi)^4 delta^((4)) -> i$ becomes
$
  scr(M) = g^2/((p_4-p_2)^2 - m_C^2 c^2)
$
There is a second diagram of order $cal(O)(g^2)$ obtained by _twisting_ the $B$ lines. This amounts to $p_3 <-> p_4$ meaning the total amplitude is
$
  cal(M) = g^2/((p_4-p_2)^2 - m_C^2 c^2) + g^2/((p_3-p_2)^2 - m_C^2 c^2)
$
We now find $dd(sigma)\/dd(Omega)$ in the CM frame. We assume $m_A = m_B = m$ and $m_C = 0$ then
$
  (p_4-p_2)^2 - m_C^2 c^2 &= p_4^2 + p_2^2 - 2 p_2 dot p_4 = - 2 p^2 (1 - cos theta) \
  (p_3-p_2)^2 - m_C^2 c^2 &= p_3^2 + p_2^2 - 2 p_3 dot p_2 = - 2 p^2 (1+ cos theta)
$
with $bold(p)$ being the incoming momentum of particle $1$. Then
$
  scr(M) = - g^2/(bold(p)^2 sin^2 theta)
$
and
$
  dv(sigma, Omega) = 1/2 ((hbar c g^2)/(16 pi E bold(p)^2 sin^2 theta))^2
$
since $S = 1\/2$.
