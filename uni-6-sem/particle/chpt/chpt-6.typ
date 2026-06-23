#import "../../temp.typ": *
#show: chpt-note.with()

= Computing $Gamma$ and $sigma$
We can probe interactions through processes like decays and scattering. Now, we would like to be able to compute stuff related to these processes. We essentially need to bridge the gap between experiment and theory.

== $Gamma$
Now, consider a decay
$
  A -> B + C + dots
$
We can define the decay rate $Gamma$ by
$
  Gamma = "probability per unit time that particle" A "decays",
$
or
$
  dv(N, t) = - Gamma N,
$
implying
$
  N = N_0 e^(-t\/tau),
$
where $tau$ is the mean lifetime
$
  tau = 1/Gamma
$
Now, consider multiple decay channels for a particle $A$. We can then define a total decay rate $Gamma_"tot"$ by
$
  Gamma_"tot" = sum_i Gamma_i,
$
where $Gamma_i$ is the decay rate of channel $i$.#footnote[Now, $tau = 1\/Gamma_"tot"$.] We also define the branching ratio by
$
  "branching ratio for the" i^"th" "channel" = Gamma_i/Gamma_"tot"
$
Then computing $tau$ and the branching ratios requires the $Gamma_i$.

== $sigma$
Now, consider a scattering process. We are interested in the cross section $sigma$.#footnote[Think of $sigma$ as the size of our target. We could again define $sigma_"tot" = sum_i sigma_i$.]

#figure(
  diagram(
    node((0, 0), radius: .2em, fill: black),
    node((0, -1), radius: .2em, fill: gray),
    node((-2, -1), radius: .2em, fill: black),
    edge((-2, -1), (-1, -1), "-|>", label: [in]),
    edge((-2, 0), (2, 0), "--", label: [target], label-side: right),
    edge((-2, 0), (-2, -1), "~", label: $b$, label-side: left),
    edge((-1, -1), (0, -1), "--"),
    edge((0, -1), (2, -2), "-|>"),
    edge(
      (0, -1),
      (2, -1),
      "--",
      label: $theta$,
      label-sep: 0em,
      label-pos: 40%,
    ),
  ),
  caption: [Example scattering.],
)<scattering>

We consider a particle scattering off a target with scattering angle $theta(b)$, as shown in @scattering. #footnote[Classical hard-sphere scattering gives
  $ b =^"hard-sphere" R cos theta/2, $
  and Rutherford scattering gives
  $ b=^"Rutherford" (q_1 q_2)/(2 E) cot theta/2. $] Now, a particle with impact parameter between $b$ and $b + dd(b)$ will scatter into an angle between $theta$ and $theta + dd(theta)$. More generally, a particle passing through $dd(sigma)$ will scatter into a solid angle $dd(Omega)$.

We relate these by
$
  dd(sigma) = D(theta) dd(Omega),
$
where $D$ being the differential cross section.

Now, we would like an expression for $D(theta)$. We can write
$
  dd(sigma) &= abs(underbracket(dd(phi.alt)/(2 pi) 2 pi (b+dd(b))^2, "outer section") - underbracket(dd(phi.alt)/(2 pi) pi b^2, "inner section")) \
  &= abs(b dd(b, phi.alt)),
$
where we consider the circle of radius $b$ expanding by $dd(b)$ and finding the change in $dd(sigma)$. Also, we consider small angles $dd(phi.alt)$ since these lie in the same plane as $b$ whereas $theta$ is "perpendicular".

Then#footnote[We have
  $ D & =^"hard-sphere" R^2/4 =^"Rutherford" ((q_1 q_2)/(4 E sin^2 (theta\/2)))^2 $]
$
  D(theta) = dv(sigma, Omega) = abs((b dd(b, phi.alt))/(sin theta dd(theta, phi.alt))) = abs(b/(sin theta) dv(b, theta))
$

The cross section is given by#footnote[We find $ sigma & =^"hard-sphere" pi R^2 =^"Rutherford" oo $
]
$
  sigma & = integral dd(sigma),
$
all particles within $sigma$ will scatter.

/*
Consider a beam of particles with uniform luminosity $L$. Then $dd(N) = L dd(sigma)$ is the number of particles per unit time passing through the area $dd(sigma)$ implying
$
  dd(N) = L dv(sigma, Omega) dd(Omega)
$
or if $L$ and $dd(Omega)$ are fixed then
$
  dv(sigma, Omega) = dd(N)/(L dd(Omega))
$
*/

== Fermi's Golden Rule
Now, consider a particle decaying
$
  1 -> 2 + 3+ dots + n.
$
Then $Gamma$ is given by Fermi's Golden Rule
$
  Gamma &= S/(2 hbar m_1) integral abs(scr(M))^2 (2 pi)^4 underbracket(delta^((4)) (p_1-p_2 - dots - p_n), "four-momentum conservation") \
  & times product_(j=2)^n 2 pi overbracket(delta(p_j^2-m_j^2 c^2), "all outgoing particles are real") underbracket(Theta(p_j^0), "positive energies") dd(p_j, 4)/(2 pi)^4,
$
where $S$ is the symmetry factor, for each group of $s$ identical particles $S$ gets a factor of $(1\/s!)$. The dynamical information of the process are contained within the amplitude $scr(M) (p_1, dots, p_n)$ everything else is referred to as phase space (the kinematical information).#footnote[The units of $scr(M)$ are $(m c)^(4-n)$ where $n = "#external lines"$.] The $delta$ and $Theta$ functions ensure the kinematics make sense. Note, all $delta$ carry $(2 pi)$ and all $dd(p_j)$ carry $(2 pi)^(-1)$.

We will write $dd(p_j, 4) = dd(p_j^0, bold(p)_j, [1,3])$ and
$
  delta(p_j^2 -m_j^2 c^2) &= delta((p_j^0)^2 - (bold(p)_j^2 + m_j^2 c^2)) \
  &= 1/(2 sqrt(bold(p)_j^2 + m_j^2 c^2)) [delta(p_j^0 - sqrt(bold(p)_j^2 + m_j^2 c^2)) + delta(p_j^0 + sqrt(bold(p)_j^2 + m_j^2 c^2))].
$
The $Theta$ function kills the spike at $p_j^0 = -sqrt(bold(p)_j^2 + m_j^2 c^2)$ implying
$
  Theta(p_j^0) delta(p_j^2 - m_j^2 c^2) = 1/(2 sqrt(bold(p)_j^2 + m_j^2 c^2)) delta(p_j^0 - sqrt(bold(p)_j^2 + m_j^2 c^2)).
$
Now, we can compute the integral over $p_j^0$ using the $delta$ function to obtain#footnote[Since $Theta$ is unity at $p_j^0 = sqrt(bold(p)_j^2 +m_j^2 c^2)$.]
$
  Gamma = S/(2 hbar m_1) integral abs(scr(M))^2 (2 pi)^4 delta^((4)) (p_1 - dots - p_n) product_(j=2)^n 1/(2 sqrt(bold(p)_j^2 + m_j^2 c^2)) dd(bold(p)_j, 3)/(2 pi)^3,
$
with
$
  p_j^0 -> sqrt(bold(p)_j^2 + m_j^2 c^2),
$
everywhere.

== Two-particle decay
Now, consider the decay
$
  1 -> 2 + 3
$
Then
$
  Gamma = S/(32 pi^2 hbar m_1) integral abs(scr(M))^2 (delta^((4)) (p_1 - p_2 - p_3))/(sqrt(bold(p)_2^2 + m_2^2 c^2) sqrt(bold(p)_3^2 + m_3^2 c^2)) dd(bold(p)_2, bold(p)_3, [3,3]).
$
We write
$
  delta^((4)) (p_1 - p_2 - p_3) = delta(p_1^0-p_2^0-p_3^0) delta^((3)) (bold(p)_1 - bold(p)_2 - bold(p)_3),
$
and assume $bold(p)_1 = 0$ and $p_1^0 = m_1 c$ giving
$
  Gamma &= S/(32 pi^2 hbar m_1) integral abs(scr(M))^2 (delta(m_1 c - sqrt(bold(p)_2^2 + m_2^2 c^2) - sqrt(bold(p)_3^2 +m_3^2 c^2)))/(sqrt(bold(p)_2^2 + m_2^2 c^2) sqrt(bold(p)_3^2 + m_3^2 c^2)) delta^((3)) (bold(p)_2 + bold(p)_3) dd(bold(p)_2, bold(p)_3, [3,3]) \
  &= S/(32 pi^2 hbar m_1) integral abs(scr(M))^2 (delta(m_1 c - sqrt(bold(p)_2^2 + m_2^2 c^2)- sqrt(bold(p)_2^2 + m_3^2 c^2)))/(sqrt(bold(p)_2^2 + m_2^2 c^2) sqrt(bold(p)_2^2 + m_3^2 c^2)) dd(bold(p)_2, 3).
$
We use spherical coordinates with $dd(bold(p)_2, 3) = r^2 sin theta dd(r, theta, phi.alt)$ and $r = abs(bold(p)_2)$
$
  Gamma &= S/(32 pi^2 hbar m_1) integral abs(scr(M))^2 (delta(m_1 c - sqrt(r^2 + m_2^2 c^2) - sqrt(r^2 + m_3^2 c^2)))/(sqrt(r^2+m_2^2 c^2) sqrt(r^2 + m_3^2 c^2)) r^2 sin theta dd(r, theta, phi.alt)
$
We had $scr(M) = scr(M)(p_1, p_2, p_3)$ initially and after everything we now have $scr(M) = scr(M)(bold(p)_2)$. However, $scr(M)$ is a scalar implying $scr(M) = scr(M) (bold(p)_2^2 = r^2)$. We can then do the $theta$ and $phi.alt$ integrals
$
  Gamma = S/(8 pi hbar m_1) integral_0^oo abs(scr(M)(r))^2 (delta(m_1 c - sqrt(r^2 +m_2^2 c^2)-sqrt(r^2+m_3^2 c^2)))/(sqrt(r^2+m_2^2 c^2) sqrt(r^2+m_3^2 c^2)) r^2 dd(r)
$
Now, we define
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
Now, consider the process
$
  1+2 -> 3+4+dots+n
$
Then $sigma$ is given by
$
  sigma &= (S hbar^2)/(4 sqrt((p_1 dot p_2)^2 - (m_1 m_2 c^2)^2)) integral abs(scr(M))^2 (2 pi)^4 delta^((4)) (p_1+p_2-p_3 - dots - p_n) \
  & times product_(j=3)^n 2 pi delta(p_j^2 - m_j^2 c^2) Theta(p_j^0) dd(p_j, 4)/(2 pi)^4
$
We can do the $p_j^0$ integrals as before
$
  sigma &= (S hbar^2)/(4 sqrt((p_1 dot p_2)^2 - (m_1 m_2 c^2)^2)) integral abs(scr(M))^2 (2 pi)^4 delta^((4)) (p_1 + p_2 - p_3 - dots - p_n) \
  & times product_(j=3)^n 1/(2 sqrt(bold(p)_j^2 + m_j^2 c^2)) dd(bold(p)_j, 3)/(2 pi)^3
$
with
$
  p_j^0 = sqrt(bold(p)_j^2 + m_j^2 c^2)
$
everywhere.

== Two-particle scattering
Now, consider the process
$
  1 + 2 -> 3 + 4
$
We will use the CM-frame where $bold(p)_2 = - bold(p)_1$ and
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
Here, we cannot do the $theta$ and $phi.alt$ integrals since $scr(M)$ depends on the direction of $bold(p)_3$ as well as its magnitude.#footnote[We have three invariants $bold(p)_1^2$, $bold(p)_3^2$ and $bold(p)_1 dot bold(p)_3 = abs(bold(p)_1)abs(bold(p)_3) cos theta$. However, $bold(p)_1$ is fixed so the only integration variables $scr(M)$ can depend on are $abs(bold(p)_3)$ and $theta$.] We instead write
$
  dd(bold(p)_3, 3) = r^2 dd(r, Omega)
$
with $r = abs(bold(p)_3)$ implying
$
  dv(sigma, Omega) &= (hbar/( 8 pi))^2 (S c)/((E_1 + E_2) abs(bold(p)_1)) underbracket(integral_0^oo abs(scr(M))^2 (delta((E_1+E_2)\/c - sqrt(r^2 + m_3^2 c^2) - sqrt(r^2+m_4^2 c^2)))/(sqrt(r^2+m_3^2 c^2) sqrt(r^2+m_4^2 c^2)) r^2 dd(r), "as before")
$
We do the integral as before to obtain
$
  dv(sigma, Omega) &= ((hbar c)/(8 pi))^2 (S abs(scr(M))^2)/(E_1+E_2)^2 abs(bold(p)_f)/abs(bold(p)_i)
$
where $abs(bold(p)_f)$ and $abs(bold(p)_i)$ are the magnitudes of either- outgoing and incoming momentum.

== ABC theory
Now, we want to compute ampltitudes $scr(M)$. This is done by writing all the relevant Feynman diagrams for a process and then using the Feynman rules. We will consider a toy model which will serve to illustrate the method.

We imagine a world with three particles: $(A,B,C)$.#footnote[We will assume $m_A > m_B + m_C$.] We assume they are spin zero and each is its own antiparticle. There is one fundamental vertex which is shown in @ABSvert. This diagram contributes to $A -> B + C$.
#figure(
  feynman.feynman(
    (
      feynman.vertex("A"),
      feynman.vertex("B"),
      feynman.vertex("C"),
      feynman.vertex("M"),
      feynman.edge("A", "M", label: $A$, type: "line"),
      feynman.edge("M", "B", label: $B$, type: "line"),
      feynman.edge("M", "C", label: $C$, type: "line"),
    ),
    orientation: "vertical",
  ),
  caption: [The fundamental vertex in ABC theory.],
)<ABSvert>
We will also consider other processes such as $A + A -> B + B$ and $A + B -> A + B$.

Now, consider a general Feynman diagram as shown in @feynmandiagramABC.
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
    feynman.edge("i2", "node", label: $p_2$),
    feynman.edge("i3", "node", label: $p_3$),
    feynman.edge("node", "f4", label: $p_4$),
    feynman.edge("node", "f5", label: $p_5$),
    feynman.edge("node", "f6", label: $p_6$),
  )),
  caption: [A general Feynman diagram.],
)<feynmandiagramABC>
The Feynman rules are:

1. Label the incoming- and outgoing momenta $p_1, dots, p_n$ as above and label internal momenta $q_1, q_2, dots$.

2. Each vertex carries a factor $-i g$. Where $g$ is referred to as the coupling.

3. Each internal line carries a propagator#footnote[We have $q_j^2 eq.not m_j^2 c^2$ since virtual particles are off-shell.] $ i/(q_j^2 - m_j^2 c^2) $

4. Each vertex carries a $delta$ function $ (2 pi)^4 delta^((4)) (k_1+k_2+k_3) $ with $k_i$ being negative for outward momenta.

5. Each internal line carries an integral $ integral dd(q_j, 4)/(2 pi)^4 $

6. The result will carry a factor $ (2 pi)^4 delta^((4)) (p_1+p_2+ dots - p_n) $ replace this by $i$.

== ABC examples
Now, consider the decay
$
  A -> B + C
$
The Feynman diagram for this process (see @ABSvert) has no internal lines so we immediately obtain
$
  scr(M) tilde (-i g) (2 pi)^4 delta^((4)) (p_1-p_2-p_3)
$
which upon $(2pi)^4 delta^((4)) -> i$ becomes
$
  scr(M) = g + cal(O)(g^2)
$
We can then compute $Gamma$
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

Now, consider the process
$
  A + A -> B + B
$
The Feynman diagram is shown in @ABCscat.
#figure(
  feynman.feynman(
    orientation: "vertical",
    (
      feynman.vertex("A1", label: $A$),
      feynman.vertex("A2", label: $A$),
      feynman.vertex("M1"),
      feynman.vertex("M2"),
      feynman.vertex("B1", label: $B$),
      feynman.vertex("B2", label: $B$),
      feynman.edge("A1", "M1", label: $p_2$, type: "fermion"),
      feynman.edge("M1", "B1", label: $p_4$, type: "fermion"),
      feynman.edge("A2", "M2", label: $p_1$, type: "fermion"),
      feynman.edge("M2", "B2", label: $p_3$, type: "fermion"),
      feynman.edge("M1", "M2", label: $C$, type: "antifermion"),
      feynman.edge("M2", "M1", type: "line", label: $q$),
    ),
  ),
  caption: [Lowest order diagram for $A + A -> B + B$.],
)<ABCscat>
With the Feynman rules we obtain
$
  scr(M)& tilde (-i g)^2 integral dd(q, 4)/(2 pi)^4 i/(q^2 - m_C^2 c^2) (2 pi)^4 delta^((4))(p_1 - q- p_3) (2 pi)^4 delta^((4))(p_2 + q- p_4) \
  &tilde -i g^2 (2 pi)^4 (delta^((4)) (p_1 + p_2 - p_3 - p_4))/((p_4 - p_2)^2 - m_C^2 c^2)
$
which upon $(2 pi)^4 delta^((4)) -> i$ becomes
$
  scr(M) = g^2/((p_4-p_2)^2 - m_C^2 c^2)
$
We can include a second diagram of order $cal(O)(g^2)$ found by twisting the $B$ lines. This amounts to $p_3 <-> p_4$ implying
$
  scr(M)_"tot" = g^2/((p_4-p_2)^2 - m_C^2 c^2) + g^2/((p_3-p_2)^2 - m_C^2 c^2)
$
Now, we can compute $D(theta)$ in the CM-frame. We assume $m_A = m_B = m$ and $m_C = 0$ implying

$
  (p_4-p_2)^2 - m_C^2 c^2 &= p_4^2 + p_2^2 - 2 p_2 dot p_4 = - 2 abs(bold(p))^2 (1 - cos theta) \
  (p_3-p_2)^2 - m_C^2 c^2 &= p_3^2 + p_2^2 - 2 p_3 dot p_2 = - 2 abs(bold(p))^2 (1+ cos theta)
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
