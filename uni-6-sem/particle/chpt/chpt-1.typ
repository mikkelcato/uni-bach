#import "../../temp.typ": *
#show: chpt-note.with()

= The Zoo
== The "classical" particles
Thomson discovered the electron $e^-$ in 1897 and gave birth to particle physics.#footnote[With $bold(F) = q bold(E) + q bold(v) times bold(B)$ we can determine sign of $q$, the magnitude of $v$, and the ratio $q\/m$.] Later in 1911, Rutherford discovered that atoms contain a massive positive nucleus.#footnote[By firing $alpha$-particles onto a gold foil.] As should be familiar the nucleus of hydrogen is the proton $p^+$. The last of the three "classical particles" was found in 1932 when Chadwick discovered the neutron $n$.

== The photon
Planck quantised electromagnetic radiation in 1900 by assuming it came in small packages or "quanta" of energy $hbar omega$.#footnote[With no physical motivation.] Later in 1905, Einstein explained the photoelectric effect by assuming the electromagnetic field was quantised. This idea was accepted in 1923 when Compton showed that light scattered from a particle with mass $m$ was shifted in wavelength
$
  lambda' = lambda+ underbracket(h/(m c), lambda_c) (1 - cos theta).
$
This formula can be derived by assuming light is a particle with mass $m_gamma = 0$ and energy $hbar omega$. We call this particle the photon $gamma$.

== Mesons
The model above has a problem, since it does not explain why the nucleus is held together. We need a new short-ranged force, which we call the strong force. Yukawa suggested in 1934 that this force occurs due to the exchange of a particle called the pion $pi$. The pion has $m_pi tilde 300m_e tilde 1/6 m_p$ and is rather heavy.

$
  "lepton" tilde & e^- \
   "meson" tilde & pi \
  "baryon" tilde & p^+", "n
$

Later in 1947, the pion $pi$ and muon#footnote["heavy electron"] $mu$ is found in cosmic rays by Powell. Now, we know there are actually three pions: $pi^0$ and $pi^plus.minus$.

== Antiparticles
Dirac predicted the existence of antiparticles in 1930. Later in 1931, the positron $e^+$ was found by Anderson.#footnote[We have also found $mu^+$, $p^-$, $overline(n)$ etc.] Some particles are their own antiparticles, with the photon being an example $overline(gamma) = gamma$. Now, consider the process
$
  A + B -> C + D,
$
which we assume is known to occur. Then by "crossing symmetry" the following processes are also allowed#footnote[e.g. $gamma + e^- -> gamma + e^-$ is "the same as" $e^-+e^+ -> 2 gamma$]
$
                          A & -> overline(B) + C + D \
            A + overline(C) & -> overline(B) + D \
  overline(C) + overline(D) & -> overline(A) + overline(B).
$
However, some could be forbidden kinematically. The "reverse" process
$
  C + D -> A + B
$
is also allowed.#footnote[This follows from the principle of "detailed balance".]

== Neutrinos
Consider $beta$-decay as in
$
  A -> B + e^-.
$
We require $Q_A = Q_B - 1$ due to charge conservation. Now, we know $n in A$ and $p^+ in B$. $beta$-decay in the above form is problematic since the electron has too little energy. Pauli resolved this in 1930 by proposing the existence of the neutrino#footnote[and antineutrino $overline(nu)$] $nu$ changing $beta$-decay as
$
  n -> p^+ + e^- + nu\/overline(nu),
$
where the neutrino is neutral and interacts weakly, thereby making detection hard. Likewise, the neutrino was proposed to be involved in other decays#footnote[When the energy of the muon or electron varies there are multiple neutrinos.]
$
  pi & -> mu + nu\/overline(nu) \
  mu & -> e^- + 2 nu\/overline(nu).
$

We determine whether a process involves neutrinos or antineutrinos by assigning particles a "lepton number"#footnote[We can use this as the defining difference between neutrinos and antineutrinos.]
$
  L(e^-, mu^-, nu) = +1",  " L(e^+,mu^+, overline(nu)) = -1",  " L("other") = 0,
$
and proposing lepton number conservation.#footnote[We also have baryon number conservation.] Then
$
     n & -> p^+ + e^- + overline(nu) \
  pi^- & -> mu^- + overline(nu) \
  pi^+ & -> mu^+ + nu \
  mu^- & -> e^- + nu + overline(nu) \
  mu^+ & -> e^+ + nu + overline(nu).
$

However, the process
$
  mu^- -> e^- + gamma
$
never occurs even though no rules seem to broken. This implies lepton number conservation of within each generation and the existence of multiple neutrinos. We call these the electron-neutrino $nu_e$ and muon-neutrino $nu_mu$. The lepton number is assigned as before, but now we have an electron number $L_e$ and muon number $L_mu$. Then
$
     n & -> p^+ + e^- + overline(nu)_e \
  pi^+ & -> mu^+ + nu_mu \
  mu^- & -> e^- + overline(nu)_e + nu_mu
$

The neutrino was found in 1955 using inverse $beta$-decay
$
  overline(nu)_e + p^+ -> n + e^+
$
Later in 1962, the existence of two neutrino types was shown using processes like#footnote[With $overline(nu)_mu$ from $pi^-$-decay.]
$
  overline(nu)_mu + p^+ & -> n + mu^+ \
  overline(nu)_mu + p^+ & -> n + e^+,
$
with the second never occuring.

== Strangeness
The kaon $K^0$ was found in 1947 decaying as
$
  K^0 -> pi^+ + pi^-.
$
Later in 1949, another kaon $K^+$ was found decaying as
$
  K^+ -> 2 pi^+ + pi^-.
$
The $K^+$ being a charged $K^0$ was realised in 1956.#footnote[We also have $overline(K)^0, K^-$ etc.] The kaons are mesons like pions.#footnote[We have since found many more mesons: $rho^0, rho^+, rho^-, eta, omega, phi.alt$ etc.]

The $Lambda$ was found in 1950 decaying as
$
  Lambda -> p^+ + pi^-.
$
$Lambda$ is a baryon by baryon number conservation.#footnote[We have since found many more baryons: $Sigma^0, Sigma^+, Sigma^-, Xi^-, Xi^0, Delta^-, Delta^0, Delta^+, Delta^(+ +)$ etc.]

These mesons and baryons were "strange" particles, since they are produced in large amounts at a time scale $t tilde 10^(-23)"s"$ and decay slowly at a time scale $t tilde 10^(-10)"s"$. Pais proposed their production and decay mechanisms differed.#footnote[and "associated production".] Now, we would say they are produced by strong interactions and decay by weak interactions. This lead to Gell-Mann and Nishijima assigning each particle a "strangeness" in 1953, which is conserved in strong and electromagnetic interactions.

$
  S(p,n) = 0",  " S(Sigma s, Lambda) = -1",  " S(Xi^-, Xi^0) = -2 \ S(eta, pi s) = 0",  " S(K^0, K^+) = 1",  " S(overline(K)^0,K^-) = -1
$
We know "strangeness" is violated in weak interactions since we have decays like
$
   Lambda & -> p^+ + pi^- \
  Sigma^+ & -> p^+ + pi^0.
$
These particles can be organised according to "The Eightfold Way" as shown in @baryon-octet, @meson-octet and @baryon-decuplet.#footnote[We have three quarks ${u,d,s}$ we have $3 times.o 3 times.o 3 = 10 times.o 8 times.o 8 times.o 1$. The octet is spin-$1/2$, while the decuplet is spin-$3/2$.]

#let baryon-octet = diagram(
  node-stroke: .1em,
  edge-stroke: .1em,
  node-fill: gray,
  node-shape: circle,
  spacing: 6em,

  node((0, 0), move(dy: -1.2em)[$n$], name: <n>, radius: .2em),
  node((1, 0), move(dy: -1.2em)[$p$], name: <p>, radius: .2em),
  node((0, 2), move(dy: .8em, dx: -.3em)[$Xi^-$], name: <Xm>, radius: .2em),
  node((1, 2), move(dy: .8em, dx: -.3em)[$Xi^0$], name: <X0>, radius: .2em),
  node(
    (-.6, 1),
    move(dy: -1em, dx: -.5em)[$Sigma^-$],
    name: <Sm>,
    radius: .2em,
  ),
  node(
    (1.6, 1),
    move(dy: -1.2em, dx: 1.2em)[$Sigma^+$],
    name: <Sp>,
    radius: .2em,
  ),
  node((.5, .95), move(dy: -1.2em)[$Sigma^0$], name: <S0>, radius: .2em),
  node((.5, 1.05), move(dy: .5em)[$Lambda$], name: <L>, radius: .2em),
  edge(<n>, <p>, "-"),
  edge(<Xm>, <X0>),
  edge(<n>, <Sm>),
  edge(<Sm>, <Xm>),
  edge(<X0>, <Sp>),
  edge(<p>, <Sp>),
  edge(
    (2.65, 2.8),
    (1.6, 1.0),
    [$Q = +1$],
    "--|>",
    mark-scale: 50%,
    label-pos: 0%,
  ),
  edge(
    (1.45, 2.8),
    (1., 2.),
    [$Q = 0$],
    "--|>",
    mark-scale: 50%,
    label-pos: 0%,
  ),
  edge(
    (0.45, 2.8),
    (0., 2.),
    [$Q = -1$],
    "--|>",
    mark-scale: 50%,
    label-pos: 0%,
  ),
  edge((-1.5, 0), (0, 0), [$S = 0$], "--|>", mark-scale: 50%, label-pos: 0%),
  edge((-1.5, 1), (-.6, 1), [$S = -1$], "--|>", mark-scale: 50%, label-pos: 0%),
  edge((-1.5, 2), (0, 2), [$S = -2$], "--|>", mark-scale: 50%, label-pos: 0%),
)

#figure(baryon-octet, caption: [The baryon octet.])<baryon-octet>

#let meson-octet = diagram(
  node-stroke: .1em,
  edge-stroke: .1em,
  node-fill: gray,
  node-shape: circle,
  spacing: 6em,

  node((0, 0), move(dy: -1.2em)[$K^0$], name: <K0>, radius: .2em),
  node((1, 0), move(dy: -1.2em)[$K^+$], name: <Kp>, radius: .2em),
  node((0, 2), move(dy: .8em, dx: -.3em)[$K^-$], name: <Km>, radius: .2em),
  node(
    (1, 2),
    move(dy: .8em, dx: -.3em)[$overline(K)^0$],
    name: <Ka0>,
    radius: .2em,
  ),
  node(
    (-.6, 1),
    move(dy: -1em, dx: -.5em)[$pi^-$],
    name: <pm>,
    radius: .2em,
  ),
  node(
    (1.6, 1),
    move(dy: -1.2em, dx: 1.2em)[$pi^+$],
    name: <pp>,
    radius: .2em,
  ),
  node((.5, .95), move(dy: -1.2em)[$pi^0$], name: <p0>, radius: .2em),
  node((.5, 1.05), move(dy: .5em)[$eta$], name: <e>, radius: .2em),
  edge(<K0>, <Kp>, "-"),
  edge(<Km>, <Ka0>),
  edge(<K0>, <pm>),
  edge(<pm>, <Km>),
  edge(<Ka0>, <pp>),
  edge(<Kp>, <pp>),
  edge(
    (2.65, 2.8),
    (1.6, 1.0),
    [$Q = +1$],
    "--|>",
    mark-scale: 50%,
    label-pos: 0%,
  ),
  edge(
    (1.45, 2.8),
    (1., 2.),
    [$Q = 0$],
    "--|>",
    mark-scale: 50%,
    label-pos: 0%,
  ),
  edge(
    (0.45, 2.8),
    (0., 2.),
    [$Q = -1$],
    "--|>",
    mark-scale: 50%,
    label-pos: 0%,
  ),
  edge((-1.5, 0), (0, 0), [$S = +1$], "--|>", mark-scale: 50%, label-pos: 0%),
  edge((-1.5, 1), (-.6, 1), [$S = 0$], "--|>", mark-scale: 50%, label-pos: 0%),
  edge((-1.5, 2), (0, 2), [$S = -1$], "--|>", mark-scale: 50%, label-pos: 0%),
)

#figure(meson-octet, caption: [The meson octet.])<meson-octet>

#let baryon-decuplet = diagram(
  node-stroke: .1em,
  edge-stroke: .1em,
  node-fill: gray,
  node-shape: circle,
  spacing: 6em,

  node((0, 0), move(dy: -1.2em)[$Delta^-$], name: <Dm>, radius: .2em),
  node((1, 0), move(dy: -1.2em)[$Delta^0$], name: <D0>, radius: .2em),
  node((2, 0), move(dy: -1.2em)[$Delta^+$], name: <Dp>, radius: .2em),
  node((3, 0), move(dy: -1.2em)[$Delta^(++)$], name: <Dpp>, radius: .2em),

  node(
    (.5, 1),
    move(dy: -1.2em, dx: 1.4em)[$Sigma^(* -)$],
    name: <Sm>,
    radius: .2em,
  ),
  node(
    (1.5, 1),
    move(dy: -1.2em, dx: .5em)[$Sigma^(* 0)$],
    name: <S0>,
    radius: .2em,
  ),
  node(
    (2.5, 1),
    move(dy: -1.2em, dx: -.5em)[$Sigma^(* +)$],
    name: <Sp>,
    radius: .2em,
  ),

  node(
    (1, 2),
    move(dy: -1.2em, dx: 1.4em)[$Xi^(* -)$],
    name: <Xm>,
    radius: .2em,
  ),
  node(
    (2, 2),
    move(dy: -1.2em, dx: -.5em)[$Xi^(* 0)$],
    name: <X0>,
    radius: .2em,
  ),

  node(
    (1.5, 3),
    move(dy: .3em, dx: 1.8em)[$Omega^-$],
    name: <Om>,
    radius: .2em,
  ),

  edge(<Dm>, <D0>),
  edge(<D0>, <Dp>),
  edge(<Dp>, <Dpp>),
  edge(<Dm>, <Sm>),
  edge(<Sm>, <Xm>),
  edge(<Xm>, <Om>),
  edge(<Om>, <X0>),
  edge(<X0>, <Sp>),
  edge(<Sp>, <Dpp>),

  edge((3.5, 1), (3, 0), [$Q = +2$], "--|>", mark-scale: 50%, label-pos: 0%),
  edge((3, 2), (2.5, 1), [$Q = +1$], "--|>", mark-scale: 50%, label-pos: 0%),
  edge((2.5, 3), (2, 2), [$Q = 0$], "--|>", mark-scale: 50%, label-pos: 0%),
  edge((2, 4), (1.5, 3), [$Q = -1$], "--|>", mark-scale: 50%, label-pos: 0%),

  edge((-1, 0), (0, 0), [$S = 0$], "--|>", mark-scale: 50%, label-pos: 0%),
  edge((-.5, 1), (.5, 1), [$S = -1$], "--|>", mark-scale: 50%, label-pos: 0%),
  edge((0, 2), (1, 2), [$S =-2$], "--|>", mark-scale: 50%, label-pos: 0%),
  edge((0.5, 3), (1.5, 3), [$S =-3$], "--|>", mark-scale: 50%, label-pos: 0%),
)

#figure(baryon-decuplet, caption: [The baryon decuplet.])<baryon-decuplet>



== Quarks
Gell-Mann and Zweig proposed hadrons were made of quarks $q$ in 1964 to explain the patterns of The Eightfold Way. We have three "flavors" which can be arranged as in @quark.#footnote[The antiquarks are "upside-down".]

#let quark = diagram(
  node-stroke: .1em,
  edge-stroke: .1em,
  node-fill: gray,
  node-shape: circle,
  spacing: 6em,

  node((0, 0), move(dy: -1.2em)[$d$], radius: .2em, name: <d>),
  node((1, 0), move(dy: -1.2em)[$u$], radius: .2em, name: <u>),
  node((0.5, 1), move(dy: .2em, dx: 1.2em)[$s$], radius: .2em, name: <s>),

  edge(<d>, <u>),
  edge(<d>, <s>),
  edge(<u>, <s>),

  edge((1.5, 1), (1, 0), [$Q=2/3$], "--|>", mark-scale: 50%, label-pos: 0%),
  edge((1, 2), (.5, 1), [$Q=-1/3$], "--|>", mark-scale: 50%, label-pos: 0%),

  edge((-1, 0), (0, 0), [$S=0$], "--|>", mark-scale: 50%, label-pos: 0%),
  edge((-.5, 1), (.5, 1), [$S=-1$], "--|>", mark-scale: 50%, label-pos: 0%),
)

#figure(quark, caption: [The quarks.])<quark>



All baryons are composed of three quarks, while all mesons are composed of a quark and an anti-quark
$
  B & = q q q \
  M & = q overline(q)
$
With three quarks we can construct ten baryons and nine mesons.


There are two problems with quarks. The first problem is that we have never observed a quark. The only experimental evidence we have is that the proton has substructure indicating the existence of quarks. The second is that they appear to violate the Pauli exclusion principle since e.g. $Delta^(++)$ contains three $u$ with spin-$1/2$. To solve this Greenberg proposed quarks carry an additional degree of freedom. We call this color and quarks can be red, green or blue. All particles we see in nature must be colorless meaning they have zero color or equal amounts of each.

== The Standard Model
As seen in the previous section a bunch of particles were found. The basic building blocks being leptons and quarks
$
   "quarks:"#h(1em) & u,d,c,s,t,b \
  "leptons:"#h(1em) & e, mu, tau, nu_e, nu_mu, nu_tau
$
With those not mentioned previously eventually being found.

The weak interaction is mediated by the intermediate vector bosons $W^+\/W^-$ and $Z$, the electromagnetic interaction is mediated by photons $gamma$, and strong interaction is mediated by eight gluons.

The final particle found in 2012 was the Higgs boson which is responsible for particles having mass!

