#import "../../temp.typ": *
#show: chpt-note.with()

= Dynamics
As far as we know nature is governed by four fundamental forces: electromagnetic, strong, weak, and gravitational. We only consider the first three.

== QED
Quantum electrodynamics is the theory that describes the electromagnetic force and is by far the most simple and best understood of the aforementioned theories.

All electromagnetic phenomena are reducible to the process in @f1. #footnote[Here, we use $e^-$ as a placeholder "charged particle". We could have $u -> u + gamma$.]
#let f1 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),
    feynman.edge("i1", "a", type: "fermion", label: $e^-$),
    feynman.edge("a", "f2", type: "fermion", label: $e^-$),
    feynman.edge("a", "f1", type: "photon", label: $gamma$),
  ),
)

#figure(
  scale(f1, 100%),
  caption: [A charged particle $e^-$ enters and emits or absorbs a photon $gamma$ before exiting. This is an example of a Feynman diagram. Here time increases towards the right, and there is a "femion flow". We see charge is conserved.],
)<f1>

This is the only allowed vertex in QED and by patching these processes together we get more complicated interactions. An example is Møller scattering ($e^- e^-$ scattering) seen in @f2.

#let f2 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("a"),
    feynman.vertex("b"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),
    feynman.edge("i1", "a", type: "fermion", label: $e^-$),
    feynman.edge("i2", "b", type: "fermion", label: $e^-$),
    feynman.edge("a", "f1", type: "fermion", label: $e^-$),
    feynman.edge("b", "f2", type: "fermion", label: $e^-$),
    feynman.edge("a", "b", type: "photon", label: $gamma$),
  ),
  orientation: "vertical",
)

#figure(
  scale(f2, 100%),
  caption: [Feynman diagram for Møller scattering. Two electrons $e^-$ scatter by exchanging a photon $gamma$. This is an example of a tree diagram, diagrams with "loops" are called loop diagrams.],
) <f2>

Another example is Bhabha scattering ($e^- e^+$ scattering) seen in @f3.

#let f3 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("a"),
    feynman.vertex("b"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),
    feynman.edge("i1", "a", type: "fermion", label: $e^-$),
    feynman.edge("i2", "a", type: "antifermion", label: $e^+$),
    feynman.edge("a", "b", type: "photon", label: $gamma$),
    feynman.edge("b", "f1", type: "fermion", label: $e^-$),
    feynman.edge("b", "f2", type: "antifermion", label: $e^+$),
  ),
)

#figure(
  scale(f3, 100%),
  caption: [Feynman diagram for Bhabha scattering. An electron $e^-$ and positron $e^+$ annihilate to form a photon $gamma$ which produces a new electron-positron pair. Here, antiparticles ($e^+$) travel "backward" in time.],
)<f3>

This is just the diagram for Møller scattering on its side! We interpret the particles travelling "backward in time" as antiparticles. A second diagram also contributes to Bhabha scattering. This diagram is seen in @f4.

#let f4 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("a"),
    feynman.vertex("b"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),
    feynman.edge("i1", "a", type: "fermion", label: $e^-$),
    feynman.edge("i2", "b", type: "antifermion", label: $e^+$),
    feynman.edge("a", "b", type: "photon", label: $gamma$),
    feynman.edge("a", "f1", type: "fermion", label: $e^-$),
    feynman.edge("b", "f2", type: "antifermion", label: $e^+$),
  ),
  orientation: "vertical",
)

#figure(
  scale(f4, 100%),
  caption: [Feynman diagram for Bhabha scattering. An electron $e^-$ and positron $e^+$ scatter by exchanging a photon $gamma$.],
)<f4>

As a final example we consider Compton scattering ($e^- gamma$ scattering) seen in @f5.

#let f5 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("a"),
    feynman.vertex("b"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),
    feynman.edge("i1", "a", type: "fermion", label: $e^-$),
    feynman.edge("i2", "a", type: "photon", label: $gamma$),
    feynman.edge("a", "b", type: "fermion", label: $e^-$),
    feynman.edge("b", "f1", type: "fermion", label: $e^-$),
    feynman.edge("b", "f2", type: "photon", label: $gamma$),
  ),
)

#figure(
  scale(f5, 100%),
  caption: [Feynman diagram for Compton scattering.],
)<f5>

Almost identical diagrams can be made for pair annihilation and pair production.

We could allow additional vertices in the above to get more complicated diagrams. However, each additional vertex comes with a factor of $alpha = e^2\/hbar c tilde 1/137$ being the fine-structure constant meaning diagrams with additional vertices are less important. The observed process is only represented by external lines. The internal lines represent "virtual particles" and describe the responsible mechanism. The Feynman diagram is a purely symbolic tool! We compute things from Feynman diagrams using Feynman rules which we will discuss later. These rules enforce energy and momentum conservation at each vertex. #footnote[This implies the QED vertex itself cannot represent a process, i.e. $e^- -> e^- + gamma$ and $e^- + e^+ -> gamma$ are unphysical.]

== QCD
Quantum chromodynamics is the theory that describes the strong force. Now, in QCD color plays the role of charge and the fundamental process is seen in @f6.

#let f6 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("a"),
    feynman.vertex("f2"),
    feynman.vertex("f1"),
    feynman.edge("i1", "a", type: "fermion", label: $q$),
    feynman.edge("a", "f1", type: "fermion", label: $q$),
    feynman.edge("a", "f2", type: "gluon", label: $g$),
  ),
)

#figure(
  scale(f6, 100%),
  caption: [Fundamental vertex for QCD. Here, flavor is conserved.],
)<f6>

We can again patch these together to form more complicated interactions. An example is seen in @f7.

#let f7 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("a"),
    feynman.vertex("b"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),
    feynman.edge("i1", "a", type: "fermion", label: $q$),
    feynman.edge("a", "f1", type: "fermion", label: $q$),
    feynman.edge("a", "b", type: "gluon", label: $g$),
    feynman.edge("i2", "b", type: "fermion", label: $q$),
    feynman.edge("b", "f2", type: "fermion", label: $q$),
  ),
  orientation: "vertical",
)

#figure(
  scale(f7, 100%),
  caption: [Feynman diagram for quark binding.],
)<f7>

We have three types of color and quarks#footnote[Note, leptons are colorless i.e. do not experience strong interactions.] can change color with the difference being carried by gluons, as seen in @f8.


#let f8 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("a"),
    feynman.vertex("f2"),
    feynman.vertex("f1"),
    feynman.edge("i1", "a", type: "fermion", label: $q(b)$),
    feynman.edge("a", "f1", type: "fermion", label: $q(r)$),
    feynman.edge("a", "f2", type: "gluon", label: $g(b, overline(r))$),
  ),
)

#figure(
  scale(f8, 100%),
  caption: [Conversion of a blue quark into a red quark with the help of a "bicolored gluon". This implies the existence of nine different gluons. However, there are only eight. We see color is conserved.],
)<f8>

Now, since gluons are colored they can couple with themselves.#footnote[Unlike the photon, which is neutral.] We get two new kinds of vertices as seen in @f10.

#let f9 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),
    feynman.edge("i1", "a", type: "gluon"),
    feynman.edge("a", "f1", type: "gluon"),
    feynman.edge("a", "f2", type: "gluon"),
  ),
)

#let f10 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),
    feynman.edge("i1", "a", type: "gluon"),
    feynman.edge("i2", "a", type: "gluon"),
    feynman.edge("a", "f1", type: "gluon"),
    feynman.edge("a", "f2", type: "gluon"),
  ),
)

#figure(
  grid(
    columns: 2,
    gutter: 12pt,
    scale(f9, 100%), scale(f10, 100%),
  ),
  caption: [Two additional kinds of fundamental vertices in QCD. Here, the color is also conserved.],
)<f10>

This additional coupling necessarily makes QCD way more complicated than QED.

With these vertices one can show $p p -> p p$ (with $p = u u d$) through the exchange of a pion $pi^0 = u overline(u)$.

== Weak interactions
All fermions (quarks and leptons) carry a "weak charge" meaning they can all interact through weak interactions. Now, we distinguish between charged interactions mediated by $W^(plus.minus)$ and neutral interactions mediated by $Z$.

The fundamental neutral vertex is seen in @f11.

#let f11 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),

    feynman.edge("i1", "a", type: "fermion", label: $f$),
    feynman.edge("a", "f2", type: "fermion", label: $f$),
    feynman.edge("a", "f1", type: "boson", label: $Z$),
  ),
)

#figure(
  scale(f11, 100%),
  caption: [The fundamental neutral vertex. Here, $f$ is a fermion.],
)<f11>

Note, that for any process mediated by $gamma$ there is an analogous process mediated by $Z$. An example of a more complicated interaction is seen in @f12.

#let f12 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("i3"),
    feynman.vertex("i4"),
    feynman.vertex("a"),
    feynman.vertex("b"),
    feynman.vertex("c"),
    feynman.vertex("d"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),
    feynman.vertex("f3"),
    feynman.vertex("f4"),

    feynman.edge("i1", "a", type: "fermion", label: $nu_mu$),
    feynman.edge("a", "f1", type: "fermion", label: $nu_mu$),
    feynman.edge("a", "b", type: "boson", label: $Z$),
    feynman.edge("i2", "b", type: "fermion", label: $u$),
    feynman.edge("b", "f2", type: "fermion", label: $u$),
    feynman.edge("i3", "c", type: "fermion", label: $u$),
    feynman.edge("c", "f3", type: "fermion", label: $u$),
    feynman.edge("i4", "d", type: "fermion", label: $d$),
    feynman.edge("d", "f4", type: "fermion", label: $d$),
  ),
  layout: "layered",
)

#figure(
  scale(f12, 100%),
  caption: [The scattering of a muon-neutrino $nu_mu$ with a proton $p$. Here, the two up-quarks $u$ are called spectator quarks.],
)<f12>

The fundamental charged vertices are seen in @f14.

#let f13 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),

    feynman.edge("i1", "a", type: "fermion", label: $l^-$),
    feynman.edge("a", "f2", type: "fermion", label: $nu_l$),
    feynman.edge("a", "f1", type: "boson", label: $W^-$),
  ),
)

#let f14 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),

    feynman.edge("i1", "a", type: "fermion", label: $q^(-1\/3)$),
    feynman.edge("a", "f2", type: "fermion", label: $q^(+2\/3)$),
    feynman.edge("a", "f1", type: "boson", label: $W^-$),
  ),
)

#figure(
  grid(
    columns: 2,
    gutter: 12pt,
    scale(f13, 100%), scale(f14, 100%),
  ),
  caption: [The fundamental charged vertices. The opposite direction is mediated by $W^+$. The quarks in the second diagram have the same color but different flavors.],
)<f14>

An example of a very important process is seen in @f15

#let f15 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("a"),
    feynman.vertex("b"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),

    feynman.edge("i1", "a", type: "fermion", label: $d$),
    feynman.edge("a", "f1", type: "fermion", label: $u$),
    feynman.edge("a", "b", type: "boson", label: $W^-$),
    feynman.edge("i2", "b", type: "fermion", label: $nu_e$),
    feynman.edge("b", "f2", type: "fermion", label: $e^-$),
  ),
  orientation: "vertical",
)

#figure(
  scale(f15, 100%),
  caption: [The semileptonic process $d + nu_e -> u + e^-$.],
)<f15>

This diagram gives us $beta$-decay $n -> p + e^- + overline(nu)_e$! Since $u d d -> u u d$ with two spectator quarks. We also have $mu^- + nu_e -> e^- + nu_mu$ implying $mu^- -> e^- + nu_mu + overline(nu)_e$.

Now, what we have done so far implies weak interactions only couple pairs $(u,d)$. However, weak interactions actually couple pairs $(u,d')$ with#footnote[Also, $(c, s')$ and $(t, b')$.]
$
  vec(d', s', b') = mat(V_(u d), V_(u s), V_(u b); V_(c d), V_(c s), V_(c b); V_(t d), V_(t s), V_(t b)) vec(d, s, b)
$
Which we call the CKM-mechanism and accounts for the decay of $Lambda -> p^+ + pi^-$ and $Omega^- -> Lambda + K^-$.#footnote[Note, this allows the construction of "penguin diagrams".]

The $W^(plus.minus)$ and $Z$ also couple directly, and since the $W^(plus.minus)$ are charged they also couples to photons $gamma$. These interactions are all shown in @f21.

#let f16 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),

    feynman.edge("i1", "a", type: "boson", label: $W$),
    feynman.edge("a", "f2", type: "boson", label: $W$),
    feynman.edge("a", "f1", type: "boson", label: $Z$),
  ),
)

#let f17 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),

    feynman.edge("i1", "a", type: "boson", label: $W$),
    feynman.edge("i2", "a", type: "boson", label: $W$),
    feynman.edge("a", "f1", type: "boson", label: $W$),
    feynman.edge("a", "f2", type: "boson", label: $W$),
  ),
)

#let f18 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),

    feynman.edge("i2", "a", type: "boson", label: $Z$),
    feynman.edge("i1", "a", type: "boson", label: $W$),
    feynman.edge("a", "f1", type: "boson", label: $Z$),
    feynman.edge("a", "f2", type: "boson", label: $W$),
  ),
)

#let f19 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),

    feynman.edge("i1", "a", type: "boson", label: $W$),
    feynman.edge("a", "f2", type: "boson", label: $W$),
    feynman.edge("a", "f1", type: "boson", label: $gamma$),
  ),
)

#let f20 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),

    feynman.edge("i2", "a", type: "boson", label: $gamma$),
    feynman.edge("i1", "a", type: "boson", label: $W$),
    feynman.edge("a", "f1", type: "boson", label: $Z$),
    feynman.edge("a", "f2", type: "boson", label: $W$),
  ),
)

#let f21 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),

    feynman.edge("i2", "a", type: "boson", label: $gamma$),
    feynman.edge("i1", "a", type: "boson", label: $W$),
    feynman.edge("a", "f1", type: "boson", label: $gamma$),
    feynman.edge("a", "f2", type: "boson", label: $W$),
  ),
)



#figure(
  grid(
    columns: 3,
    gutter: 12pt,
    scale(f16, 100%), scale(f17, 100%), scale(f18, 100%),
    scale(f19, 100%), scale(f20, 100%), scale(f21, 100%),
  ),
  caption: [The couplings between $W^(plus.minus), Z$ and $gamma$.],
)<f21>


== Decays and laws
All particles want to disintegrate into some lighter particles unless some law prevents it. The photon is stable since $m_gamma = 0$ while the electron is stable since it is the lightest charged particle so by charge conservation it cannot decay. Likewise, the proton is stable due to baryon number conservation. However, most particles want to decay and most "exotic" particles only exist in relatively small amounts.

Any decay is governed by one of the three interactions described above. Some examples are
$
  Delta^(+ +) & -> p + pi^+"   strong interaction" \
         pi^0 & -> gamma + gamma"   electromagnetic interaction" \
      Sigma^- & -> n + e^- + overline(nu)_e"   weak interaction"
$
Some interactions are less obvious
$
  Sigma^- & -> n + pi^-"   weak interaction" \
  Delta^- & -> n + pi^-"   strong interaction"
$
where in the $Sigma^-$-decay flavor changes.

We have seen all fundamental vertices and by seeing what is conserved in these we can figure out what is conserved in general. We have conservation of

1. charge.

2. color.

3. baryon number (number of quarks).

4. lepton number.

5. flavor in strong and electromagnetic.

== The Higgs particle
We mention the fundamental vertices of the Higgs particle for completeness. The Higgs particle is neutral and can couple with itself as in @f23.

#let f22 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),
    feynman.edge("i1", "a", type: "scalar", label: $h$),
    feynman.edge("a", "f1", type: "scalar", label: $h$),
    feynman.edge("a", "f2", type: "scalar", label: $h$),
  ),
)

#let f23 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),
    feynman.edge("i1", "a", type: "scalar", label: $h$),
    feynman.edge("i2", "a", type: "scalar", label: $h$),
    feynman.edge("a", "f1", type: "scalar", label: $h$),
    feynman.edge("a", "f2", type: "scalar", label: $h$),
  ),
)

#figure(
  grid(
    columns: 2,
    gutter: 12pt,
    scale(f22, 100%), scale(f23, 100%),
  ),
  caption: [Self-coupling of the Higgs particle.],
)<f23>

The Higgs also couples to $W^(plus.minus)$ and $Z$ as in @f27.

#let f24 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),

    feynman.edge("i1", "a", type: "scalar", label: $h$),
    feynman.edge("i2", "a", type: "boson", label: $W^+$),
    feynman.edge("a", "f1", type: "scalar", label: $h$),
    feynman.edge("a", "f2", type: "boson", label: $W^-$),
  ),
)

#let f25 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("i2"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),

    feynman.edge("i1", "a", type: "scalar", label: $h$),
    feynman.edge("i2", "a", type: "boson", label: $Z$),
    feynman.edge("a", "f1", type: "scalar", label: $h$),
    feynman.edge("a", "f2", type: "boson", label: $Z$),
  ),
)

#let f26 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),

    feynman.edge("i1", "a", type: "scalar", label: $h$),
    feynman.edge("a", "f1", type: "boson", label: $W^-$),
    feynman.edge("a", "f2", type: "boson", label: $W^+$),
  ),
)

#let f27 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),

    feynman.edge("i1", "a", type: "scalar", label: $h$),
    feynman.edge("a", "f2", type: "boson", label: $Z$),
    feynman.edge("a", "f1", type: "boson", label: $Z$),
  ),
)


#figure(
  grid(
    columns: 2,
    gutter: 12pt,
    scale(f24, 100%), scale(f25, 100%),
    scale(f26, 100%), scale(f27, 100%),
  ),
  caption: [The couplings between $W^(plus.minus), Z$ and $h$.],
)<f27>

We also have an interaction with fermions as in @f28.

#let f28 = feynman.feynman(
  (
    feynman.vertex("i1"),
    feynman.vertex("a"),
    feynman.vertex("f1"),
    feynman.vertex("f2"),

    feynman.edge("i1", "a", type: "scalar", label: $h$),
    feynman.edge("a", "f2", type: "fermion", label: $f$),
    feynman.edge("a", "f1", type: "antifermion", label: $overline(f)$),
  ),
)

#figure(
  scale(f28, 100%),
  caption: [Coupling between fermions and the Higgs particle.],
)<f28>

We care about the Higgs since the interaction in @f28 is dependent on $m_f$ implying the mass of a given particle is dependent on how much they interact with the Higgs.#footnote[The Higgs can be created by "gluon fusion" $2 g ->^(t "loop") h$, we can then observe decays $h ->^(W "loop")_(t "loop") 2 gamma$.]
