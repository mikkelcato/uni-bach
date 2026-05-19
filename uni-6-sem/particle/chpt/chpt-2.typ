#import "chpt-temp.typ": *
#show: chpt-note.with()

= Dynamics Overview
As far as we know nature is governed by four fundamental forces: _strong_, _electromagnetic_, _weak_, and _gravitational_. We will only consider the first three as quantum gravity is a different beast which has yet to be conquered.

== QED
_Quantum electrodynamics_ is the theory that describes the electromagnetic force and is by far the most simple and best understood of the aforementioned theories.

All electromagnetic phenomena are reducible to the process in @f1
#let f1 = feynman(
  (
    vertex("i1"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),
    edge("i1", "a", type: "fermion", label: $e^-$),
    edge("a", "f2", type: "fermion", label: $e^-$),
    edge("a", "f1", type: "photon", label: $gamma$),
  ),
)

#figure(
  scale(f1, 100%),
  caption: [A charged particle $e^-$ enters and emits or absorbs a photon $gamma$ before exiting. This is an example of a _Feynman diagram_.],
)<f1>

This is the only allowed vertex in QED and by patching these processes together we get more complicated interactions. An example is _Møller scattering_ seen in @f2

#let f2 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("a"),
    vertex("b"),
    vertex("f1"),
    vertex("f2"),
    edge("i1", "a", type: "fermion", label: $e^-$),
    edge("i2", "b", type: "fermion", label: $e^-$),
    edge("a", "f1", type: "fermion", label: $e^-$),
    edge("b", "f2", type: "fermion", label: $e^-$),
    edge("a", "b", type: "photon", label: $gamma$),
  ),
  orientation: "vertical",
)

#figure(
  scale(f2, 100%),
  caption: [Feynman diagram for Møller scattering. Two electrons $e^-$ scatter by exchanging a photon $gamma$.],
) <f2>

Another example is _Bhabha scattering_ seen in @f3

#let f3 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("a"),
    vertex("b"),
    vertex("f1"),
    vertex("f2"),
    edge("i1", "a", type: "fermion", label: $e^-$),
    edge("i2", "a", type: "antifermion", label: $e^+$),
    edge("a", "b", type: "photon", label: $gamma$),
    edge("b", "f1", type: "fermion", label: $e^-$),
    edge("b", "f2", type: "antifermion", label: $e^+$),
  ),
)

#figure(
  scale(f3, 100%),
  caption: [Feynman diagram for Bhabha scattering. An electron $e^-$ and positron $e^+$ annihilate to form a photon $gamma$ which produces a new electron-positron pair. ],
)<f3>

This is just the diagram for Møller scattering on its side! We interpret the particles going _backwards_ in time as antiparticles. A second diagram also contributes to Bhabha scattering. This diagram is seen in @f4

#let f4 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("a"),
    vertex("b"),
    vertex("f1"),
    vertex("f2"),
    edge("i1", "a", type: "fermion", label: $e^-$),
    edge("i2", "b", type: "antifermion", label: $e^+$),
    edge("a", "b", type: "photon", label: $gamma$),
    edge("a", "f1", type: "fermion", label: $e^-$),
    edge("b", "f2", type: "antifermion", label: $e^+$),
  ),
  orientation: "vertical",
)

#figure(
  scale(f4, 100%),
  caption: [Feynman diagram for Bhabha scattering. An electron $e^-$ and positron $e^+$ scatter by exchanging a photon $gamma$.],
)<f4>

As a final example we consider _Compton scattering_ seen in @f5

#let f5 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("a"),
    vertex("b"),
    vertex("f1"),
    vertex("f2"),
    edge("i1", "a", type: "fermion", label: $e^-$),
    edge("i2", "a", type: "photon", label: $gamma$),
    edge("a", "b", type: "fermion", label: $e^-$),
    edge("b", "f1", type: "fermion", label: $e^-$),
    edge("b", "f2", type: "photon", label: $gamma$),
  ),
)

#figure(
  scale(f5, 100%),
  caption: [Feynman diagram for Compton scattering.],
)<f5>

Almost identical diagrams can be made for _pair annihilation_ and _pair production_.

We could allow additional vertices in the above to get more complicated diagrams. However each additional vertex comes with a factor of $alpha$ being the _fine-structure constant_ meaning diagrams with additional vertices are less important. Since the observed process is only represented by _external lines_. The _internal lines_ represent _virtual particles_ and describes the responsible mechanism. The Feynman diagram is a _purely symbolic tool_! We compute things from Feynman diagrams using _Feynman rules_ which we will discuss later. These rules enforce energy and momentum conservation at each vertex.

== QCD
_Quantum chromodynamics_ is the theory that describes the strong force. In QCD _color_ plays the role of charge and the fundamental process is seen in @f6

#let f6 = feynman(
  (
    vertex("i1"),
    vertex("a"),
    vertex("f2"),
    vertex("f1"),
    edge("i1", "a", type: "fermion", label: $q$),
    edge("a", "f1", type: "fermion", label: $q$),
    edge("a", "f2", type: "gluon", label: $g$),
  ),
)

#figure(
  scale(f6, 100%),
  caption: [_Primitive vertex_ for QCD.],
)<f6>

We can again patch these together to form more complication interactions. An example is seen in @f7

#let f7 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("a"),
    vertex("b"),
    vertex("f1"),
    vertex("f2"),
    edge("i1", "a", type: "fermion", label: $q$),
    edge("a", "f1", type: "fermion", label: $q$),
    edge("a", "b", type: "gluon", label: $g$),
    edge("i2", "b", type: "fermion", label: $q$),
    edge("b", "f2", type: "fermion", label: $q$),
  ),
  orientation: "vertical",
)

#figure(
  scale(f7, 100%),
  caption: [Feynman diagram for quark binding.],
)<f7>

There exists three types of color and quarks can change color with the difference being carried by gluons. An example is seen in @f8


#let f8 = feynman(
  (
    vertex("i1"),
    vertex("a"),
    vertex("f2"),
    vertex("f1"),
    edge("i1", "a", type: "fermion", label: $q(b)$),
    edge("a", "f1", type: "fermion", label: $q(r)$),
    edge("a", "f2", type: "gluon", label: $g(b, overline(r))$),
  ),
)

#figure(
  scale(f8, 100%),
  caption: [Conversion of a blue quark into a red quark with the help of a _bicolored_ gluon. This implies the existence of nine different gluons. However there are only eight.],
)<f8>

Due to gluons being colored they also couple to themselves. (unlike the photon $gamma$ which is neutral!) We get two new kinds of vertices as seen in @f10

#let f9 = feynman(
  (
    vertex("i1"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),
    edge("i1", "a", type: "gluon"),
    edge("a", "f1", type: "gluon"),
    edge("a", "f2", type: "gluon"),
  ),
)

#let f10 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),
    edge("i1", "a", type: "gluon"),
    edge("i2", "a", type: "gluon"),
    edge("a", "f1", type: "gluon"),
    edge("a", "f2", type: "gluon"),
  ),
)

#figure(
  grid(
    columns: 2,
    gutter: 12pt,
    scale(f9, 100%), scale(f10, 100%),
  ),
  caption: [Two additional kinds of primitive vertices in QCD.],
)<f10>

This additional coupling necessarily makes QCD way more complicated than QED.

== Weak interactions
All quarks and leptons carry a _weak charge_ meaning they can all interact through weak interactions. We distinguish between _charged_ interactions mediated by $W^(plus.minus)$ and _neutral_ interactions mediated by $Z$.

The fundamental neutral vertex is seen in @f11

#let f11 = feynman(
  (
    vertex("i1"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),

    edge("i1", "a", type: "fermion", label: $f$),
    edge("a", "f2", type: "fermion", label: $f$),
    edge("a", "f1", type: "boson", label: $Z$),
  ),
)

#figure(
  scale(f11, 100%),
  caption: [The fundamental neutral vertex. Here $f in {q, l}$],
)<f11>

Then for any process mediated by $gamma$ there is an analogous process mediated by $Z$. An example of a more complicated interaction is seen in @f12

#let f12 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("i3"),
    vertex("i4"),
    vertex("a"),
    vertex("b"),
    vertex("c"),
    vertex("d"),
    vertex("f1"),
    vertex("f2"),
    vertex("f3"),
    vertex("f4"),

    edge("i1", "a", type: "fermion", label: $nu_mu$),
    edge("a", "f1", type: "fermion", label: $nu_mu$),
    edge("a", "b", type: "boson", label: $Z$),
    edge("i2", "b", type: "fermion", label: $u$),
    edge("b", "f2", type: "fermion", label: $u$),
    edge("i3", "c", type: "fermion", label: $u$),
    edge("c", "f3", type: "fermion", label: $u$),
    edge("i4", "d", type: "fermion", label: $d$),
    edge("d", "f4", type: "fermion", label: $d$),
  ),
  layout: "layered",
)

#figure(
  scale(f12, 100%),
  caption: [The scattering of a muon-neutrino $nu_mu$ with a proton $p$. Here the two up-quarks $u$ are called _spectator quarks_.],
)<f12>

The fundamental charged vertices are seen in @f14

#let f13 = feynman(
  (
    vertex("i1"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),

    edge("i1", "a", type: "fermion", label: $l^-$),
    edge("a", "f2", type: "fermion", label: $nu_l$),
    edge("a", "f1", type: "boson", label: $W^-$),
  ),
)

#let f14 = feynman(
  (
    vertex("i1"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),

    edge("i1", "a", type: "fermion", label: $q^(-1\/3)$),
    edge("a", "f2", type: "fermion", label: $q^(+2\/3)$),
    edge("a", "f1", type: "boson", label: $W^-$),
  ),
)

#figure(
  grid(
    columns: 2,
    gutter: 12pt,
    scale(f13, 100%), scale(f14, 100%),
  ),
  caption: [The fundamental charged vertices. The opposite direction is mediated by $W^+$. The quarks in the second diagram have the same color but different flavors since flavor is simply not conserved.],
)<f14>

An example of a very important process is seen in @f15

#let f15 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("a"),
    vertex("b"),
    vertex("f1"),
    vertex("f2"),

    edge("i1", "a", type: "fermion", label: $d$),
    edge("a", "f1", type: "fermion", label: $u$),
    edge("a", "b", type: "boson", label: $W^-$),
    edge("i2", "b", type: "fermion", label: $nu_e$),
    edge("b", "f2", type: "fermion", label: $e^-$),
  ),
  orientation: "vertical",
)

#figure(
  scale(f15, 100%),
  caption: [The _semileptonic_ process $d + nu_e -> u + e^-$.],
)<f15>

This diagram gives us $beta$-decay $n -> p + e^- + overline(nu)_e$! Since $u d d -> u u d$ with two spectator quarks.

The above picture is a bit simplistic since it seems like the weak interaction only couple pairs ${u,d}$. What actually happens is the coupling of pairs ${u,d'}$ with
$
  vec(d', s', b') = mat(V_(u d), V_(u s), V_(u b); V_(c d), V_(c s), V_(c b); V_(t d), V_(t s), V_(t b)) vec(d, s, b)
$
This is called the _CKM-mechanism_ and accounts for the decay of $Delta$ and $Omega^-$.

The $W^(plus.minus)$ and $Z$ also couple directly and since $W^(plus.minus)$ is charged it also couples to photons $gamma$. These are all seen in @f21

#let f16 = feynman(
  (
    vertex("i1"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),

    edge("i1", "a", type: "boson", label: $W$),
    edge("a", "f2", type: "boson", label: $W$),
    edge("a", "f1", type: "boson", label: $Z$),
  ),
)

#let f17 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),

    edge("i1", "a", type: "boson", label: $W$),
    edge("i2", "a", type: "boson", label: $W$),
    edge("a", "f1", type: "boson", label: $W$),
    edge("a", "f2", type: "boson", label: $W$),
  ),
)

#let f18 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),

    edge("i2", "a", type: "boson", label: $Z$),
    edge("i1", "a", type: "boson", label: $W$),
    edge("a", "f1", type: "boson", label: $Z$),
    edge("a", "f2", type: "boson", label: $W$),
  ),
)

#let f19 = feynman(
  (
    vertex("i1"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),

    edge("i1", "a", type: "boson", label: $W$),
    edge("a", "f2", type: "boson", label: $W$),
    edge("a", "f1", type: "boson", label: $gamma$),
  ),
)

#let f20 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),

    edge("i2", "a", type: "boson", label: $gamma$),
    edge("i1", "a", type: "boson", label: $W$),
    edge("a", "f1", type: "boson", label: $Z$),
    edge("a", "f2", type: "boson", label: $W$),
  ),
)

#let f21 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),

    edge("i2", "a", type: "boson", label: $gamma$),
    edge("i1", "a", type: "boson", label: $W$),
    edge("a", "f1", type: "boson", label: $gamma$),
    edge("a", "f2", type: "boson", label: $W$),
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

These are not super important for us however.

== Decays and laws
All particles want to disintegrate into some lighter particles unless some law prevents it. The $gamma$ is stable since $m_gamma = 0$ while $e^-$ is stable since it is the lightest charged particle so by charge conservation it can not decay. Similarly the $p$ is stable due to baryon number conservation. However most particles want to decay and most _exotic_ particles only exist in relatively small amounts.

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
We have seen all fundamental vertices and by seeing what is conserved in these we can figure out what is conserved in general. We have conservation of

1. charge.

2. color.

3. baryon number (number of quarks).

4. lepton number.

5. flavor in strong and electromagnetic.

== The Higgs particle
We briefly mention the fundamental vertices of the Higgs particle for completeness. The Higgs particle is neutral and can couple with itself as in @f23

#let f22 = feynman(
  (
    vertex("i1"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),
    edge("i1", "a", type: "scalar", label: $h$),
    edge("a", "f1", type: "scalar", label: $h$),
    edge("a", "f2", type: "scalar", label: $h$),
  ),
)

#let f23 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),
    edge("i1", "a", type: "scalar", label: $h$),
    edge("i2", "a", type: "scalar", label: $h$),
    edge("a", "f1", type: "scalar", label: $h$),
    edge("a", "f2", type: "scalar", label: $h$),
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

It also couples to $W^(plus.minus)$ and $Z$ as in @f27

#let f24 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),

    edge("i1", "a", type: "scalar", label: $h$),
    edge("i2", "a", type: "boson", label: $W^+$),
    edge("a", "f1", type: "scalar", label: $h$),
    edge("a", "f2", type: "boson", label: $W^-$),
  ),
)

#let f25 = feynman(
  (
    vertex("i1"),
    vertex("i2"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),

    edge("i1", "a", type: "scalar", label: $h$),
    edge("i2", "a", type: "boson", label: $Z$),
    edge("a", "f1", type: "scalar", label: $h$),
    edge("a", "f2", type: "boson", label: $Z$),
  ),
)

#let f26 = feynman(
  (
    vertex("i1"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),

    edge("i1", "a", type: "scalar", label: $h$),
    edge("a", "f1", type: "boson", label: $W^-$),
    edge("a", "f2", type: "boson", label: $W^+$),
  ),
)

#let f27 = feynman(
  (
    vertex("i1"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),

    edge("i1", "a", type: "scalar", label: $h$),
    edge("a", "f2", type: "boson", label: $Z$),
    edge("a", "f1", type: "boson", label: $Z$),
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

We also have an interaction with fermions as in @f28

#let f28 = feynman(
  (
    vertex("i1"),
    vertex("a"),
    vertex("f1"),
    vertex("f2"),

    edge("i1", "a", type: "scalar", label: $h$),
    edge("a", "f2", type: "fermion", label: $f$),
    edge("a", "f1", type: "antifermion", label: $overline(f)$),
  ),
)

#figure(
  scale(f28, 100%),
  caption: [Coupling between fermions and the Higgs particle.],
)<f28>

We care about the Higgs particle since the interaction in @f28 is dependent on $m_f$ implying the mass of a given particle is dependent on how much they interact with the Higgs particle.
