//**** init-ting
#import "temp.typ": *
#import "fig.typ": *

#show: thesis.with(
  title: [
    Entanglement, Black Holes#linebreak() and the ER = EPR Conjecture
  ],
  name: [
    Mikkel K. Holst
  ],
  prof: [
    Martin S. Sloth
  ],
)


= Introduction
The relationship between quantum mechanics and gravity remains one of the central open problems in theoretical physics. Currently, we have two theories that describe their respective domains with remarkable success, yet all attempts to reconcile them have run into difficulties. This tension becomes especially clear in the context of black holes, where seemingly simple assumptions lead to deep conceptual contradictions. Historically, however, such contradictions have often led to progress and driven major advances in our understanding of nature, including modern ideas suggesting that spacetime itself may emerge from quantum entanglement @vanraamsdonk2010commentsquantumgravityentanglement.

One example is the black hole information paradox, originating from Hawking's discovery that black holes emit thermal radiation and may delete information in the process @1975CMaPh..43..199H @HawkingBHinfo. The paradox was later made precise by Almheiri, Marolf, Polchinski, and Sully in the now famed AMPS paradox, which  highlights the crucial role of entanglement @Almheiri_2013.

This thesis attempts to motivate and introduce the $"ER"="EPR"$ conjecture proposed by Susskind and Maldacena as a possible resolution to the AMPS paradox @Maldacena_2013. The conjecture can be summarised as the proposed equivalence between entanglement (EPR) and Einstein-Rosen bridges (ER) @EPRpaper @ERbridge, offering a concrete realisation of the broader idea that spacetime may be connected to entanglement.

We proceed in roughly three parts. We discuss the $"EPR"$ and $"ER"$ aspects of $"ER"="EPR"$ in chapters $1$ and $2$ respectively, before turning to the information paradox and $"ER"="EPR"$ in chapters $3$ and $4$.

For simplicity, and due to matters of scope, we minimise the use of the AdS/CFT correspondence. However, since it plays an important role in motivating $"ER"="EPR"$, we briefly introduce the correspondence and consider an example in this context. Throughout, we use natural units with $c = hbar = k_B = 1$, while keeping $G$ explicit.

#pagebreak()
= Quantum Entanglement
Quantum mechanics usually concerns quantum states $ket(Psi)$ and their dynamics. However, such a framework is less than ideal when considering mixtures of states, or ensembles, and composite systems. This chapter introduces the nicer density matrix formalism which is well suited to handle these, especially when subsystems may be inaccessible to an observer.#footnote[i.e. the interior of a black hole.] We also discuss several constructions which will play an important role later, including the thermofield double state and the Rindler decomposition of the Minkowski vacuum.

== Quantum ensembles
The idea of an ensemble is best illustrated with an example. Consider a Stern-Gerlach-like setup in which an oven emits spin-$1/2$ particles with randomly oriented spins. Before any measurement is performed, the beam is described by the mixture or ensemble
$
  "beam" tilde {50\%: ket(Psi_+)",  " 50\%: ket(Psi_-)},
$
which cannot in general be represented by a single state $ket(Psi_"SG")$ @sakurai2020modern. One may think the beam is described by the state
$
  ket(Psi_"SG") = 1/sqrt(2) (ket(Psi_+) + ket(Psi_-)).
$
However, then all particles in the beam would be prepared in the same state. Above half the particles are prepared in the state $ket(Psi_+)$, while the other half is prepared in the state $ket(Psi_-)$.

More generally, consider an ensemble of states $ket(Psi_i)$, where the system is prepared in the state $ket(Psi_i)$ with probability $p_i$. We impose the $p_i$ sum to unity
$
  sum_i p_i = 1,
$<sum-unity>
as probabilities should. The expectation value of an operator $cal(O)$ is then given by @nielsen2010quantum @sakurai2020modern @Preskill_notes
$
  expval(cal(O)) & = sum_i p_i braket(Psi_i, cal(O), Psi_i).
$<mmment>

Now, we define the density matrix @nielsen2010quantum @sakurai2020modern @Preskill_notes
$
  rho equiv sum_i p_i ketbra(Psi_i),
$<density-mat>
which allows us to write the expectation value @mmment as
$
  expval(cal(O)) = tr [rho cal(O)],
$
implying the density matrix contains all measurable information about the system. Also, the normalisation condition @sum-unity becomes
$
  tr rho = 1.
$<normalise>
We call the state pure if the density matrix can be written as
$
  rho = ketbra(Psi),
$
implying $rho^2 = rho$ and $tr rho^2 = 1$. Otherwise, the state is called mixed.

=== Unitary evolution
Consider evolving the system with a unitary operator $U(t)$. The states $ket(Psi_i)$ evolve according to @sakurai2020modern
$
  ket(Psi_i (t)) = U(t) ket(Psi_i (0)),
$
while the probabilities $p_i$ remain unchanged. The density matrix therefore evolves as
$
  rho(t) & = sum_i p_i ketbra(Psi_i (t)) = U(t) rho(0) U^dagger (t).
$<unitary>
Now using $U^dagger (t) U(t) = bb(1)$ and the cyclicity of the trace we obtain from @unitary
$
  tr [rho(t)^2] & = tr [rho(0)^2].
$
Hence, a pure state cannot become mixed under unitary evolution.

== The entangled state
Consider a system with multiple degrees of freedom, such as the string of qubits shown in @qubit-string.

#figure(
  scale(100%, d1),
  caption: [A string of qubits.],
)<qubit-string>

Suppose the total system is in a pure state $ket(Psi)$ with density matrix
$
  rho_(A B) = ketbra(Psi).
$<rho-tot>
We imagine dividing the total system into two subsystems $A$ and $B$. This allows us to decompose the Hilbert space as
$ cal(H)_(A B) = cal(H)_A times.o cal(H)_B. $

Now, consider two observers: Alice and Bob. Alice only has access to $A$, while Bob only has access to $B$. This means Alice can perform measurements on $A$, but has no access to the degrees of freedom in $B$. To Alice the system is therefore fully described by the reduced density matrix @nielsen2010quantum @sakurai2020modern @Preskill_notes
$
  rho_A equiv tr_B [rho_(A B)],
$<reduce>
where $tr_B$ denotes the partial trace over $B$. Applying $tr_B$ removes the degrees of freedom inaccessible to Alice. To see why the definition @reduce is natural consider any operator corresponding to a measurement Alice can perform. Such an operator acts only on states in $cal(H)_A$ and therefore takes the form
$ cal(O) = cal(O)_A times.o bb(1)_B. $
Using @mmment we then obtain @nielsen2010quantum @sakurai2020modern @Preskill_notes
$
  expval(cal(O)) & = tr_(A B) [rho_(A B) ( cal(O)_A times.o bb(1)_B)] \
                 & = tr_A [tr_B (rho_(A B)) cal(O)_A] \
                 & = tr_A [rho_A cal(O)_A].
$
So the result of any measurement Alice can perform is completely determined by the reduced density matrix.

The structure of the reduced density matrix depends on the form of the total state $ket(Psi)$. We call the state $ket(Psi)$ separable if it can be written as a product
$
  ket(Psi) = ket(Psi)_A times.o ket(Psi)_B,
$
implying
$
  rho_(A B) = rho_A times.o rho_B,
$
where both $rho_A$ and $rho_B$ are pure. Conversely, if the total state $ket(Psi)$ is pure but not separable, the reduced density matrices
$ rho_A & = tr_B [rho_(A B)]",  " rho_B & = tr_A [rho_(A B)], $ are mixed. Such states are called entangled.

=== The EPR pair I
The simplest entangled states are the EPR pairs.#footnote[or Bell states.] For a system of two qubits with bases ${ket(0)_A, ket(1)_A}$ and ${ket(0)_B, ket(1)_B}$ an example of such a pair is @nielsen2010quantum @Preskill_notes @sakurai2020modern
$
  ket("EPR") &= 1/sqrt(2) (ket(0)_A times.o ket(1)_B + ket(1)_A times.o ket(0)_B) equiv 1/sqrt(2) (ket(01) + ket(10)).
$<EPR-state>

Now, imagine Alice performs a measurement on her qubit and obtains $0$. The state @EPR-state then collapses as
$
  ket("EPR") ~^"Alice measurement" ket(0 1),
$
implying that Bob's qubit must be in the state $ket(1)$. Likewise, if Alice measures $1$, Bob's qubit must be in the state $ket(0)$. Their measurement outcomes are therefore perfectly anti-correlated. These correlations persist even when Alice and Bob are separated by arbitrarily large distances. This seems highly non-local and was the basis of the Einstein-Podolsky-Rosen (EPR) paradox @EPRpaper or "spooky action at a distance". However, no violation of locality occus, since the correlations only become observable once Alice and Bob compare their measurements.

== Entanglement entropy
We would like a number quantifying the mixedness of a state. This number is the von Neumann entropy @nielsen2010quantum @sakurai2020modern @Preskill_notes
$
  S & equiv - tr [rho ln rho].
$<von-neu>
The density matrix can always be diagonalised#footnote[This follows from $rho^dagger = rho$. ] in which case the von Neumann entropy becomes
$
  S & eq - sum_i rho_(i i) ln rho_(i i),
$<von-neu-diag>
where $rho_(i i)$ are the eigenvalues of the density matrix.

To see the definition @von-neu is sensible, consider:

1. A pure state. In this case $ rho = ketbra(Psi_j), $ so $ rho_(j j) & = 1",  "
              rho_(i i) & = 0, $ for $i eq.not j$. The entropy @von-neu-diag therefore vanishes $ S = - ln 1 = 0. $

2. A completely mixed state. Suppose we have $N$ states $ket(Psi_i)$ occurring with equal probability $p_i = N^(-1)$. Then $ rho = 1/N bb(1)_N, $ or $rho_(i i) = N^(-1)$. The entropy @von-neu-diag is therefore $ S = - sum_(i=1)^N 1/N ln 1/N = ln N. $

We find the von Neumann entropy vanishes for pure states and is positive for mixed states, with the maximal value being $ln N$, as expected for a measure of mixedness.

Now, for a pure state $ket(Psi)$, we define the entanglement entropy as the von Neumann entropy of the reduced density matrix @nielsen2010quantum @Preskill_notes
$
  S_"EE" equiv S_A eq - tr [rho_A ln rho_A].
$<ent-ent>
The entanglement entropy measures how entangled the state $ket(Psi)$ is. To see this, consider taking $ket(Psi)$ to be:

1. A separable state. In this case $rho_A$ is pure, implying $ S_"EE" = 0. $

2. An entangled state. In this case $rho_A$ is mixed, implying $ S_"EE" eq.not 0. $

The difference occurs because we discard information stored in correlations between $A$ and $B$ when tracing out $B$, thereby increasing the entropy $S_"EE" eq.not 0$.

=== The EPR pair II
As an example, consider again the EPR pair. The density matrix of @EPR-state is
$
  rho_(A B) = 1/2 (ketbra(01) + ketbra(01, 10) + ketbra(10, 01) + ketbra(10) ),
$
and taking the partial trace over $B$ gives
$
  rho_A = 1/2 ketbra(0) + 1/2 ketbra(1),
$
which is maximally mixed, implying $S_A = ln 2$.

== Thermal systems
With the essential formalism in place, we now turn to systems in which ensembles show up naturally: thermal systems, as is well known from statistical mechanics @blundell.

Consider a system in thermal equilibrium described by the Hamiltonian $H$. Let $ket(E_i)$ denote a complete set of energy eigenstates satisfying
$
  H ket(E_i) = E_i ket(E_i).
$
With the basis ${ket(E_i)}$, the density matrix @density-mat becomes
$
  rho = sum_i p_i ketbra(E_i),
$
where $p_i$ is the probability that the system occupies the state $ket(E_i)$. For a canonical ensemble at temperature $T$, the probabilities are given by the Boltzmann distribution
$
  p_i prop e^(-beta E_i),
$<boltzmann>
where $beta = T^(-1)$ @blundell. Substituting @boltzmann into the density matrix gives#footnote[A less handwavy justification is given in @schmitt_tft.]
$
  rho & prop sum_i e^(-beta E_i) ketbra(E_i).
$
The $ket(E_i)$ form an eigenbasis of $H$, so we may write
$
  rho prop e^(-beta H).
$
Then using the normalisation condition @normalise we obtain the thermal density matrix
$
  rho = e^(-beta H)/Z(beta)",  " Z(beta) = tr e^(-beta H) = sum_i e^(-beta E_i),
$<thermal-den>
with $Z(beta)$ being the partition function @sakurai2020modern @schmitt_tft.

The von Neumann entropy of @thermal-den is
$
  S = beta expval(H) + ln Z,
$<neu-therm>
which coincides with the thermodynamic entropy @blundell. Now, consider:

1. The high-temperature limit $beta -> 0$. For a $N$-dimensional Hilbert space $ rho -> 1/N bb(1)_N, $ implying $ S -> ln N. $

2. The low-temperature limit $beta -> oo$. Assuming a non-degenerate ground state $ rho -> ketbra(E_0), $ implying $ S -> 0. $

Hence, the thermal density matrix is mixed for all $T > 0$, becoming pure as $T -> 0$.

=== The thermofield double
Now, since the thermal density matrix is mixed, we should be able to obtain it from a pure state living in a larger Hilbert space. The choice of purification is not unique, but a convenient one is found by introducing a second copy of the system. This leads to the thermofield double state @Cottrell_2019 @witten2025introductionblackholethermodynamics @Harlow_2016 @hartman_bhi_notes
$
  ket("TFD") = 1/sqrt(Z(beta)) sum_i e^(-beta E_i\/2) ket(E_i)_A times.o ket(E_i)_B.
$<TFD>
This state purifies the thermal density matrix. We can see this by tracing over $B$ giving
$
  rho_A &= tr_B [ketbra("TFD")] \
  &= 1/Z(beta) sum_i sum_j e^(-beta (E_i+E_j)\/2) ""_A ket(E_i) bra(E_j)_A underbracket(sum_k ""_B braket(E_k, E_i)_B ""_B braket(E_j, E_k)_B, = delta_(i j)) \
  &= 1/Z(beta) sum_i e^(-beta E_i) ""_A ket(E_i) bra(E_i)_A,
$
or
$
  rho_A & = e^(-beta H)/Z(beta),
$
which is the thermal density matrix.

Physically, the thermofield double describes an entangled pure state between two identical copies of a system. Since the reduced density matrix of either subsystem is thermal, an observer with access to only one copy sees a thermal state. We will find this state shows up everywhere.

== A detour
We take a brief detour to introduce a  method for determining the ground state of a system using path integrals. Later, we apply this method to the Minkowski vacuum and black holes. The starting point is the time-evolution operator#footnote[assuming $H$ is time-independent.] $ U(T) = e^(-i H T), $ whose transition amplitudes can be written as a path integral @peskin1995introduction @hartman_bhi_notes @witten2025introductionblackholethermodynamics
$
  braket(phi.alt_2, e^(-i H T), phi.alt_1) = integral_(phi.alt (t = 0) = phi.alt_1)^(phi.alt (t = T) = phi.alt_2) dd(phi.alt, d: D) e^(-S [phi.alt]),
$
with $S$ being the Lorentzian action. The Euclidean path integral is given by a Wick rotation $t = -i t_E$
$
  braket(phi.alt_2, e^(- beta H), phi.alt_1) = integral_(phi.alt (t_E = 0) = phi.alt_1)^(phi.alt (t_E = beta) = phi.alt_2) dd(phi.alt, d: D) e^(-S_E [phi.alt]),
$<amplitude>
with $S_E$ being the Euclidean action.

We think of $phi.alt_(1,2)$ as boundary conditions. The meaning of this depends on the space we integrate over. As an example for $RR^d$ we would integrate over the strip $RR^(d-1) times I_beta$ as shown for $d = 2$ in @geometry-1 with $phi.alt_(1,2)$ being as indicated. @hartman_bhi_notes

#figure(
  scale(100%, d2),
  caption: [Example path integral for $RR^2$. @hartman_bhi_notes ],
)<geometry-1>

While if space were a sphere $S^d$ we would integrate over the cylinder $S^(d-1) times I_beta$.

So the path integral requires boundary "data" on the cuts $t_E = 0$ and $t_E = beta$. We can think of the path integral with the cut $t_E = beta$ open as the state
$
  ket(Psi) = e^(-beta H) ket(phi.alt_1),
$
where we integrate over the geometry in @geometry-2. @hartman_bhi_notes
#figure(
  scale(100%, d3),
  caption: [Example path integral for $RR^2$ with an open cut. @hartman_bhi_notes ],
)<geometry-2>

The path integral is then a wavefunctional of the boundary field $Phi[phi.alt_2] = braket(phi.alt_2, Psi)$. This idea can be used to determine the state $ket(Psi)$, and in general any path integral with an open cut $Sigma$ defines a state on $Sigma$. @hartman_bhi_notes

=== Determining the ground state
Now, we can define the ground state of a system by a path integral extending infinitely in $t_E$ as in @geometry-3. @hartman_bhi_notes
#figure(
  scale(100%, d4),
  caption: [Example path integral to prepare the ground state. @hartman_bhi_notes ],
)<geometry-3>

To see this we expand $ket(Psi)$ in the energy eigenstates of the Hamiltonian denoted by $ket(E_i)$
$
  ket(Psi) = sum_i c_i ket(E_i).
$
We evolve in $t_E$ to obtain
$
  ket(Psi(t_E)) & = sum_i c_i e^(-E_i t_E) ket(E_i) \
                & = e^(-E_0 t_E) sum_i c_i e^(-(E_i - E_0) t_E) ket(E_i).
$
Now, taking $t_E -> oo$ we find#footnote[assuming $ket(E_0)$ is non-degenerate.]
$
  ket(Psi(t_E)) prop ket(E_0),
$
since $E_i-E_0 > 0$. So taking $t_E -> oo$ projects onto the ground state.

== The entangled vacuum
Consider a quantum field $phi.alt(x)$ in its vacuum state $ket(Omega)$ and divide space into two regions $A$ and $B$. The two-point correlator is given by @Harlow_2016 @peskin1995introduction
$
  braket(Omega, phi.alt(x) phi.alt(y), Omega) eq^"scalar field" 1/(4 pi^2) m/abs(x-y) K_1 (m abs(x-y)),
$<2pt>
where $K_1$ is the modified Bessel function of the second kind. Take $abs(x-y) << m^(-1)$ then
$
  K_1 (m abs(x-y)) tilde 1/(m abs(x-y)),
$
so @2pt becomes
$
  braket(Omega, phi.alt(x) phi.alt(y), Omega) tilde 1/abs(x-y)^2,
$
which diverges as $abs(x-y) -> 0$. Then the correlation between points $x$ and $y$ becomes formally infinite as we bring them together.#footnote[We have assumed $phi.alt$ is a scalar field, however, this behaviour is seen more generally.]
#figure(
  scale(entanglementsurf, 130%),
  caption: [Example divison of space into regions $A$ and $B$.],
)<entanglementsurface>

Now, take $x in A$ and $y in B$ as in @entanglementsurface. Then @2pt implies the regions $A$ and $B$ are entangled. To see this assume the vacuum is not entangled
$
  ket(Omega) = ket(Omega)_A times.o ket(Omega)_B.
$
Then
$
  braket(Omega, phi.alt (x) phi.alt(y), Omega) &= ""_A braket(Omega, phi.alt(x), Omega)_A ""_B braket(Omega, phi.alt(y), Omega)_B,
$
and since $braket(Omega, phi.alt(x), Omega) = 0$ we find
$
  braket(Omega, phi.alt (x) phi.alt(y), Omega) = 0 tilde "contradiction",
$
implying the vacuum is entangled.

== The Rindler decomposition
To better understand the entanglement structure of the Minkowski vacuum we can determine the ground state $ket(Psi)_"M"$. Computing $ket(Psi)_"M"$ is also a nice warm-up before diving into black holes where we will do something completely analogous. We do this by considering the Rindler decomposition.

The idea is to pick a spatial coordinate, say $x$, and split the Hilbert space as @Harlow_2016 @carroll2004spacetime @witten2025introductionblackholethermodynamics
$
  cal(H) = cal(H)_L times.o cal(H)_R,
$
with $cal(H)_R$ being acted on by fields with $x > 0$ and $cal(H)_L$ being acted on by fields with $x < 0$. We would like a simple basis for $cal(H)_L$ and $cal(H)_R$. To this end we consider the generator of Lorentz boosts. The metric in Minkowski spacetime is#footnote[We supress $dd(y^2) + dd(z^2)$.]
$
  dd(s^2) = - dd(t^2) + dd(x^2),
$<mink2d>
and the boost generator $K_x$ is given by @witten2025introductionblackholethermodynamics @peskin1995introduction
$
  K_x = x partial_t + t partial_x.
$
We want coordinates adapted to the flow of $K_x$ such that using these coordinates $K_x$ becomes the generator of time translations.

These coordinates are the Rindler coordinates $(eta, xi)$ defined for $x > abs(t)$ by @carroll2004spacetime
$
  t = xi sinh eta",  " x = xi cosh eta,
$<rindlercoord>
with $xi > 0$. We call $x > abs(t)$ the right Rindler wedge, which is shown in @rindlerwedgefig.#footnote[An analogous definition works for the left Rindler wedge $x < -abs(t)$.]

#figure(
  scale(rindlerwedge, 90%),
  caption: [The right Rindler wedge in Minkowski coordinates $(x,t)$. The dashed lines $(x = plus.minus t)$ are the Rindler horizons and the red hyperbola $(x^2-t^2 = xi^2)$ is representative of the trajectory for a Rindler observer.],
)<rindlerwedgefig>



With @rindlercoord the metric @mink2d becomes
$
  dd(s^2) = - xi^2 dd(eta^2) + dd(xi^2),
$<Rindler-metric>
which is the Rindler metric.  The surfaces $x = plus.minus t$ define the Rindler horizons and correspond to $xi = 0$ in the Rindler coordinates. A uniformly accelerated observer follows a trajectory of constant $xi > 0$, satisfying#footnote[One can show $xi = alpha^(-1)$ with $alpha$ being the proper acceleration, see e.g. @carroll2004spacetime.] @carroll2004spacetime
$
  x^2 - t^2 = xi^2,
$
which corresponds to a hyperbola in Minkowski spacetime. Such observers remain in the right Rindler wedge and can never receive signals from behind the Rindler horizons.


However, as opposed to other horizons these are a property of the observer and not the spacetime. This is obvious since the spacetime is Minkowski in disguise. We can think of these horizons arising due to the observer "accelerating away" from the spacetime behind them. @carroll2004spacetime

With the Rindler coordinates @rindlercoord the boost generator becomes
$
  partial_eta & = x partial_t + t partial_x = K_x.
$
so $K_x$ acts as the Hamiltonian generating translations in the Rindler time $eta$.

Now, to compute the path integral and determine $ket(Psi)_"M"$ we need an appropriate Euclidean geometry. This is obtained by a Wick rotation $eta -> -i eta_E$, after which the Rindler coordinates @rindlercoord become
$
  t_E = xi sin eta_E",  " x = xi cos eta_E,
$<rindlerrotate>
and the metric @Rindler-metric becomes
$
  dd(s^2) = xi^2 dd(eta_E^2) + dd(xi^2),
$<Rindler-euclid>
which is the polar metric
$
  dd(s^2) = rho^2 dd(theta^2) + dd(rho^2),
$<polar-flat>
with the identification $rho = xi$ and $theta = eta_E$.

Since spacetime is still Minkowski nothing weird should happen, so we require the metric @Rindler-euclid is regular or "smooth" at the horizon $xi = 0$. This implies $eta_E$ is periodic
$
  eta_E tilde eta_E + 2 pi.
$
The Minkowski wavefunctional on the slice $t = 0$ is now defined by evaluating the path integral over $t_E <= 0$ @Harlow_2016. From @rindlerrotate we see that $K_x$ generates rotations in the $x t_E$-plane. Euclidean time evolution generated by $K_x$ is then given by $e^(-eta_E K_x)$, where $eta_E$ is an angle in the $x t_E$-plane. The region $t_E <= 0$ therefore corresponds to evolution by an angle $pi$, as shown in @rindlergeometry.

#figure(
  scale(rindlergeometry, 100%),
  caption: [Geometry for defining the Minkowski wavefunctional using the Rindler decomposition. We use the lower-half plane as an initial condition.],
)<rindlergeometry>

With the above reasoning the path integral @amplitude can be written as#footnote[We have $phi.alt_1 = phi.alt (eta_E = 0) = phi.alt_R$ and $phi.alt_2 = phi.alt (eta_E = pi) = phi.alt_L$.] @Harlow_2016 @witten2025introductionblackholethermodynamics @hartman_bhi_notes
$
  braket(phi.alt_L phi.alt_R, Psi) prop braket(phi.alt_R, e^(- pi K_R) Theta, phi.alt_L),
$
where $K_R$ is the restriction of $K_x$ to the right Rindler wedge and $Theta$ is the anti-unitary CPT operator. The role of $Theta$ is to map states in $cal(H)_L$ to corresponding states in $cal(H)_R$ by $Theta ket(phi.alt_L) in cal(H)_R$, allowing the boost operator $K_R$ to act on them. Now, we can use completeness of the $K_R$ eigenstates $ket(i)_R$ to write @Harlow_2016 @witten2025introductionblackholethermodynamics
$
  braket(phi.alt_L phi.alt_R, Psi) &prop sum_i bra(phi.alt_R) e^(-pi K_R) ket(i)_R bra(i)_R Theta ket(phi.alt_L) \
  & prop sum_i e^(-pi E_i) bra(i)_R Theta ket(phi.alt_L) braket(phi.alt_R, i)_R \
  & prop sum_i e^(-pi E_i) braket(phi.alt_L, i^*)_L braket(phi.alt_R, i)_R,
$
where $ket(i^*)_L = Theta^dagger ket(i)_R$. Then by inspection the Minkowski vacuum is
$
  ket(Psi)_M = 1/sqrt(Z) sum_i e^(- pi E_i) ket(i^*)_L times.o ket(i)_R,
$<mm-vacuum>

which we recognise as the thermofield double state @TFD with $beta = 2 pi$. The entanglement structure is now completely manifest, with the left and right Rindler wedges clearly being entangled.

Note, the reduced density matrix of @mm-vacuum is
$
  rho_R = 1/Z sum_i e^(-2 pi E_i) ""_R ket(i) bra(i)_R,
$
which is thermal with the Unruh temperature#footnote[With the omission of a length scale.] @Unruh_1984
$
  T_"Unruh" = 1/(2 pi).
$
This is the Unruh effect where accelerating observers experience a temperature. Which we now see is related to the left Rindler wedge being inaccessible to these observers.

#pagebreak()
= Black Holes
Black holes show up as solutions to the Einstein field equations and provide the setting for Einstein-Rosen bridges. This chapter introduces classical- and quantum black holes in Minkowski spacetime, including the Hartle-Hawking state and black hole thermodynamics. We also discuss black holes in anti-de Sitter spacetime, where an $"ER"="EPR"$ like connection appears.

== The Schwarzschild metric
Consider the source-free Einstein field equations
$
  R_(mu nu) = 0.
$<EFE-vacuum>
By the Birkhoff-Jebsen theorem, the unique spherically symmetric solution of @EFE-vacuum that asymptotically approaches Minkowski spacetime is the Schwarzschild metric. The Schwarzschild metric is given by @carroll2004spacetime @harlow2023blackholesquantumgravity @Harlow_2016 @tonggr
$
  dd(s^2) = - f(r) dd(t^2) + f(r)^(-1) dd(r^2) + r^2 dd(Omega^2),
$<SCH-metric>
where
$
  f(r) equiv 1-r_s/r = (r-r_s)/r,
$<sch>
with $r_s equiv 2 G M$ being the Schwarzschild radius. The metric @SCH-metric describes the spacetime outside any spherically symmetric mass distribution of total mass $M$ centered at $r=0$.

We see that the metric components diverge at $r = 0$ and $r = r_s$. The point $r = 0$ is a true singularity, since curvature diverges as $r -> 0$. The divergence at the Schwarzschild radius, however, is caused by a bad choice of coordinates, and nothing special happens locally as $r -> r_s$.#footnote[at least classically.] We call the Schwarzschild radius "the horizon".

The horizon is significant due to the sign of $f(r)$ flipping upon crossing $r_s$. Then continuing toward smaller $r$ corresponds to moving forward in time. Any object entering $r < r_s$ is therefore destined to reach the singularity $(r=0)$ since doing otherwise corresponds to travelling backwards in time. @harlow2023blackholesquantumgravity @Harlow_2016 @carroll2004spacetime Thus crossing the horizon is a one-way trip. We can define black holes to be objects whose physical size is smaller than their corresponding Schwarzschild radius.

== The fully-extended geometry
To better understand the geometry near the horizon, we seek coordinates in which the metric @SCH-metric is non-singular at $r = r_s$. We first define the tortoise coordinate $r_*$ by @carroll2004spacetime @tonggr
$
  dd(r_*) = dd(r)/f(r),
$
which can be integrated giving
$
  r_* = r + r_s log ((r-r_s)/r_s).
$
Now, we define the Kruskal-Szerkeres coordinates by @carroll2004spacetime @tonggr

$
  U & eq - exp((r_*-t)/(2 r_s))",  " V & eq exp((r_*+t)/(2 r_s)),
$<Kruskal-coord>
which satisfy
$
  U V = (r_s-r)/r_s e^(r\/r_s).
$
With @Kruskal-coord the metric @SCH-metric becomes
$
  dd(s^2) = - (2r_s^3)/r e^(-r\/r_s) (dd(U, V) + dd(V, U)) +r^2 dd(Omega^2),
$<Kruskal-1>
where $r$ is implicitly defined by $r eq r(U,V)$. We can make @Kruskal-1 diagonal by defining
$
  T & eq (V+U)/2",  " X & eq (V-U)/2.
$
Then @Kruskal-1 takes the form
$
  dd(s^2) = (4 r_s^3)/r e^(-r\/r_s) (-dd(T^2) + dd(X^2)) + r^2 dd(Omega^2).
$<Kruskal-True>
We see @Kruskal-True is non-singular at the horizon, while also being conformally Minkowksi since#footnote[We supress $dd(Omega^2)$.]
$
  dd(s^2) prop -dd(T^2) + dd(X^2),
$
implying null radial geodesics satisfy
$
  dd(T) = plus.minus dd(X).
$
The fully-extended geometry described by @Kruskal-True is shown in @kruskal-chart

#figure(
  scale(100%, kruskal),
  caption: [The Kruskal chart. @carroll2004spacetime @tonggr The dashed lines correspond to the horizon where $U = 0$ or $V = 0$, implying $T = plus.minus X$, while the jagged lines correspond to the singularity.],
)<kruskal-chart>

== Penrose diagrams and Einstein-Rosen bridges
We are often interested in the causal structure of spacetime. This motivates the construction of Penrose diagrams, which compactify spacetime by mapping infinities to finite distances in a conformal way. Starting with @Kruskal-1 we can make the coordinates finite by applying the transformation @carroll2004spacetime
$
  tilde(U) = arctan U",  " tilde(V) = arctan V,
$
where now $tilde(U),tilde(V) in (-pi\/2, pi\/2)$. Then @Kruskal-1 becomes
$
  dd(s^2) & = Omega^(-2) (tilde(U), tilde(V)) dd(tilde(U), tilde(V)),
$<conformal-kruskal>
where
$
  Omega^(-2) (tilde(U),tilde(V)) &eq (4 r_s^3)/r e^(-r\/r_s)/(cos^2 tilde(U) cos^2 tilde(V)).
$
We again define
$
  tilde(T) & eq (tilde(V)+tilde(U))/2",  " tilde(X) & eq (tilde(V)-tilde(U))/2.
$
Then @conformal-kruskal takes the form
$
  dd(s^2) prop -dd(tilde(T)^2) + dd(tilde(X)^2),
$
which has the same causal structure as @Kruskal-True.

#figure(
  scale(100%, double-pen),
  caption: [Penrose diagram for the fully-extended Schwarzschild geometry. The dashed lines correspond to the horizons $(r=r_s)$ and the jagged lines correspond to the singularities $(r=0)$. $cal(i)^minus.plus$ are the past and future timelike infinities, $cal(i)^0$ is spatial infinity, and $cal(I)^minus.plus$ are the past and future null infinities. @Harlow_2016 @carroll2004spacetime @tonggr],
) <two-sided-minkowski>

The Penrose diagram for the fully-extended Schwarzschild geometry is shown in @two-sided-minkowski. @tonggr @Harlow_2016 @carroll2004spacetime The spacetime contains two asymptotically Minkowski exterior regions, $L$ and $R$, where $r -> oo$, together with past and future singularities at $r = 0$. The two exteriors appear connected through an Einstein-Rosen bridge or wormhole @ERbridge. However, the diagram makes clear that $L$ and $R$ are causally disconnected. An observer beginning in either $L$ or $R$ may cross the horizon and enter the black hole interior $I$. We know this is a one-way journey towards the singularity and the diagram makes this obvious.

Now, consider two observers, one originating in $L$ and the other in $R$, each crossing their respective horizons. Although the exterior regions are causally disconnected, the observers may nevertheless meet inside the interior $I$ before reaching the singularity. So the Einstein-Rosen bridge is non-traversable, yet the future interiors of $L$ and $R$ overlap causally.

$L$ and $R$ are often interpreted as distinct asymptotically Minkowski universes connected through the Einstein-Rosen bridge. Alternatively, one may attempt to interpret them as widely separated regions of the same universe, corresponding heuristically to two black holes separated by a large distance. With this interpretation the bridge appears highly non-local. Observers beginning arbitrarily far apart in the exterior geometry may encounter each other inside $I$. However, this interpretation is only approximate due to the mutual attraction between the black holes. @Maldacena_2013

Note, the fully-extended Schwarzschild geometry described above is highly idealised and does not correspond to real black holes formed by collapse. The Penrose diagram for such a black hole can be constructed by combining the Penrose diagrams for Minkowski space (see @a1) and the extended Schwarzschild geometry as shown in @b1. The resulting spacetime contains only a single interior region and does not include an Einstein-Rosen bridge or second asymptotic region. @Harlow_2016 @carroll2004spacetime

#subpar.grid(
  figure(
    scale(mink-pen, 100%),
    caption: [Penrose diagram for Minkowski space.],
  ),
  <a1>,

  figure(
    scale(pen-bh, 80%),
    caption: [Penrose diagram for classical black hole formation by collapse of a photon cloud. @Harlow_2016],
  ),
  <b1>,

  columns: (1fr, 1fr),
  caption: [Penrose diagrams.],
  label: <pen-black-hole>,
)

== The Hartle-Hawking state
We would like to determine the ground state for the fully-extended Schwarzschild geometry. This is done by computing the path integral as when we determined $ket(Psi)_"M"$.

Consider the Schwarzschild metric @SCH-metric
$
  dd(s^2) = - (r-r_s)/r dd(t^2) + r/(r-r_s) dd(r^2) + r^2 dd(Omega^2).
$
The appropriate Euclidean geometry is obtained after a Wick rotation $t = -i t_E$
$
  dd(s^2) = (r-r_s)/r dd(t_E^2) + r/(r-r_s) dd(r^2) + r^2 dd(Omega^2).
$<euclidean>
We now write
$
  r = r_s + delta,
$
with $0 < delta << r_s$. This is called the near-horizon limit and will prove useful @tonggr @Harlow_2016 @witten2025introductionblackholethermodynamics. Now, using $f(r_s) = 0$ we can expand $f(r)$ around $r_s$ as
$
  f(r) &tilde.eq^"Taylor" underbracket(f(r_s), 0) + dv(f(r_s), r) delta &tilde.eq delta/r_s.
$
Also, $r^2 tilde.eq r_s^2$ so @euclidean becomes#footnote[We suppress $dd(Omega^2)$.]
$
  dd(s^2) tilde.eq delta/r_s dd(t_E^2) + r_s/delta dd(delta^2).
$<eucl2>

We can define $rho^2 equiv 4 r_s delta$ to write @eucl2 as
$
  dd(s^2) tilde.eq (rho/(2 r_s))^2 dd(t_E^2) + dd(rho)^2.
$<rindler>
@rindler is the Rindler metric @Rindler-euclid so we say the horizon is locally Rindler. Like in the Minkowski case the metric @rindler becomes the polar metric @polar-flat upon defining $theta$ by
$
  t_E equiv 2 r_s theta.
$
We call the geometry described by @rindler the Euclidean cigar @Harlow_2016 @witten2025introductionblackholethermodynamics @hartman_bhi_notes, which is shown in @euclidcigarfig.

#figure(
  scale(image("fig/euclidcigar.png"), 25%),
  caption: [The Euclidean cigar. Adapted from @Harlow_2016.],
)<euclidcigarfig>

Now, by the equivalence principle nothing special happens at the horizon. We therefore require the Euclidean cigar is smooth at $rho = 0$. This implies $theta tilde theta + 2 pi$ or equivalently
$
  t_E tilde t_E + beta,
$
with $beta equiv 4 r_s pi$.

The Hartle-Hawking wavefunctional @Hartle:1976tp @ISRAEL1976107 on the slice $t = 0$ is then defined by evaluating the path integral over the lower-half of the Euclidean cigar. This region corresponds to the interval
$
  -beta/2 < t_E < 0,
$
and the construction is shown in @hh-geometry. Euclidean time evolution generated by $H$ is given by $e^(-t_E H)$, so evolution over the interval above is given by $e^(-beta H\/2)$.
#figure(
  scale(100%, hh-state-r),
  caption: [Geometry for defining the Hartle-Hawking wavefunctional. The upper-half is just the upper-half of the usual Penrose diagram. Adapted from @Harlow_2016.],
) <hh-geometry>

With the above reasoning the path integral @amplitude can be written as#footnote[We have $phi.alt_1 = phi.alt (t_E = 0) = phi.alt_R$ and $phi.alt_2 = phi.alt (t_E = beta/2) = phi.alt_L$.] @Harlow_2016 @hartman_bhi_notes @witten2025introductionblackholethermodynamics
$
  braket(phi.alt_L phi.alt_R, Psi) prop braket(phi.alt_R, e^(-beta H_R\/2) Theta, phi.alt_L),
$
where $H_R$ is the restriction of the Hamiltonian to $R$ and $Theta$ is the CPT-operator. As with the Minkowski vacuum we use completeness to find @Harlow_2016
$
  braket(phi.alt_L phi.alt_R, Psi) &prop sum_i e^(-beta E_i\/2) braket(phi.alt_L, i^*)_L braket(phi.alt_R, i)_R,
$
where $ket(i^*)_L = Theta^dagger ket(i)_R$. We can then read off the Hartle-Hawking state
$
  ket(Psi)_"HH" = 1/sqrt(Z(beta)) sum_i e^(-beta E_i\/2) ket(i^*)_L times.o ket(i)_R,
$<HH>
which we recognise as the thermofield double state @TFD. Therefore the regions $L$ and $R$ are quantum mechanically entangled, while geometrically they are connected by an Einstein-Rosen bridge. The fully-extended Schwarzschild geometry can, as noted, be viewed as corresponding to two black holes separated by a large distance. Then the Hartle-Hawking state tells us these will be entangled#footnote[in a heuristic manner.].

We briefly note the Hartle-Hawking state was derived using two widely accepted assumptions:

1. Validity of the equivalence principle, so the horizon is smooth.

2. Validity of semiclassical effective field theory, i.e. quantum fields propagating on a classical background.

With the second being implicit in our use of the path integral.

== Black holes as thermal systems
The reduced density matrix of the Hartle-Hawking state is thermal with the Hawking temperature @1975CMaPh..43..199H @ISRAEL1976107 @Harlow_2016
$
  T_H = beta^(-1) = 1/(4 pi r_s).
$<hawkingTemp>
This means an observer confined to either region $L$ or $R$ sees a thermal state with temperature $T_H$ because $L$ and $R$ are entangled.

The Hartle-Hawking state describes a black hole in equilibrium with its own radiation. @Harlow_2016 This equilibrium arises because we consider a spacetime which has contained a black hole for an infinite amount of time.#footnote[which is clearly non-physical.] The cause of the thermal bath can then only be the black hole with the system having plenty of time to reach equilibrium. This implies the black hole has temperature $T_H$ and emits radiation. We could now imagine placing or forming a black hole in a cold universe. This black hole still has temperature $T_H$ and will radiate. However, this black hole would not be in equilibrium and will begin evaporating.

Since a black hole has a temperature $T_"H"$ and an energy which we can identify as $E = M$ it will have some entropy. We have by definition @blundell
$
  dv(S, E) = 1/T_"H".
$
Now, using $S(E=0) = 0$ we find
$
  S = A/(4G) = 2 pi A/cal(l)_p^2,
$<bekenstein>
with $A = 4 pi r_s^2$ being the area of the horizon and $cal(l)_p = sqrt(8 pi G)$ being the Planck length.#footnote[The factor $sqrt(8 pi)$ is a matter of convention.] The quantity @bekenstein is referred to as the Bekenstein-Hawking entropy. @Bekenstein_1973 @1975CMaPh..43..199H

== in anti-de Sitter
The Schwarzschild metric @SCH-metric asymptotically approaches Minkowski spacetime which has a vanishing cosmological constant. We now want a metric that asymptotically approaches anti-de Sitter spacetime. The metric for empty#footnote[Containing only a negative cosmological constant.] anti-de Sitter spacetime has the same static, spherically symmetric form as @SCH-metric
$
  dd(s^2) = -f(r) dd(t^2) + f(r)^(-1) dd(r^2) + r^2 dd(Omega^2),
$<ads>
with @tonggr
$
  f(r) = 1+r^2/r_"ads"^2",  " r_"ads"^2 equiv -3/Lambda.
$<adsf>

where $Lambda < 0$. We combine @sch and @adsf to obtain the metric describing a Schwarzschild black hole in anti-de Sitter spacetime by using
$
  f = 1 - r_s/r + r^2/r_"ads"^2",  " r_"ads"^2 = -3/Lambda.
$<bh-ads>
The form of @bh-ads can be justified by considering:

1. For $r << r_"ads"$, the anti-de Sitter term becomes negligible and $ f tilde.eq 1-r_s/r, $ recovering the Schwarzschild metric @sch.

2. For $r -> oo$, $ f tilde.eq 1+r^2/r_"ads"^2, $ recovering the anti-de Sitter metric @adsf.

So @bh-ads reproduces the Schwarzschild geometry at small $r$ and approaches anti-de Sitter spacetime as $r -> oo$. This is exactly the behavior we wanted. The Penrose diagram describing @bh-ads is shown in @pen-AdS. @Harlow_2016
#figure(
  scale(100%, pen-ads),
  caption: [Penrose diagram for an eternal black hole in anti-de Sitter spacetime. Here, $L$ and $R$ are two timelike boundaries. @Harlow_2016 @Maldacena_2003],
) <pen-AdS>

The state describing this geometry is again the Hartle-Hawking state @HH @Maldacena_2003 @Harlow_2016.

We can determine the temperature of the black hole by defining a horizon radius $r_+$ with $f(r_+) = 0$ and expanding $f(r)$ as#footnote[Now the Schwarzschild radius $r_s$ is just a constant related to the mass $M$.]
$
  f(r) & tilde.eq^"Taylor" underbracket(f(r_+), 0) + evaluated(dv(f (r), r))_(r=r_+) delta,
$<fTaylor>
with $delta = r-r_+$ and
$
  f(r_+) = 1 - r_s/r_+ + r_+^2/r_"ads"^2 =^! 0.
$<fconstraint>
Then @fTaylor becomes
$
  f(r) & tilde.eq (r_s/r_+^2 + (2r_+)/r_"ads"^2) delta eq^(#[by @fconstraint]) ((r_"ads"^2 + 3 r_+^2)/(r_+ r_"ads"^2)) delta.
$
The Euclidean metric after a Wick rotation $t = - i t_E$ is  then
$
  dd(s^2) = ((r_"ads"^2 + 3 r_+^2)/(r_+ r_"ads"^2)) delta dd(t_E^2) + ((r_+ r_"ads"^2)/(r_"ads"^2 + 3 r_+^2)) dd(delta^2)/delta.
$<AdSBH>
We can define
$
  rho^2 equiv (4 r_+ r_"ads"^2)/(r_"ads"^2+3 r_+^2) delta,
$
and write @AdSBH as
$
  dd(s^2) = ((r_"ads"^2 + 3 r_+^2)/(2 r_+ r_"ads"^2))^2 rho^2 dd(t_E^2) + dd(rho^2).
$<AdSBH2>
@AdSBH2 is again the Rindler metric @Rindler-euclid and becomes the polar metric @polar-flat if we define $theta$ by
$
  t_E equiv (2 r_+ r_"ads"^2)/(r_"ads"^2 + 3 r_+^2) theta.
$
The horizon is smooth when $theta tilde theta + 2 pi$ implying
$
  t_E tilde t_E + beta,
$
with
$
  beta & equiv 4 pi (r_+ r_"ads"^2)/(r_"ads"^2 + 3 r_+^2).
$
The temperature of a Schwarzschild black hole in anti-de Sitter spacetime is then @witten2025introductionblackholethermodynamics
$
  T_H = beta^(-1) = 1/(4 pi) (r_"ads"^2 + 3 r_+^2)/(r_+ r_"ads"^2) = 1/(4 pi) (1/r_+ + (3 r_+)/r_"ads"^2),
$<AdStemp>
with $r_+$ defined implicitly by $f(r_+) = 0$.

Now, consider:

1. When $r_+ << r_"ads"$ we recover $r_+ tilde.eq r_s$ from @fconstraint and @AdStemp becomes$ T_H tilde.eq 1/(4 pi r_s), $ as we had before @hawkingTemp.

2. When $r_+ >> r_"ads"$ we have from @fconstraint $ r_+ - r_s + r_+^3/r_"ads"^2 =^! 0 => r_+ tilde.eq (r_s r_"ads"^2)^(1\/3), $ so @AdStemp becomes $ T_H tilde.eq 1/(4 pi) (3 r_+)/r_"ads"^2 = 3/(4 pi) (r_s/r_"ads"^4)^(1\/3), $<AdStempBH> which increases linearly with the size of the black hole.

These are referred to as small and large black holes respectively. @AdStempBH implies large black holes are stable and behave normally#footnote[i.e they have positive heat capacities.], whereas small black holes are unstable. One interesting feature of anti-de Sitter spacetime is the two timelike boundaries $L$ and $R$ in @pen-AdS. Radiation can reach these boundaries in a finite time before being reflected and returning, so the anti-de Sitter boundary effectively acts like a reflecting box. Large black holes can therefore reach stable equilibrium with their own radiation, while small black holes remain unstable and continue to evaporate. Then eternal black holes in anti-de Sitter spacetime are "meaningful", as opposed to eternal black holes in Minkowski spacetime, which are unphysical. @Maldacena_2003 @Harlow_2016

== The AdS/CFT correspondence<AdSCFTcorr>
The usefulness of black holes in anti-de Sitter spacetime may seem unfounded, however, we care because of the AdS/CFT correspondence due to Maldacena @Maldacena_1999.

The correspondence essentially states that a bulk theory of quantum gravity in anti-de Sitter spacetime is equivalent to a boundary conformal field theory. @Ramallo_2013 @Harlow_2016 This idea is illustrated in @AdSCFTcyl.

#figure(
  scale(AdSCFT, 80%),
  caption: [The AdS/CFT correspondence. A theory of gravity in the bulk AdS is equivalent to a CFT on the boundary. Adapted from @harlow2023blackholesquantumgravity.],
)<AdSCFTcyl>

Since CFTs are generally better understood than theories of quantum gravity, this duality provides a powerful framework for studying quantum gravity. Any substantial introduction to the correspondence is way beyond the scope of this project, but its key feature is that it provides an exact definition of the bulk theory in terms of the boundary CFT.

Previously, we interpreted the Hartle-Hawking state @HH as describing entanglement between the two bulk exteriors $L$ and $R$, which are connected geometrically by an Einstein-Rosen bridge. The AdS/CFT correspondence refines this picture for black holes in anti-de Sitter spacetime. Now each asymptotic region is dual to a separate boundary CFT, and the Hartle-Hawking state can more fundamentally be understood as an entangled state of these two CFTs.  @Maldacena_2003 @hartman_bhi_notes @Harlow_2016

The Hartle-Hawking state then corresponds to a thermofield double @TFD of two independent CFTs which is dual to a single connected spacetime @Maldacena_2003 @hartman_bhi_notes. The Einstein-Rosen bridge can thus be interpreted as the geometric realisation of the entanglement in the boundary theory. We say the spacetime is "emergent", and symbolically we can write:
$
  "ER-bridge" <=>^"AdS/CFT" "TFD entanglement",
$
which is the most important realisation of an $"ER"="EPR"$-like equivalence. This can be considered the canonical example of $"ER"="EPR"$ @Maldacena_2013.

#pagebreak()
= The Information Paradox
Black holes being thermal with the temperature $T_H$ leads to a tension between two otherwise sacred principles: unitarity and the equivalence principle. This chapter introduces this tension through the black hole information paradox and its entanglement-based formulation in the AMPS paradox.

== Black hole evaporation
Black holes formed by collapse emit radiation and slowly lose mass. This process reduces the horizon area, so the black hole shrinks over time. The Penrose diagram for a black hole formed by collapse then becomes as shown in @pen-bh-evap.

#figure(
  scale(80%, pen-evap),
  caption: [Penrose diagram for an evaporating black hole. @Harlow_2016],
)<pen-bh-evap>

This evaporation was shown to be problematic by Hawking @HawkingBHinfo @Mathur_2009 @Harlow_2016 @Almheiri_2021 @Guo_2021 @raju2021lessonsinformationparadox. To see why, consider a universe containing a black hole formed by collapsing a cloud of photons. When the black hole is formed and before it has had time to emit radiation the total state
$ ket(Psi) = ket(Psi)_"BH", $
is pure. During evaporation, we assume the Hilbert space may be decomposed as
$
  cal(H)_"total" = cal(H)_"BH" times.o cal(H)_"R",
$
where $cal(H)_"BH"$ describes the interior degrees of freedom in the black hole and $cal(H)_"R"$ describes the degrees of freedom in the emitted radiation. The total state is still pure
$ rho_"total" = ketbra(Psi). $
However, an exterior observer has no access to the internal degrees of freedom $cal(H)_"BH"$ so their system is described by the reduced density matrix
$
  rho_"R" = tr_"BH" rho_"total",
$
which is thermal.#footnote[assuming the horizon is smooth.] This can be quantified by saying
$ S_R > 0, $
and as the black hole evaporates $rho_"R"$ becomes increasingly mixed, so $S_R$ increases. The natural implication is that this proceeds until the black hole evaporates completely leaving us with a thermal state.#footnote[ignoring remnants.] The state has then evolved from a pure state $ket(Psi)_"BH"$ to a mixed state $rho_"R"$ violating unitarity. This leads to information loss since we lose information encoded in the matter used to form the black hole when it is emitted as radiation. This is why we refer to this problem as the information paradox @Harlow_2016.#footnote[For the sake of brevity we ignore objections to Hawking's argument.]

#figure(
  scale(pagecurve, 80%),
  caption: [An example of the Page curve. The red dotted line represents non-unitary evaporation, while the black line represents unitary evaporation. The blue dotted line corresponds to the black hole entropy. @Harlow_2016 @Almheiri_2021 @Page_1993],
)<pagecurve>

The evolution of $S_R$ is called the Page curve @Page_1993 @Harlow_2016. We trivially have $S_R = 0$ initially and as evaporation proceeds $S_R > 0$ as mentioned. When unitarity is violated $S_R$  keeps increasing until the black hole is fully evaporated and $S_R tilde.eq S_"BH"^"initial"$. However, under the assumption that unitarity is sacred there should be some turnover point when $S_R tilde.eq S_"BH"$ such that $S_R -> 0$ eventually. We refer to the time this turnover occurs as the Page time $t_"Page"$.

An example of what the curve could look like is shown in @pagecurve. The exact details of the Page curve are not useful for our purposes. Note, however, that the final state being pure implies entanglement between late- and early radiation. This is needed to purify the late radiation and ensure $S_R -> 0$.

== The AMPS paradox
Almheiri, Marolf, Polchinski, and Sully made the information paradox formulated by Hawking precise with the AMPS paradox @Almheiri_2013. We will consider a formulation of the AMPS paradox aligned with the original formulation illustrating the importance of the entanglement structure near the horizon.#footnote[originally posed to show BHC @Susskind_1993 @tHooft1990 is inconsistent.] Later, we consider a slighty different formulation in the context of $"ER"="EPR"$ where this is helpful.

Consider an old black hole evaporating, with old meaning we are past the Page time. We assume the following: @raju2021lessonsinformationparadox @Almheiri_2013
1. Validity of semiclassical effective field theory near the horizon.

2. The process is unitary.

3. Validity of the equivalence principle, so the horizon is smooth.

The AMPS paradox claims these assumptions are mutually inconsistent, thus forming a "trilemma". We further assume, motivated by the semiclassical effective field theory description, that:

4. The Hilbert space can be decomposed as $ cal(H)_A times.o cal(H)_B times.o cal(H)_C subset cal(H), $<AMPShilbert> where $cal(H)_A$ describes degrees of freedom in the early radiation, $cal(H)_B$ describes degrees of freedom in the late radiation, and $cal(H)_C$ describes degrees of freedom inside the black hole.#footnote[We follow the naming used by AMPS in @Almheiri_2013.]

This assumption seems reasonable, as such decompositions emerge naturally in local quantum field theories.#footnote[This is not even mentioned and taken as being implicitly true by AMPS in @Almheiri_2013.] However, it is expected that such a strict decomposition fails in theories of quantum gravity. The primary objection to the AMPS paradox is for this reason the validity of this "naive" decomposition. @raju2021lessonsinformationparadox @Almheiri_2021 @Harlow_2016 @Guo_2021 @Jefferson_2019 We touch more upon this later, and will find $"ER"="EPR"$ challenges this decomposition and thereby "avoids" the AMPS paradox in a clever way.

#figure(
  scale(100%, sketch-AMPS),
  caption: [Sketch of the AMPS paradox. $A$ is a mode of early radiation which has had time to propagate, $B$ is a mode of late radiation sitting outside the horizon, and $C$ is the partner mode to $B$ inside the horizon. Adapted from @Cinti_2021.],
)<sketch-AMPS>

Now, consider a newly emitted mode of late radiation $B$ sitting just outside the horizon. By unitarity, $B$ must be entangled with the early radiation, so it has an entanglement partner $A$ in the previously emitted radiation. At the same time, the horizon is assumed to be smooth, meaning spacetime is locally Minkowski near the horizon and the horizon itself is locally Rindler. Semiclassical effective field theory then implies that the near-horizon state is the Minkowski vacuum, with entanglement across the horizon as shown with the Rindler decomposition $ket(Psi)_"M"$ in @mm-vacuum. Therefore $B$ must also be entangled with an interior partner mode $C$ behind the horizon. These modes are shown in @sketch-AMPS.

The mode $B$ is then maximally entangled with two different partners. This violates monogamy of entanglement since we assume the Hilbert space can be decomposed. This is simplest to see by using the strong subadditivity of entropy @nielsen2010quantum @Mathur_2009 @Almheiri_2013
$
  S_(A B) + S_(B C) >= S_B + S_(A B C).
$<subadditivity>
Now, old black holes have decreasing entropy so $S_(A B) < S_A$ and when the horizon is smooth $S_(B C) =0$ so $S_(A B C) = S_A$. Then @subadditivity becomes
$
  S_A > S_B + S_A,
$
implying $S_B < 0$ which fails horribly since the late radiation is thermal with $S_B > 0$.

Then the AMPS paradox as presented above implies we need to sacrifice unitarity, the equivalence principle, or our semiclassical effective field theory.#footnote[or the decomposition.] This is quite disturbing since the first two are sacred principles in physics, and we generally have faith in the validity of the effective description, especially considering we could imagine large black holes where the effect of gravity is weak. AMPS took the most conservative resolution to be giving up the smooth horizon. @Almheiri_2013

This is referred to as there being a firewall since any observer crossing the horizon will be met by a wall of highly energetic particles. This is the case since we know the spacetime near the horizon is approximately Rindler. The near-horizon state should then be the Minkowski vacuum, where entanglement ensures the horizon appears smooth. When disrupting the entanglement, as AMPS concludes we should, then the near-horizon state can no longer be the Minkowski vacuum. This leads to a highly excited state and a large energy density at the horizon. @Harlow_2016 @harlow2023blackholesquantumgravity @Almheiri_2013

#pagebreak()
= $"ER" = "EPR"$
This chapter discusses $"ER"="EPR"$ and its implications for black hole evaporation, including how $"ER"="EPR"$ "solves" the AMPS paradox.

== The conjecture
We have seen both entanglement and Einstein-Rosen bridges appear non-local. We can find another similarity by considering LOCC.

We imagine two separated systems, $A$ and $B$. LO refers to any local operations Alice or Bob can perform on their respective systems, while CC refers to classical communication. By this we mean Alice and Bob are only allowed to send _classical_ information to each other. For example, Alice may perform measurements on her system and communicate the outcomes to Bob. However, LOCC prohibits Alice from sending a qubit#footnote[or any quantum system] directly to Bob. LOCC is importantly unable to create or increase entanglement between $A$ and $B$. We can only entangle $A$ and $B$ by bringing the systems together or by using a resource of already entangled states. @Chitambar_2014 @nielsen2010quantum @Maldacena_2013

Now, consider two separated black holes with no Einstein-Rosen bridge. Creating an Einstein-Rosen bridge between these is impossible. However, we can create a pair of black holes with an Einstein-Rosen bridge and then separate them. We can also create an Einstein-Rosen bridge between two separated black holes by using another pair of black holes with an Einstein-Rosen bridge. This is done by merging members of each pair. This is completely analogous to LOCC. @Maldacena_2013

We could even imagine creating many pairs of small black holes with Einstein-Rosen bridges. We give half of each pair to Alice and the other to Bob. Now, we let Alice and Bob separate after which they merge all their small black holes. We will be left with a pair of black holes with an Einstein-Rosen bridge.


#figure(
  scale(100%, EPRgedanken),
  caption: [The imagined splitting of entangled particle pairs.],
)<gedanken>

We can push this experiment further by creating many pairs of entangled particles. We give half of each pair to Alice and the other to Bob, as illustrated in @gedanken. Now, we let Alice and Bob separate after which they collapse their share of particles. We will be left with a pair of entangled black holes. Susskind and Maldacena claim these will also have an Einstein-Rosen bridge. @Maldacena_2013



This is generalised by them claiming even a single pair of entangled particles will be connected through some highly quantum Einstein-Rosen bridge. This proposed equivalence between entanglement and Einstein-Rosen bridges is summarised by the symbolic equation @Maldacena_2013
$
  "ER" = "EPR",
$
which is admittedly a wild conjecture. However, we already saw such an equivalence occur in the case of a Schwarzschild black hole in anti-de Sitter spacetime.

== The emitted radiation
To see some implications of $"ER"="EPR"$ consider an evaporating black hole. The emitted radiation is entangled with the black hole. So by $"ER"="EPR"$ each quanta of radiation will be connected with the black hole through an Einstein-Rosen bridge. We can illustrate this as in @radiationER.

#figure(
  scale(entanglementradiation, 80%),
  caption: [The Einstein-Rosen bridge between a black hole and its emitted radiation. Adapted from @Maldacena_2013.],
)<radiationER>

We can also try to illustrate the evolution of the Einstein-Rosen bridge over time, assuming it can be described geometrically, with embedding diagrams as in @b. The familiar embedding diagram for a typical Einstein-Rosen bridge is also shown in @a.
#subpar.grid(
  figure(
    scale(ERdiagram, 80%),
    caption: [The embedding diagram for an Einstein-Rosen bridge between exteriors $L$ and $R$ in the fully-extended Schwarzschild geometry.],
  ),
  <a>,
  figure(
    image("fig/radiationLIFE.png"),
    caption: [The evolution of the Einstein-Rosen bridge between a black hole and its emitted radiation @Maldacena_2013.],
  ),
  <b>,

  columns: 1fr,
  caption: [Embedding diagrams.],
  label: <radiationlife>,
)

The far left shows a newly formed black hole with the black circle representing the horizon as seen from the outside. As the black hole evaporates, it emits radiation. The red dots show where quanta of radiation are connected to the main body of the Einstein-Rosen bridge. The green circle represents a slice through the Einsten-Rosen bridge, dividing the system into two parts. Radiation to the right of the slice was emitted earlier than radiation to the left. With this in mind we can think of these slices as being proxies of time. @Maldacena_2013


The entanglement entropy of the radiation $S_R$ should follow the Page curve as to avoid violations of unitarity. Now, we claim $S_R ("a slice") prop A_"slice"$ @Maldacena_2013.#footnote[see also @verlinde2020ereprrevisited.] This kind of "area law" relating entanglement entropy to geometry was originally motivated by the Bekenstein-Hawking entropy @bekenstein where $S_"BH" prop A$. Bombelli et al. @Bombelli_1986 and Srednicki @Srednicki_1993 have explicitly shown the entanglement entropy of quantum fields obeys
$
  S_"EE" prop A/epsilon.alt^2,
$
where $epsilon.alt$ is a UV-cutoff length scale. Now, identifying the cutoff with the Planck length $cal(l)_p$ gives
$
  S_"EE" prop A/cal(l)_p^2,
$
which is $S_"BH"$ up to a numerical factor, thereby making the black hole analogy explicit. Hence, the claim $S_R prop A$ is reasonable.

We know $S_R$ is maximal at the Page time, which corresponds to the blue circle in @b. The Einstein-Rosen bridge at later times should then begin shrinking in accordance with the Page curve since $S_R prop A$. As more and more radiation is emitted the shape eventually develops into the narrow, horn-like structure shown @Maldacena_2013. With this view the Einstein-Rosen bridge can be seen as a "geometric Page curve". However, we stress @b is speculative and assumes the Einstein-Rosen bridge even has a geometric description.

== Smoothness of the horizon
Consider an observer Bob sitting outside an old black hole. Another observer Alice, of which Bob is unaware, has been collecting the emitted radiation. When she has more than half the radiation under control we imagine she collapses it to form a black hole, as shown in @Alicecollapse. This newly formed black hole will be entangled with Bob's black hole. Now, Alice can in principle operate on her black hole and form the Hartle-Hawking state @HH.#footnote[ignoring limits on computation.] @Maldacena_2013

#figure(
  scale(81%, smoothness),
  caption: [Alice collapses emitted radiation to form a black hole entangled with Bob's black hole.],
)<Alicecollapse>

With the above setup we will now consider an operational version of the AMPS paradox. This is done since reframing the paradox allows us to see what dynamically changes upon imposing $"ER"="EPR"$.

Consider a mode $B$ which is about to be emitted from Bob's black hole. Alice has at some earlier time distilled a mode from her black hole, namely the mode that $B$ will be entangled with. We denote this mode by $A$ as to be consistent with our previous notation. By distill we mean the operation
$
  "early radiation" ->^"distill" A times.o ("everything else"),
$
so Alice separates $A$ from her black hole. Such an operation is in principle possible assuming evaporation is unitary and all information is preserved.

With $A$ in her possession Alice flies over to Bob's black hole in time to meet $B$ when it is emitted. Then following our previous logic $B$ and $A$ are entangled implying $B$ cannot be entangled with the usual partner mode $C$ behind the horizon. The disruption of the $B C$ entanglement means Alice will encounter a particle upon crossing the horizon. @Maldacena_2013

So whenever Alice performs the above experiment for some $B$ she will find the corresponding $B C$ entanglement to be disrupted. We refer to this as an operational version of the argument since the problem with monogamy only becomes apparent once Alice brings $A$ to $B$. Now AMPS assumes Alice is superfluous. We take this to mean the $B C$ entanglement is disrupted even when she does not perform the experiment @Maldacena_2013. Upon performing the experiment she merely reveals the disruption which was already present. The implication is that all possible $B C$ pairs are disrupted initially creating a firewall. This assumption is quite reasonable since nothing Alice does should lead to the creation of the particle behind the horizon. This stems from AMPS assuming we can decompose the Hilbert space as in @AMPShilbert.

By $"ER"="EPR"$ this assumption is false. The particle Alice observes upon crossing the horizon was simply created as the other end of the Einstein-Rosen bridge when she distilled $A$. @Maldacena_2013 Thereby removing the necessity for more than a single $B C$ pair being disrupted. Clearly this would not be a viable explanation without imposing $"ER"="EPR"$ since there would be no such Einstein-Rosen bridge.

We now see $"ER"="EPR"$ actually implies something "stronger" about what happens to Bob when he decides to cross the horizon. We can imagine Alice performing various operations on her black hole or the radiation. Bob will then as argued be affected by these. She could in principle arrange for Bob to be met by a shockwave or bomb such that he would experience a firewall upon crossing the horizon. She could also do nothing, in which case the horizon would be smooth. Now, her actions are in principle determined by the initial state of the system, so what Bob experiences is probabilistic in the quantum sense. We conclude he may be met by a firewall, a smooth horizon, or some bumpy middle ground. One could argue the odds of him meeting a firewall are actually pretty slim, since this requires many $B C$ pairs being disrupted. Nonetheless, we can certainly conclude the presence of a firewall is uncertain. @Maldacena_2013

#pagebreak()
= Discussion and Outlook
$"ER"="EPR"$ initially appears to be a wild proposal. However, as discussed throughout the preceding sections, the idea is motivated by several deep connections and similarities between entanglement and black holes.

We have also seen one concrete realisation of $"ER"="EPR"$ arise in the context of black holes in anti-de Sitter spacetime. Here, applying the AdS/CFT correspondence led to an equivalence between the thermofield double state of two independent CFTs and a connected spacetime containing an Einstein-Rosen bridge, where the Einstein-Rosen bridge may be interpreted as the geometric realisation of the entanglement in the boundary theory, symbolically expressed as
$
  "ER-bridge" <=>^"AdS/CFT" "TFD entanglement".
$
With this example in mind, one may formulate a weak version of $"ER"="EPR"$, namely that all entangled black holes are also connected by Einstein-Rosen bridges. This is the initial conceptual jump made by Susskind and Maldacena in @Maldacena_2013, and it appears plausible since such a relation emerges naturally in anti-de Sitter spacetime. However, this weak form of the conjecture still relies on highly idealised conditions, and whether the correspondence generalises is still non-obvious.

The stronger version of $"ER"="EPR"$ would claim that all entangled particles are connected by quantum Einstein-Rosen bridges. This represents a much larger conceptual leap and, aside from the example above, is primarily motivated by observed similarities between entanglement and Einstein-Rosen bridges. The precise meaning of the proposal is also vague, and especially the notion of a "quantum Einstein-Rosen bridge" between two particles is ill-defined. This makes $"ER"="EPR"$ hard to formulate in a precise way and correspondingly difficult to test or falsify in a meaningful manner.

We have also seen how $"ER"="EPR"$ offers a possible resolution to the AMPS paradox by providing an explicit reason why the naive Hilbert space decomposition fails, while also explaining dynamically why this occurs. More recent developments related to the black hole information paradox further strengthen the connection between entanglement and spacetime geometry. As an example we have the appearance of replica wormholes and quantum extremal surfaces in the island formula @penington2020replicawormholesblackhole @Almheiri_2020 @Engelhardt_2015  which suggests entanglement may play a role in determining spacetime geometry. While these ideas lie way beyond the scope of this thesis, they still point toward the same broader picture as $"ER"="EPR"$, and some claim these ideas are concrete realisations of $"ER"="EPR"$. @Guo_2021 @Cinti_2021

Whether $"ER"="EPR"$ ultimately proves correct or false, the proposal has nevertheless influenced modern approaches#footnote[at least partly.] to quantum gravity by suggesting a deep relationship between entanglement and spacetime geometry. With this view, $"ER"="EPR"$ may be considered not only as a proposed resolution to the AMPS paradox, but also as part of a broader shift in theoretical physics toward understanding spacetime and gravity as emergent phenomena.


#pagebreak()
#bibliography("biblo.bib", style: "american-physics-society")
