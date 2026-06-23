#import "../../temp.typ": *
#show: chpt-note.with()

= Discrete Symmetries

/*
We care about symmetries in physics because of Noether's theorem which tells us that all continuous symmetries of an action $S$ imply conserved currents.

== Group theory
The notion of a symmetry is best described in the context of groups. A group is a set $G$ along with a binary operation on $G$ that combines any two elements $a, b in G$ to form an element $a dot b in G$ such that the group axioms are satisfied:

1. For all $a, b, c in G$ we have $(a dot b) dot c = a dot (b dot c)$.

2. There exists an identity element $e in G$ such that for all $a in G$ we have $e dot a = a dot e = a$.

3. For each $a in G$ there exists an element $b in G$ such that $a dot b = b dot a = e$. We call $b = a^(-1)$ the inverse of $a$.

We call commuting groups Abelian. We can have finite or infinite groups, in particular we have continuous groups.#footnote[Within physics all such continuous groups are Lie groups.]

The groups we care about in physics can typically be written as groups of matrices. We especially care about $"U"(n)$, $"SU"(n)$, $"O"(n)$, and $"SO"(n)$.#footnote[To be concrete $"SO"(3)$ would describe rotational spatial symmetry and it appears almost identical to $"SU"(2)$.] Any group $G$ can be represented by a group of matrices. This means each group element $a in G$ has a corresponding matrix $M_a$ with the correspondence respecting group multiplication,
$
  M_a M_b = M_c
$
This representation is not necessarily faithful since multiple group elements can be represented by the same matrix.#footnote[As an example take $M_a = bb(1)$ for all $a in G$.] We would say $G_M$ is homomorphic to $G$. Any $G$ which is already a group of matrices is obviously a representation of itself. We call this the fundamental representation. But, we can also have representations of other dimensionalities. We can construct new representations by combining old ones
$
  M_a = mat(M_a^((1)), 0; 0, M_a^((2)))
$
However, these are boring and we only care about irreducible representations which cannot be decomposed as above.
*/

== Angular momentum
Within quantum mechanics we care about the orbital angular momentum $bold(L)$ and the spin#footnote["intrinsic"] angular momentum $bold(S)$. These angular momenta $bold(J)$ are very similar and formally they are essentially the same. However, there are some differences. When performing measurements of $bold(L) = bold(r) times m bold(v)$ we find these are quantised#footnote[See e.g. Sakurai for details. We measure $L^2$ since the $L_i$ are non-commuting.]
$
  L^2 & tilde l(l+1) hbar^2 \
  L_z & tilde m_l hbar
$
where $l = 0, 1, dots$ and#footnote[This implies $bold(L)$ cannot be oriented along $hat(L)_z$.]
$
  m_l = -l, -l+1, dots, l-1, l
$
When performing measurements of $bold(S)$ we similarly have
$
  S^2 & tilde s(s+1)hbar^2 \
  S_z & tilde m_s hbar
$
where $s = 0, 1/2, 1, dots$ and
$
  m_s = -s,-s+1, dots, s-1,s
$
Now, all particles can be given any $l$, but each type of particle have their own spin $s$! We call particles with half-integer spin fermions and those with integer spin bosons.#footnote[All baryons, leptons, and quarks are fermions, while all mesons and mediators are bosons.]

== Addition of angular momentum
We can represent angular momentum states by $ket(l m_l)$ and $ket(s m_s)$. We would like to understand composite systems and how to handle the addition of angular momenta
$
  bold(J) = bold(J)_1 + bold(J)_2
$
This could be the total angular momentum of a single particle
$
  bold(J) = bold(L) + bold(S)
$
or the spin of a composite system
$
  bold(J) = bold(S)_1 + bold(S)_2
$
This is made difficult since
$
  [J_i, J_j] = i hbar epsilon_(i j k) J_k
$
implying we do not have access to all components at once. We are forced to work with a single component $J_z$ and the magnitude $J^2$. We would like
$
  ket(j m) tilde ket(j_1 m_1)+ket(j_2 m_2)
$
The $J_z$ add trivially so
$
  m = m_1 + m_2
$
However, the $J^2$ is more annoying and we get all
$
  j = abs(j_1-j_2), abs(j_1-j_2)+1,dots, j_1+j_2 -1, j_1+j_2
$
These explain the structure of The Eightfold Way since $B = q q q$ and $M = q overline(q)$ leading to particles with different spins.

== P
Now, we consider the parity operator $P$ which acts on a state $ket(bold(x))$ as
$
  P ket(bold(x)) = ket(-bold(x)).
$
Lee and Yang proposed#footnote[Motivated by the $theta tau$-puzzle, where two particles seemed identical, but decayed as $theta -> 2 pi$ and $tau -> 3 pi$, thereby violating parity if $theta$ and $tau$ were the same particle.] that weak interactions violate parity around $~ 1956$ which was shown experimentally by Wu.#footnote[Shown using $beta$-decay of cobalt nuclei. The direction of spin is reversed in a mirror, however, the $beta$-rays are emitted in the same direction. Then reorienting the spin (flipping the image) we find the $beta$-rays to be emitted in the opposite direction in the mirror. Thereby showing $beta$-decay (a weak interaction) is changed under parity.]

Now, consider a particle with velocity $v$ along the $z$-axis (by definition). We define the helicity $h$ as the value of
$
  m_s/s = plus.minus 1,
$
along this axis.#footnote[Here, $m_s$ is the "value" of $S_z$. An $e^-$ with $s = +1/2$ (along $S_z$) moving along $+hat(z)$ would have $h = +1$, while an $e^-$ with $s = -1/2$ moving along $+hat(z)$ would have $h = -1$. However, an $e^-$ with $s = -1/2$ moving along $-hat(z)$ would have $h = +1$.] When the spin and velocity are parallel ($h = + 1$) we call the particle right-handed. Otherwise, we call it left-handed. Under parity we have
$
  L <-->^P R
$


However, the value of $h$ is not invariant!#footnote[An observer with velocity larger than $v$ would measure $-h$.] This is nullified when $v = c$. Now, we assume $m_nu tilde 0$ implying
$
            nu_L & tilde "left-handed" \
  overline(nu)_R & tilde "right-handed"
$
This can be shown using $pi^-$-decay
$
  pi^- -> mu^- + overline(nu)_mu
$
where $mu^-$ and $overline(nu)_mu$ are emitted back-to-back. The spin is conserved so we can determine the handedness of $overline(nu)_mu$ by measuring $mu^-$.#footnote[All $mu^-$ are observed to be right-handed implying all $overline(nu)_mu$ are right-handed. Since $S_mu + S_nu =^! 0$.] Likewise, we can use $pi^+$-decay
$
  pi^+ -> mu^+ + nu_mu
$
showing $nu_mu$ is left-handed. This is unlike photons $gamma$ which have no preferred handedness, which can be seen from $pi^0$-decay.

We can characterise objects according to how they transform under parity#footnote[The "pseudo" implies they transform "unnormally".]
$
        "scalar" & tilde P(s) = s \
  "pseudoscalar" & tilde P(p) = -p \
        "vector" & tilde P(bold(v)) = - bold(v) \
  "pseudovector" & tilde P(bold(a)) = bold(a)
$
Also, we have
$
  P^2 = bb(1)
$
implying the eigenvalues of $P$ are $plus.minus 1$. Now, we will assign
$
  P (q) = +1",  " P(overline(q)) = -1
$
implying#footnote[Note, with angular momentum $l$ we include a factor $(-1)^l$.]
$
  P (B) = +1",  " P(overline(B)) = -1",  " P(M) = -1
$
To describe photons we use the vector potential $A^mu$ so
$
  P(gamma) = -1
$

== C
Now, we consider the charge conjugation operator $C$ which acts on a state $ket(p)$ as
$
  C ket(p) = ket(overline(p)).
$
The operator $C$ changes the sign of all internal quantum numbers.#footnote[This leaves mass, energy, momentum, and spin unchanged.] Also, we have
$
  C^2 = bb(1),
$
implying the eigenvalues of $C$ are $plus.minus 1$. Then an eigenstate of $C$ satisfies
$
  C ket(p) = plus.minus ket(p),
$
which holds for particles who are their own antiparticles ($gamma, pi^0, eta$ etc). We can show a $q overline(q)$ state with $l > 0$ and spin $s$ is an eigenstate of $C$ with eigenvalue $(-1)^(l+s)$. Then pseudoscalar mesons ($l = 0$ and $s = 0$) have the eigenvalue $(+1)$, while vector mesons ($l = 0$ and $s = 1$) have the eigenvalue $(-1)$.

Now, consider acting with $C$ on $nu_L$
$
  nu_L ->^C overline(nu)_L tilde "does not exist".
$
This shows weak interactions violate invariance under charge conjugation.

== CP
Consider the $pi^+$-decay
$
  pi^+ -> mu^+ + nu_(mu, L)
$
Now, we apply $C$ and $P$
$
  pi^- -> mu^- + overline(nu)_(mu,R)
$
which does occur, even though $C$ and $P$ individually are violated. However, $C P$ is also violated.

The process
$
  K^0 <--> overline(K)^0
$
is possible by weak interactions. When measuring these particles we will then find linear combinations of $K^0$ and $overline(K)^0$.

Now, we have
$
  P ket(K^0) = - ket(K^0)",  " P ket(overline(K)^0) = - ket(overline(K)^0),
$
and
$
  C ket(K^0) = ket(overline(K)^0)",  " C ket(overline(K)^0) = ket(K^0),
$
implying
$
  C P ket(K^0) = - ket(overline(K)^0)",  " C P ket(overline(K)^0) = - ket(K^0).
$
We find two eigenstates of $C P$
$
  ket(K_1) & = 1/sqrt(2) (ket(K^0) - ket(overline(K)^0)), \
  ket(K_2) & = 1/sqrt(2) (ket(K^0) + ket(overline(K)^0)),
$
where
$
  C P ket(K_1) = ket(K_1)",  " C P ket(K_2) = - ket(K_2).
$

Now, assuming $C P$ is not violated then $K_1$ can only decay into states with $C P = +1$, while $K_2$ can only decay into states with $C P = - 1$. This implies#footnote[$2pi$ has $(+1) (+1)$ and $3 pi$ has $(-1)(+1)$, with the factors being $P C$. Also, $K$ likes decaying into $pi$.]
$
  K_1 -> 2 pi",  " K_2 -> 3 pi
$
with the $K_1$-decay being much quicker.

Now, imagine a $K^0$-beam with
$
  ket(K^0) = 1/sqrt(2) (ket(K_1) + ket(K_2))
$
Since the $K_1$-decay is quick we will see many $2 pi$-decays initially and many $3 pi$-decays near the end of our beam. The $K_1$ and $K_2$ have been found and their lifetimes were determined to be
$
  tau_1 & tilde 10^(-10) "sec" \
  tau_2 & tilde 10^(-8) "sec"
$
We can then see if $C P$ is violated since we can produce an arbitrarily pure $K_2$-beam. We would then still see $2 pi$-decays if $C P$ is violated.#footnote[The "long-lived-species" would contain small amounts of $K_1$ and would therefore not be an eigenstate of $C P$ as assumed.]

== T
Now, consider the time reversal operator $T$ taking $t -> - t$. The operator $T$ is important because of the $C P T$ theorem: "any Lorentz invariant local quantum field theory with an Hermitian Hamiltonian $H^dagger = H$ preserves $C P T$".#footnote[The standard model is such a theory.] This implies $T$ is violated since the combination $C P T$ is preserved. However, this is experimentally very difficult to test.
