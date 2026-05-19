#import "chpt-temp.typ": *
#show: chpt-note.with()

= Symmetries
We care about _symmetries_ in physics because of _Noether's theorem_ which tells us that all continuous symmetries of an action $S$ imply conserved currents.

== Group theory
The notion of a symmetry is best described in the context of groups. A group is a set $G$ along with a binary operation on $G$ that combines any two elements $a, b in G$ to form an element $a dot b in G$ such that the _group axioms_ are satisfied:

1. For all $a, b, c in G$ we have $(a dot b) dot c = a dot (b dot c)$.

2. There exists an identity element $e in G$ such that for all $a in G$ we have $e dot a = a dot e = a$.

3. For each $a in G$ there exists an element $b in G$ such that $a dot b = b dot a = e$. We call $b = a^(-1)$ the inverse of $a$.

We call commuting groups _Abelian_. We can have _finite_ or _infinite_ groups, in particular we have _continuous_ groups.#footnote[Within physics all such continuous groups are _Lie_ groups.]

The groups we care about in physics can typically be written as groups of matrices. We especially care about $"U"(n)$, $"SU"(n)$, $"O"(n)$, and $"SO"(n)$.#footnote[To be concrete $"SO"(3)$ would describe rotational spatial symmetry and it appears almost identical to $"SU"(2)$.] Any group $G$ can be _represented_ by a group of matrices. This means each group element $a in G$ has a corresponding matrix $M_a$ with the correspondence respecting group multiplication,
$
  M_a M_b = M_c
$
This representation is not necessarily _faithful_ since multiple group elements can be represented by the same matrix.#footnote[As an example take $M_a = bb(1)$ for all $a in G$.] We would say $G_M$ is _homomorphic_ to $G$. Any $G$ which is already a group of matrices is obviously a representation of itself. We call this the _fundamental_ representation. But, we can also have representations of other dimensionalities. We can construct new representations by combining old ones
$
  M_a = mat(M_a^((1)), 0; 0, M_a^((2)))
$
However, these are boring and we only care about _irreducible_ representations which cannot be decomposed as above.

== Angular momentum
Within quantum mechanics we care about the orbital angular momentum $bold(L)$ and the _spin_ angular momentum $bold(S)$. These angular momenta $bold(J)$ are very similar and formally they are essentially the same. However, there are some differences. When performing measurements of $bold(L)$ we find these are quantized#footnote[See e.g. _Sakurai_ for details.]
$
  L^2 & tilde l(l+1) hbar^2 \
  L_z & tilde m_l hbar
$
where $l = 0, 1, dots$ and
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
All particles can be given any $l$, but each type of particle have their own spin $s$. We call particles with half-integer spin _fermions_ and those with integer spin _bosons_.#footnote[All baryons, leptons, and quarks are fermions, while all mesons and mediators are bosons.]

We can represent angular momentum states by $ket(l m_l)$ and $ket(s m_s)$. We would like to understand composite systems and how to handle the addition of angular momenta
$
  bold(J) = bold(J)_1 + bold(J)_2
$
This could be the _total_ angular momentum of a single particle
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
We can explicitly decompose $ket(j_1 m_1) ket(j_2 m_2)$ as
$
  ket(j_1 m_1) ket(j_2 m_2) = sum_(j=abs(j_1-j_2))^(j_1+j_2) underbracket(C_(m m_1 m_2)^(j j_1 j_2), "Clebsch-Gordan" #linebreak() "coefficients") ket(j m)
$
with $m = m_1+ m_2$.

== Discrete symmetries
Consider the _parity_ operator $pi$ which acts on a state $ket(bold(x))$ as
$
  pi ket(bold(x)) = ket(-bold(x))
$
Classically people believed physics were invariant under parity. However, weak interactions were shown to break parity in 1956 by Wu.

Consider a particle moving with velocity $v$ we define the $z$-axis to be the direction of motion. We define the _helicity_ as the value of
$
  h tilde m_s/s = plus.minus 1
$
along this axis. When $h = +1$ the spin $s$ and velocity $v$ are parallel and we call the particle _right-handed_.#footnote[Also, note that $h$ is odd under parity.] Otherwise, we call it _left-handed_. However, the value of $h$ is not invariant since an observer moving with velocity greater than $v$ would measure the opposite. The exception is particles moving with velocity $c$ since their motion cannot be reversed. We can treat neutrinos as massless and by experiment we find#footnote[This shows parity is _maximally_ violated.]
$
            nu_L & tilde "left-handed" \
  overline(nu)_R & tilde "right-handed"
$
This occurs since weak interactions break parity. We can similarly determine the _handed-ness_ of photons and we find an even split since electromagnetic interactions respect parity.

We characterize objects depending on how they transform under parity:#footnote[We see a theory with only pseudovectors and scalars is invariant under parity. However, the addition of a vector breaks this invariance. This is exactly what happens in weak interactions.]
$
        "scalar" & tilde pi(s) = s \
  "pseudoscalar" & tilde pi(p) = -p \
        "vector" & tilde pi(bold(v)) = - bold(v) \
  "pseudovector" & tilde pi(bold(a)) = bold(a)
$
Also, by definition we have
$
  pi^2 = bb(1)
$
implying the eigenvalues of $pi$ are $plus.minus 1$. We associate
$
  pi (q) = +1";  " pi(overline(q)) = -1
$
implying
$
  pi (B) = +1";  " pi(overline(B)) = -1
$
To describe photons we use the vector potential $A^mu$ so
$
  pi(gamma) = -1
$

Consider the _charge conjugation_ operator $C$ which acts on a state $ket(p)$ as
$
  C ket(p) = ket(overline(p))
$
The operator $C$ changes the sign of _all_ internal quantum numbers.#footnote[This leaves mass, energy, momentum, and spin unchanged.] Again by definition
$
  C^2 = bb(1)
$
implying the eigenvalues of $C$ are $plus.minus 1$. Then an eigenstate of $C$ satisfies
$
  C ket(p) = plus.minus ket(p)
$
which is only true for particles who are their own antiparticles. One can show a system consisting of a spin $1/2$ particle and its antiparticle in a configuration with orbital angular momentum $l$ and total spin $s$, is an eigenstate of $C$ with the eigenvalue $(-1)^(l+s)$.#footnote[Think $M = q overline(q)$.] We again have weak interactions violating invariance under $C$ since
$
  nu_L ->^C overline(nu)_L tilde "does not exist"
$

We have seen that $pi$ is not respected due to
$
  pi^+ -> mu^+ + nu_(mu, L)
$
and $C$ is not respected since acting with $C$ gives
$
  pi^- -> mu^- + overline(nu)_(mu,L)
$
which is impossible. However, we can act with both $pi$ and $C$. Then we obtain
$
  pi^- -> mu^- + overline(nu)_(mu,R)
$
which does occur! Thus we have found that physics appears to be invariant under $C P$. This has some weird consequences. We consider the process
$
  K^0 <-> overline(K)^0
$
which is possible. This means the particles we typically observe are some linear combinations of $K^0$ and $overline(K)^0$. The $K$'s are pseudoscalars so
$
  pi ket(K^0) = - ket(K^0)";  " pi ket(overline(K)^0) = - ket(overline(K)^0)
$
Also,
$
  C ket(K^0) = ket(overline(K)^0)";  " C ket(overline(K)^0) = ket(K^0)
$
implying
$
  C P ket(K^0) = - ket(overline(K)^0)";  " C P ket(overline(K)^0) = - ket(K^0)
$
The normalized eigenstates of $C P$ are then
$
  ket(K_1) = 1/sqrt(2) (ket(K^0) - ket(overline(K)^0))";  " ket(K_2) = 1/sqrt(2) (ket(K^0) + ket(overline(K)^0))
$
with
$
  C P ket(K_1) = ket(K_1)";  " C P ket(K_2) = - ket(K_2)
$
Assuming $C P$ is respected in weak interactions then $K_1$ can only decay into a state with $C P = +1$, while $K_2$ can only decay into a state with $C P = -1$. This implies#footnote[Since each $pi$ carries $pi = -1$.]
$
  K_1 -> 2 pi";  " K_2 -> 3 pi
$
with the $K_1$ decay being much faster. Then we imagine a beam of $K^0$'s
$
  ket(K^0) = 1/sqrt(2) (ket(K_1) + ket(K_2))
$
The $K_1$ component will decay quickly meaning at the end of our beam we should have a beam of pure $K_2$'s. Then we expect many $2 pi$ decays near the source and many $3 pi$ near the end of our beam. The existence of $K_1$ and $K_2$ was shown experimentally and their lifetimes were determined to be
$
  tau_1 & tilde 10^(-10) "sec" \
  tau_2 & tilde 10^(-8) "sec"
$
This allows us to test if $C P$ is actually respected since given we have a long enough beam we can produce an arbitrarily pure sample of $K_2$'s. Then if we still observed $2 pi$ decays $C P$ would be violated. This means the long lived species is not an eigenstate of $C P$ as we assumed, but instead contains a small amount of $K_1$
$
  ket(K_L) = 1/sqrt(1 + abs(epsilon)^2) (ket(K_2) + epsilon ket(K_1))
$
with $epsilon$ measureing the amount of $C P$ violation. By experiment we determine
$
  epsilon tilde 10^(-3)
$

Consider the _time reversal_ operator $T$ taking $t -> - t$. We care about this operator because of the _$C P T$ theorem_: any Lorentz invariant local quantum field theory with an Hermitian Hamiltonian $H^dagger = H$ respects $C P T$.#footnote[The standard model is such a theory.] This implies $T$ is violated since the combination $C P T$ is respected. However, this is experimentally very difficult to test.
