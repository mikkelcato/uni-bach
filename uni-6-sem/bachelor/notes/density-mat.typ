//**** init-ting
#import "temp.typ": *


#show: note.with(
  title: [
    $rho$ and Schmidt
  ],
  name: [Mikkel Kielstrup Holst],
)

= The density operator
We would like to study _open_ systems. These are simply parts of some larger system. Consider a system of two qubits $A$ and $B$. The full system $A B$ obeys the usual axioms of quantum mechanics. We take qubit $B$ as being locked away, meaning we do not have access to it. The goal is to characterize the observations which can be made on qubit $A$ alone. Such a system is analogous to a black hole in an otherwise empty universe.

As a basis we use ${ket(0)_A, ket(1)_A}$ and ${ket(0)_B, ket(1)_B}$. Consider the state
$
  ket(psi)_(A B) = a ket(0)_A times.circle ket(0)_B + b ket(1)_A times.circle ket(1)_B
$
Here qubits $A$ and $B$ are _correlated_. This should be obvious. Consider a general observable acting only on qubit $A$ by
$
  M_A times.circle underbracket(bb(1)_B, "identity")
$
Then
$
  expval(M_A) & = braket(psi, M_A times.circle bb(1)_B, psi) \
              & = abs(a)^2 braket(0, M_A, 0) + abs(b)^2 braket(1, M_A, 1) \
              & = tr(M_A rho_A)
$
Where we have defined the _density operator_
$
  rho_A = abs(a)^2 ketbra(0, 0) + abs(b)^2 ketbra(1, 1)
$
The above generalizes to any bipartite system $cal(H)_A times.circle cal(H)_B$ with basis ${ket(i)_A times.circle ket(mu)_B}$. Any pure state of this system can be written as
$
  ket(psi)_(A B) = sum_(i, mu) a_(i mu) ket(i)_A times.circle ket(mu)_B
$
Then
$
  expval(M_A) & = sum_(i, j, mu) a_(j mu)^* a_(i mu) braket(j, M_A, i) \
              & = tr(M_A rho_A)
$
Where
$
  rho_A = tr_B (ketbra(psi)) equiv sum_(i,j,mu) a_(i mu) a^*_(j mu) ketbra(i, j)
$
which is the proper definition of the density operator. We say that $rho_A$ is obtained by doing the _partial trace_ over $B$ of the density the operator for the combined system $A B$. We can immediately infer the following
$
  & rho_A "is Hermitian" rho_A = rho_A^dagger \
  & rho_A "is positive" \
  & tr(rho_A) = 1
$
Assuming the state $ket(psi)_A$ is pure (described by a ket) then $rho^2 = rho$. Any $rho$ can be written as
$
  rho_A = sum_a p_a ketbra(a, a)
$
with ${ket(a)}$ being a diagonal basis. For a _mixed_ state  this sum contains more than one term and $rho^2 eq.not rho$, and we call $rho$ _incoherent_. We may interpret $rho$ as describing an _ensemble_ of pure states. When two states $A$ and $B$ become entangled we achieve such a state. The total state is pure, but considering only $A$ we have a loss of coherence due to entanglement. Assuming $A$ were pure then it would be described by $ket(psi)_A$ and the total state would factor $ket(psi)_(A B) = ket(psi)_A times.circle ket(psi)_B$. This is not an entangled state, so entanglement leads to mixed states. We lose the information carried by correlations between $A$ and $B$!

= Schmidt decomposition

