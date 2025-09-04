//**** init-ting
#import "@preview/physica:0.9.5": *
#import "temp.typ": *


#show: thmrules.with(qed-symbol: $square$)
#show: note.with(
  title: [
    *soft matter*
  ],
  authors: (
    (
      name: "mkh",
    ),
  ),
  abstract: [
    Notes on soft matter physics loosely following Doi's book.
  ],
)

= Introduction

== Entropy review
From statistical mechanics
$
  S =^"Gibbs" - k_B sum_i^Omega P_i ln P_i =^"Boltzmann" k_B ln Omega
$
interpreting entropy as suprise or information leads to the definition of Shannon entropy,
$
  I = -k log P_i
$
where $k$ is some arbitrary constant, take e.g. $k = 1 "bit"$ then
$
  I = - 1 "bit" log_2 P_i
$
and $k = k_B$ gives the Boltzmann entropy. If we have multiple events we average this
$
  H = expval(I) = - k sum P_i log P_i
$
for $k = k_B$ this is just the Gibbs entropy. We can relate the statistical mechanical form of entropy to information entropy by
$
  S/k_B = - sum P_i ln P_i = (H ln 2)/(1 "bit")
$
to find $1 "bit" = k_B ln 2$---so the heat required to e.g. erase a bit of information would be $Q = T dot (k_B ln 2)$.

== Entropic forces
A polymer can be treated as a random walk. Given a single state has a volume $V_0 = lambda_"th"^3$ then the amount of states in a volume $V$ is just $V\/V_0$. How many polymers can we make? If for every step we have $z$ options, and our polymer is $N$ units long, i.e. we take $N$ random steps to create it, then
$
  Omega_1 = V/V_0 z^N
$
so the entropy is
$
  S_1 = k_B ln V/V_0 + k_B N ln z
$
the first part is the usual translational entropy, but now we also get an extra term do to conformation. This is essential for polymer dynamics since every system seeks maximal entropy. Consider a free polymer and a polymer next to a wall, in this case $S_"wall" < S_"free"$ since the polymer near the wall has less movement options, i.e. smaller $z$. This entropy difference creates an entropic force. This entropic force can even give rise to an attractive force between other objects if it gives the polymers more space to live in---this could e.g. make colloids clump up and then fall out of solution.
