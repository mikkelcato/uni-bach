//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Approximation Methods
== Nondegenerate time-independent pertubation theory
The central problem in quantum mechanics is solving the Schrödinger equation to get the spectrum of a given Hamiltonian, this is hard to do. To make our lives easier we write
$
  H = H_0 + V
$
where $H_0$ is a Hamiltonian for which we know the eigenkets and energies
$
  H_0 ket(n^((0))) = E_n^((0)) ket(n^((0)))
$
and $V$ is interpreted as a pertubation.

We want to solve
$
  (H_0 + V) ket(n) = E_n ket(n)
$
to track the pertubation we write
$
  (H_0 + lambda V) ket(n) = E_n ket(n)
$
where $lambda in [0,1]$ is a continuous parameter we use to keep track---in the limit $lambda -> 0$ we simply get back the free theory. What we really want is the difference in energy defined by
$
  Delta_n equiv E_n - E_n^((0))
$
so we have
$ E_n ket(n) = (Delta_n + E_n^((0))) ket(n) $
then we can rewrite the Schrödinger equation as
$
  (E_n^((0))-H_0) ket(n) = (lambda V - Delta_n) ket(n)
$
to get some insight we try to act with $bra(n^((0)))$ on this guy giving
$
  braket(n^((0)), (lambda V - Delta_n), n) = 0
$
this is nice because it tells us that $(lambda V - Delta_n) ket(n)$ has no component along $ket(n^((0)))$. We consider the operator
$
  phi_n equiv 1 - ketbra(n^((0))) = sum_(k eq.not n) ketbra(k^((0)))
$
which is just the difference between $ket(n^((0)))$ and $ket(n)$ and $phi_n ket(n^((0))) = 0$ by construction. Now we can write
$
  ket(n) = sum_k ket(k^((0))) braket(k^((0)), n) = braket(n^((0)), n) ket(n^((0))) + phi_n ket(n)
$
then
$
  (E_n^((0))-H_0) ket(n) &= underbrace((lambda V - Delta_n) ket(n^((0))), (E_n^((0))-H_0) ket(n^((0)))=0) braket(n^((0)), n) + (lambda V - Delta_n) phi_n ket(n) \
  &= (lambda V - Delta_n) phi_n ket(n)
$
we can then write
$
  ket(n) = 1/(E_n^((0))-H_0) (lambda V-Delta_n) phi_n ket(n)
$
this inversion is possible since we first act with $phi_n$ so we avoid problems like
$
  1/(E_n^((0))-E_n^((0))) = 1/0 -> "ill-defined"
$
now in the limit $lambda -> 0$ we have $Delta_n -> 0$ and $ket(n) -> ket(n^((0)))$ but in our case $ket(n) -> 0$ for this reason write
$
  ket(n) = 1/(E_n^((0))-H_0) (lambda V - Delta_n) phi_n ket(n) + c_n (lambda) ket(n^((0)))
$
notice that if we do $(E_n^((0))-H_0) ket(n)$ then we recover the previous since $(E_n^((0))-H_0) ket(n^((0))) = 0$ meaning this term acts like an arbitrary constant. In the same limit we see $c_n (lambda) -> 1$ and generally $c_n (lambda) = braket(n^((0)), n)$ which we fix by convention
$
  braket(n^((0)), n) = c_n (lambda) = 1
$
so $ket(n)$ is not normalized $braket(n, n) eq.not 1$. We finally obtain
$
  ket(n) = ket(n^((0))) + 1/(E_n^((0))-H_0) phi_n (lambda V- Delta_n) ket(n)
$
where we've switched the order of $phi_n$ and $(lambda V - Delta_n)$, since they commute. We can also easily find
$
  Delta_n = lambda braket(n^((0)), V, n)
$
this is the basic setup. So far everything is exact.

=== Pertubative solution
To solve this we assume that we can expand $ket(n)$ and $Delta_n$ in power series of $lambda$
$
  ket(n) &= ket(n^((0))) + lambda ket(n^((1))) + lambda^2 ket(n^((2))) + dots \
  Delta_n &= underbrace(Delta_n^((0)), 0) + lambda Delta_n^((1)) + lambda^2 Delta_n^((2)) + dots
$
now we can get expressions for the various $Delta_n^((i))$ by using
$
  Delta_n = lambda braket(n^((0)), V, n)
$
and comparing powers of $lambda$. Doing this we find
$
  Delta_n^((1)) & = braket(n^((0)), V, n^((0))) \
  Delta_n^((2)) & = braket(n^((0)), V, n^((1))) \
                & dots \
  Delta_n^((N)) & = braket(n^((0)), V, n^((N-1)))
$
we do the same with
$
  ket(n) = ket(n^((0))) + 1/(E_n^((0))-H_0) phi_n (lambda V-Delta_n) ket(n)
$
by substituting the expansions for $Delta_n$ and $ket(n)$ and comparing powers (and using $phi_n ket(n^((0)))=0$) we obtain $ket(n^((1)))$
$
  ket(n^((1))) = 1/(E_n^((0))-H_0) phi_n V ket(n^((0)))
$
which gives us $Delta_n^((2))$
$
  Delta_n^((2)) & = braket(n^((0)), V phi_n/(E_n^((0))-H_0) V, n^((0)))
$
which gives us $ket(n^((2)))$
$
  ket(n^((2))) = phi_n/(E_n^((0))-H_0) V phi_n/(E_n^((0))-H_0) V ket(n^((0))) - phi_n/(E_n^((0))-H_0) braket(n^((0)), V, n^((0))) phi_n/(E_n^((0))-H_0)V ket(n^((0)))
$
and we can keep going like this.

In practice these can be rewritten,
$
  Delta_n &= lambda Delta_1 + lambda^2 Delta_2 + dots \
  &= lambda V_(n n) + lambda^2 braket(n^((0)), V 1/(E_n^((0))-H_0) phi_n V, n^((0))) + dots \
  &= lambda V_(n n) + lambda^2 sum_(k eq.not n) braket(n^((0)), V 1/(E_n^((0))-H_0), k^((0))) braket(k^((0)), V, n^((0))) + dots \
  &= lambda V_(n n) + lambda^2 sum_(k eq.not n) abs(V_(n k))^2/(E_n^((0))-E_k^((0))) + dots
$
where we've defined
$
  V_(n k) equiv braket(n^((0)), V, k^((0)))
$
in a similar manner
$
  ket(n) &= ket(n^((0))) + lambda ket(n^((1))) + dots \
  &= ket(n^((0))) + lambda sum_(k eq.not n) V_(k n)/(E_n^((0))-E_k^((0))) ket(k^((0))) + dots
$
which are easy to work with.

Now it is obvious why this fails in the degenerate case, since then $E_n^((0)) = E_k^((0))$ making both of these ill-defined. This is why we want degenerate pertubation theory.

== Degenerate time-independent pertubation theory
