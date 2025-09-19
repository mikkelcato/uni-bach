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
    Notes on soft matter based on lectures given at SDU---and very loosely on _Soft Matter Physics_ by Masao Doi.
  ],
)

= Introduction

== Entropy
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

=== Polymer example
Entropy is one of the big driving forces behind the dynamics we'll analyze in this course. Take a simple polymer as an example.

A polymer can be treated as a random walk. Given a single state has a volume $V_0 = lambda_"th"^3$ then the amount of states in a volume $V$ is just $V\/V_0$. How many polymers can we make? If for every step we have $z$ options, and our polymer is $N$ units long, i.e. we take $N$ random steps to create it, then
$
  Omega_1 = V/V_0 z^N
$
so the entropy is
$
  S_1 = k_B ln V/V_0 + k_B N ln z
$
the first part is the usual translational entropy, but now we also get an extra term do to conformation. This is essential for polymer dynamics since every system seeks maximal entropy. Consider a free polymer and a polymer next to a wall, in this case $S_"wall" < S_"free"$ since the polymer near the wall has less movement options, i.e. smaller $z$. This entropy difference creates an entropic force. This entropic force can even give rise to an attractive force between other objects if it gives the polymers more space to live in---this could e.g. make colloids clump up and then fall out of solution.

#pagebreak()
= Solutions
Many soft matter systems exist as solutions---we are mainly interested in the solute, but the solvent can change how the solute interacts with itself.

== Thermodynamics
=== Free energy
The thermodynamic state of a solution (two-component) can be specified by the temperature $T$, the pressure $P$, and the numbers $N_p$ and $N_s$ which correspond to the amount of solute and solvent molecules respectively. We can then write
$
  G equiv G(N_p, N_s, T, P)
$
we define the weight concentration (of solute)
$
  c equiv (m_p N_p)/V
$
and the volume fraction
$
  phi = (v_p N_p)/(v_p N_p + v_s N_s)
$
where the specific volume is
$
  v_i = (pdv(V, N_i))_(T, P, N_j)" with "i, j = {p,s}
$
by definition these satisfy $V = v_p N_p + v_s N_s$. To see this consider increasing $N_i$ by a factor $alpha$, then by homogeneity:
$
  alpha V(N_p, N_s, T, P) &= (alpha N_p, alpha N_s, T, P) \
  V(N_p, N_s)&=dv(, alpha) V(alpha N_p, alpha N_s) \
  V(N_p, N_s) &= pdv(V, N_p) pdv((alpha N_p), alpha) + pdv(V, N_s) pdv((alpha N_s), alpha) \
  V &= v_p N_p + v_s N_s
$
we assume $v_p$ and $v_s$ are constant (incompressible)---in this case
$
  c = m_p/v_p phi = rho_p phi
$
now:
$
  V = pdv(G, P) => G(N_p,N_s,T,P)= P V + F(N_p, N_s, T)
$
now since $F(N_p,N_s,T)$ is extensive
/*
$
                          F(alpha N_p, alpha N_s, T) & = alpha F(N_p, N_s, T) \
  alpha = v_p\/V => F( (v_p N_p)/V, (v_p N_s)/V, T ) & = v_p/V F(N_p,N_s,T) \
                V/v_p F((v_p N_p)/V, (v_p N_s)/V, T) & = F(N_p, N_s, T) \
                  V/v_p F(phi, v_p/v_s (1 - phi), T) & = F(N_p, N_s, T)
$

since
$
  v_p/v_s (1 - phi) &= v_p/v_s (1 - (v_p N_p)/V) \
  &= v_p/v_s - (v_p^2 N_p)/(v_s V) = (v_p V - v_p^2 N_p)/(v_s V) \
  &= (v_p^2 N_p + v_p v_s N_s - v_p^2 N_p)/(v_s V) \
  &= (v_p N_s)/V
$
*/
we can write
$
  F(N_p, N_s, T)= V f(phi, T)
$
where $f(phi,T)$ is the Helmholtz free energy per unit volume---which is intensive. So we can write
$
  G(N_p, N_s, T, P) = V {P + f(phi,T)}
$
=== Osmotic pressure
/*
Knowing $f(phi,T)$ tells us the behaviour of mixing. If we mix $(V_1, phi_1)$ and $(V_2, phi_2)$ then $V_"final" = V_1 + V_2$ and
$
  phi_"final" = phi = (phi_1 V_1 + phi_2 V_2)/(V_1+V_2) = x phi_1 + (1-x) phi_2 "with" x = V_1/(V_1+V_2)
$
Before mixing the free energy is $V_1 f(phi_1) + V_2 f(phi_2)$ and after it is $(V_1 + V_2) f(phi)$. For the solution to be homogeneous we require
$
           (V_1 + V_2) f(phi) & < V_1 f(phi_1) + V_2 f(phi_2) \
  => f(x phi_1 + (1-x) phi_2) & < x f(phi_1) + (1-x) f(phi_2)
$
if two solutions mix at any $V_1\/V_2$ then this has to be satisfied for all $0<=x<=1$. This is equivalent to $f(phi)$ satisfying
$
  pdv(f, phi, 2) > 0", for " phi_1 < phi < phi_2
$
if for any $phi_1 < phi < phi_2$ $ pdv(f, phi, 2) <0 $ then the system can lower $f(phi)$ by seperating.
*/
The force which tends to mix solute and solvent can be described by osmotic pressure. When a solution is touching solvent across a semi-permeable membrane, solvent will mix with the solution given it lowers the free energy. If the membrane could move all the solute would move into the solution by pushing the membrane---as to create more space---until we have a homogeneous solution. To prevent this a force, the osmotic pressure $Pi$, has to be applied. Consider a system with $F_"tot" (V)$ if we change $V$ by $dd(V)$ then a work $- Pi dd(V)$ is done to the system---which is equivalent to $dd(F_"tot" (V))$, so
$
  Pi = - pdv(F_"tot" (V), V)
$
here $V$ is the amount of solution, the solvent would have volume $V_"tot"-V$. For this reason we can write
$
  F_"tot" = V f(phi) + (V_"tot"-V) f(0)
$
it follows that
$
  Pi & = - pdv(, V) (V f(phi) + (V_"tot"-V) f(0)) \
     & = - f(phi) - V pdv(f(phi), phi) pdv(phi, V) + f(0) \
     & = - f(phi) - V pdv(f(phi), phi) (-1/V phi) + f(0) \
  Pi & = - f(phi) + phi pdv(f(phi), phi) + f(0) \
     & = phi^2 pdv(, phi)(f/phi) + f(0)
$
osmotic pressure can be treated as a measure of how strongly solute and solvent mix.

The mixing force can also be interpreted in terms of the chemical potential $mu$, defined as
$
  mu_p equiv (pdv(G, N_p))_(N_s,T,P)"    " mu_s equiv (pdv(G, N_s))_(N_p,T,P)
$
using
$
  G = (N_p v_p + N_s v_s) {P + f(phi,T)}" and "phi = (N_p v_p)/(N_p v_p + N_s v_s)
$
we can find
$
  mu_p &= pdv(, N_p) (N_p v_p + N_s v_s){P + f(phi, T)} \
  &= v_p {P + f(phi,T)} + (N_p v_p+N_s v_s) pdv(, N_p) {P+f(phi,T)} \
  &= v_p {P+f(phi,T)} + (N_p v_p + N_s v_s) pdv(f, phi) pdv(phi, N_p) \
  &= v_p {P+f(phi,T)} + (N_p v_p + N_s v_s) pdv(f, phi) (v_p/(N_p v_p + N_s v_s) - (N_p v_p^2)/(N_p v_p + N_s v_s)^2) \
  &= v_p {P + f} + (N_p v_p + N_s v_s) pdv(f, phi) ((N_p v_p^2 + N_s v_s v_p - N_p v_p^2)/(N_p v_p + N_s v_s)^2) \
  &= v_p {P+f} + v_p pdv(f, phi) (N_s v_s)/V \
  mu_p &= v_p (P + f(phi, T) + (1-phi) pdv(f, phi))
$
similarly
$
  mu_s & = pdv(, N_s) (N_p v_p + N_s v_s){P + f(phi,T)} \
       & = v_s {P + f(phi,T)} + V pdv(, N_s) {P + f(phi,T)} \
       & = v_s {P + f(phi,T)} + V pdv(f, phi) pdv(phi, N_s) \
       & = v_s {P + f(phi,T)} + V pdv(f, phi) (- (N_p v_p v_s)/V^2) \
       & = v_s (P + f(phi,T) - phi pdv(f, phi))
$
or using osmotic pressure
$
  mu_s (phi, T, P) & = v_s {P-Pi(phi, T)+f(0,T)} \
                   & = v_s {P-Pi(phi, T)} + mu_s^((0)) (T)
$
we can also differentiate these expressions to get
$
  pdv(mu_p, phi) & = pdv(, phi) v_p {P + f(phi,T) + (1-phi) pdv(f, phi)} \
                 & = v_p {pdv(f, phi) + (1-phi) pdv(f, phi, 2) - pdv(f, phi) } \
                 & = v_p (1-phi) pdv(f, phi, 2)
$
and
$
  pdv(mu_s, phi) & = pdv(, phi) v_s {P + f(phi,T) -phi pdv(f, phi)} \
                 & = v_s {pdv(f, phi) - pdv(f, phi) - phi pdv(f, phi, 2)} \
                 & = - v_s phi pdv(f, phi, 2)
$
so if
$
  pdv(f, phi, 2) > 0", for " phi_1 < phi < phi_2
$
then $dd(mu_p, d: partial)\/dd(phi, d: partial)$ is always positive, and $dd(mu_s, d: partial)\/dd(phi, d: partial)$ is always negative. So if there is a potential gradient, solute molecules will move from high to low, while solvent molecules will move from low to high. This happens until $mu_p$ increases to some common $mu$ while $mu_s$ decreases---this is diffusion.


In a dilute solution we can ignore the interaction between solute molecules, and $Pi$ is given by van't Hoff's law:
$
  Pi = (N_p k_B T)/V = (phi k_B T)/v_p
$
the self-interaction leads to correction terms
$
  Pi = (phi k_B T)/v_p + A_2 phi^2 + A_3 phi^3 + dots
$
we can write
$
  Pi = phi^2 pdv(, phi) (f/phi) + f(0)
$
and integrate to find $f(phi)$
$
  (phi k_B T)/v_p + A_2 phi^2 + dots &= phi^2 pdv(, phi) (f/phi) + f(0) \
  (k_B T)/(v_p phi) + A_2 + dots &= pdv(, phi) (f/phi) + f(0)/phi^2 \
  (k_B T)/v_p ln phi + A_2 phi + dots &= f/phi - f(0)/phi - k_0 \
  => f(phi) &= f_0 + k_0 phi + (k_B T)/v_p phi ln phi + A_2 phi^2 + A_3/2 phi^3 + dots
$

using this we obtain
$
  mu_p (phi) &= v_p {P + f(phi,T) + (1- phi) pdv(f, phi)} \
  &= v_p P + v_p f_0 + v_p k_0 phi + k_B T phi ln phi + A_2 v_p phi^2 + A_3/2 v_p phi^3 + \
  &+ (1-phi) (v_p k_0 + k_B T ln phi + k_B T + 2 A_2 v_p phi + 3/2 A_3 v_p phi^2 ) + dots \
  &= v_p P + v_p f_0 + v_p k_0 phi + k_B T phi ln phi + A_2 v_p phi^2 + v_p k_0 + k_B T ln phi +k_B T + 2 A_2 v_p phi \
  &- v_p k_0 phi - k_B T phi ln phi - k_B T phi - 2 A_2 v_p phi^2 + A_3/2 v_p phi^3 + 3/2 A_3 v_p phi^2 - 3/2 A_3 v_p phi^3 \
  &= v_p P + v_p f_0 - A_2 v_p phi^2 + v_p k_0 + k_B T ln phi +k_B T + 2 A_2 v_p phi - k_B T phi + 3/2 A_3 v_p phi^2 + dots \
  &= mu_p^((0)) + P v_p + k_B T ln phi + v_p {(2 A_2-(k_B T)/v_p) phi + (3/2 A_3 -A_2) phi^2 + dots}
$
where $mu_p^((0)) = v_p (f_0 +k_0) + k_B T equiv mu_p^((0)) (T)$. Similarly
$
  mu_s (phi) &= v_s {P + f(phi,T) - phi pdv(f, phi)} \
  &= v_s P + v_s f_0 + v_s k_0 phi + v_s/v_p k_B T phi ln phi + A_2 v_s phi^2 + A_3/2 v_s phi^3 + dots \
  &- v_s k_0 phi - v_s/v_p k_B T phi ln phi - v_s/v_p k_B T phi - 2 A_2 v_s phi^2 - 3/2 A_3 v_s phi^3 - dots \
  &= v_s f_0 +v_s P - v_s/v_p k_B T phi - v_s A_2 phi^2 - v_s A_3 phi^3 \
  &= mu_s^((0)) (T) + v_s P - v_s/v_p k_B T phi - v_s (A_2 phi^2 + A_3 phi^3 + dots)
$



=== Phases
Knowing $f(phi,T)$ tells us the behaviour of mixing. If we mix $(V_1, phi_1)$ and $(V_2, phi_2)$ then $V_"final" = V_1 + V_2$ and
$
  phi_"final" = phi = (phi_1 V_1 + phi_2 V_2)/(V_1+V_2) = x phi_1 + (1-x) phi_2 "with" x = V_1/(V_1+V_2)
$
Before mixing the free energy is $V_1 f(phi_1) + V_2 f(phi_2)$ and after it is $(V_1 + V_2) f(phi)$. For the solution to be homogeneous we require
$
           (V_1 + V_2) f(phi) & < V_1 f(phi_1) + V_2 f(phi_2) \
  => f(x phi_1 + (1-x) phi_2) & < x f(phi_1) + (1-x) f(phi_2)
$
if two solutions mix at any $V_1\/V_2$ then this has to be satisfied for all $0<=x<=1$. This is equivalent to $f(phi)$ satisfying
$
  pdv(f, phi, 2) > 0", for " phi_1 < phi < phi_2
$
if for any $phi_1 < phi < phi_2$ $ pdv(f, phi, 2) <0 $ then the system can lower $f(phi)$ by seperating. We assume this happens due to a change of pressure or temperature---then our solution will seperate into a high and low concentration region. We say $(V, phi) -> (V_1, phi_1) "and" (V_2, phi_2)$, with $V = V_1 + V_2$ and $phi V = phi_1 V_1 + phi_2 V_2$ giving
$
  V_1 = (phi_2 - phi)/(phi_2 - phi_1) V "and" V_2 = (phi-phi_1)/(phi_2-phi_1) V
$
so
$
  F = V_1 f(phi_1) + V_2 f(phi_2) = V {(phi_2-phi)/(phi_2-phi_1) f(phi_1) + (phi -phi_1)/(phi_2-phi_1) f(phi_2)}
$
note the term in ${dot}$ corresponds to the line between $P_1 (phi_1, f(phi_1))$ and $P_2 (phi_2,f(phi_2))$. Therefore to minimize $F$ we want this line to be the common tangent for $f(phi)$. Let the tangent points be $P_a$ and $P_b$, then a solution with $phi_a < phi < phi_b$ is the most stable if it seperates into two solutions with $phi_a$ and $phi_b$. The condition for $P_a P_b$ being the common tangent is
$
  f'(phi_a)=f'(phi_b)",  " f(phi_a) - f'(phi_a) phi_a = f(phi_b) - f'(phi_b) phi_b
$

This is equivalent to the chemical potentials being the same. We can subdivide the region $phi_a < phi < phi_b$ into unstable and metastable regions. Let $phi^*_b$ be the point where
$
  pdv(f, phi, 2) = 0
$
then for $phi^*_b < phi < phi_b$ is called metastable, since small deviations from $phi$ in this range can keep a solution stable. Similarly for $phi_a < phi < phi^*_a$. In $phi^*_a < phi < phi^*_b$ the solution is completely unstable. All of these $phi(T)$ can be plotted together giving a phase diagram---see pg. 14-15. Notable we have the coexistence- and spinodal curve which meet at the critical point where
$
  pdv(f, phi, 2) = pdv(f, phi, 3) = 0
$

\* spinodal decomposition ($f''<0$), nucleation & growth ($f'' > 0$)---always want $mu_a = mu_b$, "trivially" done for spinodal, harder for systems in $phi_a < phi < phi^*_a ->$ nucleation & growth.

#pagebreak()
== Regular solution theory
We'll calculate $f(phi)$ using a simple model. We assume our solution is a lattice, where every point is either occupied by a solute or solvent molecule. We assume both have the same volume $v_c$. Then
$
  V = v_c N_"tot"",  " phi = N_p/N_"tot"
$
with $N_"tot" = N_p+N_s$. The energy for the configuration $i$ can be written as
$
  E_i = epsilon_"pp" N_i^(("pp")) + epsilon_"ss" N_i^(("ss")) + epsilon_"ps" N_i^(("ps"))
$
where $epsilon$ are interaction energies. The partition function is
$
  Z = sum_i exp((-E_i)/(k_B T))
$
and
$
  f(phi,T) = (- k_B T ln Z)/V
$
to make our lives easier we do---mean field approximation
$
  Z tilde.eq W exp((- expval(Delta E_c))/(k_B T))
$
with $expval(E)$ being the average energy, and $W$ being the amount of configurations, so how many ways we can place $N_p$ molecules on $N_"tot"$ cells. We start by finding $expval(Delta E_c)$
$
  expval(Delta E_c) = epsilon_"pp" expval(N_"pp") + epsilon_"ss" expval(N_"ss") + epsilon_"ps" expval(N_"ps") - E_"pure"
$
each cell has $z$ neighbors, on average $z phi$ cells are solute while $z (1-phi)$ are solvent. So we can write
$
  expval(N_"pp") = (z phi N_p)/2 = (z N_"tot" phi^2)/2
$
and similarly
$
  expval(N_"ss") & = N_s (z(1-phi))/2 = (z N_"tot" (1-phi)^2)/2 \
  expval(N_"ps") & = N_p z(1-phi)=N_s z phi = z N_"tot" phi(1-phi)
$
and
$
  E_"pure" = epsilon_"pp" expval(N^"pure"_"pp") + epsilon_"ss" expval(N^"pure"_"ss")
$
in both cases $phi = 1$ so
$
  expval(N^"pure"_"pp") = (N_p z)/2 "   " expval(N^"pure"_"ss") = (N_s z)/2
$
we get (using $phi = N_p\/N_"tot"$ and $1-phi = N_s\/N_"tot"$)
$
  expval(Delta E_c) &= epsilon_"pp" (N_p z phi )/2 + epsilon_"ps" N_p z (1-phi) + epsilon_"ss" (N_s z (1-phi))/2 - epsilon_"pp" (N_p z)/2 - epsilon_"ss" (N_s z)/2 \
  &= (N_"tot" z)/2 (epsilon_"pp" (phi^2-phi) + 2 epsilon_"ps" phi (1-phi) + epsilon_"ss" ((1-phi)^2 - 1+phi)) \
  &= (N_"tot" z)/2 (-epsilon_"pp" phi (1-phi) + 2 epsilon_"ps" phi(1-phi) - epsilon_"ss" phi (1-phi)) \
  &= N_"tot" (z )/2 (2 epsilon_"ps" - epsilon_"pp" - epsilon_"ss") phi (1-phi) \
  &= N_"tot" k_B T chi phi(1-phi)
$
where we define (which is nice and unit-less)
$
  chi equiv z/(2 k_B T) (2 epsilon_"ps"-epsilon_"pp"-epsilon_"ss")
$
now we need the amount of configurations $W$. This is just $N_"tot" = (N_p + N_s) "choose" N_p$, and since we need the logarithm we do
$
  ln (W) &= ln ((N_p + N_s)!/(N_p! N_s!)) \
  &= ln((N_p + N_s)!) - ln(N_p !) - ln(N_s !) \
  &tilde.eq^"stirling" (N_p+N_s) ln (N_p + N_s)-N_p - N_s - N_p ln N_p + N_p - N_s ln N_s + N_s \
  &= (N_p + N_s) ln (N_p+N_s) - N_p ln N_p - N_s ln N_s \
  &= N_p ln (N_p + N_s)/N_p + N_s ln (N_p + N_s)/N_s \
  &= - N_p ln N_p/(N_p + N_s) - N_s ln N_s/(N_p+N_s) \
  &= -N_"tot" ( phi ln phi + (1-phi) ln (1-phi) )
$
now we obtain:
$
  F &= - k_B T ln Z \
  &= - k_B T ln (W e^(- expval(Delta E_c)/(k_B T))) \
  &= - k_B T (ln W - expval(Delta E_c)/(k_B T)) \
  &= - k_B T ln W + expval(Delta E_c) \
  &= N_"tot" k_B T (phi ln phi + (1-phi) ln(1-phi)) + N_"tot" k_B T chi phi(1-phi) \
  &= N_"tot" k_B T (phi ln phi + (1- phi) ln(1-phi) + chi phi(1-phi))
$
we want $f(phi,T)$,
$
  f(phi,T) &= F/V = F/(N_"tot" v_c) \
  &= (N_"tot" k_B T)/(N_"tot" v_c) (phi ln phi + (1-phi) ln(1-phi) + chi phi(1-phi)) \
  &= (k_B T)/v_c (phi ln phi + (1-phi) ln(1-phi) + chi phi(1-phi))
$
where everything in the brackets is unitless, and $k_B T\/v_c$ is an energy density. The first two terms correspond to the Gibbs' entropy---and we call it the configurational entropy (or entropy of mixing), from ($H=0$)
$
  F = - k_B T sum P_i ln P_i
$
The $chi$ term is the interaction term or the enthalpy, and we can write
$
  chi = chi_H - T chi_S
$
corresponding to $F = H - T S$, evidently at high temperature entropy dominates, while at low temperature the actual interaction between molecules (enthalpy) dominates.



#pagebreak()
== Polymers
Torsion angles---every angle $phi$ has different energies, three favorable ${g^-,0,g^+}$, on the scale of $k_B T$---at sufficient temperature this becomes a random walk.

=== Flory-Huggins
The previous analysis (regular solution theory) can be extended to polymers. Instead of seperate solute molecules $p$ we assume they are connected and move together like a random walk---given there are $N$ steps per polymer we reduce the translational entropy by $N^(-1)$:
$
  f(phi) = (k_B T)/v_c (phi/N ln phi + (1-phi) ln(1 - phi) + chi phi(1-phi))
$
we essentially say we know where every connected $p$ is knowing where one is.

We have a typical size of a polymer given by
$
  R_g^2 = (b^2 N)/6 => rho^* = 1/(R_g^2)^(3\/2)
$
where $rho^*$ is the critical density. Flory-Huggins works for $rho >= rho^*$ since in this case the mean field approximation is not bad. In the extreme case the polymers essentially drown out and we return to a regular solution.

=== Random walk
We assume steps $arrow(b)_n$ are statistically independent---so
$
                 expval(arrow(b)_n) & = 0 \
  expval(arrow(b)_n dot arrow(b)_m) & = expval(b^2) delta_(n m)
$
the length of the entire polymer (contour length) is
$
  L = N b
$
the length from start to end is given by the vector $arrow(R)$, this is a stochastic quantity so
$
  expval(arrow(R)) = expval(sum_(n=1)^N arrow(b)_n) = sum_(n=1)^N expval(arrow(b)_n) = 0
$
so this is a bad quantity---what is the length
$
  expval(arrow(R)dot arrow(R)) &= expval(sum_(n=1)^N arrow(b)_n dot sum_(m=1)^N arrow(b)_m) \
  &= expval((arrow(b)_1+arrow(b)_2+dots+arrow(b)_N) dot (arrow(b)_1+arrow(b)_2+dots+arrow(b)_N)) \
  &= expval(arrow(b)_1 dot arrow(b)_1 + arrow(b)_2 dot arrow(b)_2 + dots + arrow(b)_N dot arrow(b)_N) \
  &= expval(sum_(n=1)^N b_n^2) = sum_(n=1)^N expval(b^2) = N b^2
$
so $d tilde sqrt(expval(arrow(R)dot arrow(R))) = b sqrt(N)$---size in space. So
$
  d/L = 1/(sqrt(N))
$
with $b tilde Å <->"nm"$ while $N$ varies wildly $tilde 10^2 <-> 10^6$.

Now we are interested in the probability distribution, to get this we treat every step as a spring with spring constant $k$---then
$
  H = 1/2 k (x^2 + y^2 + z^2) = 1/2 k (arrow(b)_n dot arrow(b)_n)
$
with $arrow(b)_n = vecrow(x, y, z)$. We obtain
$
  expval(H) & = 1/2 k expval(arrow(b)_n dot arrow(b)_n) \
  3/2 k_B T & = 1/2 k expval(b_n^2) => k = (3 k_B T)/expval(b^2)
$
with this spring constant we get the same step size as for the random walk---incorporating springs let us bring dynamics into a random walk, since a random walk is static. We can find
$
  expval(arrow(R) dot arrow(R)) = (3 k_B T N)/k
$
we can ask
$
  P(arrow(b)_n) & prop exp((-H)/(k_B T)) \
  & prop exp(- k/2 (arrow(b)_n dot arrow(b)_n)/(k_B T)) \
  & prop exp(-(3 x^2)/(2expval(b^2))) exp(-(3y^2)/(2expval(b^2))) exp(-(3z^2)/(2expval(b^2)))
$
For a polymer chain we have $k_"chain" = k\/N$, so one spring corresponds to our spring steps, where each step has a $k$ such that it is equivalent to the random walk---so we could write
$
  expval(arrow(R) dot arrow(R)) = (3 k_B T)/k_"chain"
$
the Hamiltonian of the chain is
$
          H_"chain" & = 1/2 k_"chain" (R_x^2 + R_y^2 + R_z^2) \
  expval(H_"chain") & = 1/2 k_"chain" expval(R_x^2+R_y^2+R_z^2) \
        (3 k_B T)/2 & = 1/2 k_"chain" expval(arrow(R) dot arrow(R)) \
                    & => expval(arrow(R) dot arrow(R)) = (3 k_B T)/k_"chain"
$
so we recover the result.

For $P(arrow(R))$ we can write
$
  P(arrow(R)) & prop exp(- H_"chain"/(k_B T)) \
  & prop exp(- k_"chain"/2 (arrow(R) dot arrow(R))/(k_B T)) \
  & prop exp(- k_"chain"/2 (R_x^2 + R_y^2 + R_z^2)/(k_B T)) \
  & prop exp(- 3/2 (R_x^2 + R_y^2 + R_z^2)/(expval(arrow(R)dot arrow(R))))
$
so everything in the chain is Gaussian, both individual steps and $arrow(R)$---here $sigma = expval(arrow(R) dot arrow(R))\/3$, as is expected since every direction is equivalent. The exact distribution is
$
  P(arrow(R)) = (3/(2 pi expval(arrow(R) dot arrow(R))))^(3\/2) exp(- 3/2 (R_x^2+R_y^2+R_z^2)/(expval(arrow(R) dot arrow(R))))
$

Consider pulling a polymer with some $f_"ext"$ then it will resist with an entropic force. In one-dimension
$
  P(R_x) prop exp(- 3/2 R_x^2/(2 expval(arrow(R) dot arrow(R)))) prop Omega(R_x)
$
with $Omega (R_x)$ being the number of microstates,
$
  S(R_x) & = k_B ln Omega(R_x) \
         & = k_B ln ("const" exp[-3/2 R_x^2/(expval(arrow(R) dot arrow(R)))] ) \
         & = - k_B (3 R_x^2)/(2 expval(arrow(R) dot arrow(R))) + S_0
$
now
$
  F_"poly" = U - T S & = -T S \
                     & = k_B T (3 R_x^2)/(2 expval(arrow(R)dot arrow(R))) + F_0 \
                     & = (3 k_B T)/(2) R_x^2/(expval(arrow(R)dot arrow(R)))
$
giving
$
  F_"system" & = (3 k_B T)/2 R_x^2/(expval(arrow(R) dot arrow(R))) - f_"ext" R_x \
           0 & =^"equil." pdv(F_"system", R_x) \
           0 & = (3 k_B T)/expval(arrow(R) dot arrow(R)) R_x - f_"ext" \
     f_"ext" & = (3 k_B T R_x)/expval(arrow(R)dot arrow(R)) \
    f_"poly" & = - f_"ext"
$
less configurations as it gets stretched means it gets lower entropy leading to a free energy and a force.

/*
so
$
  overline(E) &= epsilon_"pp" overline(N)_"pp" + epsilon_"ss" overline(N)_"ss" + epsilon_"ps" overline(N)_"ps" \
  &= 1/2 N_"tot" z Delta epsilon phi^2 + C_0 + C_1 phi
$
with $Delta epsilon = epsilon_"pp" + epsilon_"ss" - 2 epsilon_"ps"$ being the effective interaction energy, and ${C_0,C_1}$ being constants.

The free energy can be written as
$
  F = - k_B T ln Z = - k_B T ln W + overline(E)
$
here the first term represents entropy, with the mixing entropy being given by
$
  S = k_B ln W & = k_B (ln (N_p+N_s)! - ln N_p ! -ln N_s !) \
               & = k_B ((N_p+N_s) ln(N_p+N_s) - N_p ln N_p -N_s ln N_s) \
               & = k_B (N_p ln (N_p + N_s)/N_p + N_s ln (N_p + N_s)/N_s) \
               & = k_B N_"tot" (- phi ln phi - (1-phi) ln (1-phi))
$
so
$
  F = N_"tot" [k_B T (phi ln phi + (1-phi) ln (1-phi)) + z/2 Delta epsilon phi^2]
$
here we've dropped the constant and linear terms since they don't matter. So the free energy density becomes
$
  f(phi) = (k_B T)/v_c (phi ln phi + (1-phi) ln (1-phi) + chi phi (1-phi))
$
where we've added a linear term (which doesn't matter) and defined
$
  chi = - (z Delta epsilon)/(2 k_B T)
$
the osmotic pressure is
$
  Pi = (k_B T)/v_c (- ln(1-phi) - chi phi^2) =^(phi << 1) (k_B T)/v_c (phi + (1/2 - chi) phi^2)
$
here the first term is van't Hoff's law and the second corresponds to $A_2$---if $A_2 > 0$ then the interaction between solute molecules is repulsive, $chi$ represents the energetic interaction and for $chi > 0$ this is attractive, while the $1\/2$ comes from the molecules actually taking up space on our lattice.

Our $f(phi)$ has a mirror-symmetry at $phi=1\/2$. If $chi < chi_c$ then $f(phi)$ has one minimum at $phi=1\/2$, and for $chi > chi_c$ it has two minima. The value of $chi_c$ is determined by
$
  pdv(f, phi, 2)=0 "at" phi=1/2 => chi_c = 2
$
the spinodal line is
$
  pdv(f, phi, 2) = 0 => chi = 1/2 1/(phi(1-phi))
$
a similar construction as the phase diagram can be constructed by finding the common tangent, in this case it connects the two minima of $f(phi)$---so the $phi_a (T)$ and $phi_b (T)$ in the coexistence region are given by
$
  pdv(f, phi)=0 => chi = 1/(1-2 phi) ln (1-phi)/phi
$

=== Polymer solutions
Before we assumed solute and solvent molecules had the same size---for polymer solutions this is not true. We model a polymer as N segments connected by bonds---with each segment being a point on the lattice. We assume that all segments and solvent molecules have the same size.

In this case the free energy becomes
$
  f(phi) = (k_B T)/v_c (1/N phi ln phi + (1-phi)ln(1-phi) + chi phi(1-phi))
$
where the factor $N^(-1)$ appears since the segments cannot be placed independently.

correlation effect\*
*/

#pagebreak()
= Elastic soft matter
Elastic soft matter could be rubber or gel. The main property of these materials is that we can deform them, and after they'll recover they shape---these typically consist of polymers, in solution (cross-links) or not in solution.

== Basic definitions
We can define shear stress $sigma$ and shear strain $gamma$ as
$
  sigma = F/S
$
and
$
  gamma = d/h
$
for the setup on pg. 29. If the material is an ideal elastic then
$
  sigma = G gamma
$
with $G$ being the shear modulus. If the material is an ideal viscous fluid then
$
  sigma = eta dot(gamma)
$
with $eta$ being the viscocity.

If we compress a material isotropically using some pressure $Delta P$, and the response is $V arrow V + Delta V$ then we can write
$
  Delta P = - K (Delta V)/V
$
with $K$ being the bulk modulus.

If we stretch some material uniformly $L arrow L'$ then we define $lambda = L'\/L$ and $epsilon = lambda - 1$. The elongational stress $sigma$ can be written as
$
  sigma = E epsilon
$
with $E$ being Young's modulus. We can write
$
  E = (9 K G)/(3 K + G)
$

The bulk and shear moduli both characterize how strongly a material wants to preserve its shape. We'll assume that the bulk modulus $K$ is infinite, since softness is due to a low shear modulus $G$---in this case $E = 3 G$.

== Polymer chains

== Kuhn's theory

