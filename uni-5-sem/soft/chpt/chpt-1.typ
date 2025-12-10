//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Introduction

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
