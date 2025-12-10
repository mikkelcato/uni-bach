//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

#pagebreak()
= Polymers
Before discussing polymers and polymer solutions proper we discuss _normal_ solutions.

== Regular solution theory
We want to compute $f(phi)$. We model our solution as a lattice. Any point on this lattice is then either occupied by a solute or solvent molecule. Assuming both molecules have the same volume $v_c$ then
$
  V = v_c N_"tot"",  " phi = N_p/N_"tot"
$
with $N_"tot" = N_p+N_s$. The energy for some specific configuration $i$ can be written as
$
  E_i = epsilon_"pp" N_i^(("pp")) + epsilon_"ss" N_i^(("ss")) + epsilon_"ps" N_i^(("ps"))
$
where the $epsilon$ and $N_i$ have obvious meanings. The partition function is by definition
$
  Z & = sum_i exp((-E_i)/(k_B T)) \
    & tilde.eq^"mean field" W exp((- expval(Delta E_c))/(k_B T))
$
with $W$ being the total number of configurations.

We can find $expval(Delta E_c)$ by
$
  expval(Delta E_c) = epsilon_"pp" expval(N_"pp") + epsilon_"ss" expval(N_"ss") + epsilon_"ps" expval(N_"ps") - underbracket(E_"pure", "background energy")
$
Taking any cell to have $z$ neighbors then $z phi$ of these are solute while $z (1-phi)$ are solvent. So we obtain
$
  expval(N_"pp") = z phi underbracket(N_p/2, "amount of pairs") = (z N_"tot" phi^2)/2
$
similarly
$
  expval(N_"ss") = (z N_"tot" (1-phi)^2)/2";  " expval(N_"ps") = z N_"tot" phi(1-phi)
$
The pure energy is given by
$
  E_"pure" = underbrace(epsilon_"pp" expval(N^"pure"_"pp"), "all solute") + underbrace(epsilon_"ss" expval(N^"pure"_"ss"), "all solvent")
$
these are trivial
$
  expval(N^"pure"_"pp") = (N_p z)/2";  " expval(N^"pure"_"ss") = (N_s z)/2
$
We obtain
$
  expval(Delta E_c) &= epsilon_"pp" (N_p z phi )/2 + epsilon_"ps" N_p z (1-phi) + epsilon_"ss" (N_s z (1-phi))/2 - epsilon_"pp" (N_p z)/2 - epsilon_"ss" (N_s z)/2 \
  &= (N_"tot" z)/2 (2 epsilon_"ps" - epsilon_"pp" - epsilon_"ss") phi (1-phi) \ &= N_"tot" k_B T chi phi(1-phi)
$
where we define $chi$ as
$
  chi equiv z/(2 k_B T) (2 epsilon_"ps"-epsilon_"pp"-epsilon_"ss")
$
this is a unitless measure of the interaction.

We still need $W$. We compute
$
  ln W & = ln (N_p + N_s)!/(N_p! N_s!) tilde.eq^"stirling" - N_p ln N_p/(N_p + N_s) - N_s ln N_s/(N_p+N_s) \
  & = -N_"tot" ( phi ln phi + (1-phi) ln (1-phi) )
$
Finally we find
$
  F & = - k_B T ln Z = - k_B T ln W + expval(Delta E_c) \
    & = N_"tot" k_B T (phi ln phi + (1- phi) ln(1-phi) + chi phi(1-phi))
$
or as a density
$
  f(phi,T) = F/(N_"tot" v_c) = (k_B T)/v_c [underbracket(phi ln phi, "translational entropy") + underbracket((1-phi) ln(1-phi), "mixing entropy") + underbracket(chi phi(1-phi), "enthalpy")]
$
which is quite nice! The term in $[dots]$ is unitless and $k_B T v_c^(-1)$ is an energy density. The first two terms correspond to the Gibbs entropy. To see this write
$
  F =^(H = 0) - k_B T underbracket(sum P_i ln P_i, "Gibbs'")
$
The $chi$ term is the interaction or enthalpy term. We can write
$
  chi = chi_H - T chi_S tilde H - T S
$
Then at low temperatures the interaction between molecules dominates. As one would expect.

=== Flory-Huggins theory
The above derivation can be extended to polymer solutions.

We assume polymers are random walks. Then a single polymer corresponds to $p$ solute molecules connected and moving together as a random walk. Given each polymer has $N$ steps we reduce the translational entropy by $N^(-1)$ to obtain
$
  f(phi, T) = (k_B T)/v_c [phi/N ln phi + (1-phi) ln(1 - phi) + chi phi(1-phi)]
$
this is the Flory-Huggins free energy.

The typical size of a polymer is
$
  R_g^2 = (b^2 N)/6
$
We can then define a critical density $rho^* = R_g^(-3)$. For $rho >= rho^*$ the above is a fair approximation, since the mean field approximation is valid. As $rho -> oo$ the polymers drown out and we get back a regular solution.

== As a random walk
As mentioned we assume polymers are random walks. We now show this is a useful assumption.

For a random walk we assume individual steps $bold(b)_n$ are statistically independent. Meaning
$
  expval(bold(b)_n) & = 0";  " expval(bold(b)_n dot bold(b)_m) & = expval(b^2) delta_(n m)
$
We define the contour length of a polymer to be
$
  L_"contour" equiv N b
$
The end-to-end length of a polymer is given by $bold(R)$. We would like some measure of this. Consider the average
$
  expval(bold(R)) = expval(sum_(n=1)^N bold(b)_n) = sum_(n=1)^N expval(bold(b)_n) = 0
$
so this is a bad quantity. Consider instead the variance
$
  expval(bold(R)dot bold(R)) &= expval(sum_(n=1)^N bold(b)_n dot sum_(m=1)^N bold(b)_m) = sum_(n=1)^N expval(b^2) =^(expval(b^2)=b^2) N b^2
$
so the size of a polymer is $d tilde b sqrt(N)$. We can relate this to $L_"contour"$ by
$
  d/L_"contour" = 1/(sqrt(N))
$

We are interested in the probability distribution for $bold(R)$. Treating each step as a spring with spring constant $k$ the Hamiltonian of a step is
$
  H_"step" = 1/2 k (x^2 + y^2 + z^2) = 1/2 k (bold(b)_n dot bold(b)_n)
$
with $bold(b)_n = x hat(x) + y hat(y) + z hat(z)$. Then by the equipartition theorem
$
  3/2 k_B T & = 1/2 k expval(bold(b)_n dot bold(b)_n) = 1/2 k expval(b^2)
$
meaning
$
  k = (3 k_B T)/expval(b^2) = (3 N k_B T )/(expval(bold(R) dot bold(R)))
$
using this spring constant we get the same step size as for the random walk! We use springs since they incorporate dynamics and are easy to work with.

Assuming the Boltzmann distribution is valid then the probability distribution for $bold(b)_n$ is
$
  P(bold(b)_n) & prop exp((-H_"step")/(k_B T)) = exp(-3/2(bold(b)_n dot bold(b)_n)/(expval(b^2)))
$
Similarly the Hamiltonian of a full polymer is
$
  H_"chain" & = 1/2 k_"chain" (bold(R)dot bold(R)) \
$
giving
$
  k_"chain" = (3 k_B T)/expval(bold(R) dot bold(R)) = k/N
$
and by the Boltzmann distribution
$
  P(bold(R)) & prop exp(- H_"chain"/(k_B T)) = exp(- k_"chain"/2 (bold(R) dot bold(R))/(k_B T)) = exp(- 3/2 (bold(R) dot bold(R))/(expval(bold(R)dot bold(R))))
$
normalizing we find
$
  P(bold(R)) &= (3/(2 pi expval(bold(R) dot bold(R))))^(3\/2) exp(- 3/2 (bold(R) dot bold(R))/(expval(bold(R) dot bold(R))))
$
so both individual steps, and the entire chain are Gaussian!

We can now find the free energy of a polymer by
$
  F_"polymer" =^"ideal polymer" - T S
$
since the chain is fully entropic. The entropy is $prop ln W$ with
$
  W prop P(bold(R))
$
so we find
$
  S(bold(R)) &= k_B ln W(bold(R)) \
  &= -(3 k_B)/2 (bold(R) dot bold(R))/(expval(bold(R) dot bold(R))) + S_0 \
$
and we obtain
$
  F_"polymer" &= (3 k_B T)/2 (bold(R) dot bold(R))/(expval(bold(R) dot bold(R))) + F_0 \
  &= (3 k_B T)/(2 N b^2) (bold(R) dot bold(R)) + F_0
$

Consider pulling a polymer with some external force $f_"ext"$ then the total free energy of the system becomes
$
  F_"system" & = F_"polymer" - f_"ext" dot bold(R)
$
At equilibrium we have
$
  0 & = grad_bold(R) F_"system" \
  0 & = (3 k_B T)/(N b^2) bold(R)- f_"ext" => f_"ext" & = (3 k_B T)/(N b^2) bold(R)
$
so the polymer resists with an entropic force given by $f_"polymer" = - f_"ext"$.


== Elasticity
Why we care about elasticity in relation to soft matter should be fairly obvious. Typical examples of elastic soft matter are rubbers and gels. These materials recover their shape after being deformed. We specifically care about elasticity of polymer networks. To describe these materials we need to define stress and strain.

=== Stress and strain
When deforming a body we consider two types of deformations. These are uniaxial tension and shear deformation. Any general deformation can be treated as a superposition of these. To describe deformations we define the Cauchy stress tensor $sigma_(alpha beta)$
$
  sigma_(alpha beta)^((d)) = F_alpha/A_beta
$
with the meaning being obvious. As an example tension could be $sigma_(x x)$ while shear could be $sigma_(x y)$. For the full stress tensor we include ambient pressure
$
  sigma_(alpha beta) equiv - P delta_(alpha beta) + sigma_(alpha beta)^((d))
$
We would like to describe the response of a material upon applying stress. So we want some relation between some old coordinates and the deformed coordinates $r'_alpha (r_beta)$.

Under a simple stretching we would have
$
  r'_alpha = Lambda_(alpha alpha) r_alpha
$
where $Lambda_(alpha alpha)$ is some scaling. The simplest case has $Lambda_(alpha alpha) = L_alpha'\/L_alpha = lambda$ with $L_alpha$ being the length of the body we stretch. Consider stretching an incompressible body along $hat(x)$ then $L'_y =L'_z$ and $V = V'$ giving
$
  1 =^! L'_x/L_x [L'_y/L_y]^2 = lambda_(parallel) (lambda_perp)^2 => lambda_perp = lambda_parallel^(-1\/2)
$
this is is called a volume preserving uniaxial deformation. We obtain $Lambda_(x x) = lambda$ and $Lambda_(y y) = Lambda_(z z) = lambda^(-1\/2)$.

Under a simple shear we define $gamma = tan theta$ with $theta$ being the _shear angle_. Consider shearing along $x$ then $r'_z = r_z$ and $r'_y = r_y$, but $r'_x = r_x + gamma r_y$.

We generalize the above by the deformation gradient tensor $E_(alpha beta)$
$
  E_(alpha beta) equiv pdv(r'_alpha, r_beta)";  " dd(r'_alpha) = E_(alpha beta) dd(r_beta)
$
in the examples
$
  E_(alpha beta)^"stretch" = diagonalmatrix(lambda, lambda^(-1\/2), lambda^(-1\/2))";  " E_(alpha beta)^"shear" = mat(1, gamma, ; , 1, ; , , 1)
$
Consider $dd(s^2) = dd(r_alpha, r_alpha)$ after a deformation
$
  (dd(s)')^2 & = dd(r'_alpha, r'_alpha) \
             & = E_(alpha gamma) E_(alpha delta) dd(r_gamma, r_delta) \
             & equiv C_(gamma delta) dd(r_gamma, r_delta)
$
where we define the right Cauchy-Green tensor $C_(alpha beta)$
$ C_(alpha beta) equiv (E^T E)_(alpha beta) =E^T_(alpha gamma) E_(gamma beta) $
as an example the above give
$
  C_(alpha beta)^"stretch" = mat(lambda^2, 0, 0; 0, lambda, 0; 0, 0, lambda)";  " C_(alpha beta)^"shear" = mat(1, gamma, 0; gamma, 1+ gamma^2, 0; 0, 0, 1)
$

Consider the extension
$
  (dd(s'))^2-(dd(s))^2 &= (C_(alpha beta)-delta_(alpha beta)) dd(r_alpha, r_beta) \
  &equiv 2 cal(E)_(alpha beta) dd(r_alpha, r_beta)
$
where we define the Lagrangian strain tensor $cal(E)_(alpha beta)$
$
  cal(E)_(alpha beta) equiv 1/2 (E_(gamma alpha) E_(gamma beta) - delta_(alpha beta))
$
this is nice because doing nothing gives $cal(E)_(alpha beta) = 0$.

=== Hooke's law
We want a relationship between the stress we apply $sigma_(alpha beta)$ and the strain $cal(E)_(alpha beta)$. We assume they are linearly dependent giving Hooke's law
$
  sigma_(alpha beta)^((d)) = K_(alpha beta gamma delta) cal(E)_(gamma delta)
$
with $K_(alpha beta gamma delta)$ being a four-index elasticity tensor. Assuming the material is homogeneous and isotropic then we can write
$
  K_(alpha beta gamma delta) = K delta_(alpha beta) delta_(gamma delta) + G (delta_(alpha beta) delta_(beta gamma) + delta_(alpha delta) delta_(beta gamma) - 2/3 delta_(alpha beta) delta_(gamma delta))
$
with $K$ being the bulk modulus and $G$ being the shear modulus. These two moduli fully describe how a given material will deform when we apply stress. $K$ represents compressibility and is typically very large. Usually we take $K -> oo$. Any deformation is then characterized by $G$.

For a simple shear
$
  sigma_(x y)^((d)) = sigma_(x y) = G gamma
$
For a simple stretch
$
  sigma_(x x) = - P + F_x/A_x => sigma_N equiv sigma_(x x) underbracket(- 1/2 (sigma_(y y)+sigma_(z z)), P) &= G(lambda^2 - lambda^(-1)) \
  &tilde.eq^(lambda = 1 + epsilon) 3 G epsilon equiv Y epsilon
$
where $Y$ is the Young modulus.

For an isotropic compression
$
  sigma_(x x) = sigma_(y y) = sigma_(z z) = -P +F/A = -P + Delta P
$
we define
$
  Delta P = - K dd(V, d: Delta)/V = - 3 K underbracket(epsilon_V, "volumetric strain")
$
so $sigma_(alpha alpha) = - P - 3 K epsilon_V$.

=== The free energy
Consider the free energy $f(E_(alpha beta))$. This guy should be invariant under coordinate transformations. Meaning all dependence on $E_(alpha beta)$ must lie in invariant quantities.

We define the left Cauchy-Green tensor of finger tensor $B_(alpha beta) equiv (E E^T)_(alpha beta)$. Consider the determinant
$
  det (B_(alpha beta) - lambda delta_(alpha beta)) &= - lambda^3 + I_1 lambda^2 - I_2 lambda + I_3
$
with $I_i$ being the invariants
$
  I_1 &= tr B = sum_i lambda_i^2 \
  I_2 &= 1/2 (tr (B)^2 - tr(B^2)) = lambda_1^2 lambda_2^2 + lambda_2^2 lambda_3^2 + lambda_1^2 lambda_3^2 \
  I_3 &= det B = product_i lambda_i^2
$
where $lambda_i$ are the principal stretches and $lambda_i^2$ are the eigenvalues of $B_(alpha beta)$.

#proof[

  By definition $B_(alpha beta)$ is symmetric so we can diagonalize it
  $
    B = Q diag(lambda_1^2, lambda_2^2, lambda_3^2) Q^TT
  $
  with $Q Q^TT = 1$.

  Then we can find
  $
    det(B - lambda bb(1)) &= det(diag(lambda_1^2, lambda_2^2, lambda_3^2) - lambda bb(1)) \
    &= (lambda_1^2-lambda)(lambda_2^2-lambda)(lambda_3^2-lambda) \
    &= - lambda^3 + underbracket((lambda_1^2 +lambda_2^2 + lambda_3^2), I_1) lambda^2 -underbracket((lambda_1^2 lambda_2^2 + lambda_1^2 lambda_3^2+lambda_2^2 lambda_3^2), I_2) lambda+ underbracket(lambda_1^2 lambda_2^2 lambda_3^2, I_3)
  $

  Lastly we show
  $
    I_2 &= 1/2 (tr (B)^2 - tr (B^2)) \
    &= 1/2 (lambda_1^2 lambda_2^2 + lambda_1^2 lambda_3^2 + lambda_2^2 lambda_1^2 + lambda_2^2 lambda_3^2 + lambda_3^2 lambda_1^2 + lambda_3^2 lambda_2^2) \
    &= lambda_1^2 lambda_2^2 + lambda_1^2 lambda_3^2 + lambda_2^2 lambda_3^2
  $

]

So we can write
$
  f(E_(alpha beta)) equiv f(I_1,I_2,I_3)
$
For an incompressible body $I_3 = 1$ and drops out and we can Taylor expand around the undeformed state (with $I_1 = I_2 = 3$) to find
$
  f(I_1, I_2) = C_1 (I_1 - 3) + C_2 (I_2 - 3)
$
this is called a Mooney-Rivlin solid.

== Kuhn theory
We consider a solution of polymers. These will crosslink and create a network turning the liquid solution into an elastic solid. We consider deforming this network. After a deforming free strands have $R'_alpha = E_(alpha beta) R_beta$. We find the change in free energy due to this deformation by
$
  f(E_(alpha beta)) = tilde(f) (E_(alpha beta)) - tilde(f)_0
$
We denote the probability of any strand having $bold(R)$ with $N$ by $psi(bold(R), N)$. Then
$
  f(E_(alpha beta)) &= rho_s integral dd(N) integral dd(bold(R)) psi(bold(R), N) {F'_"strand"-F_"strand"}
$
with $rho_s$ being the density of strands. We assume the probability of any strand having $bold(R)$ with $N$ is the probability of having a polymer with $bold(R)$ and $N$ multiplied by the probability of then making a strand with $N$ after cross-linking $psi(bold(R), N) = P(bold(R),N) Phi_0 (N)$. We obtain
$
  f(E_(alpha beta)) &= (3 rho_s k_B T)/(2N b^2) integral dd(N) Phi_0(N) integral dd(bold(R)) P(bold(R),N) (bold(R)' dot bold(R)'- bold(R) dot bold(R)) \
  &= (3 rho_s k_B T)/(2 N b^2) integral dd(N) Phi_0 (N) integral dd(bold(R)) P(bold(R),N) (E_(alpha beta) E_(alpha gamma) - delta_(beta gamma)) R_(beta) R_gamma
$
The $bold(R)$ integral can be computed. For $beta eq.not gamma$ is vanishes since the integral becomes $expval(R)$ which vanishes. So we only get something when $beta = gamma$ where it equals $expval(R^2) = N b^2\/3$. Therefore we obtain
$
  f(E_(alpha beta)) &= (3 rho_s k_B T)/(2N b^2) integral dd(N) Phi_0 (N) [(E_(alpha beta) E_(alpha gamma) - delta_(beta gamma)) delta_(beta gamma) (b^2 N)/3] \
  &= (rho_s k_B T)/2 (E_(alpha beta) E_(alpha beta) - delta_(beta beta)) \
  &= (rho_s k_B T)/2 (E_(alpha beta) E_(alpha beta) - 3)
$
or in terms of $lambda_i$
$
  f(lambda_i) & = (rho_s k_B T)/2 (sum_i lambda_i^2 - 3)
$
For a shear we find
$
  f(E_(alpha beta)) = (k_B T rho_s gamma^2)/2
$
and when we increase the shear strain from $gamma -> gamma + dd(gamma)$ we do work $sigma_(x y) dd(gamma)$ leading to a change $dd(f)$ meaning
$
  sigma_(x y) = pdv(f, gamma) = underbracket(rho_s k_B T, G) gamma
$
this is the obvious guess one might make.

For a stretch we find
$
  f(lambda_i) = (rho_s k_B T)/2 (lambda^2 + 2/lambda - 3) = 1/2 G (lambda^2 + 2/lambda - 3)
$
and when we stretch a unit volume by $lambda$ we do work $sigma_(x x) lambda^(-1) dd(lambda)$ due to the surface area changing by $lambda^(-1)$ meaning
$
  sigma_(x x) = lambda pdv(f, lambda) = G (lambda^2 - lambda^(-1))
$
both of these are consistent with the previous expressions we found.

