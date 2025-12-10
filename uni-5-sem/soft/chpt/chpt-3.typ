//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

#pagebreak()
= Surfaces and surfactants
== Surface tension
Consider a water-air interface, or any liquid-air interface: water molecules in the bulk will be more content since it is energetically favorable---due to molecular interactions. Water molecules at the interface are less content since this is less energetically favorable, this leads to a surface tension $gamma = dd(E)\/dd(A)$. We can also consider a thin film with surface tension $gamma$, now we try to pull on this film increasing its length by $dd(x)$ using some force $F$. We obtain
$
  dd(W) = F dd(x) = 2 gamma a dd(x)
$
with $a$ being the width of our film and the factor two appearing since our film has two sides---so $gamma = F\/2 a$, and the surface tension counteracts our applied force. Similarly pressing liquid out of a syringe into air using some $P_0 + Delta P$, with $P_0$ being ambient pressure, will lead to an increased surface area, so $dd(W) = Delta P dd(V) = gamma dd(A)$. If the formed droplet is assumed spherical then
$
  dd(V) = 4 pi r^2 dd(r)",  " dd(A) = 8 pi r dd(r)
$
giving the Laplace pressure
$
  Delta P = (2 gamma)/r
$
for real droplets this is given by the more general Young-Laplace equation---but we won't cover this.

== Wetting
Consider some surface and a droplet on the surface. In this case we have three different interfaces and three corresponding surface tensions: $gamma$ (liquid-vapor), $gamma_"SV"$ (solid-vapor) and $gamma_"SL"$ (solid-liquid). We want to know whether or not this droplet spreads. If the droplet has initial area $A$ then the energy before placing the droplet is $E_"before" = gamma_"SV" A$. Assuming the curvature ($A$ with $gamma$ same as $A$ with $gamma_"SL"$) of the drop is negligible we likewise have $E_"after" = (gamma + gamma_"SL") A$. Then we define the spreading coefficient
$
  gamma_"S" = (E_"before"-E_"after")/A = gamma_"SV" - gamma - gamma_"SL"
$
which is the energy difference per unit area---if $gamma_"S">0$ then the energy after is smaller and the droplet will spread as much as possible, if $gamma_"S" < 0$ the energy increases and the droplet stays.

Another measure that describes what happens is the contact angle $theta$. We consider what happens when our drop gets extended by $dd(x)$, this lengthens the liquid-vapor interface by $cos theta dd(x)$ and changes some solid-vapor to solid-liquid, and we obtain
$
  dd(E) = cos theta dd(x) gamma + (gamma_"SL"-gamma_"SV") dd(x)
$
in equilibrium
$
  dv(E, x) = 0 => gamma cos theta + gamma_"SL" - gamma_"SV" = 0
$
giving Young's equation
$
  cos theta = (gamma_"SV"-gamma_"SL")/gamma = (gamma_"S"+gamma)/gamma
$
what happens when we introduce gravity? The energy of the surface area is given by
$
  G_A = 4 pi gamma r^2
$
and the gravitational potential energy is
$
  E_"pot" = m g h = (4 pi rho g)/3 r^4
$
obviously $r^4 > r^2$ so for larger drops gravity wins and it will spread, but for small drops surface tension wins and it will keep its shape. They are equal at the scale
$
  G_A = E_"pot" => r^* tilde sqrt((gamma)/(rho g))
$
giving something like $r^*_"water" tilde 1.4 "mm"$.

Now consider placing a drop on a surface and letting it spread until it reaches equilibrium with height $h$. We assume that the volume $V$ of the drop is known, as well as all the surface tensions. Then ignoring curvature and using $h\/2$ as the height of the center of mass we obtain
$
  E_"pot" = V rho g h/2
$
and using $A = V\/h$
$
  G_A = (gamma + gamma_"SL") V/h - gamma_"SV" V/h = - gamma_"S" V/h
$
at equilibrium
$
  dv(E_"tot", h) = 0 => (V rho g)/2 + (gamma_S V)/(h^*)^2 = 0
$
giving
$
  h^* = sqrt((-2gamma_S)/(rho g))
$
relating this to $r^*$ can be done, and it simply gives $h^* tilde.eq r^* theta$.

=== Capillary effects
What happens when inserting a tube into some liquid? We denote the diameter of our tube by $2 a$, and the liquid rises within the tube to a height $h$ above the liquid outside. The volume within the tube is then $V_"tube"=pi a^2 h$, with center of mass $h\/2$, so
$
  E_"pot" = (pi a^2 rho g h^2)/2
$
and
$
  G_A = 2 pi a h (gamma_"SL"-gamma_"SV") = - 2 pi a h gamma cos theta
$
since it happens outside and inside, this then gives
$
  h^* = (2 gamma cos theta)/(a rho g) = (2 (r^*)^2 cos theta)/a
$

== Thermodynamics of an interface
We consider a system with two phases I and II that coexist and are in equilibrium. They have the same temperature $T$ and chemical potential $mu_i$. The total grand canonical free energy of the system can then be written as
$
  G = G_I + G_(I I) + G_A
$
with $G_A$ representing the grand canonical free energy of the interface. We can use this to define the interfacial free energy.

Since volume is an extensive quantity we can write
$
  G(V,T,mu) = V g(T,mu)
$
and similarly
$
  G_A (A,T,mu) = A g_A (T,mu)
$
We can show
$
  g = G/V = - P (T,mu)
$
#proof[

  By a Legendre transform
  $
    G(V,T,mu_i) = F(V,T,N_i) - sum_i N_i mu_i
  $
  note
  $
    dd(F) = - P dd(V) - S dd(T) + sum mu_i dd(N_i)
  $
  and by the above
  $
    dd(G) & = dd(F) - sum dd(N_i) mu_i - sum_i N_i dd(mu_i) \
          & = - P dd(V) - S dd(T) - sum_i N_i dd(mu_i)
  $
  By definition we also have
  $
    dd(G) = g dd(V) + V pdv(g, T) dd(T) + V sum_i pdv(g, mu_i) dd(mu_i)
  $
  comparing these give
  $
    g & = - P";  " S & = V pdv(P, T)";  " N_i/V & = pdv(P, mu_i)
  $
]

Similarly we can show
$
  g_A = gamma(T, mu)
$
with a positive sign since surface tension is contractile.
#proof[

  We have
  $
    G_A (A, T, mu_i) = A g_A (T, mu_i)
  $
  with
  $
    dd(G_A) & = gamma dd(A) -S_A dd(T) - sum_i^A N_(A i) dd(mu_i) \
    dd(G_A) & = g_A dd(A) + A pdv(g_A, T) dd(T) + A sum_i pdv(g_A, mu_i) dd(mu_i)
  $
  everything follows as before giving
  $
    g_A & = gamma";  " S_A & = - A pdv(gamma, T)";  " N_(A i)/A & = - pdv(gamma, mu_i)
  $
]

We quickly find the Gibbs-Duhem equation. We can write
$
  dd(G) = - V dd(P) - P dd(V)
$
and
$
  dd(G) = - S dd(T) - P dd(V) - sum_i N_i dd(mu_i)
$
combing these we obtain
$
  V dd(P) = S dd(T) +sum_i N_i dd(mu_i)
$
which is the Gibbs-Duhem equation. Similarly for $G_A = A gamma(T, mu_i)$ we obtain
$
  A dd(gamma) = -S_A dd(T)-sum_i N_(A i) dd(mu_i)
$

== Surfactants
We define the surface excess
$
  Gamma_i = N_(A i)/A
$
this is simply the amount of molecules per unit area on the surface. We denote the number density of a species $i$ by $n_i$. On either side of the interface $n_i$ will approach $n_(I i)$ or $n_(I I i)$ with these being the respective bulk densitites. Around the surface $n_i$ will differ from these. As an example if the species $i$ were attracted to the surface then $n_i$ would obviously be larger at the surface. We then have an equivalent definition of the surface excess
$
  Gamma_i = integral_(-oo)^0 n_i (z) - n_(I i) dd(z) + integral_0^oo n_i (z) - n_(I I i)
$
surfactants will usually have $Gamma_i > 0$ meaning they like being at the surface. The origin is defined by $Gamma_i = 0$.



Because of entropy not all surfactants go to the surface since this not energetically favorable if the surface is _crowded_. At sufficient concentration surfactants will start forming micelles through self-assembly even though this lowers their entropy significantly because it is still energetically favorable. We will see this later.

Marangoni effect\*

=== The Langmuir isotherm
We consider adding surfactants to a solution. Assume we have a surface with surfactant density $n$. We treat the surface as a lattice where surfactants can choose to occupy states. We denote the fraction of occupied states by $theta$. Now imagine a surfactant from the bulk hitting the surface. Given the state is empty it will begin occupying the state, otherwise if the state is already occupied then it will bounce back into the bulk. We denote the rate of hitting surface by $v_a$. We write this as
$
  v_a = K_a (1- theta)
$
with $K_a$ being the number of attempts per unit time. Similarly we denote the rate of leaving the surface by $v_d$ and write this as
$
  v_d = K_d theta
$
with $K_d$ being the number of attempts per unit time. When in equilibrium these are equal meaning
$
  K_a (1-theta_"eq") = K_d theta_"eq"
$
we solve for $theta_"eq"$
$
  theta_"eq" = K/(1 + K) " with " K equiv K_a/K_d
$
We can also write
$
  theta_"eq" = Gamma/Gamma_s
$
with $Gamma_s$ being the saturation surface excess. Similarly we can write
$
  K = n/n_s
$
with $n_s$ being the bulk density at saturation.

#proof[

  We write $K_a = alpha n$ and define saturation to be when $theta_"eq" (n_s) = 1\/2$. Then
  $
    K_d = alpha n_s
  $
  so it is required for consistency.
]

Then solving for $Gamma$ we obtain
$
  Gamma = n/(n+n_s) Gamma_s
$
which is the Langmuir isotherm.

=== Surface tension
We can now derive how the surface tension changes. We want $gamma(n)$ at constant $T$. From the Gibbs-Duhem equation we have
$
  dv(gamma, mu) & = - Gamma
$
we write this as
$
  dv(gamma, n) dv(n, mu) & = - (n Gamma_s)/(n+n_s)
$
so we need some $mu(n)$. We use
$
  mu(n) = mu_0 + k_B T ln n
$
the proof is given below.

#proof[

  Consider a gas with $n = N V^(-1)$. We assume each molecule has some zero-point energy $mu_0$. The partition function is
  $
    Z = (Z_"trans" Z_"part")^N/N!
  $
  each molecule has
  $
    E = 1/2 m v^2 + mu_0
  $
  the first term enters in $Z_"trans"$ and the second enters in $Z_"part"$.

  We know
  $
    Z_"trans" = V/lambda_"th"^3
  $
  and by definition
  $
    Z_"part" = exp(- mu_0/(k_B T))
  $
  Then the free energy is
  $
    F & = - k_B T ln Z \
      & = - k_B T [- N ln N + N - (mu_0 N)/(k_B T) + N ln (V)/lambda_"th"^3]
  $
  by definition
  $
    mu = pdv(F, N) & = dots = mu_0 + k_B T ln (N lambda_"th"^3)/V \
                   & = mu_0 + k_B T ln(n lambda_"th"^3) \
                   & =^(lambda_"th" = "const") mu_0 + k_B T ln n
  $
  with the second term just being translational entropy.

]

Then
$
  dv(mu, n) = (k_B T)/n => dv(n, mu) = n/(k_B T)
$
so we obtain
$
  dv(gamma, n) = - (k_B T Gamma_s)/(n+n_s)
$
This is easily integrated giving
$
  gamma(n) = gamma_0 - Gamma_s k_B T ln(1+n/n_s)
$
at low $n$ this is heavily dominated by $gamma_0$ as we would expect and as $n$ increases $gamma$ decreases.

== Micelles
We now consider what happens if saturate the surface.

Consider $m$ surfactant molecules in a box, these will either stay apart or self-assemble to form a micelle. A single surfactant molecule has
$
  mu_1 = mu_1^0 + k_B T ln n_1
$
and a micelle has
$
  mu_m = mu_m^0 + k_B T ln n_m
$
These are equal in equilibrium $m mu_1 = mu_m$ giving
$
  k_B T ln n_m/n_1^m & = m mu_1^0-mu_m^0 \
           n_m/n_1^m & = exp((m mu_1^0-mu_m^0)/(k_B T))
$
this is just the Boltzmann distribution! We define
$
  n_c^(1-m) equiv exp((m mu_1^0-mu_m^0)/(k_B T))
$
then we can write
$
  n_m & = n_1^(m) n_c^(1-m) = (n_1/n_c)^m n_c
$
Using $n = n_1 + m n_m$ we obtain
$
  n/n_c & = n_1/n_c + m (n_1/n_c)^m
$
Consider the case $n_1 < n_c$. Since $m$ is large we get
$
  n_1 tilde.eq n
$
so below $n_c$ we have no micelles.

Consider the case $n_1 > n_c$. Here the second term dominates
$
  n/n_c & tilde.eq m(n_1/n_c)^m \
        & => n_1 tilde.eq n_c (n/(m n_c))^(1\/m) tilde.eq n_c
$
this happens since all excess surfactant molecules get used to create micelles. We can also find the micelle concentration
$
  n_m & = n_c (n_1/n_c)^m \
      & tilde.eq n/m
$
The above shows that $n_c equiv n_"cmc"$. We also see that the surface tension stops decreasing since the unimer concentration is limited by $n_1 tilde.eq n_c$.

