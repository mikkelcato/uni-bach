#import "chpt-temp.typ": *
#show: chpt-note.with()

#let d1 = diagram(
  spacing: 8em,
  node((0, 0), name: <A>),
  edge(<A>, <B>, "<-"),
  node((0, 2), name: <B>),
  edge(<B>, <C>, "->", label: $ln a$, label-pos: 95%),
  node((3.5, 2), name: <C>),
  node((0, 1.5), name: <A1>),
  node((1.5, 1.5), name: <A2>),
  edge(
    <A1>,
    <A2>,
    "dotted",
    stroke: 1.3pt,
    label: [_inflation_],
    label-sep: -5.5em,
    label-pos: 50%,
  ),
  node((2.5, 0.25), name: <B1>),
  edge(<A2>, <B1>, label: $ln H^(-1)$, label-pos: 100%, label-sep: -0.25em),
  edge(
    <A2>,
    <B1>,
    label: [$lambda_"physical"$ outside horizon],
    label-angle: 28deg,
    label-sep: 0.5em,
    label-pos: 50%,
  ),
  node((3, 0.25), name: <B2>),
  edge(
    <B>,
    <B2>,
    label: $ln lambda_"physical"$,
    label-pos: 105%,
    label-sep: -3.75em,
  ),
  node((1.5, 2), name: <C1>),
  edge(<C1>, <A2>, "--"),
)



= The Friedmann-Robertson-Walker Spacetime
By _the cosmological principle_ the Universe appears homogeneous and isotropic on large scales. We also measure the Universe as being _flat_.#footnote[We measure $k = 0$.] The solution to the Einstein field equations with these assumptions is the _Friedmann-Robertson-Walker metric_
$
  dd(s^2) = - dd(t^2) + a^2 (t) dd(bold(x)^2)
$
where $a(t)$ is the _scale factor_.

== The scalar field
We assume the scalar field $phi.alt$ contributes a negligible energy density meaning. Then we can ignore its influence and treat $g_(mu nu)$ as fixed. We consider the free scalar field $phi.alt$ as before
$
  S = - 1/2 integral dd(x, 4) sqrt(-g) [g^(mu nu) partial_mu phi.alt partial_nu phi.alt + m^2 phi.alt^2]
$
with the equation of motion
$
  square phi.alt + m^2 phi.alt = 0
$
We can simplify things by defining the _conformal time_
$
  eta (t) equiv integral_(t_0)^t dd(t)/a(t)",  " a dd(eta) = dd(t)
$
With $eta$ the metric becomes
$
  dd(s^2) = a^2 (eta) [- dd(eta^2) + dd(bold(x)^2)]
$
implying
$
  g_(mu nu) = underbracket(a^2 (eta), Omega^2) eta_(mu nu) tilde "conformally flat"
$
We also make a _field redefinition_
$
  chi equiv a(eta) phi.alt
$
and trivially $sqrt(-g) = a^4$. Then we find
$
  S = -1/2 integral dd(x, 3) dd(eta) [-(partial_eta chi)^2 + (grad chi)^2 + m_"eff"^2 (eta) chi^2]
$
with
$
  m_"eff"^2 (eta) equiv m^2 a^2 - 1/a partial_eta^2 a
$
We see that _all_ effects of gravity are within $m_"eff"$ and $S$ is time-dependent meaning the Hamiltonian will be time-dependent. This implies the energy of $chi$ is not conserved.#footnote[This leads to particle creation.]

The equation of motion for $chi$ becomes#footnote[We use $a^2 square phi.alt = - partial_eta^2 phi.alt - 2 1/a partial_eta a partial_eta phi.alt + laplacian phi.alt$.]
$
  partial_eta^2 chi - laplacian chi + (m^2 a^2 - 1/a partial_eta^2 a) chi = 0
$
The $laplacian$ is as in Minkowski space meaning we can expand in Fourier modes as
$
  chi(eta, bold(x)) = integral dd(bold(k), 3)/(2 pi)^(3\/2) chi_bold(k) (eta) e^(i bold(k) dot bold(x))
$
We then obtain
$
  partial_eta^2 chi_bold(k) + (k^2 + m^2 a^2 (eta) - 1/a partial_eta^2 a ) chi_bold(k) = 0
$
Which can be written as
$
  partial_eta^2 chi_bold(k) + omega_bold(k)^2 (eta) chi_bold(k) = 0
$
with $omega_bold(k)^2 equiv k^2 + m_"eff"^2 (eta)$. We assume we know a function $V_bold(k) (eta)$ satisfying#footnote[$omega_bold(k)^2 = omega_(-bold(k))^2$ implies $V_bold(k) = V_(-bold(k))$.]
$
  partial_eta^2 V_bold(k) + omega_bold(k)^2 V_bold(k) = 0
$
Then the general solution $chi_bold(k) (eta)$ can be written as
$
  chi_bold(k) (eta) = 1/sqrt(2) (a_bold(k)^- V_bold(k)^* + b_bold(k) V_bold(k))
$
We assume $chi$ is real implying $chi_bold(k)^* = chi_(-bold(k))$. Then
$
  (a_bold(k)^-)^* =^! b_(-bold(k)) equiv a_bold(k)^+
$
and we obtain
$
  chi_bold(k) (eta) = 1/sqrt(2) ( a_bold(k)^- V_bold(k)^* + a_(-bold(k))^+ V_bold(k))
$
Then
$
  chi(eta, bold(x)) &= integral dd(bold(k), 3)/(2 pi)^(3\/2) 1/sqrt(2) (a_bold(k)^- V_bold(k)^* + a_(-bold(k))^+ V_bold(-k)) e^(i bold(k) dot bold(x)) \
  &= integral dd(bold(k), 3)/(2 pi)^(3\/2) 1/sqrt(2) (a_bold(k)^- V_bold(k)^* e^(i bold(k) dot bold(x)) + a_bold(k)^+ V_bold(k) e^(-i bold(k) dot bold(x)))
$
We can show
$
  a_bold(k)^-/sqrt(2) = (partial_eta V_bold(k) chi_bold(k) - V_bold(k) partial_eta chi_bold(k) )/(partial_eta V_bold(k) V_bold(k)^* - V_bold(k) partial_eta V_bold(k)^*) = W[V_bold(k), chi_bold(k)]/W[V_bold(k),V_bold(k)^*]
$

== Quantization
We quantize $chi$ by promoting $a_bold(k)^plus.minus$ to operators
$
  a_bold(k)^- -> a_bold(k)";  " a_bold(k)^+ -> a_bold(k)^dagger
$
and obtain
$
  chi(eta, bold(x)) = integral dd(bold(k), 3)/(2 pi)^(3\/2) 1/sqrt(2) (a_bold(k) V_bold(k)^* e^(i bold(k) dot bold(x)) + a_bold(k)^dagger V_bold(k) e^(-i bold(k) dot bold(x)))
$
We impose the usual commutators
$
  [a_bold(k),a_bold(k')^dagger] = delta(bold(k)-bold(k)')",   others" = 0
$
implying#footnote[Where $pi(eta, bold(x)) = pdv(cal(L), (partial_eta chi)) = partial_eta chi(eta, bold(x))$.]
$
  [chi(eta, bold(x)_1), pi(eta, bold(x)_2)] = i delta(bold(x)_1-bold(x)_2)
$
when
$
  Im V_bold(k)^* partial_eta V_bold(k) = W[V_bold(k),V_bold(k)^*] = 2 i
$
This ensures the Wronskian conditions are satisfied by the modes $u_bold(k)$#footnote[This normalizes our $V_bold(k)$.]
$
  u_bold(k) (eta,bold(x)) = 1/(2 pi)^(3\/2) 1/sqrt(2) V_bold(k) (eta) e^(i bold(k) dot bold(x))
$

We now need to _choose_ a vacuum. This choice is not unique since above $a_bold(k) tilde V_bold(k)$. The vacuum is typically chosen by imposing some physical constraint on the system.#footnote[This usually depends on which observer we are trying to describe.]

== The instantaneous vacuum
We attempt to define vacuum as the smallest energy eigenstate at some time $eta = eta_0$ of the _instantaneous_ Hamiltonian. The Hamiltonian of $chi$ is
$
  H = 1/2 pi^2 + 1/2 omega^2 (eta) chi^2
$
We can obtain
$
  E(eta) &equiv braket(0, H(eta), 0) \
  &= integral dd(bold(k), 3)/(2pi)^3 (abs(partial_eta V_bold(k))^2 + omega_bold(k)^2 abs(V_bold(k))^2)/4
$
This is minimal at $eta = eta_0$ when
$
  abs(partial_eta V_bold(k))^2 + omega_bold(k)^2 abs(V_bold(k))^2
$
is minimal for all $k$-modes at $eta = eta_0$. We define
$
  s_bold(k) = V_bold(k) (eta_0)";  " t_bold(k) = partial_eta V_bold(k) (eta_0)
$
with the Wronskian condition
$
  s_bold(k)^* t_bold(k) - s_bold(k) t_bold(k)^* = 2 i
$
We now want to minimize
$
  abs(t_bold(k))^2 + omega^2 (eta_0) abs(s_bold(k))^2
$
Taking $s_bold(k)$ to be real and writing $t_bold(k) = t_(1 bold(k)) + i t_(2 bold(k))$ we obtain
$
  s_bold(k) = 1/t_(2 bold(k))
$
Then
$
  4 E (eta_0) = integral dd(bold(k), 3)/(2 pi)^3 (t_(1 bold(k))^2 + t_(2 bold(k))^2 + (omega_bold(k)^2 (eta_0))/(t_(2 bold(k))^2))
$
When $omega_bold(k)^2 (eta_0) > 0$ we have a minima for
$
  t_(1 bold(k)) = 0";  " t_(2 bold(k)) = sqrt(omega_bold(k) (eta_0))
$
When $omega_bold(k)^2 (eta_0) < 0$ no minima exists. Then for $omega_bold(k)^2 (eta_0) > 0$ the _instantaneous vacuum_ is given by
$
  V_bold(k) (eta_0) = 1/sqrt(omega_bold(k) (eta_0))";  " partial_eta V_bold(k) (eta_0) = i sqrt(omega_bold(k) (eta_0)) = i omega_bold(k) (eta_0) V_bold(k) (eta_0)
$
This agrees#footnote[Up to a phase.] with the Minkowski spacetime result where $omega_bold(k)$ is time-independent yielding
$
  V_bold(k) (eta) = 1/sqrt(omega_bold(k)) e^(-i omega_bold(k) eta)
$
When we introduce a time-dependent background this vacuum changes over time#footnote[With each instantaneous vacuum mimicking the Minkowski vacuum.]
$
  ket(0)_eta_0 eq.not ket(0)_eta_1
$
Consider two sets of modes $V_bold(k) (eta)$ and $U_bold(k) (eta)$ describing the vacuum at $eta = eta_0$ and $eta = eta_1$ respectively and assume they are both known at $eta = eta_0$. We can relate them using the Bogolubov coefficients $alpha_bold(k)$ and $beta_bold(k)$
$
  V_bold(k)^* = alpha_bold(k) U_bold(k)^* + beta_bold(k) U_bold(k) ";  " partial_eta V_bold(k)^* = alpha_bold(k) partial_eta U_bold(k)^* + beta_bold(k) partial_eta U_bold(k)
$
where
$
  alpha_bold(k) = evaluated((partial_eta U_bold(k) V_bold(k)^* - U_bold(k) partial_eta V_bold(k)^*)/(2 i))_eta_0"  " beta_bold(k)^* = evaluated((partial_eta U_bold(k) V_bold(k) - U_bold(k) partial_eta V_bold(k))/(2 i))_eta_0
$
We could have used $eta -> -oo$ instead of $eta = eta_0$. All states with finite $eta_0$ would then appear excited. We would have particle creation as the Universe evolved ex nihilo!

== The adiabatic vacuum
When $omega_bold(k)$ varies slowly with $eta$ the _adiabatic vacuum_ is useful. We still have
$
  partial_eta^2 V_bold(k) + omega_bold(k)^2 V_bold(k) = 0",  " omega_bold(k)^2 = k^2 + m_"eff"^2
$
This has a formal _WKB-type_ solution
$
  V_bold(k) = 1/sqrt(W_bold(k)) exp[-i integral^eta W_bold(k) (eta') dd(eta')]
$
with#footnote[This is our _dispersion relation_.]
$
  W_bold(k)^2 (eta) = omega_bold(k)^2 (eta) - 1/2 [1/W_bold(k) partial_eta^2 W_bold(k) - 3/2 1/W_bold(k)^2 (partial_eta W_bold(k))^2]
$
The _adiabatic limit_ has
$
  W_bold(k) tilde.eq W_bold(k)^((0)) (eta) = omega_bold(k) (eta)
$
and we obtain
$
  V_bold(k)^((0)) = 1/sqrt(omega_bold(k)) exp[-i integral_(eta_0)^eta omega_bold(k) (eta') dd(eta')]
$
This solution approximately satisfies the instantaneous minimum energy condition. We have
$
  V_bold(k) (eta_0) = V_bold(k)^((0)) (eta_0)";  " evaluated(partial_eta V_bold(k) (eta))_eta_0 = evaluated(partial_eta V_bold(k)^((0)) (eta))_eta_0
$
implying
$
  V_bold(k) (eta_0) = 1/sqrt(omega_bold(k) (eta_0))";  " evaluated(partial_eta V_bold(k) (eta))_eta_0 = [- i omega_bold(k) - 1/(2 omega_bold(k)) partial_eta omega_bold(k)] evaluated(1/sqrt(omega_bold(k)))_eta_0
$
Then the energy becomes
$
  E_bold(k) &= 1/4 (abs(partial_eta V_bold(k))^2 + omega_bold(k)^2 abs(V_bold(k))^2) \
  &= 1/2 omega_bold(k) + 1/16 1/omega_bold(k)^3 (partial_eta omega_bold(k))^2 \ &tilde.eq 1/2 omega_bold(k)
$
which is the same as the instantaneous vacuum for $eta = eta_0$.

We similarly define the second order adiabatic vacuum using
$
  (W_bold(k)^((2)))^2 = omega_bold(k)^2 - 1/2 [1/omega_bold(k) partial_eta^2 omega_bold(k) - 3/2 1/omega_bold(k)^2 (partial_eta omega_bold(k))^2]
$
This will differ from the instantaneous vacuum.

= The FRW cosmology
The FRW metric can be written as
$
  dd(s^2) = - dd(t^2) + a^2 (t) [dd(r^2)/(1-k r^2) + r^2 dd(Omega^2)]
$
and we can treat the Universe as a co-moving ideal fluid with
$
  T_(mu nu) = diag(rho, p, p, p)
$
Then from the EFE we find the _Friedmann equations_
$
              H^2 & = (8 pi G_N)/3 rho - k/a^2 \
  dot.double(a)/a & = - (4 pi G_N)/3 (rho + 3 p)
$
with $H$ being the _Hubble parameter_
$
  H equiv dot(a)/a
$
We also have
$
  D^nu T_(mu nu) = 0
$
implying
$
  dot(rho) + 3 H ( rho + p) = 0
$
which is simply the _continuity equation_.#footnote[The continuity equation is not independent and follows directly from the Friedmann equations.] We usually assume an equation of state of the form $p = w rho$ where $w$ is constant. Then
$
  rho tilde a^(-3(w + 1))
$
which in a spatially flat Universe implies
$
  a tilde t^(2\/3(1+w))
$
We care about specific values of $w$,
$
  w_Lambda = -1 & => rho_Lambda tilde "constant" \
      w_r = 1/3 & => rho_r tilde a^(-4) \
        w_m = 0 & => rho_m tilde a^(-3)
$
By comparing the scaling of the components we see the Universe will experience _eras_ of domination with certain components dominating the rest.

We can quantify the amount of curvature by defining the _critical_ energy density
$
  rho_c = (3 H^2)/(8 pi G_N)
$
This is the energy density we would measure if the Universe were spatially flat. Then we can define _density parameters_
$
  Omega = rho/rho_c
$
implying
$
  Omega - 1 = k/(a^2 H^2)
$
Then $Omega$ fully determines the spatial curvature $k$.

== The Big Bang
The Friedmann equations and measurements imply
$
  a(t=t_"BB") -> 0
$
with $t_"BB" tilde 0$ being some finite time in the past. We also have
$
  rho_r -> oo
$
as $a -> 0$, implying the Universe expanded from a hot and dense radiation-dominated state. This is called the _Big Bang model_. By assumption photons originating from the Big Bang can only have travelled a finite distance. This implies the _observable_ Universe is finite, however, the _true_ Universe may be infinite.

We can make the notion of _size_ precise by defining the _particle horizon_. We consider photons moving radially and satisfying#footnote[We have $d_H tilde H^(-1)$ when $a tilde t^beta$ for $beta < 1$.]
$
  dd(r) = dd(t)/a(t)
$
The particle horizon is then the _physical_ distance photons have travelled since $t_"BB"$ at some time $t$
$
  d_H (t) = a(t) integral_0^t (dd(t'))/a(t')
$
When the Universe was radiation-dominated we had $a tilde t^(-1\/2)$ implying
$
  d_H (t) = 1/H
$
Then any physical scale $lambda$ outside the horizon
$
  lambda >> 1/H
$
is not in _causal contact_. We call $H^(-1)$ the _Hubble radius_.

When
$
  t_"rec" tilde 380000 "years"
$
the Universe had cooled sufficiently due to the expansion that electrons and protons could form neutral hydrogen
$
  e^- + p -> H
$
Then after $t_"rec"$ the Universe becomes transparent to photons leaving us with the _cosmic microwave background radiation_. The CMB is a blackbody at the temperature
$
  T_"CMB" tilde 2.735 "K"
$
with inhomogeneities on the order
$
  (Delta T)/T tilde 10^(-4)
$
This is a problem.

== The horizon problem
To see why consider the size of the observable Universe at $t_"LS"$#footnote[This is the largest physical scale $lambda$.]
$
  lambda_H (t_"LS") = d_H (t_0) (a_"LS"/a_0) = d_H (t_0) (T_0/T_"LS")
$
We compare this with the Hubble radius. Assuming the Universe was matter-dominated at $t_"LS"$ then
$
  H^2 tilde T^3
$
implying
$
  H_"LS"^(-1) = d_H (t_"LS") = d_H (t_0) (T_"LS"/T_0)^(-3\/2)
$
Then
$
  lambda_H^3/d_H^3 tilde 10^6
$
We find there were $tilde 10^6$ causaully disconnected _Hubble volumes_ within the volume corresponding to the observable Universe at $t_"LS"$. This is weird since we observe the CMB to be homogeneous. This is called the _horizon problem_.

The reason the above is a problem is $d_H$ being finite. This happens since $dot.double(a) < 0$ implying the lower limit of $d_H$ converges. The problem would then be solved if we had $dot.double(a) > 0$ implying the lower limit of $d_H$ diverges. Then a phase just after $t_"BB"$ and before radiation-domination with $dot.double(a) > 0$ would _push_ out $d_H$ allowing the Universe to _thermalize_ thereby solving the horizon problem as in @inflation. This period is called _cosmic inflation_.

#figure(
  scale(d1, 85%),
  caption: [How inflation solves the horizon problem.],
)<inflation>

Consider
$
  dot.double(a)/a = - (4 pi G_N)/3 (rho+3p)
$
we find $dot.double(a) > 0$ implies $ p < - rho/3 $ This would occur if $p eq - rho$ implying#footnote[This could be through a cosmological constant $Lambda$.]
$
  rho eq "constant"
$
Then
$
  H^2 = (8 pi G_N)/3 rho
$
implies
$
  H tilde.eq "constant"
$
during inflation. We find
$
  a tilde e^(H_"inf" t)
$
which is called _de Sitter spacetime_.

= Inflation
To describe inflation we seek a field theoretic description of a fluid with
$
  p < - rho/3
$
This is simplest by considering a scalar field $phi.alt$ referred to as the _inflaton_. The action describing such a field is#footnote[Assuming minimal coupling.]
$
  S &= integral dd(bold(x), 4) sqrt(-g) cal(L) \
  &= integral dd(t, bold(x), [1,3]) a^3 [1/2 dot(phi.alt)^2 - 1/2 partial_i phi.alt partial^i phi.alt - V(phi.alt)]
$
Using the Euler-Lagrange equations we find
$
  dot.double(phi.alt) + 3 H dot(phi.alt) - 1/a^2 laplacian phi.alt + V_phi.alt (phi.alt) = 0
$
with
$
  V_phi.alt = dv(V, phi.alt)
$
The energy-momentum tensor is
$
  T_(mu nu) & = - 2/sqrt(-g) dv(S, g^(mu nu), d: delta) \
            & = partial_mu phi.alt partial_nu phi.alt + g_(mu nu) cal(L)
$
We compare this with the energy-momentum tensor of an ideal fluid
$
  T_(mu nu) = diagonalmatrix(rho, a^2 p, a^2 p, a^2 p)
$
Then by inspection
$
      T_(0 0) & = rho_phi.alt \
              & = 1/2 dot(phi.alt)^2 + V(phi.alt) + 1/(2 a^2) (grad phi.alt)^2 \
  T_(i i)/a^2 & = p_phi.alt \
              & = 1/2 dot(phi.alt)^2 - V(phi.alt) - 1/(6 a^2) (grad phi.alt)^2
$
We see if the $grad phi.alt$ term dominated then
$
  p_phi.alt = - rho_phi.alt/3
$
which is insufficient to drive inflation. Alternatively, if the $dot(phi.alt)^2$ term dominated then
$
  p_phi.alt = rho_phi.alt
$
which is obviously also insufficient. However, if the $V(phi.alt)$ term dominates then
$
  p_phi.alt = - rho_phi.alt
$
as required. We then want to describe inflation by a homogeneous $V(phi.alt)$ dominated scalar field. Therefore we expand $phi.alt$ as a pertubation on some homogeneous classical background
$
  phi.alt(bold(x), t) = phi.alt_c (t) + dd(phi.alt(bold(x), t), d: delta)
$
with
$
  braket(0, phi.alt (bold(x),t), 0) = phi.alt_c (t)
$
By assumption $dd(phi.alt, d: delta)$ is small so
$
  rho_phi.alt_c & = 1/2 dot(phi.alt)_c^2 + V(phi.alt_c) \
    p_phi.alt_c & = 1/2 dot(phi.alt)_c^2 - V(phi.alt_c)
$
Then if $V(phi.alt_c) >> dot(phi.alt)_c^2$ we have $p_phi.alt_c tilde.eq - rho_phi.alt_c$.

== Slow-roll inflation
The background is described by
$
  dot.double(phi.alt)_c + 3 H dot(phi.alt)_c + V_phi.alt (phi.alt_c) = 0
$
Assuming $dot(phi.alt)^2 << V(phi.alt)$ we need to also assume $dot.double(phi.alt)$ is small, otherwise inflation would stop quickly. These two assumptions are referred to as _slow-roll inflation_. Then we find the _slow-roll equations_
$
  3 H dot(phi.alt)_c & tilde.eq - V_phi.alt (phi.alt_c) \
                 H^2 & tilde.eq (8 pi G_N)/3 V(phi.alt_c)
$
These are valid when the _slow-roll conditions_
$
       dot(phi.alt)_c^2 & << V(phi.alt_c) \
  dot.double(phi.alt_c) & << 3 H dot(phi.alt)_c
$
are satisfied. These imply
$
  V_phi.alt^2/V << H^2";  " V_(phi.alt phi.alt) << H^2
$
We define the _slow-roll parameters_
$
  epsilon.alt_phi.alt & = - dot(H)/H^2 = 1/(16 pi G_N) (V_phi.alt/V)^2 \
  eta_phi.alt & = 1/(8 pi G_N) V_(phi.alt phi.alt)/V = V_(phi.alt phi.alt)/(3H^2)
$
Then the slow-roll conditions become#footnote[We formally have inflation whenever $epsilon.alt_phi.alt < 1$. ]
$
  epsilon.alt_phi.alt << 1";  " eta_phi.alt << 1
$
We quantify the amount of inflation by the $\#"e-foldings"$#footnote[From $a_f = a_i e^(N)$ since $a tilde e^(H t)$ during inflation. ]
$
  N = ln a_f/a_i = integral_(t_i)^(t_f) H dd(t)
$
Usually we need $N tilde 60$.

There are many $V(phi.alt)$ satisfying the slow-roll conditions. We consider two broad classes of inflation: _large field models_ and _small field models_. The large field models are described by a potential of the form
$
  V(phi.alt) = Lambda^4 (phi.alt/mu)^p
$
with $p > 0$. Then
$
  eta_phi.alt & = M_p^2 V_(phi.alt phi.alt)/V = p(p-1) M_p^2/phi.alt^2
$
and $eta_phi.alt << 1$ implies $phi.alt > M_p$ for $p eq.not 1$. Small field models are described by a potential of the form#footnote[We can treat large field models as being models where there is nothing preventing $phi.alt$ from reaching the minima, while small fields models can exhibit false vacua. Then the large fields models can be understood with the notion of _chaotic inflation_ where distinct pockets of inflation could emerge in the early Universe due to inhomogeneities. However, we may expect operators with $phi.alt$ to very large orders become important, as they would be unsupressed when $phi.alt > M_p$. The small field models are motivated by particle physics, but are generally more complicated and less understood.]
$
  V(phi.alt) = lambda^4 [1 - (phi.alt/mu)^p]
$
with $p > 0$. Then
$
  eta_phi.alt = M_p^2 V_(phi.alt phi.alt)/V = -p(p-1) (phi.alt/mu)^p M_p^2/phi.alt^2
$
assuming $mu << M_p$ then $eta_phi.alt << 1$ implies
$
  phi.alt << (mu/M_p)^(2\/(p-2)) mu << mu
$
for $p > 2$.

#pagebreak()
= Field Pertubations
One reason why inflation had success was due to it making the Universe $tilde$ uniform. However, due to the rapid expansion small fluctuations in the inflation field become stretched to cosmological scales. This naturally induces ripples in spacetime due to the EFE coupling fluctuations in the inflation field with the metric. These ripples eventually seed perturbations in the CMB and matter. Thereby laying the blueprint for structure formation as these over- and underdensities becomes unstable in the late Universe. The reason we exist can then be attributed to these quantum fluctuations in the very early Universe!#footnote[Assuming inflation holds water.]

== Quantum fluctuations of a scalar field in quasi de Sitter
We consider a scalar field $phi.alt(bold(x), t)$ on a time-dependent background. We split the field as above
$
  phi.alt(bold(x), t) = phi.alt_c (t) + dd(phi.alt(bold(x), t), d: delta)
$
We would now like to quantize $dd(phi.alt(bold(x), t), d: delta)$. Note,
$
  dd(phi.alt(bold(x), t), d: delta) = integral dd(k, 3)/(2 pi)^(3\/2) e^(i bold(k) dot bold(x)) dd(phi.alt_k (t), d: delta)
$
We assume $phi.alt$ is _sub-dominant_ allowing us to ignore any backreaction on the metric. The equation of motion for the field becomes to linear order
$
  dot.double(phi.alt)_c + dd(dot.double(phi.alt), d: delta) + 3 H (dot(phi.alt)_c + dd(dot(phi.alt)_c, d: delta)) - laplacian/a^2 dd(phi.alt, d: delta) + underbracket(V_phi.alt + V_(phi.alt phi.alt) dd(phi.alt, d: delta), "Taylor expansion") = 0
$
We apply the background equation of motion
$
  dot.double(phi.alt)_c + 3 H dot(phi.alt)_c + V_phi.alt = 0
$
and insert the Fourier transform to find
$
  dd(dot.double(phi.alt)_k, d: delta) + 3 H dd(dot(phi.alt)_k, d: delta) + k^2/a^2 dd(phi.alt_k, d: delta) + m^2 dd(phi.alt_k, d: delta) = 0
$
with $m^2 equiv V_(phi.alt phi.alt)$. We can rewrite this using conformal time
$
  eta = integral dd(t)/a = integral dd(a)/(a^2 H)
$
Assuming _plain_ de Sitter with $H tilde "constant"$ we have
$
  eta = -1/(a H)
$
With slow-roll inflation we can consider time scales of order one _Hubble time_ $dd(t, d: Delta) = H^(-1)$ and perform a slow-roll expansion
$
  eta & = integral e^(-H t) dd(t) \
      & = integral exp[-(H_* + dot(H)_* /H_*) t] dd(t) \
      & = integral exp[-H_* (1 + dot(H)_* / H_*^2) t] dd(t) \
      & = integral e^(-H_* (1 - epsilon.alt_*) t) dd(t) \
      & = - 1/(a H_* (1-epsilon.alt_*))
$
Then to leading order in the slow-roll expansion
$
  eta = - 1/(a H (1- epsilon.alt))
$
We also introduce the canonical field $chi = a phi.alt$ giving
$
  pdv(dd(chi_k, d: delta), eta, 2) + (k^2 + m^2 a^2 - 1/a pdv(a, eta, 2)) dd(chi_k, d: delta) = 0
$
This is the equation for a SHO with
$
  M^2 (eta) equiv m^2 a^2 - 1/a pdv(a, eta, 2)
$
Note,
$
  dot.double(a)/a = dot(H) + H^2";  " epsilon.alt = - dot(H)/H^2
$
Then
$
  pdv(dd(chi_k, d: delta), eta, 2) + (k^2 + M^2) dd(chi_k, d: delta) = 0
$
with
$
  M^2(eta) = 1/eta^2 (m^2/H^2 - 2 - 3 epsilon.alt) + cal(O)(epsilon.alt^2)
$
We define $nu_phi.alt$ by
$
  M^2 = - 1/eta^2 (nu_phi.alt^2 - 1/4) => nu_phi.alt^2 = 9/4 + 3 epsilon.alt - m^2/H^2
$
Then the equation of motion becomes a _Bessel equation_
$
  pdv(dd(chi_k, d: delta), eta, 2) + [k^2 - 1/eta^2 (nu_phi.alt^2-1/4)] dd(chi_k, d: delta) = 0
$
with solutions in terms of the _Hankel functions_ for $nu_phi.alt in RR$
$
  dd(chi_k, d: delta) = sqrt(-eta) [A_k H_nu^((1)) (-k eta) + B_k H_nu^((2)) (-k eta)]
$
The Hankel functions are defined in terms of the Bessel functions of the first and second kind as
$
  H_nu^((1)) (x) = J_nu (x) + i Y_nu (x)";  " H_nu^((2)) (x) = J_nu (x) - i Y_nu (x)
$
We can then quantize the fluctuation field by promoting $dd(chi_k, d: delta)$ to an operator and expanding as
$
  dd(hat(chi)_k, d: delta) = dd(chi_k, d: delta) a_k + dd(chi^*_k, d: delta) a_(-k)^dagger
$

== The $m=0$ scalar field in pure de Sitter
Assuming the scalar field is massless and we are in pure de Sitter then $v_phi.alt^2 = 9/4$ and
$
  dd(chi_k, d: delta) = A_k (1 + i/(k eta)) e^(-i k eta) + B_k (1-i/(k eta)) e^(i k eta)
$
On sub-horizon scales $k >> a H$ so $- eta k >> 1$ and
$
  dd(chi_k, d: delta) tilde.eq A_k e^(-i k eta) + B_k e^(i k eta)
$
We could take
$
  A_k = 1/sqrt(2k)";  " B_k = 0
$
and recover the Minkowski vacuum in the infinite past
$
  dd(chi_k, d: delta) = 1/sqrt(2 k) e^(-i k eta)";  " eta -> - oo
$
This is also the $0^"th"$ order adiabatic vacuum, and we recover the usual dispersion relation
$
  omega_k^2 (eta) = k^2 + M^2 = k^2 - 2/eta^2 -> k^2
$
Note, the instantaneous vacuum only exists for $omega_k^2 > 0$
$
  k^2 - 2/eta^2 > 0
$
meaning only on sub-horizon scales with $(k eta)^2 gt.tilde 2$ can the energy of a mode be normalized.#footnote[Then for $eta_0 -> -oo$ all modes can be taken in the instantaneous minimal configuration.] The instantaneous vacuum was determined by
$
  V_k (eta_0) = 1/sqrt(2 omega_k (eta_0))";  " 1/(V_k (eta_0)) pdv(V_k (eta_0), eta) = -i omega_k (eta_0)
$
Using $omega_k^2 (eta_0) -> k^2$ as $eta_0 -> -oo$ we immediately find
$
  dd(chi_k, d: delta) = 1/sqrt(2 k) e^(-i k eta)";  " eta -> -oo
$
We call this the _Bunch-Davies vacuum_.#footnote[This is simply the instantaneous vacuum in de Sitter for $eta -> -oo$.] On super-horizon scales $-k eta << 1$ we find
$
  dd(chi_k, d: delta) tilde.eq i/sqrt(2k) (a H)/k
$
implying
$
  abs(dd(phi.alt_k, d: delta)) tilde.eq H/sqrt(2 k^3)
$


== The scalar field in quasi de Sitter
We consider the solutions valid for slow-roll inflation
$
  dd(chi_k, d: delta) = sqrt(-eta) [A_k H_nu^((1)) (-k eta) + B_k H_nu^((2)) (-k eta)]
$
We consider the sub-horizon limit $-k eta >> 1$. Note, in this limit
$
  H_nu^((1)) tilde sqrt(2/(-pi k eta)) exp[i(-k eta - (pi nu_phi.alt)/2 - pi/4)]";  " H_nu^((2)) tilde sqrt(2/(-pi k eta)) exp[-i(-k eta - (pi nu_phi.alt)/2 - pi/4)]
$
We normalize to the Bunch-Davies vacuum by requiring
$
  dd(chi_k, d: delta) = 1/sqrt(2 k) e^(-i k eta)";  " -k eta >> 1
$
implying
$
  A_k = sqrt(pi)/2 exp[i(nu_phi.alt +1/2) pi/2]";  " B_k = 0
$
Then we consider the super-horizon limit $-k eta << 1$. Note,
$
  H_nu^((1)) tilde sqrt(2/pi) e^(-i pi \/2) 2^(nu_phi.alt-3\/2) Gamma(nu_phi.alt)/Gamma(3\/2) (-k eta)^(-nu_phi.alt)
$
implying
$
  dd(chi_k, d: delta) = exp[i(nu_phi.alt-1/2)pi/2] 2^(nu_phi.alt - 3\/2) Gamma(nu_phi.alt)/Gamma(3\/2) 1/sqrt(2k) (-k eta)^(1\/2-nu_phi.alt)
$
and#footnote[We ignore some superfluous prefactors.]
$
  abs(dd(phi.alt_k, d: delta)) tilde.eq H/sqrt(2 k^3) (k/(a H))^(3\/2-nu_phi.alt)";  " -k eta << 1 => k/(a H) << 1
$
We see this reduces to the $m=0$ case when $nu_phi.alt = 3\/2$.

== The power spectrum
Consider a generic $g(bold(x),t)$ with
$
  g(bold(x),t) = integral dd(k, 3)/(2 pi)^(3\/2) e^(i bold(k) dot bold(x)) g_k (t)
$
We define the _power spectrum_ through
$
  braket(0, g_(k_1)^* g_(k_2), 0) equiv delta^((3)) (bold(k)_1-bold(k)_2) (2 pi^2)/k^3 P_g (k_1)
$
Then
$
  braket(0, abs(g^2 (bold(x),t)), 0) &= integral dd(k_1, k_2, [3,3])/(2pi)^3 e^(i (bold(k)_1-bold(k)_2) dot bold(x)) braket(0, g_(k_1)^* g_(k_2), 0) \
  &= integral dd(k_1, k_2, [3,3])/(2 pi)^3 e^(i (bold(k)_1-bold(k)_2) dot bold(x)) delta^((3)) (bold(k)_1-bold(k)_2) (2 pi^2)/k^3 P_g (k_1) \
  &= integral dd(k_2, 3)/(2 pi)^3 (2 pi^2)/k_2^3 P_g (k_2) \
  &=^"spherical" integral dd(k)/k P_g (k)
$
Then the _variance_ of $g(bold(x),t)$ is related to the power spectrum by
$
  braket(0, abs(g^2 (bold(x),t)), 0) = integral dd(k)/k P_g (k)
$
We can compute this for the scalar field
$
  braket(0, abs(dd(phi.alt, d: delta))^2, 0) &= braket(0, abs(integral dd(k, 3)/(2 pi)^(3\/2) (dd(phi.alt_k, d: delta) a_k e^(i bold(k) dot bold(x)) + dd(phi.alt_k^*, d: delta) a_k^dagger e^(-i bold(k) dot bold(x))))^2, 0) \
  &= braket(0, integral dd(k_1, k_2, [3,3])/(2pi)^3 abs(dd(phi.alt_k, d: delta))^2 a_(k_1) a_(k_2)^dagger, 0) \
  &= integral dd(k_1, k_2, 3)/(2pi)^3 abs(dd(phi.alt_k, d: delta))^2 delta^((3)) (bold(k)_1-bold(k)_2) \
  &= integral dd(k, 3)/(2 pi)^3 abs(dd(phi.alt_k, d: delta))^2 \
  &= integral dd(k)/k k^3/(2 pi^2) abs(dd(phi.alt_k, d: delta))^2 \
  &=^"by definition" integral dd(k)/k P_(dd(phi.alt, d: delta)) (k)
$
implying
$
  P_(dd(phi.alt, d: delta)) (k) = k^3/(2 pi^2) abs(dd(phi.alt_k, d: delta))^2
$
With the above we find
$
  P_(dd(phi.alt, d: delta)) (k) = (H/(2 pi))^2 (k/(a H))^(3-2 nu_phi.alt)
$
The quantity $P_(dd(phi.alt, d: delta))$ gives information on $V(phi.alt)$ and slow-roll inflation through the _spectral index_#footnote[Commonly referred to as _tilt_. This is one of the primary values we wish to determine in order to test different inflationary models.]
$
  n_(dd(phi.alt, d: delta)) -1 equiv dv(ln P_(dd(phi.alt, d: delta)), ln k) = 3-2 nu_phi.alt
$
Also,
$
  nu_phi.alt^2 & = 9/4 + 3 epsilon.alt - m^2/H^2 \
  & =^(m^2 = V_(phi.alt phi.alt)) 9/4 + 3 epsilon.alt - 3 eta_phi.alt \
  nu_phi.alt & = 3/2 sqrt(1+4/3 epsilon.alt - 4/3 eta_phi.alt) \
  &tilde.eq 3/2 + epsilon.alt - eta_phi.alt
$
implying
$
  n_(dd(phi.alt, d: delta)) - 1 = 2 eta_phi.alt - 2 epsilon.alt
$
which is quite nice! The analysis when we include fluctuations of the metric is more complicated, however, the effect only changes the value above slightly.#footnote[We find $n_(dd(phi.alt, d: delta)) -1 = 2 eta_phi.alt - 6 epsilon.alt$. See below.]

#pagebreak()
= Metric Perturbations#footnote[See https://arxiv.org/pdf/hep-ph/0210162 for a review.]
Above we considered fluctuations of the scalar field $dd(phi.alt, d: delta)$. We ignored their effect on $T_(mu nu)$, however, they will induce fluctuations $dd(T_(mu nu), d: delta)$. Assuming $phi.alt$ dominates#footnote[Which it would for inflation.] the EFE then imply we will have perturbations $dd(g_(mu nu), d: delta)$ since
$
  dd(cal(R)_(mu nu), d: delta) - 1/2 dd((g_(mu nu) cal(R)), d: delta) = 8 pi G_N dd(T_(mu nu), d: delta)
$
But, this also goes the other way, so
$
  dd(phi.alt, d: delta) <=> dd(g_(mu nu), d: delta)
$

== The linear metric
We can expand $g_(mu nu)$ as
$
  g_(mu nu) = g_(mu nu)^((0)) + dd(g_(mu nu) (bold(x),t), d: delta)
$
with $g_(mu nu)^((0))$ being the FRW metric. We write
$
  g_(mu nu) (bold(x),eta) = a^2 (eta) (eta_(mu nu) + h_(mu nu))
$
with
$
  g^(mu nu) (bold(x),eta) = a^(-2) (eta) (eta^(mu nu) + h^(mu nu))
$
where $h^(mu nu) = eta^(mu rho) eta^(nu alpha) h_(rho alpha)$. We parametrize $h_(mu nu)$ as#footnote[This parametrizes the metric degrees of freedom.]
$
  h_(mu nu) equiv mat(-2 phi.alt, B_i; B_i, -2 psi delta_(i j)+ E_(i j))
$
with $phi.alt$ and $psi$ being scalars, $B_i$ being a vector and $E_(i j)$ being a traceless tensor.

We can decompose the vectors using the _Helmholtz decomposition_
$
  u_i = partial_i V + V_i
$
where $partial_i V$ has no curl and $V_i$ has no divergence. Then
$
  B_i = B_i^"scalar" + B_i^"vector"
$
with
$
  B_i^"scalar" = partial_i B";  " partial^i B_i^"vector" = 0
$
Likewise
$
  E_(i j) = E_(i j)^"scalar" + E_(i j)^"vector" + E_(i j)^"tensor"
$
with
$
  E_(i j)^"scalar" & = [partial_i partial_j - delta_(i j)/3 laplacian] E \
  E_(i j)^"vector" & = -1/2 (partial_j E_i + partial_i E_j)
$
We want $E_(i j)$ to be traceless. This implies $E_i$ has no divergence $partial^i E_i = 0$. Also,
$
  partial_i E_(i j)^"tensor" = 0";  " tensor(E^"tensor", i, -i) = 0
$
We can now determine the total degrees of freedom
$
  "dof" & = "scalars" + underbracket("vector", B_i^"vector" "and" E_i) times 2 + underbracket("tensor", E_(i j)^"tensor") = 10
$
Note, at linear order these evolve independently.

== Gauge transformations
Consider the change of coordinates due to an infinitesimal gauge transformation
$
  tilde(x)^mu = x^mu + dd(x^mu, d: delta)
$
We decompose $dd(x^mu, d: delta) = xi^mu (x^nu)$ as
$
  dd(x^0, d: delta) = xi^0 (x^mu)";  " dd(x^i, d: delta) = xi^i (x^mu) + partial^i xi (x^mu)";  " partial_i xi^i = 0
$
Then we find two scalar gauge degrees of freedom $xi^0$ and $xi$. Which we can use to eliminate two scalar degrees of freedom above.

Consider a scalar of the form
$
  q (eta,bold(x)) = q_0 (eta) + dd(q(eta,bold(x)), d: delta)
$
and apply the gauge transformation $eta -> tilde(eta) = eta + xi^0 (eta,bold(x))$#footnote[Comes from $tilde(q) (x) = q(x - xi) = q(x) - xi^mu partial_mu q$]
$
  dd(q(tilde(eta),bold(x)), d: delta) = dd(q(eta,bold(x)), d: delta) + pdv(q_0 (eta), eta) xi^0
$
Then $dd(s^2) =^! dd(tilde(s)^2)$ implies
$
  a^2 (tilde(x)^0) (1+ 2 tilde(phi.alt)) (dd(tilde(x)^0))^2 =^! a^2 (x^0) (1+2 phi.alt) (dd(x^0))^2
$
with
$
  a^2 (tilde(x)^0)^2 & = a^2 (x^0) + 2 a partial_eta a xi^0 \
  dd(tilde(x)^0) & =^"by chain-rule" (1 + partial_eta xi^0 ) dd(x^0) + partial_i xi^0 dd(x^i)
$
implying
$
  tilde(phi.alt) = phi.alt - pdv(xi^0, eta) - 1/a pdv(a, eta) xi^0
$
Likewise,
$
    tilde(B) & = B + xi^0 + pdv(xi, eta) \
  tilde(psi) & = psi - 1/3 laplacian xi + 1/a pdv(a, eta) xi^0 \
    tilde(E) & = E + 2 xi
$

== Gauge invariants

We would like gauge invariant quantities since these are nice to work with. An example are the _Bardeen potentials_
$
  Phi & = - phi.alt + 1/a pdv(, eta) [(-B + 1/2 pdv(E, eta)) a] \
  Psi & = -psi - 1/6 laplacian E + 1/a pdv(a, eta) (B-1/2 pdv(E, eta))
$
Also, consider#footnote[Note, we assume $k = 0$.]
$
  ""^((3)) cal(R) = 4/a^2 laplacian psi
$
which is the _spatial curvature_ on hypersurfaces of $eta = "constant"$. However, the _curvature perturbation_ $psi$ does depend on the gauge. We instead consider#footnote[We use $phi.alt$ to mean the scalar field.]
$
  cal(R) = psi + 1/a pdv(a, eta) (pdv(phi.alt, eta))^(-1) dd(phi.alt, d: delta) = psi + H dd(phi.alt, d: delta)/dot(phi.alt)
$
which is gauge invariant. We can pick the _comoving gauge_ $dd(phi.alt, d: delta) = 0$ where
$
  cal(R) = evaluated(psi)_(dd(phi.alt, d: delta) = 0)
$
Then $cal(R)$ is actually the curvature perturbation on a comoving slice! Likewise we can define
$
  zeta = psi + 1/a pdv(a, eta) (pdv(rho, eta))^(-1) dd(rho, d: delta) = psi + H dd(rho, d: delta)/dot(rho)
$
and pick the _uniform density gauge_ $dd(rho, d: delta) = 0$ where
$
  zeta = evaluated(psi)_(dd(rho, d: delta) = 0)
$
Then $zeta$ is the curvature perturbation on spatial slices of uniform $rho$.

We would like to understand the construction of these quantities. We do this by finding the gauge transformation bringing us to the comoving gauge by picking $xi^0$ as to cancel $dd(phi.alt, d: delta)$. We have
$
  psi -> psi + 1/a pdv(a, eta) dd(tau, d: delta)
$
under $t -> t + dd(tau, d: delta)$. Where $dd(tau, d: delta) equiv xi^0$ and $xi = 0$. Also,
$
  dd(phi.alt, d: delta) &-> dd(phi.alt, d: delta) - pdv(phi.alt, eta) dd(tau, d: delta) = 0 => dd(tau, d: delta) = (pdv(phi.alt, eta))^(-1) dd(phi.alt, d: delta) \
  dd(rho, d: delta) &-> dd(rho, d: delta) - pdv(rho, eta) dd(tau, d: delta) = 0 => dd(tau, d: delta) = (pdv(rho, eta))^(-1) dd(rho, d: delta)
$
This implies $psi$ transforms as
$
  psi & -> psi + 1/a pdv(a, eta) (pdv(phi.alt, eta))^(-1) dd(phi.alt, d: delta) \
  psi & -> psi + 1/a pdv(a, eta) (pdv(rho, eta))^(-1) dd(rho, d: delta)
$
as we had above.

We can also pick the _flat gauge_ $psi_"flat" = 0$. Then
$
  dd(tau, d: delta) = - psi/scr(H)
$
where $scr(H)$ is a _Hubble-like_ parameter
$
  scr(H) equiv 1/a pdv(a, eta)
$
implying
$
  dd(phi.alt, d: delta) -> dd(phi.alt, d: delta) + 1/scr(H) pdv(phi.alt, eta) psi equiv Q
$
with $Q$ being the _Mukhanov-Sasaki variable_
$
  Q & = dd(phi.alt, d: delta) + 1/scr(H) pdv(phi.alt, eta) psi \
    & = dd(phi.alt, d: delta) + dot(phi.alt)/H psi \
    & equiv dot(phi.alt)/H cal(R)
$
which is gauge invariant and
$
  Q = evaluated(dd(phi.alt, d: delta))_(psi = 0)
$

== The longitudinal gauge
We choose $B = E = 0$. This is called the _longitudinal gauge_. The Bardeen potentials become
$
  Phi =- phi.alt";  " Psi =- psi
$
so $phi.alt$ and $psi$ become gauge invariant. Also, assuming $T_(i j) = 0$ for $i eq.not j$ the EFE imply
$
  psi = phi.alt
$
We also find from $T_(0 i)$#footnote[We now use $phi.alt$ to mean the scalar field.]
$
  dot(psi) + H psi & = 4 pi G_N dot(phi.alt) dd(phi.alt, d: delta) \
                   & = epsilon.alt H^2 dd(phi.alt, d: delta)/dot(phi.alt)
$
and from $T_00$
$
  -3 H (dot(psi) + H psi) + (laplacian psi)/a^2 = 4 pi G_N [dot(phi.alt) dd(dot(phi.alt), d: delta) - dot(phi.alt)^2 psi - V_phi.alt dd(phi.alt, d: delta)]
$
and from $T_(i i)$
$
  -[2 dot.double(a)/a + dot(a)^2/a^2] psi - 3 H dot(psi)- dot.double(psi) = - [dot(phi.alt) dd(dot(phi.alt), d: delta) - dot(phi.alt)^2 psi -V_phi.alt dd(phi.alt, d: delta)]
$
Then using
$
  dot(H) = 4 pi G_N dot(phi.alt)^2";  " dot.double(phi.alt) + 3 H dot(phi.alt) + V_phi.alt = 0
$
we can obtain the trace of the EFE
$
  dot.double(psi)_k + [H+2 dot.double(phi.alt)/phi.alt] dot(psi)_k + 2 [dot(H) - H dot.double(phi.alt)/dot(phi.alt)] psi_k + k^2/a^2 psi_k = 0
$
or
$
  pdv(psi_k, eta, 2) + 2 scr(H) (eta_phi.alt - epsilon.alt) pdv(psi_k, eta) + 2 scr(H)^2 (eta_phi.alt - 2 epsilon.alt) psi_k + k^2 psi_k = 0
$
Assuming super-horizon scales and $epsilon.alt, eta_phi.alt -> 0$ we have
$
  pdv(psi_k, eta, 2) = 0 => a^2 [H dot(psi)_k + dot.double(psi)_k] = 0
$
implying $dot(psi)_k tilde.eq dot.double(psi)_k tilde.eq 0$. Then
$
  psi_k = epsilon.alt H dd(phi.alt_k, d: delta)/dot(phi.alt)
$
and
$
  cal(R)_k tilde.eq H dd(phi.alt_k, d: delta)/dot(phi.alt)
$
Then
$
  P_cal(R) (k) = k^3/(2 pi^2) H^2/dot(phi.alt)^2 abs(dd(phi.alt_k, d: delta))^2 = k^3/(4 pi^2) 1/M_p 1/epsilon.alt abs(dd(phi.alt_k, d: delta))^2
$
Also, the Klein-Gordon equation becomes
$
  dd(dot.double(phi.alt)_k, d: delta) + 3 H dd(dot(phi.alt)_k, d: delta) + k^2/a^2 dd(phi.alt_k) + V_(phi.alt phi.alt) dd(phi.alt_k, d: delta) = underbracket(-2 psi_k V_phi.alt + 4 dot(psi)_k dot(phi.alt)_c, "terms from "dd(g, d: delta))
$
with
$
  abs(4 dot(psi)_k dot(phi.alt)_c) << abs(psi_k V_phi.alt)";  " psi_k = epsilon.alt H dd(phi.alt_k, d: delta)/dot(phi.alt)_c";  " V_phi.alt tilde.eq - 3 H dot(phi.alt)_c
$
implying
$
  dd(dot.double(phi.alt)_k, d: delta) + 3 H dd(dot(phi.alt)_k, d: delta) + (V_(phi.alt phi.alt) + 6 epsilon.alt H^2) dd(phi.alt_k, d: delta) = 0
$
Then with $chi_k equiv dd(phi.alt_k, d: delta)\/a$ we find
$
  pdv(dd(chi_k, d: delta), eta, 2) - 1/eta^2 (nu^2 - 1/4) dd(chi_k, d: delta) = 0
$
with
$
  nu^2 = 9/4 + 9 epsilon.alt - 3 eta_phi.alt
$
Then on super-horizon scales
$
  abs(dd(phi.alt, d: delta)_k) tilde.eq H/sqrt(2 k^3) (k/(a H))^(3/2 - nu)
$
as before. This implies
$
  P_cal(R) (k) = 1/(2 M_p^2 epsilon.alt) (H/(2 pi))^2 (k/(a H))^(n_cal(R)-1)
$
where
$
  n_cal(R) - 1 = dv(ln P_cal(R), ln k) = 3 - 2 nu = 2 eta_phi.alt - 6 epsilon.alt
$
Also, we can see
$
  dot(cal(R))_k tilde.eq dot(zeta)_k -> 0
$
implying the comoving perturbation is conserved on super-horizon scales. Note, we can compute
$
  P_Phi & =^"radiation dominated" (2/3)^2 P_cal(R) \
        & =^"matter dominated" (3/5)^2 P_cal(R)
$
and
$
  dd(T, d: delta)/T = 1/3 Phi = 1/5 cal(R)
$
