//**** init-ting
#import "@preview/physica:0.9.5": *
#import "chpt-temp.typ": *

#show: thmrules.with(qed-symbol: $square$)
#show: chpt-note.with()

= Structures and formation
= Baryons and photons
Baryons can interact with photons through emission, absorption, scattering etc. which complicates structure formation---especially before and after decoupling.

== Overdensities
Today we observe
$
  rho_("bary",0) tilde 4.2 times 10^(-28)" kg"" m"^(-3)
$
but in our solar neighborhood we find an overdensity
$
  delta_"sn" = (rho_"sn"-rho_("bary",0))/rho_("bary",0) tilde 2 times 10^7
$
and the Sun has $delta_dot.circle tilde 3 times 10^30$. Most of the observed $rho_("bary",0)$ is in intergalactic gas ($tilde 85%$), about half is diffuse with $T < 10^5 "K"$ and $delta <= 0$ while the rest is warm-hot with $10^5 "K" < T < 10^7 "K"$ and $3 < delta < 300$---the other $tilde 15%$ is in stars, galaxies, clusters and the gas between these. Another characteristic of the intergalactic gas is that $X tilde 1$, this is unexpected since in the early universe we saw that $X$ fell to $tilde 0.1$ around $z tilde 1200$. So something must have reionized the hydrogen in the universe---to know what this something is then we need to know a more precise time for when this something happened.

== Reionization
We turn to the CMB for information. The reionized gas, and the free electrons within, acts to scatter CMB photons. This happens at a rate
$
  Gamma = n_e sigma_e c
$
if reionization started at $t_*$ then
$
  tau_* = integral_(t_*)^t_0 Gamma (t) dd(t) = c sigma_e integral_(t_*)^(t_0) n_e (t) dd(t)
$
if $tau_* >> 1$ all CMB photons would be scattered many times, and thus give no information---the CMB would be smeared. We observe great detail in the CMB so it seems like $tau_* < 1$---in fact we do observe smearing, and we see slight suppression at high $l$. This is consistent with $tau_* = 0.066 plus.minus 0.016$.

We assume all baryonic matter is hydrogen $n_H + n_p = n_"bary" = n_("bary",0) a^(-3)$---if we assume all hydrogen reionizes at $t=t_*$ instantly then
$
  n_e = n_p = n_("bary",0) a^(-3)
$
for $t > t_*$ while $n_e=0$ for $t < t_*$---so
$
  tau_* = underbracket(Gamma_0, c sigma_e n_("bary",0)) integral_(t_*)^(t_0) dd(t)/(a(t)^3) = Gamma_0 integral_(a(t_*))^1 dd(a)/(dot(a)a^3) = Gamma_0 integral_(a(t_*))^1 dd(a)/(H(a) a^4) = Gamma_0 integral_0^z_* ((1+z)^2 dd(z))/(H(z))
$
in a lambda-matter universe
$
  H(z) = H_0 [Omega_(m,0) (1+z)^3 + Omega_(Lambda,0)]^(1\/2)
$
giving
$
  tau_* = 2/(3 Omega_(m,0)) Gamma_0/H_0 ([Omega_(m,0) (1+z_*)^3 + Omega_(Lambda,0)]^(1\/2)-1)
$
and we obtain $z_* = 7.8 plus.minus 1.3 => t_* tilde 650 "Myr"$---a very short time.

One way to ionize hydrogen is with photons of energy $h f > 13.6 "eV"$. The most apparent sources for these photons are galaxies and active galactic nuclei---both of which would have been present at $t_*$.

== First structures
The time before formation of the first stars and AGN is collectively known as the dark ages (absence of starlight). It is estimated that the first stars began to form around $z tilde 50 => t tilde 50 "Myr"$, and in turn they began to reionize gas---which at $z tilde 8$ began to merge and form a single ionized intergalactic medium, mostly due to luminous stars (since quasars and lower-luminosity AGN are insufficient).

=== Galaxies
Stars tend to be contained within galaxies, which are simply a relatively small concentration of stars and interstellar gas within a larger halo, consisting of mainly dark matter.

One can define the luminosity function
$
  Phi(L) dd(L) = Phi^* (L/L^*)^alpha exp(- L/L^*) dd(L)/L^*
$
telling us the number density of galaxies with luminosity in the range $L -> L + dd(L)$. It shows that galaxies with $L > L^*$ are rare---since they need to be large which usually required cannibalizing other galaxies---this also _limits_ the baryonic mass.

To see why this is the case we consider a simple model of galaxy formation: Take a spherical overdense region at $t_"rm"$, this will eventually become a galaxy. The mass of the sphere is
$
  M = (4 pi)/3 rho_m (t) [1 + delta(t)]R(t)^3
$
with $rho_m (t) = rho_(m,0) (1+z)^3$. Initially $delta(t_"rm") = delta_"rm" << 1$, and the sphere's expansion is dominated by Hubble flow. At some point it will begin to collapse when $delta(t_"coll") tilde 1$. Given $delta prop a prop t^(2\/3)$ we can approximate $t_"coll" tilde delta_"rm"^(-3\/2) t_"rm"$ or $ 1+z_"coll" tilde delta_"rm" (1+z_"rm") $
when the sphere starts collapsing
$
  overline(rho) (t_"coll") tilde 2 rho_m (t_"coll") tilde 2 rho_(m,0) (1+z_"coll")^3
$
after collapsing it will reach equilibrium with $R_"halo" tilde R(t_"coll")\/2$---this process is called _virilization_---note $overline(rho)_"halo" = 8 overline(rho)(t_"halo")$, so from the density you can determine the time of collapse. If we treat the gas as ideal we have
$
  P_"gas" = (rho_"gas" k_B T_"gas")/(mu)
$
and for it to be in equilibrium
$
  dv(P_"gas", r) = - (G M(r) rho_"gas" (r))/r^2
$
combining these immediately give
$
  M(r) = (k_B T_"gas" (r) r)/(G mu) [- dv(ln rho_"gas", ln r) - dv(ln T_"gas", ln r)]
$
assuming a simple power law $rho_"gas" prop r^(-beta)$ and a uniform temperatur $T_"gas"$ with $r = R_"halo"$ and $M(r) = M_"tot"$ gives
$
  k_B T_"gas" = (G M_"tot" mu)/(beta R_"halo")
$
with
$
  R_"halo" = ((3 M_"tot")/(4 pi overline(rho)_"halo"))^(1\/3)
$
we can write
$
  k T_"gas" = 4/beta (pi/3)^(1\/3) G mu rho_(m,0)^(1\/3) M_"tot"^(2\/3) (1+z_"coll")
$
importantly this depends on $M_"tot"$ which is problematic for massive halos. To form a dense galaxy at the center of a halo the gas must be able to cool by radiating light---else it will stay supported. In general halos with $T_"gas" < 10^6 "K"$ can cool quickly, since their hydrogen and helium is not fully ionized. While halos with $T_"gas" > 10^6 "K"$ are fully ionized, meaning they primarily cool by bremsstrahlung (free-free emission)---acceleration of electrons by free positive ions. For halos with $M_"tot" > 10^12 M_dot.circle$ it is statistically unlikely that they have collapsed---since they cool too slowly. However, this model is obviously simplistic and in reality not all gas is heated to $T_"gas"$ leading to _cold flows_ of gas which collapse quicker leading to the formation of galaxies at higher redshifts. But this model does explain why the hot gas between clusters doesn't form a single massive galaxy.

=== Stars
What happens now to the cooled gas as it collapses? By observation we know it doesn't form one giant object or many, many tiny objects---instead it forms stars with $M_* tilde 0.1 M_dot.circle$. This happens as cores of molecular clouds (regions of dense and cold gas) collapse. Related to this are multiple scales:
$
  t_"dyn" &= 1/(4 pi G rho_"core")^(1\/2) \
  c_s &= ((k_B T_"core")/mu)^(1\/2) \
  lambda_J &= 2 pi c_s t_"dyn" => M_J = (4 pi rho_"core")/3 lambda_J^3 prop rho^(-1\/2) T^(3\/2)
$
with objects smaller than the Jeans mass being pressure-supported---naively this indicates that molecular clouds with $M < 15 M_dot.circle$ cannot collapse and form stars, which is nonsense. To see where this fails we consider a cloud collapsing---as this happens it's thermal energy increases as
$
  dv(E_"core", t) = - P_"core" dv(V_"core", t) = - 3 N k_B T_"core" (1/R_"core" dv(R_"core", t)) = - 2 E_"core" (1/R_"core" dv(R_"core", t))
$
if $T_"core"$ is constant then it must radiate (at $T_"core" = 20"K"$)
$
  L_"core"^"iso" = - dv(E_"core", t) tilde (2 E_"core")/t_"dyn" tilde 0.015 L_dot.circle (M_"core"/(15 M_dot.circle))(rho_"core"/(10^(-15) "kg"" m"^(-3)))^(1\/2)
$
assuming the cloud acts like a blackbody the luminosity is
$
  L_"core" = 4 pi R_"core"^2 underbracket(f_e, "correction") sigma_"SB" T_"core"^4
$
then (at $T_"core" = 20 "K"$)
$
  L_"core" tilde 1100 L_dot.circle f_e (M_"core"/(15 M_dot.circle))^(2\/3) (rho_"core"/(10^(-15) "kg"" m"^(-3)))^(-2\/3)
$
for $f_e > 10^(-5)$ this is larger than $L_"core"^"iso"$ meaning the temperature is held constant. This means that the Jeans mass $M_J$ decreases as the core collapses and $rho$ increases while $T$ stays constant. Eventually the core becomes unstable and splits to form a pair of fragments which continue to collapse until they split etc. This hierarchical fragmentation naturally produces a power-law distribution of $M_*$ (given it is a bit inefficient).

After $n$ rounds of fragmentation each fragment will have $M_"frag" = 2^(-n) M_"core"$ and $rho_"frag" = 2^(2n) rho_"core"$. This relation ensures that $L_"core"^"iso"$ is the same for every fragment, so at some point this process stops since luminosity depends on size. The maximal is found to be $n = 9$ (where $f_e tilde 3$ which is unphysical). These fragments end up becoming protostars if they are massive enough $M > 0.08 M_dot.circle$.

The above process relies on $f_e$ and the molecular clouds being able to cool through fragmentation. If $f_e$ is very small due to e.g. small amounts of dust (as in the very early universe) then the process would stop much earlier. For this reason we expect early stars to have been much larger and comparable to $M_J$. Again this model is very simplistic.

= Active galaxies
== Basic observations
Observations show activity in the centers of young, remote galaxies that isn't found in nearer galactic nuclei---these are broadly called active galaxies.

A notably feature of many active galactic nuclei (AGN) the heart of an active galaxy, is that their spectral energy density is very persistent over many frequencies---unlike a blackbody. Generally their spectra seem like they can be decomposed into thermal- (big blue bump and infrared bump) and nonthermal sources. The bumps are believed to be caused by an optically thick accretion disk and dust grains respectively. The spectra also generally show a turnover at lower frequencies which is steep for radio quite AGN and less steep for radio loud AGN---which could be due to synchrotron radiation or other thermal- and nonthermal emission.

Typically broad-line emission is variable, while narrow-line emission is non-variable.

=== Seyferts
The first active galaxies discovered appeared to have very bright nuclei and broad emission lines. Today these are called Seyfert galaxies---with Seyfert 1 having broad and narrow emission lines, while Seyfert 2 only have narrow emission lines. All Seyfert galaxies spectra also show a continuum originating from a central source---in Seyfert 1 galaxies the luminosity due to this typically overwhelms the light from the entire galaxy.

Generally these are quite variable---e.g. the Seyfert 1 X-ray emission changes hourly---and most are spiral galaxies.

=== Radio galaxies
As their name implies radio galaxies are very bright at radio wavelengths. We again get two types: broad-line radio galaxies ($tilde$ Seyfert 1) and narrow-line radio galaxies ($tilde$ Seyfert 2). BLRGs have bright, starlike nuclei surrounded by a faint envelope, while NLRGs are giant elliptical galaxies. These are similar to Seyferts, but Seyferts are relatively quiet at radio wavelengths, and nearly all Seyferts are spiral galaxies---all strong radio galaxies are elliptical.

Radio galaxies sometimes display very large extended radio lobes---the optical galaxy is flanked by these---it may radiate its energy from the nucleus and from a halo with a size similar to the optical galaxy. These are sometimes connected to the galaxy through a jet---strong galaxies typically have one-sided jets (with a dim counterjet), while weaker galaxies have two-sided jets.

=== Quasars
Quasars or quasi-stellar radio sources were initially discovered as the optical counterpart to strong radio sources (either from the core or from lobes). Quasars typically have very high redshift, meaning optically they appear as very bright, starlike nuclei with fuzzy halos (the parent galaxy)---meaning they are extremely luminous.

They emit an excess of UV, meaning they are relatively blue (big blue bump). Their spectra also have absorption lines which can be used as probes (see later). They have broad and narrow emission lines. Note that quasi-stellar objects or radio-quiet quasars also exist, and are in abundance.

== Unified model of AGN
We'd like a unified model able to describe all of the above. The general consensus is that all active galaxies are fueled by accretion onto a central supermassive black hole---with the observational differences being orientation dependent.

=== Black hole & accretion
One can compute that due to high variability the size of any AGN must be $R tilde 10 "AU"$ which is very small. Further since we have found very luminous objects we can try to estimate the size, this is done using the Eddington limit
$
  L_"ed" = (4 pi G c)/overline(kappa) M tilde 1.5 times 10^31 "W" (M/M_dot.circle)
$
which is the maximal luminosity a spherically symmetric object can have and still be stable. For a typical quasar with $L =5 times 10^(39) "W"$ this gives a minimal mass
$
  M = 3.3 times 10^8 M_dot.circle
$
which is enourmous---and suggests that AGN are supermassive black holes.

The luminosity generated by a black hole can be written as
$
  L_"disk" = eta dot(M) c^2
$
the disk temperature can be written as ($R = G M\/c^2$)
$
  T_"disk" = ((3 c^6 dot(M))/(8 pi sigma G^2 M^2))^(1\/4)
$
defining $f_"ed" equiv L_"disk"\/L_"ed"$ we can find
$
  T_"disk" = ((3 c^5 f_"ed")/(2 overline(kappa) sigma G M eta))^(1\/4)
$
so $T_"disk" prop M^(-1\/4)$---it is believed that this causes the big blue bump. Making a model for the accretion disk is difficult, one suggestion is a three-part disk. In the first part, closest to the black hole, radiation pressure exceeds gas pressure leading to a thick, hot disk. Outside this is a thinner disk supported by gas pressure, which becomes thicker at larger distances---this flaring means it can be irradiated by the thick disk. Beyond this the disk breaks up into smaller clouds.

=== Spectra
The continuous synchrotron spectra of AGNs could be described by the magnetic field of the black hole---and the ionized diskrotating could induce a large electric field, leading to radiation (using accretion energy). Blandford-Znajek suggests that the black hole itself acts like a spinning conductor, creating an emf by taking rotational energy from the black hole---this energy could be in the form of radiation. X-rays can also be created by low-energy photons inverse Compton scattering with relativistic electrons---this can also produce gamma rays. The generation of charged particles, their acceleration and subsequent ejection is thought to be the origin of jets and accompanying radio lobes.

The accretion disk, obscuring torus and jets all contribute to the continuum spectra of AGNs---essentially a power law $F_nu prop nu^alpha$.

The broad- and narrow emission lines of AGN spectra are the result of photoionization by the continuum radiation---the nature of these lines give clues about where they arise. They are formed in the broad-line and narrow-line regions respectively. The broad-line region is close to the center (from observation)---and there is general agreement that it consists of clumpy partially ionized clouds of gas, these orbit the black hole very quickly leading to doppler broadening. The unified model postulates that a large optically thick obscuring torus of gas and dust surrounds the cloud of gas in the broad-line region. This torus acts to conceal the broad-line region and the central source form direct view (e.g. when observing Seyfert 2). Outside this torus is the narrow-line region, and is probably composed of a spherical-ish distribution of clouds---observation indicates these are moving radially outward from the center. Some of these can be ionized by the center, others can't due to the torus.

So in a typical spectrum we'd have a continuum with the addition of the broad- and narrow emission lines.

== Quasars as probes
As mentioned quasars have very high redshift, and are thus very far away. We can then use them as probes since their light is affected by everything between us and them.

One way this is done is through gravitational lensing.

Another way we'll focus on is the Lyman-$alpha$ forest. High-redshift quasars show many narrow absorption lines superimposed on the quasar's continuum spectrum. These come when light from a quasar passed through some material along the line of sight. If this material is far from us then these lines appear redshifted. So if light passes through multiple things we'll see multiple lines redshifted by different amounts. The Lyman-$alpha$ forest is a dense collection of hydrogen absorption lines. Ionized metals also lead to absorption lines (from galaxies, star formation, etc.). And both types are generally found in the UV, since the material is moving at very high speeds. What one observes is a very strong Lyman-$alpha$ emission line by the quasar itself and then many, many Lyman-$alpha$ absorption lines at shorter wavelengths (smaller redshifts). One can combine this with lensing to deduce things like size of intergalactic clouds, etc.
