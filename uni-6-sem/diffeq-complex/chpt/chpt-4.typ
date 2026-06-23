#import "../../temp.typ": *
#show: chpt-note.with()

= Complex Analysis
We extend $RR$ to a _field_ wherein the solution to
$
  x^2 + 1 = 0
$
exists. To this end we define $i^2 equiv -1$ and form expressions $z = a + i b in CC$ with $a, b in RR$.#footnote[Artihmetic operations are defined as usual and it is trivial to show that $CC$ forms a field.] We will assume basic knowledge of complex numbers.

== Complex functions
Let $z = x + i y in D subset.eq CC$. Then $f(z)$ is a _complex function_ of the _complex variable_ $z$ if each $z in D$ corresponds to one or more values $f(z)$. We write
$
  f(z) = underbracket(u(x,y), "real") + i underbracket(v(x,y), "imaginary")
$
and will consider _single-valued_ functions unless specified.

Let $D subset.eq CC$ be an open set with $z_0 in D$ and $f:D-> CC$ be a function. We say $f$ is _differentiable_ in $z_0$ if and only if the limit
$
  f'(z_0) equiv lim_(z->z_0) (f(z)-f(z_0))/(z-z_0) tilde "complex derivative"
$
exists and is unique.#footnote[_Unique_ meaning independent of how we take $Delta z -> 0$.] When $f$ is differentiable for all $z_0 in D$ we call $f$ _holomorphic_ on $D$. A single-valued and differentiable $f$ on $D subset.eq CC$ is called _analytic_ on $D$. When $f$ is analytic on $D$ expect at a finite number of points in $D$ then these are called _singularities_ of $f$.

#theorem(name: "Cauchy-Riemann")[
  Let $D subset.eq CC$ be an open interval, $z_0 in D$ and $f : D -> CC$ be a function. Let $z = x + i y$ and $z_0 = x_0 + i y_0$. Writing $f(x,y) = u(x,y) + i v(x,y)$ the following are equivalent:

  1. $f$ is complex differentiable in $z_0$.

  2. $f$ is _totally real differentiable_ in $z_0$ and the _Cauchy-Riemann equations_ hold
  $
    pdv(u(x_0,y_0), x) = pdv(v(x_0,y_0), y)";  " pdv(u(x_0,y_0), y) = - pdv(v(x_0,y_0), x)
  $
]

Assuming $f$ is differentiable in $z_0$ we immediately find
$
    f'(z_0) & = pdv(u, x) + i pdv(v, x) = pdv(f, x) tilde "along" x \
  i f'(z_0) & = i pdv(v, y) + pdv(u, y) = pdv(f, y) tilde "along" y
$
Where we use Cauchy-Riemann and uniqueness.

Consider the infinite series
$
  sum_(n=0)^oo a_n (z-z_0)^n
$
with $a_n in CC$. This is called a _power series_ with coefficients $a_n$ around $z_0$. We call the set of all $z$ that converge the _domain of convergence_.

#theorem[
  To all $sum_(n=0)^oo a_n (z-z_0)^n$ there exists _exactly_ one $r in [0,oo) union {oo}$ called the _radius of convergence_ with the following properties:

  1. in every _closed circle_ $abs(z-z_0) <= rho < r$ the series is _uniformly absolutely convergent_.

  2. for $abs(z-z_0) > r$ the series is _divergent_.

  3. the _Cauchy-Hadamard formula_ holds
  $
    1/r equiv lim_(n -> oo) sup (root(n, abs(a_n)))
  $
]

#theorem[
  Let $sum_(n=0)^oo a_n (z-z_0)^n$ be a power series with radius of convergence $r > 0$. Let the function $f(z)$ be specified by this series on $B_r (z_0)$#footnote[$B_r (z_0)$ is the _ball_ of radius $r$ around $z_0$.] meaning for all $z in B_r (z_0)$ we have
  $
    f(z) = sum_(n=0)^oo a_n (z-z_0)^n
  $
  Then $f$ is holomorphic on $B_r (z_0)$ and
  $
    f'(z) = sum_(n=1)^oo n a_n (z-z_0)^(n-1)
  $
  for all $z in B_r (z_0)$.
]

By induction we obtain for all $k in NN$
$
  f^((k)) (z) = sum_(n=0)^oo (n+k)!/n! a_(n+k) (z-z_0)^n
$
for all $z in B_r (z_0)$. When $z = z_0$ we find
$
  f^((k)) (z_0) = k! a_k => a_k = (f^((k)) (z_0))/k!
$
Which upon substituting gives us the _Taylor series_
$
  f(z) = sum_(n=0)^oo (f^((k)) (z_0))/k! (z-z_0)^n
$
for all $z in B_r (z_0)$.

== Elementary functions
We can now redefine familiar functions in terms of complex variables.

Consider $z in CC$ then we can write
$
  e^z = sum_(n=0)^oo z^n/n!
$
which is the _complex exponential function_. We have the following properties:

1. The radius of convergence is $r = oo$.

2. The derivative of $e^z$ is $e^z$.

3. $e^(z+w) = e^z e^w$ for all $z, w in CC$.

4. With $z = x + i y$ we have $e^z = e^x e^(i y)$ with $abs(e^z) = abs(e^x)$.

We also seek solutions $w in CC$ to
$
  e^w = z
$
The solution $w$ is only defined modulo $2 pi i$. To see this we can write $ z = r e^(i theta) $ for $- pi < theta <= pi$. Then for any $k in ZZ$ there is an equivalent expression $ z = r e^(i(theta + 2pi k)) $ The _naive logarithm_ is then
$
  w = ln z = ln r + i (theta + 2 pi k)
$
which is then infinitely multi-valued. The _principal value_ of the logarithm is
$
  "Ln"(z) & = ln r + i"Arg"(z) \
          & = ln r + i theta
$
with $-pi <= theta < pi$. The above corresponds to introducing a _branch cut_ on $RR^-$. We can now define the general _complex power_
$
  t^z = e^(z ln t)
$
for $t in CC\\{0}$. When $z = n^(-1)$ for $n in NN$ we have $n$ distinct roots of $t$. This follows from
$
  t^(1\/n) = r^(1\/n) exp[i (theta + 2 pi k)/n]
$
and the $2pi$-periodicity of $e^z$.

== Singularities and zeros
Consider $D subset.eq CC$ open. When $f$ is singular at $z = z_0$ and analytic for all $z in D \\ {z_0}$ then $z_0$ is called an _isolated singularity_ of $f$.

When $f$ has the form
$
  f(z) = g(z)/(z-z_0)^n
$
with $n in NN$ and $g(z)$ is analytic on a neighbourhood $D$ with $z_0 in D$ and $g(z_0) eq.not 0$ then $f$ has a _pole of order_ $n$ at $z = z_0$. When $n = 1$ we call the pole _simple_. Alternatively, we have
$
  lim_(z->z_0) (z-z_0)^n f(z) = a
$
with $a in CC\\{0}$. This definition has some consequences:

1. When $a$ vanishes then $z_0$ is either a pole of order less than $n$ or $f$ is analytic at $z_0$.

2. When $a -> oo$ the pole is of order higher than $n$.

3. When $f$ has a pole at $z_0$ then $abs(f) -> oo$ as $z_0$ is approached from any direction.

4. When no finite $n$ yields a value $a$ then $z_0$ is called an _essential singularity_.

To define the behaviour of $f(z)$ at _infinity_ we need to specify what _infinity_ means. We do this by considering
$
  f(z) = f(xi^(-1))
$
at $xi = 0$.

#example[
  Consider $f(z) = e^z$ then in terms of $xi$ $ f(xi^(-1)) = sum_(n=0)^oo xi^(-n)/n! $ at $xi = 0$ we have an essential singularity and therefore at $z = oo$.
]

When $f(z_0) = 0$ we call $z_0$ a _zero_ of $f$. These are classified analogously to poles, so if we can write
$
  f(z) = (z-z_0)^n g(z)
$
with $n in NN$ and $g(z_0) eq.not 0$ then $z_0$ is called a _zero of order_ $n$ of $f$. When $n = 1$ we call the zero _simple_. When $z_0$ is a zero of order $n$ of $f$ then by definition it is a pole of order $n$ of $f^(-1)$.

== Complex integrals
When $f(z)$ is a complex single-valued continuous function in $D subset.eq CC$#footnote[With $D$ open.] then we can define the _complex integral_ of $f(z)$ between two points $A$ and $B$ along some curve $Gamma subset D$. We consider a path $Gamma$ described by the parametrisation $z(t) = x(t) + i y(t)$ where $t in [alpha,beta] in RR$. Then along $Gamma$ we have
$
  cal(I) = integral_Gamma f(z) dd(z)
$
or
$
  cal(I) = integral_alpha^beta u dv(x, t) dd(t) - integral_alpha^beta v dv(y, t) dd(t) + i integral_alpha^beta u dv(y, t) dd(t) + i integral_alpha^beta v dv(x, t) dd(t)
$
implying $cal(I)$ exists when
$
  dv(x, t) " and " dv(y, t)
$
are continuous.

We note the _standard estimate_. Assume $abs(f(z)) <= M in RR^+$ along $Gamma$ then
$
  abs(integral_Gamma f(z) dd(z)) <= integral_Gamma abs(f(z))abs(dd(z)) <= M integral_Gamma dd(Gamma) = M L
$
with the length of $Gamma$ being $L$.

== Cauchy's theorem
#theorem(name: "Cuachy's integral theorem")[
  Let $D subset.eq CC$ open, $f : D -> CC$ holomorphic and $f'(z)$ continuous at all points in and on a closed contour $Gamma subset D$. Then
  $
    integral.cont_Gamma f(z) dd(z) = 0
  $
]

#proof[
  Let $p, q : RR^2 -> RR$ be functions with continuous first derivatives in and on a closed contour $Gamma subset D$ bounding a domain $R subset D$. Then
  $
    integral.double_R (pdv(p, x) + pdv(q, y)) dd(x, y) = integral.cont_Gamma (p dd(y) - q dd(x))
  $
  This is _Green's theorem_ and follows from Gauss' thoerem with $hat(n) dd(s) = hat(x) dd(y) - hat(y) dd(y)$. We also have
  $
    cal(I) &= integral.cont_Gamma (u dd(x) - v dd(y)) + i integral.cont_Gamma (v dd(x) + u dd(y)) \
    &= -integral.double_R (pdv(v, x) + pdv(u, y)) dd(x, y) + i integral.double_R (pdv(u, x) - pdv(v, y)) dd(x, y) \
    &=^"by Cauchy-Riemann" 0
  $
  Where the Cauchy-Riemann equations hold since $f$ is assumed holomorphic.
]

Consider now two closed contours in $CC$ denoted by $Gamma$ and $gamma$ where $gamma$ lies inside $Gamma$. Then if $f$ is analytic in the region between $Gamma$ and $gamma$ we have
$
  integral.cont_Gamma f(z) dd(z) = integral.cont_gamma f(z) dd(z)
$

#proof[
  We imagine cutting open $Gamma$ and $gamma$ introducing two contours $cal(C)_1$ and $cal(C)_2$ connecting $Gamma$ and $gamma$. We now have a closed _keyhole_ contour $cal(C)$ consisting of four pieces. By assumption $f$ is holomorphic on $cal(C)$ and the region it encloses to by Cauchy's theorem
  $
    integral.cont_cal(C) f dd(z) = integral_Gamma f dd(z) + integral_cal(C)_1 f dd(z) + integral_cal(C)_2 f dd(z) + integral_gamma f dd(z) = 0
  $
  We integrate along $cal(C)_1$ and $cal(C)_2$ in opposite directions so in the limit where the introduced cut vanishes these contributions cancel. Likewise, we integrate along $Gamma$ and $gamma$ in opposite directions so
  $
    integral.cont_Gamma f dd(z) = integral.cont_gamma f dd(z)
  $
]

== Cauchy's integral formula
#theorem(name: "Cauchy's integral formula")[
  Let $D subset.eq CC$ open and $f : D -> CC$ holomorphic. Let $Gamma subset D$ be a closed contour and $z_0$ a point within $Gamma$. Then
  $
    f(z_0) = 1/(2 pi i) integral.cont_Gamma f(z)/(z-z_0) dd(z)
  $
]

#proof[
  Let $gamma$ be a circle of radius $rho$ around the point $z_0$ such that $gamma$ is completely inside $Gamma$. Since $f$ is holomorphic inside $Gamma$ the integrand
  $
    f(z)/(z-z_0)
  $
  is analytic in the region between $Gamma$ and $gamma$. Then by the above corollary $ integral.cont_Gamma = integral.cont_gamma $ Any point $z$ on $gamma$ is parametrised by
  $
    z = z_0 + rho e^(i theta)
  $
  so
  $
    cal(I) & = integral.cont_gamma f(z)/(z-z_0) dd(z) \
           & = i integral_0^(2 pi) f(z_0 + rho e^(i theta))
  $
  We can take the limit $rho -> 0$ so
  $
    cal(I) = 2 pi i f(z_0)
  $
  and the result follows.
]

Then the value of any analytic function inside $Gamma$ is fully determined by the values on $Gamma$!

#theorem[
  Let $D subset.eq CC$ open and $f : D -> CC$ holomorphic. Let $Gamma subset D$ as above. When $f$ is arbitrarily often complex differentiable on $D$ then
  $
    f^((n)) (z_0) = n!/(2 pi i) integral.cont_Gamma f(z)/(z-z_0)^(n+1) dd(z)
  $
]

#proof[
  Trivial by induction.
]

This allows us to determine the value of all derivatives within $Gamma$ by knowing the values on $Gamma$.

#theorem(name: "Cauchy's inequality")[
  Let $f$ be analytic inside and on a circle $gamma$ of radius $R$ centered on the point $z_0 in CC$. If $abs(f) <= M$ on $gamma$ with $M in RR^+$ then
  $
    abs(f^((n)) (z_0)) <= (M n!)/R^n
  $
]

#proof[
  We have
  $
    abs(f^((n)) (z_0)) = n!/(2 pi) abs(integral.cont_gamma f(z)/(z-z_0)^(n+1) dd(z))
  $
  Using $abs(z-z_0) = R$ and $L_gamma = 2 pi R$ the standard estimate gives
  $
    abs(f^((n)) (z_0)) <= (M n!)/R^n
  $
]

We call a function $f$ _entire_ when it is holomorphic on $CC$.

#theorem(name: "Liouville's theorem")[
  An entire function which is bounded has to be a constant.
]

#proof[
  We consider arbitrary $z_0 in CC$. Then by Cauchy's inequality
  $
    abs(f^((1)) (z_0)) <= M/R
  $
  as $R -> oo$ we have $abs(f^((1)) (z_0)) = 0$ implying $f^((1)) (z_0) = 0$ since $M$ exists. Given $f$ is an entire function and $z_0$ is arbitrary then $f^((1)) (z) = 0$ for all $z in CC$ implying $f(z) tilde$ constant.

]

== Taylor- and Laurent series
#theorem(name: "Taylor's theorem")[
  Let $f$ be analytic inside and on a circle $gamma$ of radius $R$ centered at the point $z = z_0$ with $z$ in $gamma$. Then
  $
    f(z) = sum_(n=0)^oo a_n (z-z_0)^n";  " a_n = (f^((n)) (z_0))/n!
  $
  The _Taylor series_ above is valid inside the region of analyticity and for any $z_0$ can be shown to be unique.
]

#proof[
  Because $f$ is analytic inside $gamma$ the Cauchy integral formula is valid
  $
    f(z) = 1/(2 pi i) integral.cont_gamma f(xi)/(xi-z) dd(xi)
  $
  Recall
  $
    sum_(k=0)^oo a r^k = a/(1-r)
  $
  when $abs(r) < 1$. Then
  $
    (xi-z_0)/(xi-z) & = 1/(1-(z-z_0)/(xi-z_0)) \
                    & = sum_(n=0)^oo ((z-z_0)/(xi-z_0))^n
  $
  implying
  $
    1/(xi-z) = 1/(xi-z_0) sum_(n=0)^oo ((z-z_0)/(xi-z_0))^n
  $
  Then
  $
    f(z) &= 1/(2 pi i) integral.cont_gamma f(xi)/(xi-z_0) sum_(n=0)^oo ((z-z_0)/(xi-z_0))^n dd(xi) \
    &= 1/(2 pi i) sum_(n=0)^oo (z-z_0)^n integral.cont_gamma f(xi)/(xi-z_0)^(n+1) dd(xi) \
    &= 1/(2 pi i) sum_(n=0)^oo (z-z_0)^n 2 pi i (f^((n)) (z_0))/n! \
    &= sum_(n=0)^oo (f^((n)) (z_0))/n! (z-z_0)^n
  $
  and we are done.
]

#theorem(name: "Identity theorem")[
  Let $f$ and $g$ be analytic in some region $R$ and $f = g$ in $S subset R$ containing an _accumulation point_ $z_0$. Then $f = g$ in $R$.#footnote[The identity theorem implies _analytic continuations_ are unique.]
]

#proof[
  Consider the difference $h = f - g$. We want to show that if $h$ vanishes on $S$ then it vanishes on $R$. Let $z_0 in S$ be an accumulation point. Then
  $
    h(z) = h(z_0) + h^((1)) (z_0) (z-z_0) + 1/2 h^((2)) (z_0) (z-z_0)^2 + dots
  $
  with $h^((n)) (z_0) = f^((n)) (z_0) - g^((n)) (z_0)$. This converges on a circle $gamma$ extending to the boundary of $R$ since $h$ is analytic in $R$. Since $z_0 in S$ is an accumulation point we have
  $
    h(z_0) = h^((1)) (z_0) = h^((2)) (z_0) = dots = 0
  $
  implying $h(z) = 0$ inside $gamma$. We can repeat this by choosing a new arbitrary $z_0$ inside $gamma$. We eventually find $h(z) = 0$ in $R$.

]

When $f$ has a singularity inside $gamma$ at some $z = z_0$ we cannot write a Taylor series.

#theorem(name: "Laurent series")[
  Let $f$ have a pole of order $p$ at $z = z_0$ and be analytic at all other points inside and on $gamma$. Then
  $
    f(z) = underbracket(a_(-p)/(z-z_0)^p + dots + a_(-1)/(z-z_0), "principal part") + underbracket(a_0 + a_1 (z-z_0) + a_2 (z-z_0)^2 + dots, "analytic part")
  $
  where $a_(-p) eq.not 0$. This extension of the Taylor series is called a _Laurent series_.
]

#proof[
  By definition we can write $g(z) = (z-z_0)^p f(z)$ with $g(z_0) = b_0 eq.not 0$. This function is analytic at $z = z_0$ so we can write a Taylor series
  $
    g(z) = sum_(n=0)^oo b_n (z-z_0)^n
  $
  implying
  $
    f(z) & = 1/(z-z_0)^p sum_(n=0)^oo b_n (z-z_0)^n \
         & = sum_(n=0)^oo b_n (z-z_0)^(n-p) \
         & = sum_(n=-p)^oo a_n (z-z_0)^n
  $
  where $b_n equiv a_(n-p)$ and $a_(-p) = b_0 eq.not 0$.
]

We have
$
  b_n & = (g^((n)) (z_0))/n! \
      & = 1/(2 pi i) integral.cont g(z)/(z-z_0)^(n+1) dd(z)
$
implying
$
  a_n & = 1/(2 pi i) integral.cont g(z)/(z-z_0)^(n+1+p) dd(z) \
      & = 1/(2 pi i) integral.cont f(z)/(z-z_0)^(n+1) dd(z)
$
We can have cases where the series around $z_0$ has infinite terms
$
  f(z) = sum_(n=-oo)^oo a_n (z-z_0)^n
$
The principal part converges for $CC \\ C_1$ defined by
$ abs(z-z_0) > R_("pp") $
while the analytic part converges for $C_2$ defined by
$
  abs(z-z_0) < R_"ap"
$
The Laurent series then converges within a punctured disk $D = CC\\C_1 inter C_2$ when $R_"ap" > R_"pp"$. Then any $f$ which is analytic on $D$ can be expressed a Laurent series about $z_0$. When the principal part is finite $C_1$ contains a single point $z_0$ and $C_2$ could extend to infinity.

The Laurent series can be used to _classify_ $z_0$. This is done as:

1. If $f$ is analytic in $z = z_0$ then $a_n = 0$ for $n < 0$.

2. Let $a_n = 0$ for $n < 0$ and $a_0 = a_1 = dots = a_(m-1) = 0$. Then if for $m > 0$ we have $a_m (z-z_0)^m eq.not 0$ then $f$ has a zero of order $m$ at $z_0$.

3. Let $f$ be non-analytic at $z_0$ and assume there exists some $p in NN$ such that $a_(-p) eq.not 0$ and $a_(-p-k) = 0$ for all $k in NN$. Then $f(z) = a_(-p) (z-z_0)^(-p) + dots$ has a pole of order $p$ at $z =z_0$.

4. Let $f$ be non-analytic at $z_0$ and assume there does not exist some $p in NN$ such that $a_(-p) eq.not 0$ and $a_(-p-k) = 0$ for all $k in NN$. Then there are infinitely many negative powers of $(z-z_0)$ and $z_0$ is an essential singularity.

We call $a_(-1)$ the _residue_ of $f$ at the pole $z_0$. We denote this by $Res_(z=z_0) f(z)$ or $Res(z_0)$.

#example[
  Consider
  $
    f(z) = 1/(z(z-2)^3)
  $
  which has singularities at $z = 0$ and $z = 2$. We consider $z = 0$. Then
  $
    f(z) & = - 1/(8 z) 1/(1-z/2)^3 \
    &=^("analytic at" z = 0) -1/(8z) [1 + (-3) (-z/2) + ((-3)(-4))/2! (-z/2)^2 + dots] \
    &= -1/(8z) - 3/16 - 3/16 z - 5/32 z^2 + dots
  $
  We see $z = 0$ is a pole of order $p = 1$ with residue $-1/8$.
]

#theorem[
  When $f$ has a pole of order $m$ then
  $
    Res_(z=z_0) f(z) = lim_(z-> z_0) [1/(m-1)! dv(, z, m-1) ((z-z_0)^m f(z))]
  $
]

#proof[
  Consider the Laurent series around $z_0$
  $
    f(z) = a_(-m)/(z-z_0)^m + dots + a_(-1)/(z-z_0) + a_0 + a_1 (z-z_0) + dots
  $
  Then
  $
    (z-z_0)^m f(z) = a_(-m) + a_(-m+1) (z-z_0) + dots + a_(-1) (z-z_0)^(m-1) + a_0 (z-z_0)^m + dots
  $
  and
  $
    dv(, z, m-1) ((z-z_0)^m f(z)) & = (m-1)! a_(-1) + sum_(n=1)^oo b_n (z-z_0)^n
  $
  yielding the result upon taking $z -> z_0$.
]

When the pole is simple $m = 1$ the above reduces to
$
  Res_(z=z_0) f(z) = lim_(z->z_0) ((z-z_0) f(z))
$
We can simplify this further by assuming
$
  f(z) = g(z)/h(z)
$
Then by _l'Hôpital's rule_
$
  Res_(z=z_0) f(z) = g(z_0)/(h^((1)) (z_0))
$

== The residue theorem
We now consider integrals on a closed contour $Gamma$ which encircles points where $f$ is non-analytic. Consider $f$ has a pole of order $m$ at $z = z_0$. Then
$
  f(z) = sum_(n=-m)^oo a_n (z-z_0)^n
$
We seek to compute
$
  cal(I) = integral.cont_Gamma f(z) dd(z)
$
with $z_0$ inside $Gamma$ and no other singularities of $f$. By Cauchy's integral theorem
$
  cal(I) = integral.cont_gamma f(z) dd(z)
$
with $gamma$ being a circle inside $Gamma$ centered at $z_0$. We parametrise $gamma$ by $z = z_0 + rho e^(i theta)$ giving
$
  cal(I) = sum_(n=-m)^oo a_n integral_0^(2 pi) i rho^(n+1) e^(i (n+1) theta) dd(theta)
$
We see all terms with $n eq.not -1$ vanish so we obtain
$
  cal(I) = 2pi i a_(-1)
$
implying the value of $cal(I)$ is fully determined by the residue of $f$ at $z_0$. The generalisation of the above is the _residue theorem_
#theorem(name: "Residue theorem")[
  Let $f$ be continuous within and on a closed contour $Gamma$ and analytic, except for a finite number of poles within $Gamma$. Then
  $
    integral.cont_Gamma f(z) dd(z) = 2 pi i sum_j Res_(z=z_j) f(z)
  $
]
#proof[
  The proof follows by the above observation. We can encircle each singularity and connect these in such a way that the internal contours cancel.
]

#example(name: "Green's function for the wave equation")[
  Consider
  $
    (-1/c^2 pdv(, t, 2) + laplacian) u(bold(r),t) = f(bold(r),t)
  $
  We demand
  $
    (-1/c^2 pdv(, t, 2) + laplacian) G(bold(r),t;bold(r)',t') = delta(bold(r)-bold(r)') delta(t-t')
  $
  Then
  $
    u(bold(r),t) = integral G(bold(r),t; bold(r)',t') f(bold(r)',t') dd(bold(r)', t, [3,1])
  $
  We will now find $G$. We are interested in $t >= t'$ and require for $t < t'$ then $G = 0$. This is referred to as the _retarded Green's function_. Consider
  $
    delta(bold(r)-bold(r)')delta(t-t') = 1/(2pi)^4 integral dd(k, omega, [3]) e^(i bold(k) dot (bold(r)-bold(r)')) e^(-i omega (t-t'))
  $
  and
  $
    G(bold(r),t ; bold(r)',t') = integral dd(k, 3)/(2pi)^(3\/2) integral dd(omega)/(2 pi)^(1\/2) hat(G)(bold(k),omega) e^(i bold(k) dot (bold(r)-bold(r)')) e^(-i omega (t-t'))
  $
  These reduce the problem to
  $
    hat(G)(bold(k),omega) (omega^2/c^2-bold(k)^2) = 1/(4 pi^2)
  $
  implying
  $
    hat(G) (bold(k),omega) = - c^2/(4 pi^2) 1/(k^2 c^2 - omega^2)
  $
  Then
  $
    G(bold(r),t ; bold(r)',t') = - c^2/(4 pi^2) integral dd(k, 3)/(2 pi)^(3\/2) integral dd(omega)/(2 pi)^(1\/2) 1/(k^2 c^2 - omega^2) e^(i bold(k) dot (bold(r)-bold(r)')) e^(-i omega (t-t'))
  $
  We need to compute
  $
    cal(I) = integral_(-oo)^oo dd(omega) e^(-i omega tau)/(k^2 c^2 - omega^2)
  $
  where $tau equiv t-t'$. The integral $cal(I)$ is over $RR$
  and has poles at $omega = plus.minus k c$ where $k = abs(bold(k))$. We wish to use the residue theorem so we apply the _$i epsilon$-prescription_ shifting the poles to
  $
    plus.minus k c -> plus.minus k c minus i epsilon
  $
  with $epsilon$ small. This allows us to perform the integral along $RR$. We then enclose the poles by a semi-circle contour $Gamma$ with radius $rho > k c$. The corresponding integral is
  $
    cal(I)_epsilon = integral.cont_Gamma dd(omega) e^(-i omega tau)/(k^2 c^2 - (omega + i epsilon)^2)
  $
  We would like $cal(I)_epsilon -> cal(I)$ as $rho -> oo$.#footnote[This is somewhat trivial since $omega -> - i oo$ so we have exponential supression.] This happens if the contribution from the semi-circle vanishes. The semi-circle is parametrised by $omega = rho e^(i phi)$ with $pi <= phi < 2 pi$ so in the $epsilon -> 0$ limit for simplicity we have
  $
    abs(e^(-i omega tau)/(k^2 c^2 - omega^2)) &= abs(e^(-i rho (cos phi + i sin phi) tau)/(k^2 c^2 - rho^2 (cos 2 phi + i sin 2 phi))) \
    &= e^(rho sin phi tau)/abs(k^2 c^2 - rho^2 (cos 2 phi + i sin 2 phi))
  $
  and $abs(omega^2) = rho^2$. Since $sin phi <= 0$ for $pi <= phi < 2 pi$ and demanding $tau >= 0$ we find for $rho >> k c$
  $
    abs(e^(-i omega tau)/(k^2 c^2 - omega^2)) tilde^(rho -> oo) e^(rho sin phi tau)/rho^2
  $
  implying
  $
    abs(integral_gamma dd(omega) e^(-i omega tau)/(k^2 c^2- omega^2)) & <= integral_pi^(2 pi) dd(phi) rho e^(rho sin phi tau)/abs(k^2 c^2 - rho^2) \ &->^(rho -> oo) 0
  $
  Then $cal(I)_epsilon = cal(I)$. Then by the residue theorem#footnote[We have an extra $-$ sign since we close the contour in the lower-half plane.]
  $
    cal(I) = - 2 pi i (Res_(omega = k c) + Res_(omega = - k c))
  $
  We have
  $
    e^(-i omega tau)/(k^2 c^2 - omega^2) = e^(-i omega tau)/((k c - omega)(k c + omega))
  $
  implying
  $
    Res_(omega=k c) &= lim_(omega -> k c) e^(-i omega tau)/((k c -omega)(k c+ omega)) ( omega-k c) \
    &= - e^(-i k c tau)/(2 k c) \
    Res_(omega=-k c) &= lim_(omega -> - k c) e^(-i omega tau)/((k c - omega)(k c + omega)) (omega+ k c)
    \ &= e^(i k c tau)/(2 k c)
  $
  Then
  $
    cal(I) & = - (pi i)/(k c) (e^(i k c tau) - e^(-i k c tau)) \
           & = - (pi i)/(k c) 2 i sin k c tau \
           & = (2 pi sin k c tau)/(k c)
  $
  The last thing required to determine $G$ is then computing the $k$-space integral. This is annoying so we will not do it.
]
