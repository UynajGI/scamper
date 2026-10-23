# Percolation thresholds and fractal dimensions for square/cubic lattices with correlated disorder

> Phys. Rev. E 97, 062141 (2018) · [arXiv](https://arxiv.org/abs/1708.02296)
> Local copy: `docs/ref/Pruessner_2018_correlated_thresholds.pdf` (12 pages)

<!-- page 1 -->

Percolation thresholds and fractal dimensions for square and cubic lattices with
long-range correlated defects
Johannes Zierenberg,1, 2, 3, 4 Niklas Fricke,1, 2 Martin Marenz,1, 2
F. P. Spitzner,1 Viktoria Blavatska,2, 5 and Wolfhard Janke1, 2
1Institut f¨ur Theoretische Physik, Universit¨at Leipzig Postfach 100 920, 04009 Leipzig, Germany
2Doctoral College for the Statistical Physics of Complex Systems, Leipzig-Lorraine-Lviv-Coventry (L4),
Postfach 100 920, 04009 Leipzig, Germany
3Max Planck Institute for Dynamics and Self-Organization, Am Fassberg 17, 37077 G¨ottingen, Germany
4Bernstein Center for Computational Neuroscience, Am Fassberg 17, 37077 G¨ottingen, Germany
5Institute for Condensed Matter Physics of the National Academy of Sciences of Ukraine, 79011 Lviv, Ukraine
We study long-range power-law correlated disorder on square and cubic lattices. In particular, we
present high-precision results for the percolation thresholds and the fractal dimension of the largest
clusters as function of the correlation strength.
The correlations are generated using a discrete
version of the Fourier ﬁltering method. We consider two diﬀerent metrics to set the length scales
over which the correlations decay, showing that the percolation thresholds are highly sensitive to
such system details. By contrast, we verify that the fractal dimension df is a universal quantity
and unaﬀected by the choice of metric. We also show that for weak correlations, its value coincides
with that for the uncorrelated system. In two dimensions we observe a clear increase of the fractal
dimension with increasing correlation strength, approaching df →2. The onset of this change does
not seem to be determined by the extended Harris criterion.
I.
INTRODUCTION
Structural obstacles (impurities) play an important
role for a wide range of physical processes as most sub-
strates and surfaces in nature are rough and inhomo-
geneous [1, 2].
For example, the properties of mag-
netic crystals are often altered by the presence of ex-
tended defects in the form of linear dislocations or re-
gions of diﬀerent phases [3, 4]. Another important class
of such disordered media are porous materials, which of-
ten exhibit large spatial inhomogeneities of a fractal na-
ture. Such fractal disorder aﬀects a medium’s conductiv-
ity, and diﬀusive transport can become anomalous [5–8].
This aspect is relevant, for instance, for the recovery of
oil through porous rocks [9, 10], for the dynamics of ﬂu-
ids in disordered media [11, 12], or for our understanding
of transport processes in biological cells [13, 14].
Disordered systems are conveniently studied in the
framework of lattice models with randomly positioned
defects (or empty sites).
Of particular interest is the
situation where the concentration of occupied (i.e., non-
defect) lattice sites is near the percolation threshold and
clusters of connected occupied sites become fractal. The
case where defects are uncorrelated is a classic text-
book model, whose properties have been studied exten-
sively [15]. In nature, however, inhomogeneities are often
not distributed completely at random but tend to be cor-
related over large distances. To understand the impact
of this, it is useful to consider the limiting case where
correlations asymptotically decay by a power law rather
than exponentially with distance:
C(r) ∼|r|−a.
(1)
An illustration of such power-law correlations for con-
tinuous and discrete site variables on a square lattice is
shown in Fig. 1. If the correlation parameter a is smaller
than the spatial dimension D, the correlations are con-
sidered long-range or “inﬁnite”.
The problem of power-law correlated disorder has ﬁrst
been investigated in the context of spin systems and later
for percolation [16, 17].
The relevance of the disorder
was shown to be characterized by an extension of the
Harris criterion for uncorrelated defects [18]: the crit-
ical behavior of the system deviates from the uncorre-
lated case if the minimum of D and a is smaller than
2/ν (where ν denotes the correlation-length exponent for
the ordered system). It was furthermore argued that in
the regime of strong correlations, the critical correlation-
length exponent for strong disorder is universally given
by 2/a. Since D is always larger than 2/ν for percolation,
the correlation-length exponent for long-range correlated
percolation is given by
-4
-3
-2
-1
0
1
2
3
4
FIG. 1. Illustration of long-range correlated defects on a 20482
lattice with correlation parameter a = 0.5. (left) Continu-
ous correlated Gaussian random variables. The color reﬂects
the value at the respective lattice site. (right) Correspond-
ing lattice of discrete variables at the ﬁnite-size percolation
threshold pL
c = 0.522 with defects shown in black.
arXiv:1708.02296v2  [cond-mat.stat-mech]  2 Jan 2018

<!-- page 2 -->

2
νa =
(
2/a
for a < aH = 2/ν
ν
for a ≥aH = 2/ν.
(2)
The extended Harris criterion is still slightly controver-
sial [19, 20], but it has to some extent been supported by
numerical investigations [21–23]. These studies made use
of the Fourier ﬁltering method (FFM) [21–31] to generate
power-law correlated disorder and have yielded estimates
for critical exponents and fractal dimensions characteriz-
ing the system in 2D. However, they in part used semi-
analytical implementations of the FFM, involving various
approximations and free parameters. In this work we use
a numerical version without free parameters, and whose
errors are fully controlled.
The remainder of the article is organized as follows:
Section II gives a detailed description of the FFM, so that
our implementation is easily reproducible [32].
There-
after, in Sec. III, we specify how the mapping to discrete
site variables is carried out.
In the following Sec. IV,
we present our results for the percolation thresholds on
square and cubic lattices. Our main ﬁndings, regarding
the fractal dimension df for long-range correlated perco-
lation clusters in 2D and 3D, are discussed in Sec. V.
Finally, our results and conclusions are summarized in
Sec. VI.
II.
GENERATING LONG-RANGE
CORRELATED DISORDER
We start with the more general problem of how to ob-
tain a hyper-cubic lattice LD of identically distributed
random variables τx ∈R that exhibit correlations of the
form
⟨τxτx+r⟩= Cr,
(3)
where ⟨. . .⟩denotes the expectation value and Cr is a
(discrete) correlation function. Cr should be symmetric
around zero and periodic along all spatial dimensions,
i.e., Cr+Lei = Cr for all unit vectors e1, . . . , eD. It is
furthermore convenient to choose τx as Gaussian random
variables with mean ⟨τx⟩= ⟨τ⟩= 0 and variance σ2
τ = 1.
Otherwise, we consider Cr to be an arbitrary function for
now. (Note that we use the index notation for explicitly
discrete functions).
We use a variant of the Fourier ﬁltering method that
employs discrete Fourier transforms (DFT) and is similar
to that from Ref. [30, 31]. The key idea of the FFM is to
correlate random variables in Fourier space. The result of
the inverse transform will in general be complex numbers,
ϕx = ϕre
x + iϕim
x . To explain how the method works, let
us now assume that we already have a lattice of complex
random variables ϕx. Let us further assume that {ϕre
x }
and {ϕim
x } are independent sets of random variables, each
spatially correlated according to Eq. (3), i.e.,

ϕre
x ϕre
x+r

=

ϕim
x ϕim
x+r

= Cr,

ϕre
x ϕim
x+r

=

ϕim
x ϕre
x+r

= 0,
(4)
and see what that implies for the distributions of Fourier
coeﬃcients.
As we are interested in a discrete lattice with periodic
boundary conditions of linear size L and volume N = LD,
we consider a DFT of the form
Φk =
X
x
ϕxe
2πikx
L
,
(5)
ϕx = 1
N
X
k
Φke−2πikx
L
,
(6)
where P
x denotes the D-dimensional sum over possible
realizations of the vector x on the hypercubic lattice. In
practice, we employ a numerical fast Fourier transform
(FFT) [33] and follow the convention that xi ∈[0, L) and
ki ∈[0, L).
As shown in Appendix A, the correlation function is
connected to the Fourier coeﬃcients via
2Cr =

ϕ∗
xϕx+r

=
1
N 2
X
k

|Φk|2
e−2πikr
L
.
(7)
The discrete spectral density
Sk =
X
r
Cre
2πikr
L
(8)
can thus be written as
Sk =
X
r
1
2N 2
X
k′

|Φk′|2
e−2πik′r
L
e
2πikr
L
=
1
2N
X
k′

|Φk′|2
δk′,k =
1
2N

|Φk|2
=
1
2N
D
Φre
k
2E
+
D
Φim
k
2E
.
(9)
In return, this means we can generate complex real-
space random variables with the desired correlation from
Fourier-space random variables that satisfy Eq. (9). It
is convenient to consider distributions of Φk with zero
mean, so that Eq. (9) can be expressed in terms of the
variance:
2NSk = σ2
Φk = σ2
Φre
k + σ2
Φim
k .
(10)
Hence, we can simply draw real and imaginary parts of
Φk independently from identical distributions (for each
frequency k):
Φre/im
k
=
p
SkU,
(11)
where U is a random variable with mean ⟨U⟩= 0 and
variance σ2
U = N. Transforming Φk back to x-space, we
get two sets of variables, {ϕre
x } and {ϕim
x }, each with zero
mean and spatial correlations Cr. Thanks to the orthog-
onality of the Fourier transform, the two sets are statisti-
cally independent. Each can be associated with the real
random site variables τx in Eq. (3) and used for further
analysis. We draw U from a Gaussian distribution, and

<!-- page 3 -->

3
so the resulting distributions will also be Gaussian. (In
fact, they would be Gaussian anyway for large systems
due to the central limit theorem.)
The derivation above did not use any assumptions re-
garding the correlation function Cr.
However, we see
from Eq. (9) that its Fourier transform Sk needs to be
positive. Any Cr that is symmetric (around zero) will
give rise to real Sk, but the positivity constraint is some-
what problematic. For the continuum Fourier transform,
it is in fact also implied by the symmetry [28], but for
discrete systems, some values of Sk can become negative.
This has to do with the restricted frequency range, lead-
ing to an aliasing eﬀect that causes periodic modulations
on the signal. Note, however, that this is not just an arti-
fact of the method, but rather implies that some correla-
tions are fundamentally not possible on a ﬁnite discrete
lattice. In practice, we can simply ﬁx this problem by
setting all negative values of Sk to zero (“zero-cutoﬀ”).
While this will inevitably modify the resulting correla-
tions, the eﬀect is usually negligible and vanishes rapidly
with increasing system size, see Appendix B.
In short, our version of the FFM can be summarized
as follows:
1. Choose a discrete correlation function Cr that is
symmetric around zero. For optimal performance,
the linear size of the lattice should be L = 2l with
integer l.
2. Perform a DFT, Cr →Sk, and set Sk = 0 for all
Sk < 0 (zero-cutoﬀ). This step only needs to be
done once for the whole disorder ensemble.
3. Construct real and imaginary parts of each com-
ponent independently, Φre/im
k
= √SkU, where U
is drawn from a Gaussian distribution with mean
⟨U⟩= 0 and variance σ2
U = N.
4. Perform an inverse DFT, Φk →ϕx, to obtain two
independent sets of long-range correlated variables
{ϕre
x } and {ϕim
x }. Each can be associated with a
set of real random variables {τx}.
No free parameter is involved in the process. The only
minor issue is a potential zero-cutoﬀ(only for strong cor-
relations), but the practical impact of this intervention
is small and can be assessed a priori (see Appendix B).
Here we are interested in long-range power-law cor-
related, Gaussian random variables with the following
properties:
⟨τxτx+r⟩∼|r|−a,
for
|r| ≫1,
σ2
τ = C0 = 1.
We follow the suggestion by Makse et al. [28] and consider
the correlation function
C(r) =
 1 + r2−a/2 ,
(12)
which satisﬁes the above conditions. More generally, cor-
relations of the form C(r, α) = (1 + |r|α)−a/α with α > 0
10−3
10−2
10−1
100
100
101
C(r)
|r|
a = 0.50
a = 1.00
a = 1.50
a = 2.00
FIG. 2.
Correlation function C(r) compared to the measured
site-site correlation function ⟨Cr⟩R of continuous variables
along the x-direction on a 162 lattice with R = 106 disorder
replicas. The continuous random variables are obtained via a
discrete Fourier transform of Cr and satisfy σ2
τ = ⟨C0⟩R = 1.
are all suitable and may be chosen depending on the de-
sired behavior of convergence to the asymptotic limit.
To verify the correlations numerically, we measure the
site-site correlation function along the “x-direction” (unit
vector e1) with periodic boundary conditions,
⟨Cr⟩R =
*
1
N
X
x
(τx −⟨τ⟩)
 τx+|r|e1 −⟨τ⟩

+
R
.
(13)
Here ⟨. . .⟩R denotes the disorder average over R replicas,
and the expectation value ⟨τ⟩is zero, which we veriﬁed
numerically.
With increasing sample size R, the mea-
sured correlation function rapidly converges to the envis-
aged Cr. As can be seen in Fig. 2 for a two-dimensional
lattice, the agreement is striking even for very small sys-
tems (162), despite the zero-cutoﬀ. This is one of the ben-
eﬁts of a fully discrete implementation of the FFM over
semi-analytical techniques, which often cannot faithfully
reproduce the desired distributions for small systems.
For a short review of other variants to generate long-
range power-law correlations and a discussion of some of
the diﬃculties, see Appendix C.
III.
MAPPING TO LONG-RANGE
CORRELATED DEFECTS
To study percolation, we have to map the correlated
continuous variables τx to correlated discrete values tx ∈
{0, 1}. For this, we need to specify the mean density of
available sites p (considering defects as tx = 0). Here,
we use a global or grand-canonical [34] approach and ﬁx
the expectation value ⟨P
x tx/N⟩= ⟨t⟩= p. Therefore,
we introduce a threshold θ such that sites are considered
defects if τx < θ.
In the disorder average the τx are

<!-- page 4 -->

4
10−5
10−4
10−3
10−2
10−1
100
100
101
102
103
⟨Cr⟩R / ⟨C0⟩R
|r|
global
local
a = 0.50
a = 1.50
FIG. 3.
Normalized correlation function of discrete ran-
dom variables with long-range power-law correlation at the
percolation threshold pc(a) (see Table II) on a 10242 lattice
averaged over 104 disorder replicas.
Continuous correlated
random variables were obtained as described in Sec. II. The
mapping to discrete variables is performed via a global ap-
proach (open symbols), i.e., on the level of the disorder aver-
age [see Eq. (14)], and via a local approach, i.e., on the level of
each disorder realization (ﬁlled symbols). The measured site-
site correlation function ⟨Cr⟩R along the x-axis is normalized
with respect to the variance ⟨C0⟩R for discrete site variables.
normally distributed, such that the threshold is tied to p
via
p = p(θ) =
Z ∞
θ
P(τ)dτ = 1
2erfc
 
θ
p
2σ2τ
!
,
(14)
where erfc denotes the standard complementary error
function and σ2
τ = 1 by construction.
Note that for
strong correlations, the densities on individual replica
ﬂuctuate signiﬁcantly. If we measure the site-site corre-
lation function of discrete variables according to Eq. (13)
(where we replace ⟨τ⟩with ⟨t⟩= p), we observe ⟨C0⟩R =
σ2
t < 1. The variance of discrete site variables is no longer
unity but is instead connected to the variance of uncor-
related random lattices, σ2
t = p(1 −p). Figure 3 (open
symbols) shows the discrete site-site correlation function
averaged over 104 lattices of size 10242. It can be seen
that the average site-site correlations on discrete lattices
mapped via Eq. (14) decay according to Eq. (12) over a
long range, though the amplitudes are somewhat dimin-
ished.
Alternatively, one might use a local or canonical [34]
approach, adjusting P
x tx/N = p for each replica by
sorting the continuous correlated variables and adjusting
θ until P
x Θ(τx −θ)/N = p, where Θ is the unit step
function [21]. However, ﬁxing p on every lattice tends to
suppress correlations on a macroscopic scale. As can be
seen in Fig. 3 (ﬁlled symbols), this results in a decay rate
of the correlation function that is faster than polynomial.
This eﬀect is most signiﬁcant for strong correlations and
small systems and can be expected to vanish in the limit
of inﬁnite system size. By contrast, the global approach
Eq. (14) described above works reliably for any lattice
size and appears thus generally preferable.
IV.
PERCOLATION THRESHOLD
The value of the percolation threshold pc is not a uni-
versal quantity. It may not only depend on the type of
lattice but also on local aspects of the correlation func-
tion Cr and hence on the implementation of the FFM.
Numerical results given in this section therefore only ap-
ply for the speciﬁc settings we used and cannot be quan-
titatively compared to those from previous studies, e.g.,
Ref. [22]. We were careful to be explicit about these set-
tings to ensure that our results for the fractal dimensions
are reproducible, and so that future studies may use our
estimates for pc.
We use the correlation function Eq. (12) and perform
a discrete numerical Fourier transform as discussed in
the previous section. The radial distance |r| is usually
considered in the Euclidean metric, but here we also use
the Manhattan metric, i.e., the minimum number of steps
on the lattice. This is done to demonstrate the sensitivity
of pc to changes of the correlation function that are not
captured in the correlation parameter a. Later, we also
use the Manhattan metric to test the robustness of our
estimates for the fractal dimensions, which should be the
same for both variants.
To deﬁne percolation on a ﬁnite lattice, we apply the
horizontal wrapping criterion: a cluster percolates if it
closes back on itself across one speciﬁc lattice boundary.
This choice has the beneﬁt of being translationally in-
variant and is known to give relatively small ﬁnite-size
errors [35]. The percolation threshold pL
c for the ﬁnite
system of extension L is then deﬁned as the average oc-
cupation density at which a percolating cluster emerges.
We estimate this value by determining the maximum
threshold θc for each replica of continuous variables at
which a percolating cluster exists for the subset of sites
with τx ≥θc. We then take the average of the mapped
values,
pL
c = ⟨p (θc)⟩,
(15)
where the mapping is carried out according to Eq. (14).
A.
Square lattice
In 2D, we extrapolate to the percolation threshold for
the inﬁnite system, pc := p∞
c , via the standard ﬁnite-
size scaling approach [15] without higher-order correction
terms
|pc −pL
c | ∼L−1/νa.
(16)
Here νa denotes the critical exponent of the correlation
length. The value of νa is determined by Eq. (2) with the

<!-- page 5 -->

5
TABLE I. Estimates of the percolation threshold for square
lattices with correlated disorder. For consistency all ﬁts in-
clude L ≥64.The extended Harris criterion, Eq. (2), modiﬁes
νa̸ = ν for a < 2/ν = 1.5.
a
νa
pc (Euclid.)
χ2/dof pc (Manh.) χ2/dof
∞
4/3 [36] 0.592746 [37]
3
4/3
0.561406(4)
0.92
2.5
4/3
0.556214(4)
0.87
2
4/3
0.550143(5)
0.90
0.528397(5) 1.9
1.75 4/3
0.546717(7)
0.41
1.5
4/3
0.54299(1)
3.5
0.519991(8) 4.0
1.25 8/5
0.53895(2)
1.4
1
2
0.53452(4)
0.87
0.51226(4)
2.28
0.75 8/3
0.5296(1)
0.63
0.5
4
0.5239(3)
0.53
0.5054(3)
0.66
0.25 8
0.516(1)
0.38
0.1
20
0.508(4)
1.2
uncorrelated correlation-length exponent ν = 4/3 [15].
This assumed behavior of νa has been numerically sup-
ported for percolation on a 2D triangular lattice [21].
To obtain our numerical estimates, we randomly gener-
ated 105 replicas for each size L = 2l where l = 6, . . . , 13
(L = 64–8192). Some of the results for pL
c are shown in
Fig. 4 (top) together with least-squares ﬁts of Eq. (16)
over the range L ≥64.
The estimates for pc are the
y-intercepts of the ﬁt curves.
The values are listed in
Table I, where we also give the reduced χ2-values per
degree-of-freedom (dof) of the respective ﬁts. The last
columns show our results for systems where the Manhat-
tan metric is used to set the distance |r| for the correla-
tions. Here the estimates for pc are considerably smaller
than for the Euclidean case, underlining the strong de-
pendence of pc on the details of the correlation function.
In both cases the χ2-values are mostly close to one, justi-
fying the simple scaling ansatz. However, they are quite
large at the “crossover” value of aH = 1.5, where the be-
havior is supposed to change according to the extended
Harris criterion, Eq. (2). This suggests the presence of
additional correction terms in the vicinity of aH, possibly
of logarithmic nature.
An overview of the results for the percolation thresh-
olds as a function of a is shown in the bottom plot of
Fig. 4.
As can be seen, correlations tend to lower pc,
which is intuitive as they promote the emergence of larger
clusters.
As noted in Ref. [22] the value of pc for the
square lattice must eventually approach 1/2. This bound
can be understood considering that a cluster of occu-
pied sites that wraps the system in one direction exists
if and only if no cluster of defects wraps the system in
the orthogonal direction, where the defects are allowed
to connect via next-nearest neighbors (diagonally). For
a −→0, the relevance of these next-nearest neighbor con-
nections becomes negligible, and the resulting symmetry
between clusters of defects and occupied sites demands
0.51
0.52
0.53
0.54
0.55
0.56
0
0.05
0.1
0.15
0.2
0.25
0.3
0.35
0.4
pL
c
L−1/νa
a = 2.00
a = 1.50
a = 1.00
a = 0.50
0.5
0.51
0.52
0.53
0.54
0.55
0.56
0.57
0.58
0.59
0.6
0
0.5
1
1.5
2
2.5
3
uncorrelated
pc
a
Euclidean
Manhattan
FIG. 4. Top: Measured percolation thresholds pL
c for varying
lattice size L and diﬀerent values of the correlation parame-
ter a in 2D (Euclidean metric). Lines show the best ﬁts of
Eq. (16) to the data, whose intercepts represent our estimates
for the inﬁnite-system value pc. Bottom: Estimates for pc as
function of the correlation strength for the square lattice with
distances measured in the Euclidean metric (squares) and the
Manhattan metric (crosses). The horizontal line shows the
value for the system without correlations.
pc = 1/2. Note that for the Manhattan metric, diagonal
correlations are weaker to begin with. The strong devi-
ations do not only depend on the chosen metric but are
already aﬀected by the details of the employed method,
as can be seen by comparing to results we obtained with
the continuous FFM on a square lattice [38], which qual-
itatively look similar but do not agree within error bars.
When a is increased, i.e., when the correlation strength
is diminished, pc must converge towards the value for the
uncorrelated system as long as C(r)/C(0) −→0 for all
|r| > 0. Note, however, that the uncorrelated value is
only reached in the limit a →∞and not at a = D, where
the correlations become eﬀectively short range. This is
contrary to the results from previous studies due to dif-
fering deﬁnitions of the correlation function C(r), which
at a = 2 has a vanishing amplitude in Ref. [22] and a
divergent variance C(0) in Ref. [27].

<!-- page 6 -->

6
TABLE II. Estimates of the percolation threshold for cubic
lattices with correlated disorder. For consistency all ﬁts in-
clude L ≥32. The extended Harris criterion, Eq. (2), modiﬁes
the correlation-length exponent νa̸ = ν for a < 2/ν ≈2.28.
a
νa
pc (Euclid.) χ2/dof pc (Manh.) χ2/dof
∞0.8762 0.311610(2) 0.44
4
0.8762 0.238778(4) 0.10
3
0.8762 0.208438(5) 0.83
0.209315(4) 0.75
2.5 0.8762 0.188289(7) 1.9
0.189801(5) 6.0
2
1
0.16302(2)
3.0
0.16514(1)
0.54
1.5 4/3
0.13022(5)
0.51
0.13251(4)
0.37
1
2
0.0863(3)
1.1
0.0878(2)
0.86
0.5 4
0.025(3)
1.4
0.030(2)
0.25
B.
Cubic lattice
The version of the FFM described in Sect. II can di-
rectly be applied in three (or more) dimensions as well,
which allowed us to study percolation with long-range
correlated disorder on the cubic lattice. We looked at
systems with linear extensions in the range L = 8–512,
and we again generated 105 random replicas for each size.
Unlike in 2D, however, the simple ﬁnite-size scaling ap-
proach to estimate the percolation threshold pc, Eq. (16),
proved unsuccessful, suggesting the need of higher-order
terms (see Ref. [39] for a discussion of ﬁnite-size scaling
for uncorrelated systems):
|pc −pL
c | ∼L−1/νa(A + BL−ω + CL−1/νa + . . .), (17)
where νa is given by Eq. (2) with ν the correlation-length
exponent for uncorrelated percolation (0.8764(12) [40],
0.8762(12) [41], 0.8751(11) [42]).
In practice, the correction to Eq. (16) seems to be de-
scribed well by the latter (quadratic) term alone, suggest-
ing that the correction-to-scaling exponent ω is relatively
large. This is in fact the case for the uncorrelated system,
where previous estimates locate the correction-to-scaling
exponent between ω ≈1.62 [29] and ω ≈1.2 [40]. We
thus used the ansatz
|pc −pL
c | ∼L−1/νa(A + CL−1/νa),
(18)
which we ﬁtted to the data for L ≥32. The correspond-
ing ﬁt curves and our results for pL
c (a) are shown for
selected correlations in Fig. 5 (top), and the resulting es-
timates for pc = p∞
c
are listed in Table II. Again, we see
that the changing behavior predicted by the extended
Harris criterion (at aH = 2/ν ≈2.28) manifests itself
in a poorer quality of the ﬁts for nearby values (a = 2
and a = 2.5). Our estimate for the uncorrelated case,
pc(∞) = 0.311 610(2), is in decent agreement with pre-
vious estimates (0.311 608 1(11) [29], 0.311 607 7(2) [40],
0.311 607 68(15) [41]).
In contrast to the 2D situation, using the Manhattan
metric in place of the Euclidean metric to measure the
0
0.1
0.2
0.3
0.4
0.5
0.6
0.7
0
0.05
0.1
0.15
0.2
0.25
0.3
0.35
0.4
pL
c
L−1/νa
a = ∞
a = 3.00
a = 2.00
a = 1.00
0
0.05
0.1
0.15
0.2
0.25
0.3
0.35
0
0.5
1
1.5
2
2.5
3
3.5
4
uncorrelated
pc
a
Euclidean
Manhattan
FIG. 5. Top: Measured percolation thresholds pL
c for varying
lattice size L and diﬀerent values of the correlation parame-
ter a in 3D (Euclidean metric). Lines show the best ﬁts of
Eq. (18) to the data, whose intercepts represent our estimates
for the inﬁnite-system value pc. Bottom: Estimates for pc as
function of the correlation strength for the cubic lattice with
distances measured in the Euclidean metric (squares) and the
Manhattan metric (crosses). The horizontal line shows the
value for the system without correlations.
distance for the correlation function does not signiﬁcantly
lower the percolation threshold. As can be seen in Ta-
ble II and Fig. 5 (bottom), the values are even slightly
larger.
That is plausible since the argument why the
Manhattan metric should lower pc in 2D does not apply
in 3D, where wrapping clusters of defects and occupied
sites can coexist. This also means that there is no obvi-
ous lower bound for pc in 3D other than zero. Indeed, our
estimates for strong correlations are very small, and the
overview shown in Fig. 5 (bottom) even seems to suggest
the extrapolation pc →0 for a →0.
We should note, however, that the scaling ansatz
Eq. (18) is mainly motivated empirically. Especially for
small a, some of the ﬁnite-size corrections have a diﬀerent
origin as in the uncorrelated system, namely that smaller
systems are not self-averaging: For small a and small L,
the continuous site variables τx within each individual

<!-- page 7 -->

7
TABLE III. Fractal dimension for square lattices with corre-
lated disorder. For consistency all ﬁts include sizes L ≥128.
In two dimensions we consider only the leading-order behav-
ior at the ﬁnite-size percolation transition. Results from Eu-
clidean and Manhattan metric are in good agreement.
a
df (Euclid.) χ2/dof df (Manh.) χ2/dof
∞
91/48 ≈1.89583...
3
1.8961(2)
0.74
2.5
1.8962(2)
1.2
2
1.8966(2)
4.5
1.8964(2)
1.2
1.75 1.8964(2)
2.8
1.5
1.8965(3)
1.6
1.8956(3)
2.3
1.25 1.8950(3)
1.2
1
1.8961(3)
1.2
1.8952(3)
0.29
0.75 1.9006(4)
1.2
0.5
1.9128(5)
0.47
1.9126(6)
0.17
0.25 1.9360(6)
0.085
0.1
1.9602(8)
0.39
replica tend to be very similar, and about half the en-
semble has mostly negative values, while the other has
mostly positive values. In the limit a →0 (at ﬁxed L)
the values τx across each replica become constant, so that
a wrapping cluster in the discrete system emerges when a
threshold θc = τx is used for the mapping. Since the over-
all distribution of the τ-values is symmetric (Gaussian)
and we deﬁne pL
c according to Eq. (14) and Eq. (15), the
a →0 limit at ﬁxed L is pL
c (0) = 1/2. This “segregation”
ﬁnite-size eﬀect might play a signiﬁcant role for the most
strongly correlated cases (a = {0.5, 1}), and our respec-
tive estimates should therefore be taken with a pinch of
salt.
V.
FRACTAL DIMENSION
The fractal dimension df describes how the volume of
a critical percolation cluster increases with its linear size.
It can conveniently be estimated via
⟨V ⟩∼Ldf,
(19)
where L is the lattice extension and ⟨V ⟩denotes the aver-
age number of sites in the largest cluster [15]. Note that
for correlated systems, it is important to include replicas
with no percolating cluster. It is possible to either con-
sider all systems at the same, asymptotic concentration
pc or to take size-dependent values, pL
c . We opted for
the latter approach, so we would not have to rely on the
ﬁtting ansatz for pc.
A.
Square lattice
In 2D, ﬁnite-size corrections again turned out to be
small, so that ﬁtting Eq. (19) without any higher-order
0.2
0.3
0.4
0.5
10
100
1000
10000
⟨V ⟩/N
L
a = 0.50
a = 1.00
a = 1.50
a = 2.00
uncorrelated
1.88
1.9
1.92
1.94
1.96
1.98
2
0
0.5
1
1.5
2
2.5
3
df
a
Euclidean
Manhattan
uncorrelated
FIG. 6. Top: Average volume fraction of the largest cluster vs.
lattice extension in 2D plotted for diﬀerent correlations a on
a double-logarithmic scale (Euclidean metric). Colored lines
are least-squares ﬁts to Eq. (19). The black line represents the
behavior of the uncorrelated system with a slope of 91/48−2.
Bottom: Overview of our estimates for the fractal dimensions
as function of a with the horizontal line again corresponding
to the uncorrelated system.
correction terms worked well. Figure 6 (top) shows the
average volume of the largest cluster relative to the total
number of sites, ⟨V ⟩/LD, for several diﬀerent values of
a plotted on a double-logarithmic scale. The lines cor-
respond to least-squares ﬁts of Eq. (19) over the range
L ≥128, and their slopes show the diﬀerences to the
Euclidean dimension, df −2. Our resulting estimates for
df can be found in Table III together with the reduced
χ2-values of the ﬁts. Also listed are estimates obtained
using the Manhattan instead of the Euclidean metric.
Here, the ﬁts yielded smaller amplitudes, but the expo-
nents resulted very similar. This can be seen in Fig. 6
(bottom), which shows an overview of the estimates for
df. The data veriﬁes that df is universal, i.e., indepen-
dent of system details. For weak correlations the uncor-
related value, df = 91/48 [36], seems to be recovered in
accordance with the extended Harris criterion [Eq. (2)]
and earlier numerical ﬁndings [21–23, 43]. Interestingly

<!-- page 8 -->

8
TABLE IV. Fractal dimension for cubic lattices with corre-
lated disorder. For consistency all ﬁts include sizes L ≥16. In
three dimensions, we require higher-order corrections of the
form Eq. (20). Assuming universality, we perform simultane-
ous ﬁts including both Euclidean and Manhattan metric.
a
df
w
χ2/dof
∞2.52295(15) [40] 1.2(2) [40]
3
2.524(2)
1.2(2)
0.70
2.5 2.522(2)
0.9(2)
3.3
2
2.512(3)
0.9(2)
1.9
1.5 2.507(6)
0.6(1)
0.55
1
2.6(2)
0.20(1)
2.1
though, there seems to be no increase of df directly be-
low aH = 3/2, the crossover threshold set by the ex-
tended Harris criterion. For the Manhattan metric, the
ﬁt quality is still diminished around aH, suggesting that
the threshold may still aﬀect correction terms. However,
it is yet unclear why χ2 is largest at a = 2 for the Eu-
clidean case.
The fact that df does not increase directly below aH
was already noted in Ref. [21], where a crossover thresh-
old of ax = 2/3 (or in terms of the Hurst exponent
Hx = −ax/2 = −1/3) was suggested instead. However,
that value is not quite consistent with our ﬁndings, which
show a signiﬁcant increase already at a > 3/4. Another
disagreement regards the behavior in the correlated limit,
a →0 (H →0): our results are consistent with the idea
that the fractal dimension converges to the “Euclidean”
value of D = 2 as clusters get more and more compact,
while according to Ref. [21] the value stays well below
2. This discrepancy may be owed to the use of diﬀerent
mapping rules as discussed at the end of Sec. III.
It is interesting to compare the results for df with the
Ising model at criticality, which exhibits spin-spin corre-
lations of the form ⟨SiSj⟩∼r−(d−2+η). In two dimen-
sions η = 1/4 and the fractal dimension of the geometri-
cal Ising clusters is df = 187/96 = 1.9479 . . . [44, 45],
which is indeed quite similar to our result of df
=
1.9360(6) for a = 1/4. As already noted [22], the results
could not be expected to agree perfectly. In fact, it is in-
tuitive that df should be slightly larger for Ising clusters,
where the spin-spin correlation function is essentially the
probability that two spins belong to the same cluster. In
our system, by contrast, spins from unconnected clus-
ters still contribute to the correlation function, so that
connected clusters may be “thinner” for the same decay
exponent.
B.
Cubic lattice
The situation in 3D turned out to be more diﬃcult.
As with the percolation threshold, the scaling behavior
seems to involve strong ﬁnite-size corrections, so that
0.02
0.1
0.3
1
10
100
1000
⟨V ⟩/N
L
a = 2, Euclidean
a = 2, Manhattan
uncorrelated
2.5
2.6
2.7
2.8
2.9
3
0
0.5
1
1.5
2
2.5
3
df
a
simultaneous ﬁt
uncorrelated
FIG. 7. Top: Average volume fraction of the largest cluster
vs.
lattice extension in 3D plotted for the exemplary case
a = 2 on a double-logarithmic scale. Colored lines are the
simultaneous ﬁts according to Eq. (20); the black line rep-
resents the behavior of the uncorrelated system with a slope
of 2.52295 −3. Bottom: Overview of our estimates for the
fractal dimensions as function of a with the horizontal line
again corresponding to the uncorrelated system. Our ﬁtting
did not work properly for very strong correlations (a ≤1).
simply ﬁtting Eq. (19) would not work for the system
sizes that we considered.
Including a correction term
also failed as the ﬁt could not handle two additional pa-
rameters. What did work reasonably well, at least for
a > 1, was ﬁtting our data for the Euclidean and the
Manhattan versions simultaneously, while assuming the
exponents of the leading term and the correction to be
equal for both cases:
⟨V ⟩Euclid. = A1Ldf(1 + B1L−w),
⟨V ⟩Manh. = A2Ldf(1 + B2L−w).
(20)
This approach was motivated by general universality ar-
guments [46, 47] and our previous observation that the
fractal dimensions in 2D are the same for both versions.
We assume that equality also holds for the correction
exponents ω, which seems reasonable since ω is also
strongly believed to be universal for percolation without

<!-- page 9 -->

9
correlations, see for instance Ref. [48]. Figure 7 shows
ﬁnite-size scaling data and ﬁts for the case a = 2 as an
example (top) and an overview of the obtained estimates
for df (bottom). The values of our estimates can be found
in Table IV together with the correction exponents and
the χ2-values of the ﬁts.
Unfortunately, the data for
a ≤1 could not be convincingly ﬁtted by this approach.
For these strongly correlated cases, one would probably
need to investigate systems still much larger than 5123.
For a ≥1.5, the value for df seems to be very similar to
the one without correlations (df(∞) = 2.522 95(15) [40]).
As in 2D, the Harris threshold, aH = 2/ν ≈2.28, does
hence not determine the onset of a sudden increase in
the fractal dimension. Surprisingly, the value even seems
to decrease slightly below aH. At close inspection, this
can also be observed in 2D for a = 1.25 < aH, compare
Table III. A diminishing fractal dimension does not seem
plausible as stronger correlations should make the clus-
ters more compact. We suspect that a correction term
comes into play at aH which is not captured by our ﬁtting
approach.
VI.
CONCLUSIONS
We presented high-precision results for the percolation
thresholds on square and cubic lattices with long-range
power-law correlated disorder as well as estimates for
the fractal dimensions of the critical percolation clusters.
The correlations were generated using the Fourier ﬁlter-
ing method (FFM) based on the discrete Fourier trans-
form. We speciﬁed the details of our implementation, so
that it may easily be reproduced [32] and discussed the
diﬀerences to previous approaches regarding, e.g., how
the continuous site variables are mapped to discrete dis-
order.
The percolation threshold is dependent on the em-
ployed method and moreover on short-range details of the
model. We demonstrated this by using both the standard
Euclidean metric and the discrete Manhattan metric to
deﬁne the correlation. The eﬀect of this choice on the
percolation threshold is particularly strong for the square
lattice. This is because diagonal correlations are weaker
for the Manhattan metric, bringing the system closer to
pc = 0.5 where the percolation thresholds for occupied
sites and defects connected via next-nearest (diagonal)
neighbors coincide. In general, correlations were shown
to lower pc, and in three dimensions the value even be-
comes very close to zero for small a, i.e., strong correla-
tions.
The fractal dimension, by contrast, is a universal quan-
tity and does not depend on details of the model. We
veriﬁed that for large a (weak correlations) the fractal di-
mension of the uncorrelated model is reproduced showing
df ≈91/48 (2D) and df ≈2.52 (3D). This was expected
above the bound from the extended Harris criterion, i.e.,
for all a ≥aH with aH = 1.5 (2D) and aH ≈2.28 (3D).
However, as was previously noticed for the triangular lat-
tice [21], df ≈duncorr
f
seems to remain true also well below
the Harris bound. In two dimensions, our data suggests
that the value of df starts to rise below ax ≈1, approach-
ing df →2 as a →0. Diﬀerences to previous ﬁndings
may be attributed to diﬀerent mapping prescriptions em-
ployed. To obtain estimates for df in three dimensions,
we simultaneously ﬁtted our data for correlations with
Euclidean and Manhattan metric using a polynomial ﬁt
with a correction term. Unfortunately, this approach did
not work for very strong correlations, i.e., for a ≤1. In
the accessible range, the values were found very close to
the uncorrelated value. Below aH ≈2.28, they even re-
sulted slightly smaller, which we suspect is due to chang-
ing corrections to scaling. We conclude that while the
bound from the Harris criterion does not seem to deter-
mine a change in the leading exponent df, it does aﬀect
the system’s sub-leading behavior.
ACKNOWLEDGMENTS
We thank Martin Weigel and Martin Treﬀkorn for
helpful discussions. This work has been supported by an
Institute Partnership Grant Leipzig-Lviv of the Alexan-
der von Humboldt Foundation (AvH). Further ﬁnan-
cial support from the Deutsche Forschungsgemeinschaft
(DFG) via the Sonderforschungsbereich SFB/TRR 102
(project B04), the Leipzig Graduate School of Natural
Sciences “BuildMoNa”, as well as from the Deutsch-
Franz¨osische Hochschule (DFH-UFA) through the Doc-
toral College “L4” under Grant No. CDFA-02-07 and the
EU through the Marie Curie IRSES network DIONICOS
under Contract No. PIRSES-GA- 2013-612707 (FP7-
PEOPLE-2013-IRSES) is gratefully acknowledged. J. Z.
received ﬁnancial support from the German Ministry of
Education and Research (BMBF) via the Bernstein Cen-
ter for Computational Neuroscience (BCCN) G¨ottingen
under Grant No. 01GQ1005B.
Appendix A: Discrete Wiener-Khinchin theorem
We require ϕx to be complex random variables with in-
dependent real and imaginary contributions. For a given
disorder realization the lattice average of ϕ∗
xϕx+r can be

<!-- page 10 -->

10
written as
1
N
X
x
ϕ∗
xϕx+r
= 1
N
X
x
 
1
N
X
k
Φ∗
ke
2πikx
L
! 
1
N
X
l
Φle−2πil(x+r)
L
!
=
1
N 2
X
k
Φ∗
k
X
l
Φle−2πilr
L
1
N
X
x
e
2πi(k−l)x
L
=
1
N 2
X
k
Φ∗
k
X
l
Φle−2πilr
L δl,k
=
1
N 2
X
k
|Φk|2e−2πikr
L
.
(A1)
Here,
we
used
the
notation
of
a
D-
dimensional
Kronecker-Delta
function
δl,k = Πiδli,ki =
1
LD
P
x e2πi(k−l)x/L.
The
result
is
essentially the discrete Wiener-Khinchin theorem,
a
special case of the cross-correlation theorem.
Taking the disorder average on both sides of Eq. (A1)
and exploiting translational invariance on the left, we
thus obtain

ϕ∗
xϕx+r

=
1
N 2
X
k

|Φk|2
e−2πikr
L
.
(A2)
Appendix B: Eﬀect of zero-cutoﬀin Sk on Cr
As mentioned in Sec. II, particular choices of C(r) eval-
uated on a ﬁnite lattice may lead to unphysical negative
values of the discrete spectral density Sk. This seems to
occur only for strong correlations (small a) and becomes
more noticeable with increasing dimension. Numerically,
we deal with this issue by a zero-cutoﬀ, i.e., by using a
modiﬁed spectral density
eSk =
(
Sk,
if Sk ≥0
0,
else.
(B1)
This inevitably aﬀects the resulting correlation function.
We can directly predict the eﬀect from the inverse dis-
crete Fourier transform of eSk since
eCr = 1
N
X
k
eSke−2πikr
L
(B2)
is of course the asymptotic limit of the measured site-site
correlation function ⟨Cr⟩R for large sample size R.
It turns out that the eﬀect of this zero-cutoﬀis very
small, if present at all, with deviations mainly occur-
ring for small lattices. With increasing lattice size, the
predicted (and measured) deviations quickly converge to-
wards the desired correlation function. To demonstrate
this, we consider the example of three-dimensional lat-
tices with (strong) correlations a = 0.5 and linear exten-
sions L = {8, 16, 32} in Fig. 8. The measured correlation
0.3
0.4
0.5
0.6
0.7
0.8
0.9
1
1.1
0
1
2
3
4
5
6
7
8
9
3D
⟨Cr⟩R
|r|
target C(r) = (1 + |r|2)−0.5/2
predictions ˜
Cr
L =
8
L = 16
L = 32
10−4
10−3
10−2
10−1
100
0
1
2
3
4
5
6
7
8
9
3D
⟨Cr⟩R −C(r)
|r|
predictions ˜
Cr −C(r)
L =
8
L = 16
L = 32
FIG. 8. (Top) Demonstration of the eﬀect of the zero-cutoﬀ
Sk →eSk on the measured correlation function ⟨Cr⟩R in three
dimensions for a = 0.5. Also shown is the prediction eCr for
the asymptotic limit from Eq. (B2).
(Bottom) Plot of the
deviations from the desired function C(r).
function of continuous site variables along the x-direction
⟨Cr⟩R [Eq. (13)], is evaluated with data from Sec. IV. The
eﬀect of the zero-cutoﬀon the correlation function from
Eq. (12) is indeed perfectly predicted by Eq. (B2).
Appendix C: Diﬀerent versions of the FFM
Many diﬀerent variants of the FFM can be found in the
literature [22, 25–31]. We want to give a brief overview
of the diﬀerences and discuss the eﬀects of some of the
implied approximations.
In early works the spectral density is approximated as
S(q) = |q|−(D−a) [22]. The resulting non-trivial ampli-
tude in the correlation function C(r) = f(D −a)|r|−a
was shown to vanish for a →D, in accordance with the
picture that the uncorrelated case should be recovered
for short-range correlations (a > D). Still, the desired
correlation function could only be produced in a small
region of the system with this approach.

<!-- page 11 -->

11
Reference [27] follows a similar idea, but directly uses
C(r) = |r|−a.
This function diverges at |r| = 0, and
hence the authors interpolate by integrating the function
in the corresponding discrete bin around zero. This works
reliably in one dimension but becomes cumbersome in
more dimensions.
Moreover, this assigns a non-trivial
value to C(0) thus modifying the variance of the desired
Gaussian random variables.
The most inﬂuential works are by Makse et al. [26, 28].
They introduced the correlation function from Eq. (12),
allowing them to approach the problem both numeri-
cally [26] and (partially) analytically [28]. Their numer-
ical approach is quite similar to ours, but the analytical
one has received far more attention.
We have in fact
tried it ourselves [38], but found that it has many pit-
falls, which we want to brieﬂy discuss here. The idea is
to discretize the Fourier transform of C(r) = (1+r2)−a/2
for the inﬁnite continuum, which can be calculated ana-
lytically:
S(q) =
Z ∞
−∞
C(r)eiqrdr
= 2πD/2
Γ(a/2)
|q|
2
β
Kβ(|q|),
(C1)
where Γ is Euler’s gamma function and Kβ(|q|) is the
modiﬁed Bessel function of order β = (a−D)/2 [49]. The
variance is recovered by integrating over full continuous
space σ2
τ = C(0) =
1
2π
R ∞
−∞S(q)dq = 1.
The next step is to identify q =
2π
L k and map the
continuous result to a discrete lattice by evaluating the
function S( 2π
L k) at each lattice site k. The ﬁrst problem
here is that S(0) diverges.
This can be circumvented
by evaluating the zero-signal at a shifted frequency, i.e.,
S0 = S( 2π
L |k|0) with |k|0 ∈(0, 1) chosen “appropri-
ately” [28]. With increasing system size the choice be-
comes less relevant, and the diﬀerences can be expected
to vanish in the inﬁnite-system limit. For ﬁnite systems,
however, the eﬀect of the parameter |k|0 depends on the
dimension and the strength of the correlations. In ad-
dition, |k|0 has to be adjusted iteratively, rendering the
application of the method rather tedious.
There is another, more severe problem with the dis-
cretization, which is relevant for the mapping to discrete
site variables [38] (see Sec. III). As we are interested in
the asymptotic long-range scaling behavior, we typically
use a ﬁxed lattice spacing of unit length ∆xi = 1 and
consider the limit to inﬁnite system size rather than to
the continuum.
Thus, the frequency space is conﬁned
to qi ∈[−π, π), while the resolution increases with in-
creasing system size. As a consequence, the variance σ2
τ
will deviate from one, complicating the mapping pro-
cedure.
In fact, we can estimate the deviations via
σ2
τ =
1
2π
R π
−π S(q)dq ≤1 = C(0). We numerically veri-
ﬁed this but also found additional ﬁnite-size scaling cor-
rections of the form σ2
τ,L = σ2
τ + O(L−1).
These cor-
rections are inconvenient for ﬁnite-size scaling, e.g., for
ﬁnding the percolation threshold, because one needs to
evaluate the variance σ2
τ,L in addition to the parameter
|k|0 for each value of the correlation parameter a and
each system size. By contrast, the method sketched in
Sec. II is parameter free and always yields the correct
variance σ2
τ up to negligible eﬀects from the zero-cutoﬀ.
[1] P. Pfeifer and D. Avnir, “Chemistry in noninteger dimen-
sions between two and three. I. Fractal theory of hetero-
geneous surfaces,” J. Chem. Phys. 79, 3558 (1983).
[2] D. Avnir, D. Farin,
and P. Pfeifer, “Molecular fractal
surfaces,” Nature 308, 261 (1984).
[3] S. N. Dorogovtsev, “Critical exponents of magnets with
lengthy defects,” Phys. Lett. A 76, 169 (1980).
[4] Y. Yamazaki, A. Holz, M. Ochiai,
and Y. Fukuda,
“Static and dynamic critical behavior of extended-defect
N-component systems in cubic anisotropic crystals,”
Phys. Rev. B 33, 3460 (1986).
[5] J.-P. Bouchaud and A. Georges, “Anomalous diﬀusion
in disordered media: Statistical mechanisms, models and
physical applications,” Phys. Rep. 195, 127 (1990).
[6] K. Malek and M.-O. Coppens, “Eﬀects of surface rough-
ness on self- and transport diﬀusion in porous media
in the Knudsen regime,” Phys. Rev. Lett. 87, 125505
(2001).
[7] M. E. Foulaadvand and M. Sadrara, “Dynamics of a rigid
rod in a disordered medium with long-range spatial cor-
relation,” Phys. Rev. E 91, 012122 (2015).
[8] I. Goychuk, V. O. Kharchenko,
and R. Metzler, “Per-
sistent Sinai type diﬀusion in Gaussian random poten-
tials with decaying spatial correlations,” Phys. Rev. E
96, 052134 (2017).
[9] F. A. L. Dullien, Porous Media: Fluid Transport and
Pore Structure (Academic Press, New York, 1979).
[10] M. Sahimi, Flow and Transport in Porous Media and
Fractured Rock (VCH, Weinheim, 1995).
[11] T. O. E. Skinner, S. K. Schnyder, D. G. A. L. Aarts,
J. Horbach, and R. P. A. Dullens, “Localization dynam-
ics of ﬂuids in random conﬁnement,” Phys. Rev. Lett.
111, 128301 (2013).
[12] M. Spanner, F. H¨oﬂing, S. C. Kapfer, K. R. Mecke, G. E.
Schr¨oder-Turk,
and T. Franosch, “Splitting of the uni-
versality class of anomalous transport in crowded media,”
Phys. Rev. Lett. 116, 060601 (2016).
[13] A. Bancaud, C. Lavelle, S. Huet, and J. Ellenberg, “A
fractal model for nuclear organization: Current evidence
and biological implications,” Nucleic Acids Res. 40, 8783
(2012).
[14] F. H¨oﬂing and T. Franosch, “Anomalous transport in the
crowded world of biological cells,” Rep. Prog. Phys. 76,
046602 (2013).
[15] D. Stauﬀer and A. Aharony, Introduction to Percolation
Theory (Taylor and Francis, London, 1992).
[16] A. Weinrib and B. I. Halperin, “Critical phenomena in
systems with long-range-correlated quenched disorder,”
Phys. Rev. B 27, 413 (1983).
[17] A. Weinrib, “Long-range correlated percolation,” Phys.

<!-- page 12 -->

12
Rev. B 29, 387 (1984).
[18] A. B. Harris, “Eﬀect of random defects on the critical
behaviour of Ising models,” J. Phys. C 7, 1671 (1974).
[19] V.
V.
Prudnikov
and
A.
A.
Fedorenko,
“Critical
behaviour of 3D systems with long-range-correlated
quenched disorder,” J. Phys. A: Math. Gen. 32, L399
(1999).
[20] V. V. Prudnikov, P. V. Prudnikov, and A. A. Fedorenko,
“Field-theory approach to critical behavior of systems
with long-range-correlated defects,” Phys. Rev. B 62,
8777 (2000).
[21] K. J. Schrenk, N. Pos´e, J. J. Kranz, L. V. M. van
Kessenich, N. A. M. Ara´ujo, and H. J. Herrmann, “Per-
colation with long-range correlated disorder,” Phys. Rev.
E 88, 052102 (2013).
[22] S. Prakash, S. Havlin, M. Schwartz,
and H. E. Stan-
ley, “Structural and dynamical properties of long-range
correlated percolation,” Phys. Rev. A 46, R1724 (1992).
[23] H. A. Makse, J. S. Andrade Jr., M. Batty, S. Havlin,
and H. E. Stanley, “Modeling urban growth patterns with
correlated percolation,” Phys. Rev. E 58, 7054 (1998).
[24] D. Saupe, “Algorithms for random fractals,” in The Sci-
ence of Fractal Images, edited by H.-O. Peitgen and
D. Saupe (Springer, New York, 1988) pp. 71–136.
[25] C.-K. Peng, S. Havlin, M. Schwartz, and H. E. Stanley,
“Directed-polymer and ballistic-deposition growth with
correlated noise,” Phys. Rev. A 44, R2239 (1991).
[26] H. A. Makse, S. Havlin, H. E. Stanley, and M. Schwartz,
“Novel method for generating long-range correlations,”
Chaos, Solitons & Fractals 6, 295 (1995).
[27] N.-N. Pang, Y.-K. Yu, and T. Halpin-Healy, “Interfacial
kinetic roughening with correlated noise,” Phys. Rev. E
52, 3224 (1995).
[28] H. A. Makse, S. Havlin, M. Schwartz, and H. E. Stanley,
“Method for generating long-range correlations for large
systems,” Phys. Rev. E 53, 5445 (1996).
[29] H. G. Ballesteros and G. Parisi, “Site-diluted three-
dimensional Ising model with long-range correlated dis-
order,” Phys. Rev. B 60, 12912 (1999).
[30] B. Ahrens and A. K. Hartmann, “Critical behavior of
the random-ﬁeld Ising magnet with long-range correlated
disorder,” Phys. Rev. B 84, 144202 (2011).
[31] M. S. Simon, J. M. Sancho,
and A. M. Lacasta, “On
generating random potentials,” Fluct. Noise Lett. 11,
1250026 (2012).
[32] Our
code
(C++)
is
available
at
github.com/
CQT-Leipzig/correlated_disorder.
[33] W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P.
Flannery, Numerical Recipes 3rd edition: The Art of Sci-
entiﬁc Computing (Cambridge University Press, Cam-
bridge, 2007).
[34] S. Wiseman and E. Domany, “Self-averaging, distribu-
tion of pseudocritical temperatures, and ﬁnite size scal-
ing in critical disordered systems,” Phys. Rev. E 58, 2938
(1998).
[35] M. E. J. Newman and R. M. Ziﬀ, “Fast Monte Carlo
algorithm for site or bond percolation,” Phys. Rev. E
64, 016706 (2001).
[36] B. Nienhuis, “Critical behavior of two-dimensional spin
models and charge asymmetry in the Coulomb gas,” J.
Stat. Phys. 34, 731 (1984).
[37] R. M. Ziﬀ, “Spanning probability in 2D percolation,”
Phys. Rev. Lett. 69, 2670 (1992).
[38] N. Fricke, J. Zierenberg, M. Marenz, F. P. Spitzner,
V. Blavatska, and W. Janke, “Scaling laws for random
walks in long-range correlated disordered media,” Con-
dens. Matter Phys. 20, 13004 (2017).
[39] H. J. Herrmann and D. Stauﬀer, “Corrections to scaling
and ﬁnite size eﬀects,” Phys. Lett. A 100, 366 (1983).
[40] J. Wang, Z. Zhou, W. Zhang, T. M. Garoni,
and
Y. Deng, “Bond and site percolation in three dimen-
sions,” Phys. Rev. E 87, 052107 (2013).
[41] X Xu, H Wang, J.-P. Lv, and Y. Deng, “Simultaneous
analysis of three-dimensional percolation models,” Front.
Phys. 9, 113 (2014).
[42] H. Hu, H. W. J. Bl¨ote, R. M. Ziﬀ, and Y. Deng, “Short-
range correlations in percolation at criticality,” Phys.
Rev. E 90, 042106 (2014).
[43] H. A. Makse, S. Havlin, and H. E. Stanley, “Modelling
urban growth patterns,” Nature 377, 608 (1995).
[44] A. L. Stella and C. Vanderzande, “Scaling and fractal
dimension of Ising clusters at the d = 2 critical point,”
Phys. Rev. Lett. 62, 1067 (1989).
[45] B. Duplantier and H. Saleur, “Exact fractal dimension of
2D Ising clusters,” Phys. Rev. Lett. 63, 2536 (1989).
[46] J. Cardy, Scaling and Renormalization in Statistical
Physics (Cambridge University Press, Cambridge, 1996).
[47] J. Zinn-Justin, Phase Transitions and Renormalization
Group (Oxford University Press, New York, 2007).
[48] H. G. Ballesteros, L. A. Fern`andez, V. Mart´ın-Mayor,
A. Mu˜noz Sudupe, G. Parisi,
and J. J. Ruiz-Lorenzo,
“Scaling corrections: site percolation and Ising model in
three dimensions,” J. Phys. A: Math. Gen. 32, 1 (1999).
[49] Note that there is a typo in the argument of the Euler
gamma function in Ref. [28].
