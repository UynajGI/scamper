# Percolation with long-range correlated disorder

> Phys. Rev. E 88, 052102 (2013) · [arXiv](https://arxiv.org/abs/1309.2994)
> Local copy: `docs/ref/Schrenk_2013_correlated_disorder_percolation.pdf` (13 pages)

<!-- page 1 -->

Percolation with long-range correlated disorder
K. J. Schrenk,1, ∗N. Pos´e,1, † J. J. Kranz,1, ‡ L. V. M. van Kessenich,1, § N. A. M. Ara´ujo,1, ¶ and H. J. Herrmann1, 2, ∗∗
1Computational Physics for Engineering Materials, IfB, ETH Zurich,
Wolfgang-Pauli-Strasse 27, CH-8093 Zurich, Switzerland
2Departamento de F´ısica, Universidade Federal do Cear´a, 60451-970 Fortaleza, Cear´a, Brazil
Long-range power-law correlated percolation is investigated using Monte Carlo simulations. We
obtain several static and dynamic critical exponents as function of the Hurst exponent H which
characterizes the degree of spatial correlation among the occupation of sites.
In particular, we
study the fractal dimension of the largest cluster and the scaling behavior of the second moment
of the cluster size distribution, as well as the complete and accessible perimeters of the largest
cluster. Concerning the inner structure and transport properties of the largest cluster, we analyze
its shortest path, backbone, red sites, and conductivity. Finally, bridge site growth is also considered.
We propose expressions for the functional dependence of the critical exponents on H.
PACS numbers: 64.60.ah, 64.60.al, 89.75.Da
I.
INTRODUCTION
In percolation on a lattice, each lattice element (site or
bond) is occupied with probability p or empty with prob-
ability 1 −p. Occupied sites are connected to their near-
est neighbors and form clusters, the properties of which
depend on p [1, 2]. There is a threshold value pc, such
that for p > pc there exists a cluster spanning between
two opposite sides of the lattice. At p = pc, a continuous
transition occurs between this connected state and the
state for p < pc, where there is no spanning cluster. The
spanning cluster is only fractal at p = pc.
Percolation theory and related models have been ap-
plied to study transport and geometrical properties of
disordered systems [3, 4]. Frequently the disorder in the
system under study exhibits power-law long-ranged spa-
tial correlations. This fact has motivated some studies
of percolation models where the sites of the lattice are
not occupied independently, but instead with long-range
spatial correlation, in a process named correlated per-
colation [3–16].
The qualitative picture that emerged
from those works is that, in the presence of long-range
correlations, percolation clusters become more compact
and their transport properties change accordingly. These
ﬁndings have also been conﬁrmed by experimental stud-
ies of the transport properties of clusters in correlated
invasion percolation [17, 18].
The critical exponents of the uncorrelated percola-
tion transition in two dimensions are known rigorously
for the triangular lattice [19]. In addition, at the criti-
cal point, the correlation-length diverges and universal-
ity holds, i.e., critical exponents and amplitude ratios
do not depend on short-range details, such as lattice
∗jschrenk@ethz.ch
† posen@ifb.baug.ethz.ch
‡ jkranz@ethz.ch
§ laurensm@ethz.ch
¶ nuno@ethz.ch
∗∗hans@ifb.baug.ethz.ch
speciﬁcs [1–4, 20]. This statement has been made pre-
cise by renormalization group theory, which predicts that
the scaling functions within a universality class are the
same, while the lattice structure only inﬂuences the non-
universal metric factors [21, 22]. If, by contrast, inﬁnite-
range, power-law correlations are present, according to
the extended Harris criterion, the critical exponents can
change, depending on how the correlations decay with
spatial distance [5, 7, 16, 23, 24].
Here, we investigate a two-dimensional percolation
model where the sites of a lattice are occupied based on
power-law correlated disorder generated with the Fourier
ﬁltering method [6, 25–32]. The Hurst exponent H of the
disorder is related to the exponent of the power-law de-
cay of spatial correlations with the distance; we ﬁnd that
the fractal dimension of the largest cluster, its perimeter,
as well as the dimension of its shortest path, backbone,
and red sites depend on H [33]. Strong dependence on
H is also found for the electrical conductivity exponent
of the largest cluster and the growth of bridge sites in
the correlated percolation model. For two-dimensional
critical phenomena, conformal ﬁeld theory has been used
to obtain exact values of critical exponents in the form
of simple rational numbers [34–36]. Therefore, we make
proposals for the functional dependence of all measured
exponents on the Hurst exponent H, as being the sim-
plest rational expressions that ﬁt the numerical data.
This work is organized as follows. Section II deﬁnes
the method of generating long-range correlations and the
corresponding correlated percolation model. In Sec. III,
we consider the percolation threshold of the used lattice.
This result is applied in Sec. IV to measure the fractal
dimension of the largest cluster and the scaling behav-
ior of the second moment of the cluster size distribution
at the percolation threshold. The complete and acces-
sible perimeters of the largest cluster are investigated in
Sec. V. Section VI discusses shortest path, backbone, and
red sites of the largest cluster at the threshold. The con-
ductivity of the largest cluster is analyzed in Sec. VII. In
Sec. VIII, we discuss the growth exponent of bridge sites
in the correlated percolation model. Finally, in Sec. IX,
arXiv:1309.2994v2  [cond-mat.stat-mech]  4 Nov 2013

<!-- page 2 -->

2
we present some concluding remarks.
II.
CORRELATED PERCOLATION
To study correlated percolation on a lattice, it is con-
venient to work with a landscape of random heights h,
where h(x) is the height of the landscape at the lat-
tice site at position x [3, 5–7, 37, 38]. Recently, ranked
surfaces have been introduced, providing the adequate
framework to tackle this problem [39]. The ranked sur-
face of a discrete landscape is constructed as follows: One
ﬁrst ranks all sites in the landscapes according to their
height, from the smallest to the largest value. Then, a
ranked surface is constructed where each site has a num-
ber corresponding to its position in the rank. The fol-
lowing percolation model can then be deﬁned: Initially,
all sites of the ranked surface are unoccupied. The sites
are occupied one-by-one, following the ranking. At each
step, the fraction of occupied sites p increases by the in-
verse of the total number of sites in the surface. By this
procedure, a conﬁguration of occupied sites is obtained,
the properties of which depend on the landscape. For
example, if the heights are distributed uniformly at ran-
dom, classical percolation with fraction of occupied sites
p is obtained [40–42].
Here, we study the case where the heights h have
long-range spatial correlations. Such a power-law corre-
lated disorder can be generated using the Fourier ﬁltering
method (Ffm) [6, 16, 25–32, 43], which is based on the
Wiener-Khintchine theorem (WKt) [25, 44]. The WKt
states that the auto-correlation of a time series equals
the Fourier transform of its power spectrum, i.e. of the
absolute squares of the Fourier coeﬃcients. This fact is
exploited in the Ffm by imposing the following power-law
form of the power spectrum S(f) of the disorder:
S(f) ∼|f|−βc =
q
f 2
1 + f 2
2
−βc
,
(1)
where
βc
deﬁnes
the
Hurst
exponent
H
via
βc = 2(H + 1).
By the WKt, this gives the follow-
ing correlation function c(r) of the heights h,
c(r) = ⟨h(x)h(x + r)⟩x ∼|r|2H,
(2)
where the power-law decay of the spatial correlation is
described by the Hurst exponent H. For correlated per-
colation, one considers the range −1 ≤H ≤0 [3, 5–7].
H = −1 corresponds to βc = 0, such that the power spec-
trum in Eq. (1) is independent on the frequency, and the
landscape proﬁle is white noise. This limit recovers un-
correlated percolation. Since H ≤0, as H increases to-
wards zero, the correlation function decays more slowly.
In simulations, for a desired value of H one can generate
random Fourier coeﬃcients of the heights h with ampli-
tudes according to the power spectrum in Eq. (1) and
then apply an inverse fast Fourier transform to obtain
h(x) [6, 25–32, 43].
FIG. 1.
(Color online) Triangular lattice stripe of size L = 4
and aspect ratio A = 2.
The extended Harris criterion,
as formulated in
Refs. [5, 7, 16, 23, 24],
states that for the range
−d/2 < H < 0 the correlations do not aﬀect the crit-
ical exponents of the percolation transition if H
≤
−1/νuncorr, where νuncorr is the correlation-length critical
exponent and for d = 2, νuncorr
2D
= 4/3 [1, 19]. Whereas
for −1/νuncorr
2D
< H < 0 the critical exponents are ex-
pected to depend on the value of H. The quantitative de-
pendence of the critical exponents on H, in this regime, is
not yet entirely clear. Concerning the correlation-length
critical exponent for the correlated case νH, the analyt-
ical works in Refs. [5, 7, 23] predict that νH = −1/H.
In the case of Weinrib and Halperin [5, 23] this is a
conjecture based on renormalization group calculations;
Schmittbuhl et al. [7] found the same result by analyz-
ing hierarchical networks. Therefore, in both analytical
approaches, it is not certain that νH actually behaves
as conjectured and there is some controversy regarding
this question, as discussed, e.g., in the ﬁeld-theoretical
work of Prudnikov et al. [45, 46]. For correlated perco-
lation, the relation νH = −1/H has been supported by
the numerical work in Refs. [16, 47, 48]. Agreement has
also been reported by Prakash et al. [6], however only
approximately for the range −1/νuncorr
2D
≤H ≤−0.5. Fi-
nally, for H > 0 there is no percolation transition [7, 49].
In the following, we consider values of the Hurst exponent
in the range −1 ≤H ≤0.
III.
PERCOLATION THRESHOLD
We consider the correlated percolation model deﬁned
in Sec. II on triangular lattice stripes of length L and as-
pect ratio A, consisting of N = AL2 sites (see Fig. 1). To
investigate critical correlated percolation, one ﬁrst needs
to determine the percolation threshold pc of this lattice.
For site percolation on the triangular lattice, it is possible
to show that pc = 1/2 [1]. The argument of Sykes and
Essam [50, 51] is as follows: For certain lattices, one can
ﬁnd their corresponding matching lattice. In the context
of Refs. [50, 51], this is related to matching expansions of
the mean number of clusters for high and low p. A more
visual explanation of the concept of matching lattice is

<!-- page 3 -->

3
the following [3]: Suppose that for a lattice G1 there ex-
ists a diﬀerent lattice G2, such that each site in lattice G1
is uniquely related to one site in G2 and the other way
around. Also, assume that if a site is occupied in one
of the lattices, its partner in the other one can not be
occupied. Now, if the presence of a cluster spanning G2
in one direction prevents any cluster spanning G1 in the
perpendicular direction and, conversely, there can only
be a percolating cluster in G1 if there is no percolation
in G2, then G1 and G2 are matching lattices. For ex-
ample, the triangular lattice is its own matching lattice,
called self-matching, while the square lattice is matched
by the star lattice [51]. Sykes and Essam argued, based
on the uniqueness of the threshold pc [50–52], that for
any lattice G1 and its matching one G2, the sum of the
thresholds of both equals unity:
pG1
c
+ pG2
c
= 1.
(3)
Then, since the triangular lattice is self-matching, one
has pG1
c
= pG2
c
and it follows that pc = 1/2. The question
of which pairs of lattices match each other is independent
on the statistical properties of the heights h that deter-
mine the cluster properties. Therefore, the site percola-
tion threshold of the triangular lattice is pc = 1/2, also
for correlated percolation. We also checked this state-
ment numerically by measuring pc for diﬀerent values of
the Hurst exponent H, ﬁnding that it is compatible with
1/2, within error bars [53].
As a ﬁrst check of the theory presented in Refs. [5, 7,
23] regarding the dependence of νH on H, we consider
here the convergence of a threshold estimator, namely
the value pc,J at which the maximum change in the size
of the largest cluster smax occurs [54–61]. The expected
scaling behavior [56, 62] is
|pc,J(L) −pc| ∼L−1/νH,
(4)
where pc = 1/2. Figure 2 shows |pc,J(L) −pc| as function
of the lattice size L for diﬀerent values of H.
Within
error bars, the data is compatible with 1/νH = −H for
the considered values of H.
IV.
MAXIMUM CLUSTER SIZE AND SECOND
MOMENT
At the threshold, p = pc, the largest cluster is a fractal
of fractal dimension df, i.e., its size smax scales with the
lattice size L as
smax ∼Ldf .
(5)
This is also related to the order parameter P∞of the
percolation transition, which is deﬁned as the fraction of
sites in the largest cluster,
P∞= smax/N,
(6)
and is expected to scale at p = pc as
P∞∼L−β/ν = Ldf −d,
(7)
10-4
10-3
10-2
10-1
101
102
103
104
pc,J−pc
Lattice size L
0.75±0.02
0.68±0.05
0.41±0.05
0.13±0.05
-1
-0.7
-0.4
-0.1
FIG. 2.
(Color online) Convergence of the percolation
threshold estimator pc,J. The diﬀerence between the estima-
tor and the threshold |pc,J −1/2| is shown as function of the
lattice size L for H = −1, −0.7, −0.4, and −0.1. The data
is shifted vertically to improve visibility. Results are averages
over 105 samples. We keep track of the cluster properties with
the labeling method proposed by Newman and Ziﬀ[40, 41],
as in Ref. [63].
where β is the order parameter critical exponent and
d = 2 is the spatial dimension [1].
For uncorrelated
percolation, β = 5/36 and ν = νuncorr
2D
= 4/3, such that
df = 91/48 ≈1.8958 [1]. To measure df as function of H,
we considered the scaling of the size of the largest clus-
ter smax with the lattice size [see Fig. 3(a) and Eq. (5)].
For diﬀerent values of H, we measured smax(L) and
calculated the local slopes df(L) of the data (see e.g.
Ref. [64]),
df(L) = log[smax(2L)/smax(L/2)]/ log(4).
(8)
Finally, df(L) is extrapolated to the thermodynamic
limit, L →∞, to obtain df(H), see Fig. 4(a). The frac-
tal dimension is, within error bars, independent on H,
for H ≲−1/3. For H approaching zero, the value of df
does increase. While this behavior is in agreement with
Ref. [6], it is in strong contrast to the behavior of all
other fractal dimensions considered in this work, whose
values depend strongly on H. Based on the data, we pro-
pose the following dependence of df on H (in the range
−1/3 ≤H ≤0) as being the simplest rational expression
that ﬁts the numerical data:
df(H) = 91
48 + 13
80
1
3 + H

.
(9)
The hyperscaling,
d = γ
ν + 2β
ν = γ
ν + 2(d −df),
(10)
relates the fractal dimension df to the susceptibility crit-
ical exponent γ and the correlation-length critical expo-
nent ν [1]. For uncorrelated percolation, γ = 43/18 and

<!-- page 4 -->

4
(a)
0.2
0.5
101
102
103
104
Largest cluster fraction smax/N
Lattice size L
-1
-0.85
-0.7
-0.55
-0.4
-0.25
-0.1
0
(b)
100
101
102
103
104
105
106
101
102
103
104
Second moment M2’
Lattice size L
-1
-0.85
-0.7
-0.55
-0.4
-0.25
-0.1
0
FIG. 3. (Color online) (a) Fraction of sites in the largest clus-
ter smax/N as function of the lattice size L for diﬀerent values
of H. (b) Second moment of the cluster size distribution M ′
2
as function of L for the same values of H as in (a). The data
is shifted vertically to improve visibility.
Solid black lines
are guides to the eye. Results are averages over at least 104
samples.
therefore γ/νuncorr
2D
= 43/24 ≈1.7917 [1]. To test the va-
lidity of Eq. (10) for diﬀerent values of H, we measure
γH/νH, where γH and νH are the susceptibility and cor-
relation length critical exponents for a certain H, by con-
sidering the scaling behavior of the second moment M ′
2,
deﬁned as,
M ′
2 = M2 −s2
max/N,
(11)
where
M2 =
X
k
s2
k/N,
(12)
and the sum goes over all clusters with sk being the num-
ber of sites in cluster k. At p = pc, the following scaling
with the lattice size L is expected [1, 65]:
M ′
2 ∼LγH/νH.
(13)
(a)
1.75
1.80
1.85
1.90
1.95
2.00
-1.0
-0.8
-0.6
-0.4
-0.2
0.0
Hurst exponent H
df
γH/νH
(b)
1.94
1.96
1.98
2.00
2.02
2.04
2.06
2.08
2.10
-1.0
-0.8
-0.6
-0.4
-0.2
0.0
γH/νH+2βH/νH
Hurst exponent H
Data
Hyperscaling
FIG. 4.
(Color online) (a) Fractal dimension of the largest
cluster df and critical exponent ratio γH/νH as function of
the Hurst exponent H, where γH is the susceptibility critical
exponent and νH is the correlation-length critical exponent.
For H > −1/3, the solid lines show the expressions of Eqs. (8)
and (14). (b) With df = d −βH/νH, where βH is the order
parameter critical exponent, the hyperscaling relation reads:
2 = d = γH/νH + 2βH/νH [1].
One observes that the data
agrees, within error bars, with the hyperscaling relation.
In Fig. 3(b), one sees M ′
2 as function of L for diﬀerent
values of H. Figure 4(a) shows γH/νH, while γH/νH +
2βH/νH is plotted in Fig. 4(b) for diﬀerent values of H.
One observes that the hyperscaling relation, Eq. (10),
is fulﬁlled, within error bars. Based on this result, we
propose for the functional dependence of γH/νH on the
Hurst exponent H, in the range −1/3 ≤H ≤0, as being
the simplest rational expression that ﬁts the numerical
data:
γH
νH
= (76 + 13H)/40.
(14)

<!-- page 5 -->

5
FIG. 5.
(Color online) Complete and accessible perime-
ter. The blue (ﬁlled) sites of the triangular lattice are part
of the largest cluster, while the white (empty) sites are unoc-
cupied. Bonds of the dual lattice are shown as dashed lines.
Assume that the largest cluster percolates in the vertical di-
rection and does not touch the left or right boundaries of the
lattice. Consider a walker starting on the left-bottom side of
the lattice, which never visits a bond twice and traces out
the complete perimeter, turning left or right depending on
which of the two available bonds separates an occupied from
an empty site. The complete perimeter is fully determined
when the top side of the lattice is reached. Performing the
same walk, but with the additional constraint that fjords with
diameter ≤
√
3/3 (in lattice units) are not accessible, yields
the accessible perimeter. The solid green (thick) lines on the
honeycomb lattice form the accessible perimeter, while dashed
green (thick) lines indicate bonds that are part of the com-
plete perimeter but not of the accessible one. A similar walk
yields the two perimeters on the right hand side of the cluster.
V.
CLUSTER PERIMETERS
Here, we consider triangular lattice stripes of aspect
ratio A = 8 (Fig. 1). For every largest cluster that spans
the lattice vertically (between the long sides of the lat-
tice, Fig. 1) and does not touch its vertical boundaries,
there are two contours that can be deﬁned: the complete
and accessible perimeters [37, 66–76]. Figure 5 shows the
deﬁnition of the two perimeters, which live on the hon-
eycomb lattice, in the case of the triangular lattice: The
complete perimeter consists of all bonds of the honey-
comb lattice that separate sites belonging to the span-
ning cluster from unoccupied sites that can be reached
from the vertical boundaries of the lattice without cross-
ing sites belonging to the largest cluster. If, in addition,
fjords of the perimeter with diameter less than
√
3/3 (lat-
tice units) are inaccessible, the accessible perimeter is ob-
tained. Figure 6 shows the left hand side complete and
accessible perimeters of a percolating cluster on a lattice
of size L = 128. In the upper inset of Fig. 7, the length
of the complete perimeter Mcp is observed to scale with
the lattice size L as,
Mcp ∼Ldcp,
(15)
where for the uncorrelated case, given by H = −1, it is
known that dcp = 7/4 [66, 67, 71, 72].
In addition to
considering the scaling of Mcp with L, we also determined
the fractal dimension dcp using the yardstick method [82,
83]. There, one measures the number of sticks S(m) of
size m needed to follow the perimeter from one end to
the other.
Figure 8 shows that, for intermediate stick
lengths, S(m) scales as
S ∼m−dcp.
(16)
We measured the value of the fractal dimension with this
method for diﬀerent lattice sizes L (see Fig. 8) and then
extrapolated the results to L →∞to obtain dcp. The
fractal dimension dcp(H) determined by this method is
compatible with the one obtained from the scaling of the
length of the perimeter [see Eq. (15)] and we combined
both measurements for the ﬁnal estimates: In Fig. 7,
one sees the fractal dimension of the complete perimeter
as function of the H. For H approaching zero, dcp de-
creases and ﬁnally converges towards 3/2, in agreement
with previous results [37, 75, 76].
The fractal dimension of the accessible perimeter dap
is deﬁned by the scaling of the length of the accessible
perimeter Map with L (see lower inset of Fig. 7),
Map ∼Ldap.
(17)
For uncorrelated percolation the fractal dimension of the
accessible perimeter is known to be dap = 4/3 [71, 72, 77,
78]. Figure 7 shows dap(H), determined using the scaling
of Map and the yardstick method.
For the critical Q-state Potts model [84], Duplantier
[80, 85] established the following duality relation between
the fractal dimension of the complete perimeter dcp and
of the accessible perimeter dap:
(dap −1)(dcp −1) = 1/4.
(18)
The case Q = 1 corresponds to uncorrelated percolation
[86]. Having measured dcp and dap as functions of H, we
see in Fig. 9 that the duality relation of Eq. (18) holds,
within error bars, for −1 ≤H ≤0.
Therefore, taking
the known results for H = −1 and H = 0 into account,
we propose the following functional dependence of the
complete perimeter fractal dimension on H (in the range
−1/νuncorr
2D
≤H ≤0, see Ref. [76]):
dcp = 3
2 −H
3 ,
(19)
which, assuming the validity of the duality relation also
for correlated percolation, implies the following form of
the accessible perimeter fractal dimension:
dap = 9 −4H
6 −4H .
(20)

<!-- page 6 -->

6
(a) H = −1
(b) H = −0.5
(c) H = −0.25
(d) H = 0
FIG. 6.
(Color online) Snapshots of typical complete and accessible perimeters. The accessible perimeter is shown in bold
solid blue lines. In addition, the parts of the complete perimeter that do not belong to the accessible perimeter are drawn
with thin black lines. The snapshots are taken for (a) H = −1, (b) −0.5, (c) −0.25, and (d) 0, on a lattice of (vertical) length
L = 128.
VI.
SHORTEST PATH, BACKBONE, AND RED
SITES
For uncorrelated percolation, the shortest path be-
tween two sites in the largest cluster is a fractal of di-
mension dsp ≈1.131 [87–90]. For a given conﬁguration,
it can be identiﬁed using the burning method [87]: On
the cluster spanning the lattice vertically [with aspect ra-
tio A = 1 (see Fig. 1)], we select one cluster site in the top
row and one in the bottom row, such that their Euclidean
distance is minimized and ﬁnd the number of sites Msp
in the shortest path between them. The following scaling
of the length with the lattice size L is observed:
Msp ∼Ldsp,
(21)
which can be used to determine the fractal dimension
dsp(H) using the local slopes [see Eq. (8)], shown in
Fig. 10. These results are also compatible with the ones
obtained using the yardstick method (not shown). For in-
creasing correlation, dsp deceases and is compatible with
unity for H = 0, as also reported in Ref. [91].
Using
this observation and the literature results for uncorre-
lated percolation [87–90], we propose the following de-
pendence of dsp on the Hurst exponent H (in the range
−1/νuncorr
2D
≤H ≤0):
dsp(H) = 147
130 −
3/4 + H
195/34 + H .
(22)
In addition to measuring the length of the shortest

<!-- page 7 -->

7
1.3
1.4
1.5
1.6
1.7
1.8
-1.0
-0.8
-0.6
-0.4
-0.2
0.0
Perimeter fractal dimension
Hurst exponent H
Complete
Accessible
Eqs. (19), (20)
101
102
101
102
103
Accessible/L
Lattice size L
101
102
101
102
103
Complete/L
Lattice size L
-1
-0.85
-0.7
-0.55
-0.4
-0.25
-0.10
FIG. 7.
(Color online) Main plot: Fractal dimension of the
complete perimeter dcp and of the accessible perimeter dap
as function of the Hurst exponent H. For H = −1 (uncorre-
lated), our results, dcp = 1.75 ± 0.02 and dap = 1.34 ± 0.02,
are in agreement with values previously reported [66, 67,
70, 77–80].
With increasing H, both fractal dimensions
seem to approach 3/2, compatible with the data of Kalda
et al.
[37, 75, 76, 81].
In the range −1/νuncorr
2D
≤H ≤0,
the solid lines show the expressions dcp = 3/2 −H/3 and
dap = (9 −4H)/(6 −4H). Insets: Length of the complete and
of the accessible perimeter as function of the lattice size L for
the values of H shown in the main plot.
100
101
102
103
104
105
100
101
102
103
Number of sticks S
Stick length m
H=0
−1.49±0.03
16
32
64
128
256
512
1024
FIG. 8.
(Color online) Yardstick method to measure the
fractal dimension of the complete perimeter.
The number
of sticks needed to follow the perimeter S is shown as func-
tion of the stick length m, for diﬀerent lattice sizes L, and
H = 0. The numerical value of the complete perimeter frac-
tal dimension dcp(H) obtained with the yardstick method,
dcp(0) = 1.49 ± 0.03, agrees, within error bars, with the re-
sults of the analysis of the local slopes of the perimeter length
(see Fig. 7), as well as with the literature [37, 75, 76].
0.20
0.25
0.30
0.35
-1.0
-0.8
-0.6
-0.4
-0.2
0.0
(dap-1)(dcp-1)
Hurst exponent H
Data
Duality
FIG. 9.
(Color online) Left hand side of the duality relation
for cluster perimeters, (dap −1)(dcp −1) = 1/4 [80, 85], as
function of the Hurst exponent H.
0.95
1.00
1.05
1.10
1.15
1.20
1.25
-1.0
-0.8
-0.6
-0.4
-0.2
0.0
Shortest path dsp
Hurst exponent H
Data
Eq. (22)
101
102
103
104
101
102
103
Path length
Lattice size L
-1
-0.85
-0.7
-0.55
-0.4
-0.25
-0.10
FIG. 10.
(Color online) Fractal dimension of the short-
est path dsp of the largest cluster as function of the Hurst
exponent H.
The inset shows the number of sites in the
shortest path as function of the lattice size L for the same
value of H as in the main plot.
For uncorrelated disor-
der, i.e.
H = −1, we ﬁnd dsp = 1.130 ± 0.005, in agree-
ment with the literature [87–90].
With increasing Hurst
exponent, dsp approaches unity [91].
This behavior is
due to the backbone becoming increasingly compact as H
approaches 0, see Fig. 11.
The solid line is the graph
of the proposed behavior of the shortest path fractal di-
mension:
dsp(H) = 147/130 −(3/4 + H)/(195/34 + H), for
−3/4 ≤H ≤0,
and
dsp(−1
≤
H
≤
−1/νuncorr
2D
)
=
dsp(−1/νuncorr
2D
).
path between two sites in the largest cluster, one can
also ask which sites would carry non-zero current if the
occupied sites would be resistors and a potential diﬀer-
ence were applied between these two sites. This subset
of sites of the largest cluster is called backbone and it is
the union of all non-self-crossing paths between these two

<!-- page 8 -->

8
1.6
1.7
1.8
1.9
2.0
2.1
-1.0
-0.8
-0.6
-0.4
-0.2
0.0
Backbone dbb
Hurst exponent H
Data
Eq. (24)
101
102
103
104
105
106
107
101
102
103
Backbone size
Lattice size L
-1
-0.85
-0.7
-0.55
-0.4
-0.25
-0.1
0
FIG. 11.
(Color online) Fractal dimension of the back-
bone dbb as function of the Hurst exponent H.
With in-
creasing H, the backbone becomes more compact and, con-
sequently, dbb increases, while the fractal dimension of the
shortest path (see Fig. 10) decreases [6].
For uncorrelated
disorder, H = −1, we measure dbb = 1.64 ± 0.02, compatible
with the results reported in Refs. [64, 87, 89, 92, 93, 95].
The solid line is the graph of the following interpolation:
dbb(H) = 39/20(1 + H) −166/101H. Inset: Backbone size as
function of the lattice size L for the same values of H as in
the main plot.
sites [1, 64, 87, 89, 92–95]. Some sites of the backbone are
singly connected, i.e., the connectivity between the two
ends of the backbone is broken if any one of these sites
is removed. These sites are called red sites [39, 96, 97].
Algorithmically, for a given cluster, the backbone and its
red sites can be found with the burning method [87]. The
total number of sites in the backbone Mbb scales with the
lattice size L,
Mbb ∼Ldbb,
(23)
where dbb is the backbone fractal dimension, see inset of
Fig. 11. With increasing H, dbb increases and is compati-
ble with the fractal dimension of the largest cluster for H
approaching zero. Similarly to Ref. [6], for the functional
dependence of dbb on H, we propose to interpolate lin-
early between the best known value for uncorrelated per-
colation, dbb(−1) = 1.6434 ± 0.0002 [95] and the fractal
dimension of the largest cluster for H = 0 [see Eq. (9)]:
dbb(H) = 39
20(1 + H) −166
101H.
(24)
The backbone becomes more compact with increasing
correlation, which is also compatible with the fact that
the shortest path fractal dimension is decreasing in this
-0.2
0.0
0.2
0.4
0.6
0.8
-1.0
-0.8
-0.6
-0.4
-0.2
0.0
Red sites drs
Hurst exponent H
Data
Theory
101
102
103
104
101
102
103
Red sites
Lattice size L
-1
-0.85
-0.7
-0.55
-0.4
-0.25
-0.1
0
FIG. 12.
(Color online) Fractal dimension of the red
sites drs as function of the Hurst exponent H.
In the
main plot, based on the result by Coniglio [39, 96, 97],
the data (squares) is compared to the theoretical pre-
diction for 1/νH
as function of H,
where νH
is the
correlation-length critical exponent of two-dimensional perco-
lation: for H < −1/νuncorr
2D
, 1/νH = 1/νuncorr
2D
= 3/4 and for
−1/νuncorr
2D
≤H < 0, 1/νH = −H [5, 7, 23].
We note that
these results are similar to measurements in Refs. [6, 7, 16,
48, 98]. Inset: Number of red sites as function of the lattice
size L for the values of H shown in the main plot.
limit (see Fig. 10). For the same reason, one would expect
the fractal dimension of the set of red sites drs to decrease
with increasing H. Coniglio [96] has shown that the red
site fractal dimension is related to the correlation-length
critical exponent νuncorr
2D
by drs = 1/νuncorr
2D
. To test the
theoretical predictions in Refs. [5, 7, 23] for 1/νH, we
measured the red site fractal dimension drs as function
of H, see Fig. 12. Although for H approaching zero the
ﬁnite size eﬀects become more severe, see inset of Fig. 12,
the relation seems to be compatible with the data, in
agreement with the results in Refs. [6, 7, 16, 47, 48]. This
is consistent with the ﬁnite-size scaling in the percolation
threshold estimation (see Sec. III).
VII.
CLUSTER CONDUCTIVITY
At the percolation threshold, the backbone of the
largest cluster is a fractal and the conductivity C be-
tween its ends has a power-law dependence on the Eu-
clidean distance r of the end sites,
C(r) ∼r−tH/νH,
(25)
where tH is the conductivity exponent and we call
tH/νH the reduced conductivity exponent [6, 64, 99–
105].
For uncorrelated percolation, tuncorr
2D
/νuncorr
2D
=
0.9826 ± 0.0008 [64].
As the backbone becomes more
compact with increasing correlation (see Sec. VI), one
might expect the conductivity to decay more slowly with

<!-- page 9 -->

9
0.3
0.5
0.7
0.9
-1.0
-0.8
-0.6
-0.4
-0.2
0.0
Reduced exponent tH/νH
Hurst exponent H
Data
Eq. (27)
10-3
10-2
10-1
100
101
101
102
103
Conductivity C
Lattice size L
-1
-0.85
-0.7
-0.55
-0.4
-0.25
-0.1
0
FIG. 13.
(Color online) Reduced conductivity exponent
tH/νH as function of the Hurst exponent H.
For increas-
ing value of H, as the backbone becomes more compact (see
Fig. 11), tH/νH decreases. For uncorrelated disorder, we ﬁnd
tH/νH(−1) = −0.992 ± 0.027 in agreement with Ref. [64].
The solid line corresponds to the expression tH/νH = 16/41−
H −7H2/25 in the range −1/νuncorr
2D
≤H ≤0 and tH/νH =
t/νuncorr
2D
for −1 ≤H ≤−1/νuncorr
2D
. Inset: Conductivity C as
function of the lattice size L, for the same values of the Hurst
exponent H as in the main plot.
the spatial separation, and, consequently, that tH/νH de-
creases [6, 106].
To measure the conductivity C of the backbone, we
solved Kirchhoﬀ’s laws and obtained for every site i in
the backbone:
X
k
(Vi −Vk) = 0,
(26)
where the sum runs over the nearest neighbors k belong-
ing to the backbone of site i, and the conductivity is unity
between neighboring sites. The boundary conditions are
chosen such that V = N on the top end of the backbone
and V = 0 on its bottom end. Solving the sparse linear
system of equations one obtains the conductivity and the
value of the potential at each site of the backbone (for
details, see e.g. Ref. [105]). The inset of Fig. 13 shows
the conductivity C as function of the lattice size L, for
diﬀerent values of H.
Since in our setup the distance
between the end points r ∼L, we use this scaling to de-
termine the reduced conductivity exponent tH/νH, see
Fig. 13. Our result for uncorrelated percolation agrees
with the literature and one observes tH/νH to decrease
with increasing H. We propose the following functional
dependence of the reduced conductivity exponent on H
(in the range −1/νuncorr
2D
≤H ≤0):
tH
νH
= 16
41 −H −7H2
25 .
(27)
 0
 4000
 8000
 12000
 16000
 20000
0.0
0.2
0.4
0.6
0.8
1.0
Number of bridges Mbr
Control parameter p
L=4096
-1
-0.85
-0.7
-0.55
-0.4
-0.25
-0.1
0
FIG. 14.
(Color online) Number of bridge sites Mbr as
function of the control parameter p, for diﬀerent values of the
Hurst exponent H, on a lattice of size L = 4096. Results are
averages over 104 samples.
VIII.
BRIDGE SITE GROWTH
To explore further the impact of correlations on the
structure of percolation clusters, we analyze the bridge
sites, which are related to red sites, at the percolation
threshold [39, 96, 107]. Consider the following modiﬁ-
cation of the percolation model: While the sites are se-
quentially occupied, starting from the empty lattice, if a
site would lead to the emergence of a spanning cluster
between the top and bottom sides of the lattice, this site
does not become occupied and is labeled as a bridge site
[39, 62, 108]. While the fraction of occupied sites p is
lower than the percolation threshold pc, the set of bridge
sites is empty, since there would be no percolating cluster
in classical percolation for p < pc [1]. At the threshold,
the number of bridge sites Mbr behaves identically to the
number of red sites and diverges with the lattice size as
Mbr ∼L1/ν,
(28)
where ν is the correlation-length critical exponent of per-
colation [39, 96, 97]. For uncorrelated disorder, at p > pc,
the number of bridge sites grows as a power law with the
distance from the threshold,
Mbr ∼(p −pc)ζ,
(29)
where ζ = 0.50 ± 0.03 [39] is called the bridge growth ex-
ponent (see also Fig. 14). When p goes to unity, the set
of bridge sites merge to a singly connected line, spanning
the lattice horizontally, which is the watershed of the
landscape of considered heights h, if the top and bottom
sides of the lattice would be connected to water outlets
[109–114]. For uncorrelated landscapes, this watershed is
a fractal path of dimension dbr = 1.2168 ± 0.0005 [107].
To determine how the bridge site growth depends on
H, we measured the number of bridge sites Mbr as func-
tion of p, for diﬀerent values of H, see Fig. 14. For values

<!-- page 10 -->

10
0.0
0.2
0.4
0.6
0.8
-0.25
0.00
0.25
0.50
Bridges Mbr / Ldbr
Distance p−pc
H=−0.85
1024
2048
4096
100
101
102
100
101
102
MbrL−1/ν
(p−pc)Lθ
ζ
FIG. 15.
(Color online) Main plot: Rescaled number of
bridge sites Mbr/Ldbr as function of the distance to the per-
colation threshold p −pc, with H = −0.85, for diﬀerent lattice
sizes L. Here, we use dbr(−0.85) = 1.211 [30, 43, 107]. Inset:
Rescaled number of bridge sites MbrL−1/νuncorr
2D
as function of
the scaling variable (p −pc)Lθ, with θ = 0.72. The solid line
is a guide to the eye with slope 0.64.
of H ≤−1/νuncorr
2D
, we observe that the data for diﬀer-
ent lattice sizes collapses, when rescaled by Ldbr(H), for
all values of p > pc (see Fig. 15). This suggests that the
same crossover scaling as in the uncorrelated case [39]
can be applied to extract the growth exponent ζ:
Mbr(p, L) = L1/νuncorr
2D
F[(p −pc)Lθ],
(30)
where the scaling function F[x] ∼xζ for large x and the
power-law behavior of Mbr in the lattice size L and p
yields
θ = (dbr −1/νuncorr
2D
)/ζ.
(31)
For
H = −0.85,
the
rescaled
data
is
shown
in
the inset of Fig. 15 and the growth exponent is
ζ(−0.85) = 0.64 ± 0.06, larger than for H = −1.
The
corresponding value of θ yielding the best collapse of the
data is θ = 0.72 ± 0.08, in agreement with the scaling re-
lation of Eq. (31), given the known dependence of the
watershed fractal dimension dbr on H [30, 43, 107].
For H ≥−1/νuncorr
2D
, the behavior of bridge sites is
qualitatively diﬀerent from the uncorrelated case. The
rescaled number of bridge sites Mbr(p)/Ldbr(H) does not
overlap for diﬀerent lattice sizes L for any value of p > pc,
except when the complete fractal line has emerged, i.e.
for p →1. An example of this behavior, for H = −0.1,
is shown in Fig. 16. To analyze this size eﬀect in more
detail, we plot in Fig. 17 the number of bridges Mbr as
function of the lattice size L, for diﬀerent values of p ≥pc.
One observes that, in contrast to the uncorrelated case
[39], for p > pc, there is no crossover to the fractal di-
mension of the continuous bridge line dbr. Precisely at
the critical point, the expected behavior Mbr ∼L1/νH is
still observed.
0.0
0.3
0.6
0.9
-0.25
0.00
0.25
0.50
Bridges Mbr / Ldbr
Distance p−pc
H=−0.1
4096
2048
1024
512
256
128
100
101
102
103
10-3
10-1
101
Mbr / L0.3
(p−pc)L0.3
FIG. 16.
(Color online) Main plot: Rescaled number of
bridge sites Mbr/Ldbr as function of p −pc, with H = −0.1,
for diﬀerent lattice sizes L. Inset: Data for the largest three
L, with Mbr/L0.3 as function of (p −pc)L0.3.
100
101
102
103
104
101
102
103
104
Number of bridges Mbr
Lattice size L
0.4±0.3
0.6±0.3
0.8±0.2
0.9±0.2
1.1±0.2
1.18±0.07
0.5
0.6
0.7
0.8
0.9
1.0
FIG. 17.
(Color online) Number of bridge sites Mbr, for
H = −0.1, as function of the lattice size, for diﬀerent values
of the fraction of occupied sites p = pc = 0.5, 0.6, 0.7, 0.8,
0.9, and unity. The solid lines are guides to the eye. The
estimated slopes are indicated on the right had side of the
ﬁgure.
IX.
FINAL REMARKS
Concluding, we studied percolation with long-range
correlation in the site occupation probabilities, as char-
acterized by the Hurst exponent H. The site percolation
threshold of the triangular lattice was argued to be 1/2,
independent of H. For H approaching zero the fractal
dimension of the largest cluster, as well as the exponent
ratio γH/νH where found to increase in accordance with
the hyperscaling relation. The fractal dimensions of the
complete and the accessible perimeter were observed to

<!-- page 11 -->

11
approach 3/2 for H →0, while the duality relation be-
tween both exponents seems to hold independently on the
value of H. As H increased, the backbone of the largest
cluster was observed to become more compact, consistent
with the scaling behavior of shortest path, red sites, and
conductivity. Finally, we found the bridge growth expo-
nent to increase with increasing H. While the qualitative
picture is consistent with previous studies in the litera-
ture, we proposed quantitative relations for the depen-
dence of the critical exponents of the percolation tran-
sition on H, as being the simplest rational expressions
that ﬁt the numerical data.
ACKNOWLEDGMENTS
We acknowledge ﬁnancial support from the ETH Risk
Center, the Brazilian institute INCT-SC, and grant num-
ber FP7-319968 of the European Research Council.
[1] D. Stauﬀer and A. Aharony, Introduction to Percolation
Theory, 2nd ed. (Taylor and Francis, London, 1994).
[2] M. Sahimi, Applications of Percolation Theory (Taylor
and Francis, London, 1994).
[3] M. B. Isichenko, Rev. Mod. Phys. 64, 961 (1992).
[4] M. Sahimi, Rev. Mod. Phys. 65, 1393 (1993).
[5] A. Weinrib, Phys. Rev. B 29, 387 (1984).
[6] S. Prakash, S. Havlin, M. Schwartz, and H. E. Stanley,
Phys. Rev. A 46, R1724 (1992).
[7] J. Schmittbuhl, J.-P. Vilotte, and S. Roux, J. Phys. A
26, 6115 (1993).
[8] M. Sahimi, J. Phys. I France 4, 1263 (1994).
[9] C. Du, C. Satik, and Y. C. Yortsos, AIChE Journal 42,
2392 (1996).
[10] H. A. Makse, S. Havlin, and H. E. Stanley, Nature 377,
608 (1995).
[11] H. A. Makse, J. S. Andrade Jr., M. Batty, S. Havlin,
and H. E. Stanley, Phys. Rev. E 58, 7054 (1998).
[12] M. Sahimi, Phys. Rep. 306, 213 (1998).
[13] H. A. Makse, J. S. Andrade Jr.,
and H. E. Stanley,
Phys. Rev. E 61, 583 (2000).
[14] A. D. Ara´ujo, A. A. Moreira, H. A. Makse, H. E. Stan-
ley,
and J. S. Andrade Jr., Phys. Rev. E 66, 046304
(2002).
[15] A. D. Ara´ujo, A. A. Moreira, R. N. Costa Filho,
and
J. S. Andrade Jr., Phys. Rev. E 67, 027102 (2003).
[16] N. Sandler, H. R. Maei, and J. Kondev, Phys. Rev. B
70, 045309 (2004).
[17] H. Auradou, K. J. M˚aløy, J. Schmittbuhl, A. Hansen,
and D. Bideau, Phys. Rev. E 60, 7224 (1999).
[18] J. Schmittbuhl, A. Hansen, H. Auradou,
and K. J.
M˚aløy, Phys. Rev. E 61, 3985 (2000).
[19] S. Smirnov and W. Werner, Math. Res. Lett. 8, 729
(2001).
[20] J. W. Essam, Rep. Prog. Phys. 43, 833 (1980).
[21] V. Privman and M. E. Fisher, Phys. Rev. B 30, 322
(1984).
[22] J. J. Binney, N. J. Dowrick, A. J. Fisher, and M. E. J.
Newman, The Theory of Critical Phenomena (Oxford
University Press, New York, 1992).
[23] A. Weinrib and B. I. Halperin, Phys. Rev. B 27, 413
(1983).
[24] W. Janke and M. Weigel, Phys. Rev. B 69, 144208
(2004).
[25] H.-O. Peitgen and D. Saupe, eds., The Science of Fractal
Images (Springer, New York, 1988).
[26] K. B. Lauritsen, M. Sahimi, and H. J. Herrmann, Phys.
Rev. E 48, 1272 (1993).
[27] H. A. Makse, S. Havlin, M. Schwartz, and H. E. Stanley,
Phys. Rev. E 53, 5445 (1996).
[28] H. G. Ballesteros and G. Parisi, Phys. Rev. B 60, 12912
(1999).
[29] B. D. Malamud and D. L. Turcotte, J. Stat. Plan. Infer.
80, 173 (1999).
[30] E. A. Oliveira, K. J. Schrenk, N. A. M. Ara´ujo, H. J.
Herrmann,
and J. S. Andrade Jr., Phys. Rev. E 83,
046113 (2011).
[31] B. Ahrens and A. K. Hartmann, Phys. Rev. B 84,
144202 (2011).
[32] P. A. Morais, E. A. Oliveira, N. A. M. Ara´ujo, H. J.
Herrmann,
and J. S. Andrade Jr., Phys. Rev. E 84,
016102 (2011).
[33] We note that the correlation parameter λ in Ref. [6]
is related to the Hurst exponent H used here by
λ = 2(H + 1). For the analogous parameter a of Ref. [5],
one has a = −2H.
[34] C. Domb and J. L. Lebowitz, eds., Phase Transitions
and Critical Phenomena, Vol. 11 (Academic, London,
1987).
[35] P. Di Francesco, P. Mathieu, and D. S´en´echal, Confor-
mal Field Theory (Springer, New York, 1997).
[36] M. Henkel and D. Karevski, eds., Conformal invariance:
An introduction to loops, interfaces and Stochastic
Loewner Evolution, Lecture Notes in Physics, Vol. 853
(Springer, Heidelberg, 2012).
[37] J. Kalda, EPL 84, 46003 (2008).
[38] J. Kondev, C. L. Henley, and D. G. Salinas, Phys. Rev.
E 61, 104 (2000).
[39] K. J. Schrenk, N. A. M. Ara´ujo, J. S. Andrade Jr., and
H. J. Herrmann, Sci. Rep. 2, 348 (2012).
[40] M. E. J. Newman and R. M. Ziﬀ, Phys. Rev. Lett. 85,
4104 (2000).
[41] M. E. J. Newman and R. M. Ziﬀ, Phys. Rev. E 64,
016706 (2001).
[42] H. Hu, H. W. J. Bl¨ote,
and Y. Deng, J. Phys. A 45,
494006 (2012).
[43] E. Fehr, D. Kadau, N. A. M. Ara´ujo, J. S. Andrade Jr.,
and H. J. Herrmann, Phys. Rev. E 84, 036116 (2011).
[44] D. K. C. MacDonald, Noise and Fluctuations (Wiley,
New York, 1962).
[45] V. V. Prudnikov and A. A. Fedorenko, J. Phys. A 32,
L399 (1999).
[46] V. V. Prudnikov, P. V. Prudnikov,
and A. A. Fe-
dorenko, Phys. Rev. B 62, 8777 (2000).
[47] T. Abete, A. de Candia, D. Lairez,
and A. Coniglio,
Phys. Rev. Lett. 93, 228301 (2004).

<!-- page 12 -->

12
[48] V. I. Marinov and J. L. Lebowitz, Phys. Rev. E 74,
031120 (2006).
[49] Z. Olami and R. Zeitak, Phys. Rev. Lett. 76, 247 (1996).
[50] M. F. Sykes and J. W. Essam, Phys. Rev. Lett. 10, 3
(1963).
[51] M. F. Sykes and J. W. Essam, J. Math. Phys. 5, 1117
(1964).
[52] H. A. Kramers and G. H. Wannier, Phys. Rev. 60, 252
(1941).
[53] We also determined psquare
c
of the square and pstar
c
of
the star lattice for various values of H and found that,
in contrast to the behavior of the triangular lattice, the
percolation threshold does depend on H. Our results
for psquare
c
(H) of the square lattice agree, within error
bars, with the ones reported in Ref. [6]. Also, we have,
within error bars, psquare
c
(H) + pstar
c
(H) = 1, consistent
with the matching property. Besides pc = 1/2, an ad-
ditional advantage of the triangular lattice is that the
cluster perimeters (see Sec. V) are well-deﬁned, avoid-
ing common problems encountered on the square lattice
[74, 77, 115, 116].
[54] J. Nagler, A. Levina, and M. Timme, Nat. Phys. 7, 265
(2011).
[55] W. Chen and R. M. D’Souza, Phys. Rev. Lett. 106,
115701 (2011).
[56] S. S. Manna, Physica A 391, 2833 (2012).
[57] K. J. Schrenk, A. Felder, S. Deﬂorin, N. A. M. Ara´ujo,
R. M. D’Souza, and H. J. Herrmann, Phys. Rev. E 85,
031103 (2012).
[58] S. D. S. Reis, A. A. Moreira,
and J. S. Andrade Jr.,
Phys. Rev. E 85, 041112 (2012).
[59] Y. Zhang, W. Wei, B. Guo, R. Zhang,
and Z. Zheng,
Phys. Rev. E 86, 051103 (2012).
[60] Y. S. Cho, S. Hwang, H. J. Herrmann, and B. Kahng,
Science 339, 1185 (2013).
[61] W. Chen, J. Nagler, X. Cheng, X. Jin, H. Shen,
Z. Zheng, and R. M. D’Souza, Phys. Rev. E 87, 052130
(2013).
[62] K. J. Schrenk, N. A. M. Ara´ujo, and H. J. Herrmann,
Phys. Rev. E 87, 032123 (2013).
[63] J. Hoshen and R. Kopelman, Phys. Rev. B 14, 3438
(1976).
[64] P. Grassberger, Physica A 262, 251 (1999).
[65] R. M. Ziﬀ, Phys. Rev. E 82, 051105 (2010).
[66] R. F. Voss, J. Phys. A 17, L373 (1984).
[67] B. Sapoval, M. Rosso, and J. F. Gouyet, J. Phys. Lett.
46, L149 (1985).
[68] R. M. Ziﬀ, P. T. Cummings, and G. Stell, J. Phys. A
17, 3009 (1984).
[69] A. Weinrib and S. A. Trugman, Phys. Rev. B 31, 2993
(1985).
[70] P. Grassberger, J. Phys. A 19, 2675 (1986).
[71] G. F. Lawler, O. Schramm, and W. Werner, Acta Math.
187, 237 (2001).
[72] G. F. Lawler, O. Schramm, and W. Werner, Acta Math.
187, 275 (2001).
[73] J. Asikainen, A. Aharony, B. Mandelbrot, E. Rausch,
and J.-P. Hovi, Eur. Phys. J. B 34, 479 (2003).
[74] D. A. Adams, L. M. Sander,
and R. M. Ziﬀ, J. Stat.
Mech. , P03004 (2010).
[75] J. Kalda, Phys. Rev. E 64, 020101(R) (2001).
[76] I. Mandre and J. Kalda, Eur. Phys. J. B 83, 107 (2011).
[77] T. Grossman and A. Aharony, J. Phys. A 20, L1193
(1987).
[78] H. Saleur and B. Duplantier, Phys. Rev. Lett. 58, 2325
(1987).
[79] M. Aizenman, B. Duplantier,
and A. Aharony, Phys.
Rev. Lett. 83, 1359 (1999).
[80] B. Duplantier, Phys. Rev. Lett. 84, 1363 (2000).
[81] J. Kondev and C. L. Henley, Phys. Rev. Lett. 74, 4580
(1995).
[82] B. B. Mandelbrot, The Fractal Geometry of Nature
(Freeman, New York, 1983).
[83] C. Tricot, J. F. Quiniou, D. Wehbi, C. Roques-Carmes,
and B. Dubuc, Revue Phys. Appl. 23, 111 (1988).
[84] F. Y. Wu, Rev. Mod. Phys. 54, 235 (1982).
[85] W. Janke and A. M. J. Schakel, Nucl. Phys. B 700, 385
(2004).
[86] C. M. Fortuin and P. W. Kasteleyn, Physica 57, 536
(1972).
[87] H. J. Herrmann, D. C. Hong,
and H. E. Stanley, J.
Phys. A 17, L261 (1984).
[88] H. J. Herrmann and H. E. Stanley, J. Phys. A 21, L829
(1988).
[89] P. Grassberger, J. Phys. A 25, 5475 (1992).
[90] Z. Zhou, J. Yang, Y. Deng, and R. M. Ziﬀ, Phys. Rev.
E 86, 061101 (2012).
[91] H. A. Makse, S. Havlin, P. C. Ivanov, P. R. King,
S. Prakash,
and H. E. Stanley, Physica A 233, 587
(1996).
[92] H. J. Herrmann and H. E. Stanley, Phys. Rev. Lett. 53,
1121 (1984).
[93] M. D. Rintoul and H. Nakanishi, J. Phys. A 25, L945
(1992).
[94] G. F. Lawler, O. Schramm, and W. Werner, Electron.
J. Probab. 7, 1 (2002).
[95] Y. Deng, H. W. J. Bl¨ote, and B. Nienhuis, Phys. Rev.
E 69, 026114 (2004).
[96] A. Coniglio, Phys. Rev. Lett. 62, 3054 (1989).
[97] O. Scholder, Int. J. Mod. Phys. C 20, 267 (2009).
[98] S. Roux, in Non-linearity and breakdown in soft con-
densed matter, Lecture Notes in Physics, Vol. 437,
edited by K. K. Bardhan, B. K. Chakrabarti,
and
A. Hansen (Springer, Heidelberg, 1994) p. 235.
[99] D. C. Hong, S. Havlin, H. J. Herrmann,
and H. E.
Stanley, Phys. Rev. B 30, 4083 (1984).
[100] B. Derrida and J. Vannimenus, J. Phys. A 15, L557
(1982).
[101] J. G. Zabolitzky, Phys. Rev. B 30, 4077 (1984).
[102] J.-M. Normand, H. J. Herrmann,
and M. Hajjar, J.
Stat. Phys. 52, 441 (1988).
[103] C. J. Lobb and D. J. Frank, Phys. Rev. B 30, 4090
(1984).
[104] D. J. Frank and C. J. Lobb, Phys. Rev. B 37, 302 (1988).
[105] N. Pos´e, N. A. M. Ara´ujo, and H. J. Herrmann, Phys.
Rev. E 86, 051140 (2012).
[106] J. S. Andrade Jr., H. J. Herrmann, A. A. Moreira, and
C. L. N. Oliveira, Phys. Rev. E 83, 031133 (2011).
[107] E. Fehr, K. J. Schrenk, N. A. M. Ara´ujo, D. Kadau,
P. Grassberger, J. S. Andrade Jr., and H. J. Herrmann,
Phys. Rev. E 86, 011117 (2012).
[108] K. J. Schrenk, N. A. M. Ara´ujo, and H. J. Herrmann,
Sci. Rep. 2, 751 (2012).
[109] E. Fehr, J. S. Andrade Jr., S. D. da Cunha, L. R. da
Silva, H. J. Herrmann, D. Kadau, C. F. Moukarzel, and
E. A. Oliveira, J. Stat. Mech. , P09007 (2009).
[110] E. Fehr, D. Kadau, J. S. Andrade Jr., and H. J. Her-
rmann, Phys. Rev. Lett. 106, 048501 (2011).

<!-- page 13 -->

13
[111] E. Daryaei, N. A. M. Ara´ujo, K. J. Schrenk, S. Rouhani,
and H. J. Herrmann, Phys. Rev. Lett. 109, 218701
(2012).
[112] R. F. S. Andrade and H. J. Herrmann, Phys. Rev. E 87,
042113 (2013).
[113] E. Daryaei and S. Rouhani, arXiv:1308.5692.
[114] E. Daryaei, arXiv:1309.1207.
[115] A. A. Saberi, J. Stat. Mech. , P07030 (2009).
[116] C. Chatelain, J. Stat. Mech. , P08004 (2010).
