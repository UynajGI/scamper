# Continuum Percolation Thresholds in Two Dimensions

> Phys. Rev. E 86, 061109 (2012) · [arXiv](https://arxiv.org/abs/1209.4936)
> Local copy: `docs/ref/Mertens_2012_PRE_continuum_percolation_thresholds_2d.pdf` (6 pages)

<!-- page 1 -->

arXiv:1209.4936v2  [cond-mat.stat-mech]  8 Dec 2012
Continuum Percolation Thresholds in Two Dimensions
Stephan Mertens1, 2, ∗and Cristopher Moore1, †
1Santa Fe Institute, 1399 Hyde Park Rd., Santa Fe, NM 87501, USA
2Institut f¨ur Theoretische Physik, Universit¨at Magdeburg, Universit¨atsplatz 2, 39016 Magdeburg, Germany
(Dated: November 27, 2024)
A wide variety of methods have been used to compute percolation thresholds. In lattice percolation,
the most powerful of these methods consists of microcanonical simulations using the union-ﬁnd
algorithm to eﬃciently determine the connected clusters, and (in two dimensions) using exact values
from conformal ﬁeld theory for the probability, at the phase transition, that various kinds of wrapping
clusters exist on the torus. We apply this approach to percolation in continuum models, ﬁnding
overlaps between objects with real-valued positions and orientations. In particular, we ﬁnd precise
values of the percolation transition for disks, squares, rotated squares, and rotated sticks in two
dimensions, and conﬁrm that these transitions behave as conformal ﬁeld theory predicts. The running
time and memory use of our algorithm are essentially linear as a function of the number of objects at
criticality.
I.
INTRODUCTION
For more than 50 years, percolation theory has been
used to model static and dynamic properties of porous
media and other disordered physical systems [1–3].
Most natural systems correspond to continuum perco-
lation, yet most analytical and numerical work has fo-
cused on lattice percolation.
This is reasonable since
continuum and lattice percolation lie in the same uni-
versality class.
For properties that are non-universal,
however, such as the location of the threshold, one has
to study discrete and continuum models individually,
and it is also satisfying to conﬁrm universality experi-
mentally by measuring critical exponents and crossing
probabilities.
In this contribution we discuss an algorithm to com-
pute the location of the transition in continuum perco-
lation models. The algorithm works in arbitrary dimen-
sions, and for arbitrarily shaped objects; here we focus
on two-dimensional percolation with disks, squares that
are aligned or randomly rotated, and randomly rotated
sticks (see Figure 1). Our algorithm is an adaption of the
union-ﬁnd algorithm of Newman and Ziﬀ[4], the fastest
known algorithm for lattice percolation. We show that it
can be adapted to continuum percolation with the aid of
some simple additional data structures, and we back up
our claim by computing numerical values of the transi-
tion points that extend the accuracy of previously known
values by several orders of magnitude.
In two-dimensional continuum percolation, a number
n of penetrable objects are thrown at random in a square
of size L2. If the mean density ρ = n/L2 is ﬁnite as n and
L go to inﬁnity, the spatial distribution of the objects’
centers is a Poisson point process with density ρ. The
system percolates if there exists a cluster of overlapping
objects that spans the square.
We follow [4] in using
∗mertens@ovgu.de
† moore@santafe.edu
Figure 1. (Color online) Continuum percolation with disks,
randomly rotated sticks, and aligned or rotated squares. In
each example, the wrapping cluster is marked by color.
periodic boundary conditions, and focusing on clusters
that wrap around horizontally, vertically, or both. These
wrapping clusters display better ﬁnite-size eﬀects than
crossing clusters on open boundary conditions.
If each object has area a, then the probability that a
percolating cluster exists in the limit L →∞clearly de-
pends only on the product η = ρa. This dimensionless
quantity is called the ﬁlling factor. It also gives the total
fraction φ of the plane covered by the objects,
φ = 1 −e−η .
(1)
While we expect continuum percolation to be in the same
universality class for any ﬁxed shape, the location of the

<!-- page 2 -->

2
transition, i.e., the critical ﬁlling factor ηc, depends on the
shape of the objects. We write η
c , η□
c , η^
c , and η×
c for the
percolation of disks, aligned squares, randomly rotated
squares, and randomly rotated sticks. In deﬁning η, we
treat sticks of length ℓas if they have area a = ℓ2.
Table I lists the most accurate numerical values for
ηc from previous work and the work presented here.
The best previous results on disk percolation are due to
Quintanilla, Torquato, and Ziﬀ[5] who varied the den-
sity of the Poisson process as a function of position and
kept track of the front of the connected cluster. The best
previous results on aligned squares are due to Torquato
and Jaio [7], who rescale an initial set of particles so that
its density is close to rigorous bounds. The best previ-
ous results on rotated squares are due to Baker et al. [6].
The best previous results on sticks are due to Li and
Zhang [8], who used an approach similar to ours but
with open boundary conditions.
Our results are consistent with the rigorous bounds
1.127 ≤η
c ≤1.12875
1.098 ≤η□
c ≤1.0995 ,
(2)
computed with 99.99% conﬁdence by Balister, Bollob´as
and Walters [9] using a Monte Carlo estimate of a high-
dimensional integral. On the other hand, it is a little sad
to dash the hope—which one might have entertained
after reading [6, 7], and which is just barely consistent
with (2)—that φ□
c is exactly 2/3.
In the following sections, we review the union-ﬁnd
algorithm of [4], how it ﬁnds wrapping clusters in peri-
odic boundary conditions, and how we extend it to the
continuous case. We show that the running time of our
algorithm is essentially linear in the number of objects,
i.e., linear in L2. In addition to estimating the threshold,
we also measure the ﬁnite-size exponent ν, giving strong
evidence that these continuum models are in the same
universality class as lattice percolation. Finally, we ﬁnd
that the probability of a wrapping cluster at criticality is
precisely that predicted by conformal ﬁeld theory.
II.
THE ALGORITHM
We will simulate percolation in the microcanonical en-
semble, i.e., where the number n of objects in the square is
ﬁxed. In each trial, we add one object at a time, stopping
as soon as a percolating cluster appears. Following [4],
we keep track of the connected components at each step
using the union-ﬁnd data structure. In union-ﬁnd, each
cluster is represented uniquely by one of its members.
We have access to two functions: find(i), which ﬁnds the
representative r(i) of the cluster to which object i belongs,
and merge(i, j), which merges i’s cluster and j’s cluster
together into a single one with the same representative.
Internally, union-ﬁnd works in a very simple way.
Each object i is linked to a unique “parent” p(i) in the
same cluster, except for the representative which has no
parent. When we call find(i), it follows the links from i
1
2
3
4
5
1
2
3
4
5
Figure 2. When we call find(i), we split and shorten the path
from i to its representative r(i) by setting the parent of each
object along the path to be its grandparent. This turns a path
of length ℓinto two paths of length ℓ/2.
to its parent p(i), its grandparent p(p(i)), and so on, until
it reaches i’s representativer(i). Similarly, merge(i, j) uses
find(i) and find(j) to obtain r(i) and r(j), and declares
one of them to be the parent of the other, unless r(i) = r(j)
and they are already in the same cluster.
The running time of find(i) is proportional to the
length of the path from i to r(i).
If merge(i, j) sensi-
bly links the smaller cluster to the larger one, setting
p(r(i)) = r(j) whenever i’s cluster is smaller than j’s, a
simple inductive argument shows that these paths never
exceed log2 n in length. However, we can make these
paths even shorter using a trick called path compression.
Since r(i) is the representative of every object j along the
path from i to r(i), we can set p(j) = r(i) for all of them,
linking them directly to their representative so that find
will work in a single step the next time we call it.
As a result, the amortized cost of the find and merge
operations—that is, the average cost per operation over
the course of many operations—is nearly constant.
Speciﬁcally, it is proportional to α(n), when α is the in-
verse of the Ackermann function [10]. The Ackermann
function grows faster than any primitive recursive func-
tion, i.e., any function that can be computed with a ﬁxed
number of for-loops: faster than an exponential, an iter-
ated tower of exponentials, and so on [11]. As a conse-
quence, α(n) grows incredibly slowly, and the smallest
value of n such that α(n) > 4 is so large that it can only
be written with exotic notation. Thus the total running
time for n objects is essentially O(n).
In our implementation, we employ a form of path
compression that is faster and almost as eﬀective: we
link each object j on the path to its grandparent, setting
p(j) = p(p(j)). This is known as path splitting, since it
turns a path of length ℓinto two paths of length ℓ/2, or
(ℓ+ 1)/2 and (ℓ−1)/2 if ℓis odd, as shown in Figure 2. It
has the advantage of requiring only one pass along the
path, and it takes just one line of code (e.g. [4, Appendix
A]). Like path compression, it guarantees an amortized
running time of O(α(n)) [12].
For lattice percolation as in [4], each time we add a
new occupied site, we can check which of its neighbors
are occupied, and merge them together with the new

<!-- page 3 -->

3
η
c
η□
c
η^
c
η×
c
previous 1.128085(2)
1.0982(3)
0.9819(6)
5.63726(2)
our work 1.12808737(6) 1.09884280(9) 0.9822723(1) 5.6372858(6)
φ
c
φ□
c
φ^
c
φ×
c
previous 0.6763475(6)
0.6665(1)
0.6254(2)
0.99643738(7)
our work 0.67634831(2) 0.66674349(3) 0.62554075(4) 0.996437475(2)
Table I. Numerical values of critical ﬁlling factors ηc and area factors φc = 1 −e−ηc in continuum percolation for disks, aligned
squares, randomly rotated squares, and randomly rotated sticks. Previous estimates are from [5–8].
Figure 3. (Color online) We divide the plane into square bins
whose width equals the diameter of the disks. Each disk in a
given bin (dashed) can only intersect with other disks in the
same bin, or in the eight neighboring bins.
site. In the continuous case, we have more work to do: if
we add a new disk (say), we have to ﬁnd which nearby
disks it intersects. To do this eﬃciently, we divide the
plane into square bins as shown in Figure 3. Each disk
belongs to whichever bin its center lies in. The width of
each bin is the diameter of the disks, so that a disk in a
given bin can only intersect with other disks in that bin
or the eight bins in its neighborhood.
On average, the number of disks in each bin is a con-
stant proportional to ρ, so we can ﬁnd all the disks inter-
secting with each new one in constant time. We use the
same approach for the other shapes; for rotated squares
of width ℓ, the bins need to have width
√
2ℓ. A similar
approach for rotated sticks was used in [8].
If we wished to detect crossing clusters—those that
connect, say, the top and bottom edges of the square—
we could add two special objects to the union-ﬁnd data
structure, which are connected by ﬁat to all the disks in
the bins along the top or bottom edge. We would then
check, at each step, whether these two objects are in the
same cluster. However, as discussed below and in [4],
the ﬁnite-size scaling is much better if we use periodic
boundary conditions instead, and look for clusters that
wrap around the torus horizontally or vertically.
We detect these wrapping clusters using a technique
originally used for detecting crossing clusters in the Potts
model [13]. We associate a vector with each object in the
union-ﬁnd data structure, recording the displacement
between it and its parent. In principle this displacement
is real-valued, but it suﬃces to record an integer vector
giving the displacement between their respective bins.
When we compress and splint a path, we sum these
vectors to get the total displacement between each object
on the path and its new parent.
Now suppose that merge(i, j) ﬁnds that two overlap-
ping disks i and j are already in the same cluster. Object
i now has two paths to its representative; one that goes
through its own parent, and another that consists of hop-
ping to j and then going through j’s parents. We sum
the displacement vectors along both these paths. If these
sums are the same, then the cluster is simply-connected.
But if they diﬀer by ±L in either coordinate, then the
cluster has a nontrivial winding number around one or
both directions on the torus.
Like the union-ﬁnd algorithm itself, the time it takes
to sum these vectors is proportional to the length of the
paths from i and j to their representative. As Figure 4
shows, the total running time of our entire algorithm—
the time it takes to carry out a trial on a lattice of size
L, adding objects one at a time until a wrapping cluster
appears—is essentially linear in the number n of objects
at criticality, or equivalently linear in L2. It slows down
somewhat when the computer is forced to switch to parts
of its memory with slower access, but this only aﬀects
the leading constant.
III.
ANALYSIS AND RESULTS
If in each trial we stop at the ﬁrst n where a wrapping
cluster appears, then the estimated probability PL(a, n)
that a wrapping cluster exists in the microcanonical en-
semble with n objects of area a is the fraction of trials that
stop on or before the nth step. To obtain the probabil-
ity RL(η) of percolation in the grand canonical ensemble
with ﬁlling fraction η, we convolve PL with the Poisson
distribution with mean λ = ρL2 = ηL2/a:
RL(η) = e−λ
∞
X
n=0
λn
n! PL(a, n) .
(3)

<!-- page 4 -->

4
10
1
10
2
10
3
10
4
10
5
10
6
10
-5
10
-4
10
-3
10
-2
10
-1
10
0
L2
CPU time in seconds
squares
disks
sticks
rot. squares
Figure 4. Average CPU time T for a single realization of an
L × L lattice up to the percolation transition. T is measured on
a laptop with Intel Core 2 Duo 2.53 GHz CPU with 3 MB cache.
The dashed lines T ∝L2 are guides for the eye. The slope
increases when the cache memory is exhausted, forcing the
computer to switch to regions of memory with slower access;
however, the running time remains linear in n ∝L2.
To avoid numerical diﬃculties where the numerator and
denominator are both very large, we compute Poisson
weights wn ∝λn/n! inductively in two sequences w¯n−k
and w¯n+k to the left and right of the peak at ¯n = ⌊λ⌋,
where we deﬁne w¯n = 1:
w¯n−k =

1
for k = 0
¯n−(k−1)
λ
w¯n−(k−1)
for k = 1, 2, . . .
and
w¯n+k =

1
for k = 0
λ
¯n+k w¯n+k−1
for k = 1, 2, . . .
The sum (3) only needs to be computed for a ﬁnite num-
ber of terms. In one direction, we only need to sum down
to the smallest n where PL(a, n) is nonzero, i.e., the small-
est value of n where we observed a wrapping cluster in
at least one trial. In the other direction, once we pass
the largest n where a wrapping cluster ﬁrst appeared,
then PL(a, n) = 1. At that point, we sum the remaining
terms until they are zero to within the numerical preci-
sion of the computer. We then normalize the entire sum
by dividing by P wn.
Equipped with the data from the microcanonical sim-
ulations and this convolution routine, we compute the
wrapping probability functions RL(η) for various system
sizes L and shapes. Like [4], we look for several kinds of
wrapping in particular. Speciﬁcally:
• Re
L(η) is the probability of any kind of wrapping
cluster.
This is indicated by a winding number
that is nonzero in either coordinate.
0.8
0.9
1.0
1.1
1.2
1.3
0
0.2
0.4
0.6
0.8
1
η
Re
L(η)
Figure5. Wrapping probabilities Re
L(η) for disk percolation and
L = 16, 32, 64, 128, 256, 512. The dashed line is the exact value of
the critical wrapping probability Re
∞(ηc) from conformal ﬁeld
theory.
• Rh
L(η) is the probability of a cluster that wraps hor-
izontally. This is indicated by a winding number
that is nonzero in the ﬁrst coordinate.
• Rb
L(η) is the probability of a cluster that wraps both
horizontally and vertically. This is indicated by a
single winding number that is nonzero in both co-
ordinates, or a pair of winding numbers that are
nonzero in the ﬁrst and second coordinates respec-
tively.
• R1
L(η) is the probability of a cluster that wraps hor-
izontally, but not vertically. This is indicated by a
winding number that is nonzero in only the ﬁrst
coordinate.
For any L and any η, these probabilities obey
Re
L = 2Rh
L −Rb
L = 2R1
L + Rb
L .
We assume here that the torus is square, so that horizon-
tal and vertical wrapping probabilities are equal.
Note that if the ﬁrst nonzero winding number ob-
served in a given trial is nonzero in both coordinates,
then a cluster of type 1 (horizontal but not vertical) does
not occur at all in that trial. Thus R1
L(η) does not tend to
1 as η increases.
In practice, we focused on Re
L and Rb
L. In each run, we
recorded the number of objects nh at which horizontal
wrapping ﬁrst occurred, and the number nv at which
vertical wrapping ﬁrst occurred. Then ne = min(nh, nv)
and nb = max(nh, nv) are our estimates, in that run, of the
values of n at which Re
L and Rb
L jump from 0 to 1.

<!-- page 5 -->

5
10
1
10
2
10
3
10
0
10
1
10
2
L
Re
L
′(ηL)
Figure6. Slope of Re
L(η) at ηL, the estimated critical ﬁlling factor.
The line is 0.361L3/4, conﬁrming the universal critical exponent
ν = 4/3 for ﬁnite-size scaling.
A beautiful fact is that, even though the percolation
threshold ηc is not known for any of our models, confor-
mal ﬁeld theory implies exact values for these probabil-
ities at the transition in the limit L →∞[4, 14]. Speciﬁ-
cally,
Rh
∞= 0.521 058 289 248 821 787 848...
Re
∞= 0.690 473 724 570 168 677 230...
Rb
∞= 0.351 642 853 927 474 898 465...
R1
∞= 0.169 415 435 321 346 889 383...
(4)
For each L, and each type of wrapping cluster, we can
estimate the critical ﬁlling factor ηL as the solution of the
equation
RL(ηL) = R∞.
(5)
For instance, Figure 5 shows Re
L(η) for disks for L ranging
up to 512. The ﬁlling factors ηL where these curves cross
Re
∞rapidly converge to ηc.
The rate of convergence is determined by two factors.
The ﬁrst comes from the fact that the width of the tran-
sition window from RL ≈0 to RL ≈1 scales as L−1/ν
where ν = 4/3 is a universal critical exponent for two-
dimensional percolation.
This scaling holds even for
small systems, as can be seen in Figure 6, where we plot
the slope of RL at the estimated critical ﬁlling factor ηL.
The slope scales perfectly like L3/4.
The second factor comes from the fact that RL(η) not
only becomes steeper but also moves upward in the crit-
ical region (see the inset in Figure 5). To measure the
contribution from this eﬀect, we computed the diﬀer-
ence RL(ηc) −Re
∞using the previously best known value
for ηc from Table I. This diﬀerence scales like L−2, as can
be seen from Figure 7. The exponent −2 correponds to
10
100
10
-5
10
-4
10
-3
L
Re
L(η!
c ) −Re
∞(η!
c )
Figure 7. Convergence with increasing system size L of Re
L(η
c )
to its known value at L = ∞for disk percolation. The line is
proportional to L−2, the conjectured convergence of Re
L(ηc).
0
2!10
-5
4!10
-5
6!10
-5
8!10
-5
1!10
-4
1.09884
1.09885
1.09886
1.09887
1.09888
1.09889
1.09890
1.09891
L−11/4
ηL
Figure 8. Estimated critical ﬁlling factors for continuum per-
colation of squares, derived from (5) with Re
L (top) and Rb
L
(bottom). Note the resolution of the η-axis.
the leading irrelevant renormalization exponent yi in the
Kac table [15]. Note that the periodic boundary condi-
tions are responsible for this decay. With open boundary
conditions, this factor scales as L−1 [16], leading to more
severe ﬁnite-size eﬀects.
These two factors combine to give
ηL −ηc ∼L−2−1/ν = L−11/4
(6)
for the rate of convergence. Hence we expect a straight
line if we plot ηL vs. L−11/4, and this is exactly what we
observe in Figure 8. Extrapolating this line to zero then
gives our estimates of ηc shown in Table I.

<!-- page 6 -->

6
How do we compute the error bars in our estimates of
ηc? First consider the ﬂuctuations in RL(η). Each of our
microcanonical experiments contributes to our estimate
of RL(η) for all η through the convolution (3). We can
imagine this as choosing n from the Poisson distribution,
adding n objects, and returning an estimate of RL(η) = 1
or 0 depending on whether they percolate or not.
If
we perform N trials, the number of trials that return 1 is
binomially distributed with mean RL(η)N, and averaging
gives an estimate of RL(η) with standard deviation
σRL =
s
RL(η)

1 −RL(η)

N
.
(7)
Depending on which kind of wrapping cluster we are
looking for, this is roughly 0.4N−1/2.
When we look for the ηL where RL(η) crosses R∞, the
error on ηL is given by
σηL =
σRL
R′
L(ηL) .
Since the slope R′
L(ηL) grows as 0.361L3/4 (see Figure 6)
this gives
σηL ≈N−1/2 L−3/4 .
These are the error bars shown in Figure 8.
The extrapolated value for ηc is computed from simu-
lations for D diﬀerent system sizes L, which in a weighted
linear regression as in Figure 8 yields an error roughly
√
D times smaller than the error bars of the underlying
data points.
Finally, we average our estimates of ηc from Re
L and
Rb
L. Assuming that these estimates are only weakly cor-
related reduces the error bars by another factor of
√
2.
The error bars shown in Table I are the result of sim-
ulating roughly D = 50 system sizes ranging from L = 8
to L = 2048, with sample sizes N ranging from 1010 for
the systems with L ≤100, to 109 for 100 < L ≤500, to 106
for 500 < L ≤2048.
We ran these simulations in parallel on several
computer clusters with greatly varying computational
power. In total, our simulations would have taken about
400 years if done only on the laptop on which this paper
was written.
IV.
CONCLUSIONS
We have shown that the union-ﬁnd approach to esti-
mating percolation thresholds introduced by Newman
and Ziﬀ[4] can be applied in the continuous case. With
the help of an algorithm for estimating ηc that runs in
essentially linear time as a function of the number of ob-
jects at criticality, we have obtained new estimates for
ηc in a variety of continuum percolation models that are
several orders of magnitude more accurate than previous
results. In the process, we have conﬁrmed the predic-
tions of conformal ﬁeld theory for these models, both for
the ﬁnite-size scaling exponent ν and the probabilities
that various kinds of wrapping clusters exist at ηc on
periodic boundary conditions.
ACKNOWLEDGMENTS
S.M. thanks the Santa Fe Institute for their hospitality.
C.M. is supported by the National Science Foundation
through grant CCF-1219117 and by the Air Force Oﬃce
of Scientiﬁc Research and the Defense Advanced Re-
search Projects Agency through grant FA9550-12-1-0432.
We are grateful to Robert Ziﬀand Mark Newman for
helpful conversations.
[1] E. N. Gilbert, Journal of the Society for Industrial and
Applied Mathematics 9, 533 (1961).
[2] D. Stauﬀer and A. Aharony, Introduction to Percolation The-
ory (Taylor & Francis, London, 1992).
[3] M. Sahimi, Applications of Percolation Theory (Taylor &
Francis, London, 1994).
[4] M. E. J. Newman and R. M. Ziﬀ, Physical Review E 64,
016706 (2001).
[5] J. A. Quintanilla and R. M. Ziﬀ, Physical Review E 76,
051115 (2007).
[6] D. R. Baker, G. Paul, S. Sreenivasan, and H. E. Stanley,
Physical Review E 66, 046136 (2002).
[7] S. Torquato and Y. Jiao, Journal of Chemical Physics 137,
074106 (2012).
[8] J. Li and S.-L. Zhang, Physical Review E 80, 040104(R)
(2009).
[9] P. Balister, B. Bollob´as, and M. Walters, Random Struc-
tures & Algorithms 26, 392 (2005).
[10] R. E. Tarjan, Journal of the Association for Computing
Machinery 22, 215 (1975).
[11] C.
Moore
and
S.
Mertens,
The
Nature
of
Computation
(Oxford
University
Press,
2011)
www.nature-of-computation.org.
[12] R. E. Tarjan and J. van Leeuwen, Journal of the Association
for Computing Machinery 31, 245 (1984).
[13] J. Machta, Y. S. Choi, A. Lucke, T. Schweizer, and L. M.
Chayes, Physical Review E 54, 1332 (96).
[14] H. T. Pinson, Journal of Statistical Physics 75, 1167 (1994).
[15] H. Hu, H. W. Bl¨ote, and Y. Deng, “Percolationinthe canon-
ical ensemble,” http://arxiv.org/abs/1210.3463 (2012).
[16] J.-P. Hovi and A. Aharony, Physical Review E 53, 235
(1996).
