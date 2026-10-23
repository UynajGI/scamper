# On the scaling of the chemical distance in long-range percolation

> Ann. Probab. 32, 2938 (2004) · [arXiv](https://arxiv.org/abs/math/0304418)
> Local copy: `docs/ref/Biskup_2003_long_range_chemical_distance.pdf` (42 pages)

<!-- page 1 -->

arXiv:math/0304418v4  [math.PR]  6 Apr 2005
The Annals of Probability
2004, Vol. 32, No. 4, 2938–2977
DOI: 10.1214/009117904000000577
c
⃝Institute of Mathematical Statistics, 2004
ON THE SCALING OF THE CHEMICAL DISTANCE IN
LONG-RANGE PERCOLATION MODELS
By Marek Biskup
University of California, Los Angeles
We consider the (unoriented) long-range percolation on Zd in
dimensions d ≥1, where distinct sites x,y ∈Zd get connected with
probability pxy ∈[0,1]. Assuming pxy = |x−y|−s+o(1) as |x−y| →∞,
where s > 0 and | · | is a norm distance on Zd, and supposing that
the resulting random graph contains an inﬁnite connected compo-
nent C∞, we let D(x,y) be the graph distance between x and y mea-
sured on C∞. Our main result is that, for s ∈(d,2d),
D(x,y) = (log |x −y|)∆+o(1),
x,y ∈C∞, |x −y| →∞,
where ∆−1 is the binary logarithm of 2d/s and o(1) is a quantity
tending to zero in probability as |x −y| →∞. Besides its interest for
general percolation theory, this result sheds some light on a question
that has recently surfaced in the context of “small-world” phenomena.
As part of the proof we also establish tight bounds on the probabil-
ity that the largest connected component in a ﬁnite box contains a
positive fraction of all sites in the box.
1. Introduction.
1.1. Motivation.
Percolation is a simple but versatile model with appli-
cations ranging from the study of phase transitions in mathematical physics
to opinion spreading in social sciences. The most well-understood questions
of percolation theory are those concerning the appearance and uniqueness
of the inﬁnite component [10], uniqueness of the critical point [1, 16, 20],
decay of connectivity functions [11, 12], and the scaling properties at the
critical point in dimensions d = 2 [24, 25] and d large enough [17, 18]. Less
well understood remain natural questions about the qualitative structural
and geometrical properties of the inﬁnite connected component, especially
Received April 2003; revised September 2003.
AMS 2000 subject classiﬁcations. Primary 60K35; secondary 82B43, 82B28.
Key words and phrases. Long-range percolation, chemical distance, renormalization,
small-world phenomena.
This is an electronic reprint of the original article published by the
Institute of Mathematical Statistics in The Annals of Probability,
2004, Vol. 32, No. 4, 2938–2977. This reprint diﬀers from the original in
pagination and typographic detail.
1

<!-- page 2 -->

2
M. BISKUP
below the upper critical dimension. In particular, this includes the tantaliz-
ing open problem concerning the absence of percolation at the percolation
threshold.
Long-range versions of the percolation model have initially been intro-
duced in order to study the eﬀect of long-range interaction on the onset
of phase transition in one-dimensional systems. On Z, the most common
setup is that, in addition to random nearest-neighbor connections with prob-
ability p ∈(0,1), a bond between x,y ∈Z is added with probability 1 −
exp{−β|x −y|−s}, where β ∈(0,∞) and s > 0. In dimension 1, the interest-
ing ranges of values of s are s < 1, where the resulting graph is almost surely
connected [23], 1 < s < 2, where an inﬁnite component appears once p is large
enough [22], and the critical case s = 2, where the inﬁnite component ap-
pears “discontinuously” for some p < 1 suﬃciently large if and only if β > 1
[2] and where the truncated connectivity function decays with a β-dependent
exponent [19] for β in the interval (1,2). The cases s > 2 are qualitatively
very much like the nearest-neighbor case (in particular, there is no percola-
tion for p < 1 and β < ∞). In dimensions d > 1, the insertion of long-range
connections is not essential for the very existence of percolation—the main
problem of interest there is to quantify the eﬀect of such connections on the
critical behavior.
In this paper we study the global scaling properties of the inﬁnite com-
ponent in long-range percolation models on Zd for arbitrary d. We focus on
the scaling of the graph distance (aka chemical distance) in the cases when
the probability that a bond is occupied falls oﬀwith exponent s ∈(d,2d).
More precisely, we let distinct x,y ∈Zd be connected independently with
probability pxy that has the asymptotics pxy = 1 −exp{−|x −y|−s+o(1)}
as |x −y| →∞. Assuming that there is a unique inﬁnite connected com-
ponent C∞almost surely, we let D(x,y) be the distance between the sites
x and y measured on C∞. Then we prove that D(x,y) scales with the Eu-
clidean distance |x −y| as
D(x,y) = (log |x −y|)∆+o(1),
x,y ∈C∞, |x −y| →∞,
(1.1)
where ∆= ∆(s,d) is given by
∆(s,d) =
log 2
log(2d/s).
(1.2)
This result should be contrasted with those of [[5]–[7, 14]] (see also [4]), where
various (other) regimes of decay of long-range bond probabilities have been
addressed. We refer to Section 1.3 for further discussion of related work and
an account of the current state of knowledge about the asymptotic behaviors
of D(x,y).

<!-- page 3 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
3
The nonlinear scaling (1.1) is a manifestation of the fact that adding
sparse (but dense-enough) long edges to a “Euclidean” graph may substan-
tially alter the graph geometry and, in particular, its scaling properties. This
is exactly what has recently brought long-range percolation renewed atten-
tion in the studies of the so-called “small-world” phenomena; see [26] for an
initial work on these problems. This connection was the motivation of the
work by Benjamini and Berger [6], who studied how the (graph) diameter
of a ﬁnite ring of N sites changes when long connections are added in. On
the basis of a polylogarithmic upper bound, the authors of [6] conjectured
(cf. Conjecture 3.2 in [6]) that in the regime when s ∈(d,2d), the diam-
eter scales as (log N)γ, where γ = γ(s) > 1. The present paper provides a
polylogarithmic lower bound in this conjecture. However, at present it is
not clear whether the exponent for the diameter growth matches that for
the typical distance between two remote points. We refer to Section 1.4 for
further discussion of “small-world” phenomena.
The remainder of this paper is organized as follows. In Section 1.2 we
deﬁne precisely the long-range percolation model and state our main theo-
rem (Theorem 1.1). In Section 1.3 we proceed by summarizing the previous
results concerning the behavior of D(x,y)—and graph diameter—for var-
ious regimes of s. In Section 1.4 we discuss the relation to “small-world”
phenomena. Section 2 is devoted to a heuristic explanation of the proof of
Theorem 1.1. The proof requires some preparations, particularly an estimate
on the size of the largest connected component in large but ﬁnite boxes. This
is the content of Theorem 3.2 in Section 3. The actual proof of our main
result comes in Sections 4.1 (upper bound) and 4.2 (lower bound).
1.2. The model and main result.
Consider the d-dimensional hypercubic
lattice Zd and let (x,y) 7→|x−y| denote a norm distance on Zd. For deﬁnite-
ness, we can take | · | to be the usual ℓ2-norm; however, any other equivalent
norm will do. Let q :Zd →[0,∞) be a function satisfying
lim
|x|→∞
log q(x)
log |x| = −s,
(1.3)
where s ≥0. (Here we set log 0 = −∞.) For each (unordered) pair of distinct
sites x,y ∈Zd, we introduce an independent random variable ωxy ∈{0,1}
with probability distribution given by P(ωxy = 1) = pxy, where
pxy = 1 −e−q(x−y).
(1.4)
Note that pxy = px−y,0 so the distribution of (ωxy) is translation invariant.
Let G be the random graph with vertices on Zd and a bond between
any pair of distinct sites x and y, where ωxy = 1. Given a realization of G ,
let us call π = (z0,z1,...,zn) a path, provided zi are all distinct sites in Zd

<!-- page 4 -->

4
M. BISKUP
and ωzi−1zi = 1 for each i ∈{1,2,... ,n}. Deﬁne the length |π| of π to be the
number of bonds constituting π (i.e., the number n above). Using Π(x,y) to
denote the (random) set of all paths π with z0 = x and z|π| = y, we let
D(x,y) = inf{|π|:π ∈Π(x,y)},
x,y ∈Zd.
(1.5)
[In particular, we have D(x,y) = ∞if Π(x,y) = ∅.] The random variable D(x,y)
is the chemical distance between x and y, that is, the distance measured on
the graph G .
Throughout the rest of the paper, it will be assumed that the random
graph G contains an inﬁnite connected component. We will focus on the
cases when s ∈(d,2d) in (1.3), in which percolation can be guaranteed, for
instance, by requiring that the minimal probability of a nearest-neighbor
connection, p, is suﬃciently close to 1. (Indeed, in d ≥2, it suﬃces that p
exceeds the percolation threshold for bond percolation on Zd, while in d =
1, this follows by the classic result of [22].) Moreover, by an extension of
Burton–Keane’s uniqueness argument due to [15], the inﬁnite component is
unique almost surely. We will use C∞to denote the set of sites in the inﬁnite
component of G .
Our main result is as follows:
Theorem 1.1.
Suppose that (1.3) holds with an s ∈(d,2d) and assume
that, P-almost surely, the random graph G contains a unique inﬁnite com-
ponent C∞. Then for all ε > 0,
lim
|x|→∞P

∆−ε ≤log D(0,x)
log log |x| ≤∆+ ε
0,x ∈C∞

= 1,
(1.6)
where ∆= ∆(s,d) is as in (1.2).
Formula (1.6) is a precise form of the asymptotic expression (1.1). The
fact that ∆−1 is the binary logarithm of 2d/s is a consequence of the fact
that the longest bonds in the shortest path(s) between two distant sites
of C∞exhibit a natural binary hierarchical structure; see Section 2 for
more explanation. Note that s 7→∆(s,d) is increasing throughout (d,2d)
and, in particular, ∆(s,d) > 1 for all s ∈(d,2d) with lims↓d ∆(s,d) = 1
and lims↑2d ∆(s,d) = ∞.
Remark 1.1.
The requirement of translation invariance is presumably
not crucial for (the essence of ) the above result. Indeed, most of our proofs
should carry through under the weaker assumption of approximate homo-
geneity on large spatial scales. Notwithstanding, some of our arguments in
Section 3 are based on previous results that require translation invariance
and so we stick with the present setting for the rest of this paper.

<!-- page 5 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
5
1.3. Discussion.
As already alluded to, several diﬀerent asymptotic be-
haviors are possible in the above problem depending on the value of the
exponent s. We proceed by reviewing the known (and conjectured) results.
Throughout, we will focus on the speciﬁc distribution
pxy = 1 −exp{−β(1 + |x −y|)−s},
(1.7)
where β ∈[0,∞). (Some of the results also required that all nearest-neighbor
connection be a priori present.) We will concentrate on the asymptotic of two
quantities: The typical graph distance D(x,y)—the focus of this paper—and
the diameter DN of the graph obtained by “decorating” a box of N ×···×N
sites in Zd by the bonds in G with both endpoints therein. There are ﬁve
distinct regimes marked by the position of s relative to the numbers d and 2d.
The cases of s < d fall into the category of problems that can be analyzed
using the concept of stochastic dimension introduced in [7]. The result is the
almost-sure equality
sup
x,y∈Zd D(x,y) =

d
d −s

;
(1.8)
see Example 6.1 in [7]. A similar asymptotic statement holds for the N →∞
limit of DN; see Theorem 4.1 of [6].
For s = d, Coppersmith, Gamarnik and Sviridenko [14] study the asymp-
totic of DN. The resulting scaling is expressed by the formula
DN = Θ(1)
log N
log log N ,
N →∞,
(1.9)
where Θ(1) is a quantity bounded away from 0 and ∞. Since the typical
distance is always less than the diameter, this shows that D(x,y) will grow
at most logarithmically with |x −y|. However, at present the appropriate
lower bound on D(x,y) is missing.
In the cases d < s < 2d, Benjamini and Berger [6] and Coppersmith,
Gamarnik and Sviridenko [14] proved polylogarithmic upper bounds on DN
[and hence on D(x,y) for |x −y| ≈N]. However, the best lower bound these
references gave was proportional to log N. The present paper provides a
sharp leading-order asymptotic for D(x,y) which constrains DN to grow
at least as fast as (log N)∆+o(1). Unfortunately, neither the bounds from
[6] and [14] nor those derived for D(x,y) in the present paper are sharp
enough to make any deﬁnitive asymptotic statements about DN. We hope
to return to this question in a future publication.
The critical cases s = 2d are at present not very well understood. Here
Benjamini and Berger [6] conjectured that
DN = N θ(β)+o(1),
N →∞,
(1.10)

<!-- page 6 -->

6
M. BISKUP
with θ(β) ∈(0,1), and we expect a similar asymptotic to be valid for the typ-
ical distance D(x,y). A general upper bound on the above θ(β) was derived
in [14]. The corresponding—but not sharp—lower bounds were derived in
[6] and [14] under the restriction to the “nonpercolative” regime d = 1 and
β < 1. (This restriction appears because the proof relies heavily on the no-
tion of a cut-point; see [22].) Surprisingly, similarly open are the cases s > 2d
where we expect that DN scales linearly with N. The latter seems to have
been proved only in d = 1 [6], or for the case of the supercritical nearest-
neighbor percolation in d ≥2 [5]; see also [4].
An important technical resource for this paper has been the recent work of
Berger [8] on long-range percolation with exponents d < s < 2d. Employing
a variant of the renormalization scheme of [22], Berger proved among other
things the absence of critical percolation and, whenever there is percolation,
the existence of a cluster of at least N d−o(1) sites in any box of volume N d
(see Theorem 3.1). An extension of this result (see Theorem 3.2) establishing
tight bounds on the probability that the largest connected component in a
ﬁnite box contains a positive fraction of all sites is essential for the proof of
the upper bound in (1.6).
1.4. Relation to “small-world” phenomena.
As already mentioned, long-
range percolation has been recently used in the study of “small-world” phe-
nomena. The catchy term “small worlds” originates in the old but still-
fun-to-read article by Milgram [21], who observed through an ingenious
experiment that two typical Americans are just six acquaintances (or six
“handshakes”) away from each other. With the rise of the overall world
connectivity in recent years due to the massive expansion of air traﬃc, elec-
tronic communications and particularly the internet, and so on, the “small-
world” phenomena experienced a fair amount of new interest. Novel exam-
ples emerged in physics and biology, particularly after the publication of [26].
Several mathematical models were devised and studied using both rigorous
and nonrigorous techniques. A brief overview of the situation from the per-
spective of the theory of random graphs (and additional references) can be
found in Section 10.6 of [9].
While we will not attempt to summarize the content of the publication
boom following the appearance of [26], let us mention that a major deﬁciency
of many models introduced so far seems to be—at least to the author of the
present paper—the unclear status of their relevance to the actual (physical,
biological or sociological) systems of interest. In particular, a large fraction
of studied models seem to unjustly ignore the underlying spatial structure
present in the practical problem of interest. (The reason for that is most
likely the reduced complexity—as in statistical mechanics, models without
underlying geometry, the so-called mean-ﬁeld models, are often exactly solv-
able.) With this problem in mind, Benjamini and Berger [6] proposed a new

<!-- page 7 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
7
class of “small-world” models based on long-range percolation on Euclidean
graphs. More precisely, as an underlying graph they consider an a priori
connected ring of N sites to which long edges are added with probability as
described in (1.7).
One of the questions discussed by Benjamini and Berger was how the
diameter of the resulting random graph depends on N for various ranges
of values of s. As detailed in Section 1.3, this behavior depends rather sen-
sitively on the value of the exponent s. In particular, “phase transitions”
occur at s = d, which is the borderline of the region with ﬁnite diameters,
and s = 2d, which separates the regions of linear and sublinear scaling. Each
of the resulting behaviors may be useful in diﬀerent contexts. For instance, if
we believe Milgram’s assertion that six is the typical graph distance between
two average Americans regardless of the population size, the exponent s
should be within the regime described by (1.8).
2. Main ideas of the proof.
The proof of Theorem 1.1 consists of two
parts where we separately prove the upper and lower bounds in (1.6). Both
parts will be based on the concept of certain hierarchies of sites whose
deﬁnition is given below. In this deﬁnition—and elsewhere in this paper—
the symbol σ denotes a hierarchical index, σ ∈{0,1}k, which can be viewed
as a parametrization of the leaves of a binary tree of depth k. Thus, for
instance, σ = 01101 means that, starting from the root, we “go” left, right,
right, left and right to reach the leaf represented by σ. Adding digits behind σ
denotes index concatenation.
Definition 2.1.
Given an integer n ≥1 and distinct sites x,y ∈Zd, we
say that the collection
Hn(x,y) = {(zσ):σ ∈{0,1}k, k = 1,2,...,n; zσ ∈Zd}
(2.1)
is a hierarchy of depth n connecting x and y if:
1. z0 = x and z1 = y.
2. zσ00 = zσ0 and zσ11 = zσ1 for all k = 0,1,...,n −2 and all σ ∈{0,1}k.
3. For all k = 0,1,... ,n −2 and all σ ∈{0,1}k such that zσ01̸ = zσ10, the
bond between zσ01 and zσ10 is occupied, that is, (zσ01,zσ10) ∈G .
4. Each bond (zσ01,zσ10) as speciﬁed in part 3 appears only once in Hk(x,y).
In the following, the pairs of sites (zσ00,zσ01) and (zσ10,zσ11) will be referred
to as “gaps.”
Remark 2.1.
By assumption 2, a hierarchy of depth n is uniquely spec-
iﬁed by its nth level. Note that we do not require the sites of the hierarchy
to be distinct and, if two points of the form zσ10 and zσ11 coincide, we do not

<!-- page 8 -->

8
M. BISKUP
Fig. 1.
A schematic picture of a hierarchy of depth 5 connecting x = z0 and y = z1. The
straight line represents a path between x and y and the arches stand for the bonds between
pairs of sites (zσ01,zσ10). The arrows indicate the sites on levels 1–4 of the hierarchy;
the ﬁfth level has been omitted for brevity. Note that, by part 2 of Deﬁnition 2.1, we
have z100 = z10, and so on.
insist on having a bond between them. The phrase “connecting x and y”
in the deﬁnition of Hn(x,y) is not to imply that Hn(x,y) is an occupied
path from x to y. Instead, Hn(x,y) should be thought of as a framework of
large-scale bonds which can be turned into a path by connecting the “gaps”
in an appropriate way; see Figure 1.
Our strategy in both the upper and lower bound will be to identify a
hierarchy of suﬃcient depth from within a given path. In an idealized situ-
ation, this hierarchy between sites at Euclidean distance N would be such
that the primary bond (z01,z10) has length (approximately) N s/(2d), the
secondary bonds have length N (s/(2d))2, and so on. The principal diﬃculty
is to “make” the hierarchy deep enough so that it already contains “most
of” the bonds in the underlying path. In particular, we will have to guaran-
tee that the “gaps”—which may still be rather spread out in the Euclidean
distance—can be spanned without substantially aﬀecting the overall length.
2.1. Upper bound.
To outline the proof of the upper bound on the graph
distance, it is convenient to start by analyzing the cases when all pairs of
nearest neighbors on Zd are a priori connected. In these situations one can
(essentially) construct a path connecting two distant sites which uses about
the optimal number of distinct occupied bonds.
Let γ ∈(s/(2d),1). The construction is based on the following observa-
tion: If x and y are two sites at distance |x −y| = N ≫1, and if B0, re-
spectively, B1, are boxes of side N γ centered at x, respectively, y, then B0
and B1 are with overwhelming probability connected by a bond in G . In-
deed, there are N dγ ways to place each endpoint of such a bond while its

<!-- page 9 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
9
Euclidean length has to be essentially equal to N. Hence, the probability
that B0 and B1 are not directly connected by an occupied bond is
P(B0 ↮B1) = exp
(
−
X
z∈B0
X
z′∈B1
q(z −z′)
)
= exp{−N 2dγ−s+o(1)}.
(2.2)
Since 2dγ > s, the right-hand side tends rapidly to zero as N →∞.
Once the bond between B0 and B1 has been selected, a similar argument
shows that the boxes B00 and B01 of side N γ2, centered at x and the “nearer”
endpoint of the primary bond, respectively, will typically be connected by
an occupied bond. Continuing the process for a while, the endpoints of
the family of bonds thus identiﬁed give rise to a hierarchy of sites in the
above sense: First we let z0 = x and z1 = y, then we let z01 and z10 be the
endpoints of the primary bond connecting B0 and B1. Next, the endpoints
of the secondary bonds connecting the boxes B00 and B01, respectively,
B10 and B11, will be denoted by z001 and z010, respectively, z101 and z110.
The higher levels will be denoted similarly. Note that, in order to have
each level of the hierarchy completely deﬁned, we need to use part 2 of
Deﬁnition 2.1 to identify z00 with z0, and so on.
Of course, the most pertinent question is now for how long we can continue
developing such a hierarchy. Proceeding as in our previous estimates, the
probability that not all pairs of boxes Bσ0 and Bσ1 with σ ∈{0,1}k will be
connected by a bond in G is bounded by
2k exp{−N γk(2dγ−s+o(1))},
(2.3)
where 2k counts the number of bonds we are trying to control at this step
and the factor N γk(2dγ−s) in the exponent originates from the fact that we
are connecting boxes of side N γk+1 which are at Euclidean distance N γk from
each other. This estimate shows that, as long as N γk(2dγ−s) ≫log log N, the
probability that the identiﬁcation procedure fails is negligible. However, this
allows us to reach the level when the pairs of sites constituting the “gaps” are
no farther than N γk = (log N)o(1) from each other. This happens for k ≈K,
where
K = log log N
log(1/γ) .
(2.4)
Now, a hierarchy of depth K consists of roughly 2K bonds and 2K “gaps.”
Using nearest-neighbor paths to span each “gap,” the total number of bonds
needed to connect all “gaps” will thus be at most 2K(log N)o(1). Hence,
the graph distance between x and y cannot exceed 2K(log N)o(1). Plugging
the value of K and passing to the limit γ ↓s/(2d), the latter is no more
than (log N)∆+o(1).

<!-- page 10 -->

10
M. BISKUP
Performing the above argument without the luxury of an a priori con-
nected graph involves quite some extra work. Indeed, we need to ensure
that the sites identiﬁed in the process are connected to x and y (and, there-
fore, to C∞) and that the bonds lie in a “backbone”—rather than a “dead-
end”—of the connection between x and y. Our solution to this nonlocal
optimization problem is to construct the hierarchy so that each site zσ for
σ ∈{0,1}k is connected to a positive fraction of all sites in the N γk neigh-
borhood of zσ. Since the distance between the endpoints of the “gaps” in
such a hierarchy is at most of the order N γk−1, the connected components
of these endpoints are still with a large probability connected by a bond
from G . Now, if k ≈K, we have N dγk = (log N)o(1) and we need no more
than (log N)o(1) steps to connect the endpoints of each “gap.” This allows
us to proceed as before.
To ensure the connectivity property, we will introduce the concept of a
dense site which is a site x that is connected to at least a (prescribed) fraction
of all sites in a suﬃciently large box centered at x. Then we need to establish
two additional facts: First, any site x ∈C∞is with overwhelming probabil-
ity dense. Second, any suﬃciently large box contains a positive fraction of
dense sites. These statements—which come as Corollaries 3.3 and 3.4—will
allow us to look for hierarchies containing only dense sites, for which the
above argument easily carries through. The proof of the two corollaries in
turn requires showing that the largest connected component in any box
contains a (uniformly) positive fraction of all sites. To maintain generality,
this statement—which comes as Theorem 3.2—has to be proved under very
modest assumptions; essentially, we only assume the asymptotic (1.3) and
the fact that there is percolation.
2.2. Lower bound.
The argument for the upper bound shows that there
exists a path that connects x to y in about (log |x −y|)∆steps. The goal
of the lower bound is to show that, among the multitude of paths possibly
connecting x and y, no path will be substantially shorter.
In an idealized situation, our reasoning would go as follows: We set N =
|x−y| and pick a path from Π(x,y) that connects x with y in less than (log N)O(1)
steps. [Here O(1) represents a ﬁxed number whose value is irrelevant in
the following.] Next we will attempt to identify a hierarchy from π. The
primary bond (z01,z10) is chosen simply as (one of ) the longest bonds
in π. Since |x −y| = N but |π| = (log N)O(1), this bond must be longer
than N/(log N)O(1). But in order for this bond to exist with a reasonable
probability, a similar argument as used in the upper bound shows that the
distances N0 = |x −z01| and N1 = |z10 −y| must be such that
N d
0 N d
1 ≳N s+o(1).
(2.5)

<!-- page 11 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
11
Supposing (without any good reason) that N0 is comparable with N1, the
removal of (z01,z10) from π would leave us with two paths that connect sites
at distance N s/(2d)+o(1) in a polylogarithmic number of steps. The argument
could then be iterated which would eventually allow us to categorize the
whole path into a hierarchical structure, with one bond of length N, two
bonds of length N s/(2d), four bonds of length N (s/(2d))2, and so on.
It is easy to check that the hierarchy thus identiﬁed would involve roughly 2K
bonds, where K is as in (2.4) with γ = s/(2d), and |π| would thus have to be
at least (log N)∆+o(1). Of course, the main problem with the above argument
is that the assumption N0 ≈N1 is not justiﬁed and presumably fails in a large
number of places. Extreme ways to violate the condition N0 ≈N1 are not so
hard to dismiss. For instance, in the case of a “gap” collapse, for example,
when N0 = N o(1), the bound (2.5) forces that N1 ≳N s/d+o(1) ≫N, implying
that (z01,z10) was not the longest bond after all. But, since we are dealing
with an exponentially growing number of bonds, even “soft” violations of
this condition could make the whole argument crumble to pieces. As we will
describe below, the solution is to work with (2.5)—and its generalizations—
the way it stands without trying to extract any information about the par-
ticular N0 and N1.
Here is what we do. We pick a γ satisfying γ < s/(2d) and show that every
path connecting x and y in less than (log N)O(1) steps, where N = |x −y|,
contains a hierarchy Hn(x,y) of depth n ≲K—with K as in (2.4)—such
that the following holds with overwhelming probability: The length of the
“gaps” is comparable with the length of the bonds that “ﬁll” them; that is,
for all k = 1,...,n −1, we have
|zσ01 −zσ10| = |zσ0 −zσ1|1−o(1),
σ ∈{0,1}k−1.
(2.6)
Moreover, the average size of the “gaps” on the kth level is about N γk; that
is, the quantities Nσ = |zσ0 −zσ1| satisfy
Y
σ∈{0,1}k
Nσ ≥N (2γ)k
(2.7)
for all k = 1,...,n −1. Obviously, for k = 1 this is a more precise form
of (2.5). Part (2.6) is a consequence of the fact that, in order to connect two
sites at distance N γk in less than (log N)O(1) steps, at least one bond in the
path must be longer than N γk/(log N)O(1). This equals N γk(1−o(1)) as long
as k ≲K. As to the proof of (2.7), let En be the event that the inequality in
(2.7) holds for k = 1,2...,n −1. We will sketch the derivation of an upper
bound on P(Ec
k+1 ∩Ek) which can then be iterated into a bound on P(Ec
n).
Fix a collection of numbers (Nσ) representing the distances between var-
ious “gaps” in the hierarchy, and let us estimate the probability that a
hierarchy with these (Nσ) occurs. In light of (2.6), the primary bond will

<!-- page 12 -->

12
M. BISKUP
cost N −s+o(1) of probability, but there are of the order (N0N1)d−1 ways to
choose the endpoints. (Remember that N0 and N1 are ﬁxed.) Similarly, the
two secondary bonds cost (N0N1)−s+o(1) of probability and their endpoints
contribute of the order (N00N01N10N11)d−1 of entropy. Applying this to the
collections (Nσ) compatible with Ec
n+1 ∩En, we get
P(Ec
n+1 ∩En)
(2.8)
≤
X
(Nσ)
(N0N1)d−1
N s−o(1)
(N00N01N10N11)d−1
(N0N1)s−o(1)
...
Q
σ∈{0,1}n N d−1
σ
Q
σ∈{0,1}n−1 N s−o(1)
σ
,
where the sum goes over all (Nσ) for which Ec
n+1 ∩En holds.
To evaluate the right-hand side of (2.8), we need to observe that the nu-
merator (N0N1)d−1 can be combined with the denominator of the next quo-
tient into a term which by s > d is summable on both N0 and N1; using (2.7),
the resulting sum over N0 and N1 is bounded by N −(s−d)(2γ)+o(1). The other
numerators will be handled analogously; the upshot is that, for all k ≤n−1,
the sum over all Nσ with σ ∈{0,1}k is bounded by N −(s−d)(2γ)k+o(1). The
last numerator has no denominator to be matched with, but here we can
use that, since (2.7) for k = n fails on Ec
n+1 ∩En, the product of Nσ’s is
now bounded from above! Consequently, the relevant sum does not ex-
ceed N d(2γ)n+o(1). Putting all these estimates together, and applying the
inequality
s −d(2γ)n + (s −d)
n−1
X
k=1
(2γ)k ≥(s −2dγ)(2γ)n−1,
(2.9)
the right-hand side of (2.8) is bounded by N −(s−2dγ)(2γ)n−1. [The inequal-
ity in (2.9) can be derived either by direct computation or by a repeated
application of the inequality s + (s −d)(2γ) ≥s(2γ) to the ﬁrst and the
third term on the left-hand side.] Summing the obtained bound over n, the
probability P(En) is shown to be essentially one as long as n ≲K.
Once we have established that (2.6) and (2.7) hold, we will use a similar
estimate to ﬁnd a lower bound on D(x,y). Here we simply have to prove
that, even though the hierarchy is already rather large, the lower bound (2.7)
requires that at least as many bonds be used to connect all of the “gaps.”
To avoid some unpleasant combinatorial estimates, we will continue under
the simplifying assumption that all of the 2n−1 “gaps” of the hierarchy are
nontrivial.
Let Fn be the event that every hierarchy of depth n satisfying (2.6) and (2.7)
requires more than 2n−1 extra steps to connect all of its “gaps.” In light of
our bound on P(Ec
n), it suﬃces to estimate the probability of Fc
n ∩En. Since
all “gaps” are nontrivial, the only way Fc
n can occur is that each “gap”

<!-- page 13 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
13
is spanned by a single bond. Now the bond spanning the “gap” (zσ0,zσ1)
costs N −s+o(1)
σ
amount of probability and so P(Fc
n ∩En) can be bounded by
P(Fc
n ∩En)
(2.10)
≤
X
(Nσ)
(N0N1)d−1
N s−o(1)
...
Q
σ∈{0,1}n−1 N d−1
σ
Q
σ∈{0,1}n−2 N s−o(1)
σ
Y
σ∈{0,1}n−1
1
N s−o(1)
σ
,
where the (Nσ)’s now obey (2.7) for all k = 0,1,...,n −1. The last prod-
uct on the right-hand side makes the entire sum convergent and (for N ≫
1) small. Thus, with overwhelming probability, Fn occurs for all n ≲K,
which means that the shortest path(s) between x and y must contain at
least 2n−1 = 2K(1−o(1)) distinct bonds. For K as in (2.4), we have 2K =
(log N)∆′, where 1/∆′ = log2(1/γ). From here the lower bound in (1.6) fol-
lows by letting γ ↑s/(2d).
3. Percolation in ﬁnite boxes.
3.1. Size of the largest connected component.
In this section we will prove
an estimate showing that the largest connected component in large but ﬁnite
boxes contains a positive fraction of all sites whenever there is percolation.
This estimate will be essential for the proof of the upper bound in (1.6).
Throughout this section, the original meaning of the quantity from (1.3)
will be substituted by a weaker form (3.1) below. We will return to the
original deﬁnition in Section 4.
We begin by quoting a result from [8]. Let us say that the collection of
probabilities (pxy)x,y∈Zd is percolating, if the associated i.i.d. measure has
an inﬁnite cluster almost surely. Let L ≥1 be an integer and let ΛL be
a box in Zd of side L containing Ld sites. Consider the percolation prob-
lem restricted to the sites of ΛL (and, of course, only the bonds with both
endpoints in ΛL) and let |CL| denote the size of the largest connected com-
ponent in ΛL. In [8], Berger proved that once (pxy)x,y∈Zd are percolating, ΛL
contains a large cluster. The precise formulation is as follows:
Theorem 3.1 ([8], Lemma 2.3).
Let d ≥1 and suppose that the col-
lection of probabilities (pxy)x,y∈Zd, where pxy = 1 −e−q(x−y), is percolating.
Suppose that, for some s ∈(d,2d),
liminf
|x|→∞|x|s(1 −e−q(x)) > 0.
(3.1)
Then for each ε > 0 and each ζ ∈(0,∞), there exists an L such that
P(|CL| < ζLs/2) ≤ε.
(3.2)

<!-- page 14 -->

14
M. BISKUP
We note that once (3.1) holds for some s, then it holds also for any s′ >
s. Therefore, Theorem 3.1 actually guarantees that the largest connected
component in ΛL will contain at least Ld−o(1) sites. (Note that the statement
forbids us to take s = 2d, and an inspection of Lemma 2.3 in [8] reveals that
this is nontrivially rooted in the proof.) However, for our purposes we need
to work with the event that |CL| is proportional to Ld and, in addition, we
also need a more explicit estimate on the probability of such an event. Our
extension of Berger’s result comes in the following theorem:
Theorem 3.2.
Let d ≥1 and consider the probabilities (pxy)x,y∈Zd such
that (3.1) holds for some s ∈(d,2d). Suppose that (pxy)x,y∈Zd are percolating.
For each s′ ∈(s,2d) there exist numbers ρ > 0 and L0 < ∞such that for
each L ≥L0,
P(|CL| < ρ|ΛL|) ≤e−ρL2d−s′
.
(3.3)
In particular, once L is suﬃciently large, the largest connected component
in ΛL typically contains a positive fraction of all sites in ΛL.
Theorem 3.2 alone would allow us to establish the existence of a hierar-
chy between two sites, but it would not ensure that the “gaps” are prop-
erly connected (which is what we need to turn the hierarchy into a path).
Fortunately, the structure of the proof of Theorem 3.2 allows us to make
this conclusion anyway. To state the relevant mathematical claims, for each
x ∈Zd and any odd integer L ≥1, let ΛL(x) be the box of side L centered
at x and let CL(x) be the set of sites in ΛL(x) that are connected to x by
an occupied path in ΛL(x). Then we have:
Corollary 3.3.
Under the conditions of Theorem 3.2, there exists a
constant ρ > 0 such that
lim
L→∞P(|CL(x)| < ρ|ΛL(x)|, x ∈C∞) = 0
(3.4)
holds for each x ∈Zd.
Corollary 3.4.
Given ℓ< L, let D(ρ,ℓ)
L
be the set of sites x ∈ΛL such
that |Cℓ(x)| ≥ρ|Λℓ(x)|. Under the conditions of Theorem 3.2, for each s′ ∈
(s,2d) there are constants ℓ0 < ∞and ρ > 0 such that
P(|D(ρ,ℓ)
L
| < ρ|ΛL|) ≤e−ρL2d−s′
(3.5)
holds for any ℓwith ℓ0 ≤ℓ≤L/ℓ0.

<!-- page 15 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
15
Theorem 3.2 and the two corollaries are what this section contributes
to the proof of the main result of this paper. An impatient (or otherwise
uninterested) reader should feel free to skip the rest of this section on a ﬁrst
reading and pass directly to Sections 4.1 and 4.2. For those who stay put,
we proceed by discussing the main ideas of the proof and a breakdown of its
steps into the corresponding technical lemmas. The actual proofs appear in
Section 3.5.
3.2. Outline of the proof.
Our strategy of the proof of Theorem 3.2 is
as follows. First we combine a one-step renormalization with Theorem 3.1
to convert the problem into a similar question for site-bond percolation. An
important feature of this reformulation is that the occupation probability of
sites and bonds can be made as close to 1 as we wish.
Given an odd integer K ≥1, let C (x)
K
denote the largest connected compo-
nent in ΛK(x); in the case of a tie we will choose the component containing
the site that is minimal in the standard lexicographic order on Zd. For any
two distinct x,y ∈Zd, we will say that ΛK(Kx) and ΛK(Ky) are directly
connected if there is an occupied bond connecting a site from C (Kx)
K
to a
site from C (Ky)
K
. We will use {ΛK(Kx) ↮ΛK(Ky)} to denote the event
that ΛK(Kx) and ΛK(Ky) are not directly connected. Then we have:
Lemma 3.5.
Under the assumptions of Theorem 3.1, for any s ∈(d,2d)
the following is true: For each β < ∞and r < 1 there exist a number δ > 0
and an odd integer K < ∞, such that
P(|C (Kx)
K
| < δ |ΛK(Kx)|) ≤1 −r
(3.6)
and
P(ΛK(Kx) ↮ΛK(Ky)) ≤e−β|x−y|−s
(3.7)
hold for all distinct x,y ∈Zd.
Regarding boxes of side K as new sites and the pairs of maximal-connected
components connected by a bond from G as new bonds, Lemma 3.5 allows
us to set up a renormalization scheme of [22]. Clearly, by (3.6) and (3.7), a
site is occupied with probability at least r and two occupied sites x and y
are connected with probability at least 1−e−β|x−y|−s. (For the sites that are
not occupied such a connection will not be relevant, so we will often assume
that the latter holds for all sites.) This puts us into a position where we can
apply the following “idealized” version of the desired claim:
Lemma 3.6.
Let d ≥1 and consider the site-bond percolation model
on Zd with sites being occupied with probability r ∈[0,1] and the bond be-
tween sites x and y being occupied with probability
pxy = 1 −exp{−β|x −y|−s},
(3.8)

<!-- page 16 -->

16
M. BISKUP
where s ∈(d,2d) and β ≥0. Let |CN| denote the size of the largest connected
component of occupied sites and occupied bonds in ΛN. For each s′ ∈(s,2d)
there exist numbers N0 < ∞, ϑ > 0 and β0 < ∞such that
Pβ,r(|CN| < ϑ|ΛN|) ≤e−ϑβN2d−s′
(3.9)
holds true for all N ≥N0 whenever β ≥β0 and r ≥1 −e−ϑβ.
Remark 3.1.
The fact that the exponent in (3.9) is proportional to β
will not be needed for the proofs in this paper. The addition of β represents
only a minor modiﬁcation of the proof and the stronger result will (hopefully)
facilitate later reference.
Once Lemma 3.6 is inserted into the game, the proof of Theorem 3.2
will be easily concluded. To prove Lemma 3.6, we will invoke a combina-
tion of coarse-graining, stochastic domination and a corresponding estimate
for the complete graph. (Let us recall that a complete graph of n vertices
is a graph containing a bond for each unordered pair of distinct numbers
from {1,2,... ,n}.) The relevant complete-graph statement is extracted into
the following lemma:
Lemma 3.7.
Consider a complete graph of n vertices. Let each site be
occupied with probability r and each bond be occupied with probability p.
Let Pp,r
n
be the resulting i.i.d. measure and let |Cn| denote the number of sites
in the largest connected component of occupied sites and occupied bonds. For
each q,q′ ∈[0,1] with q′ < q, there exists a number ψ(q′,q) > 0 such that for
each r′ ∈[0,r), each p′ ∈[0,p) and all n ≥1,
Pp,r
n (|Cn| ≤p′r′n) ≤e−nψ(r′,r) + e−(1/2)(n2(r′)2−n)ψ(p′,p).
(3.10)
Moreover, for each α ∈[0,1), there exists a constant C = C(α) < ∞such
that
ψ(q′,q) ≥(1 −α)(1 −q′)

log
1
1 −q −C

(3.11)
holds true for all q,q′ ∈[0,1] satisfying the bound (1 −q′) ≥(1 −q)α.
Remark 3.2.
While Lemma 3.6 can presumably be proved without ref-
erence to the complete graph, in our case the passage through Lemma 3.7
has the advantage of easily obtained quantitative estimates. As mentioned
before, the present proof of Theorem 3.2 invokes a renormalization scheme
for which Theorem 3.1—whose proof is based on a similar renormalization
scheme—serves as a starting point. It would be of some conceptual inter-
est to see whether a more direct proof of Theorem 3.2 based on a single
renomalization is possible.

<!-- page 17 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
17
Remark 3.3.
In d ≥2, the decay rate in (3.9)—and, consequently, in
(3.3)—is not always optimal. The reason is that, for β ≫1 and 1−r ≪1, the
site-bond percolation problem dominates the nearest-neighbor percolation
on Zd for which it is expected (and essentially proved, see [3, 13]) that the
probability in (3.9) should decay exponentially with N d−1. [For s ∈(d,2d),
this is sometimes better and sometimes worse than N 2d−s.] This alternative
decay rate is not reﬂected in our proofs because, to apply equally well in all
dimensions d ≥1, they consistently rely only on long-range connections.
Having outlined the main strategy and stated the principal lemmas, we
can plunge into the proofs. First we will prove technical Lemmas 3.5 and 3.7.
Then we will establish the site-bond percolation Lemma 3.6. Once all prepara-
tory statements have been dispensed with, we will assemble the ingredients
into the proofs of Theorem 3.2 and Corollaries 3.3 and 3.4.
3.3. Preparations.
Here we will prove Lemmas 3.5 and 3.7. First we will
attend to the one-step renormalization from Lemma 3.5 whose purpose is to
wash out the short-range irregularities of the pxy’s and to ensure that the
constants β and r in Lemma 3.6 can be chosen as large as required.
Proof of Lemma 3.5.
The principal task before us is to choose K so
large that both bounds follow from Theorem 3.1 and the assumption (3.1).
Let d ≥1 and let β < ∞and r < 1. Let ε = 1−r and pick a d′ ∈(s/2,d) such
that d′ −(d −1) > d −d′. By Theorem 3.1 and the paragraph that follows,
for each ζ > 0, there exists a K = K(ε,d′,ζ) such that |CK| ≥ζKd′ occurs
with probability exceeding 1 −ε. Since Theorem 3.1 allows us to make K
arbitrarily large (indeed, the constraint |CK| ≤Kd forces K to exceed a
positive power of ζ), we can also assume that, for some constant α > 0,
pxy ≥1 −e−α|x−y|−s
once |x −y| ≥Kd−d′.
(3.12)
Here we rely on (3.1).
Let b > 0 be a constant such that bK exceeds the diameter of ΛK in the
| · |-norm for all K. We will show that (3.6) and (3.7) hold once K is large
enough and, in particular, so large that
1
4
αζ2
(2b)s K2d′−s ≥β.
(3.13)
(Note that 2d′ > s so the left-hand side increases with K.) Consider a par-
titioning of Zd into disjoint boxes of side K; that is, let us write Zd =
S
x∈Zd ΛK(Kx). We will call a box ΛK(Kx) occupied if the bond conﬁgura-
tion restricted to ΛK(Kx) contains a connected component larger than ζKd′.
By the choice of K, each ΛK(Kx) is occupied independently with probability
exceeding 1 −ε = r. This proves (3.6) with δ = ζKd′−d.

<!-- page 18 -->

18
M. BISKUP
To prove also (3.7), we need to ensure that suﬃcient portions of the com-
ponents in the two K-blocks are so far from each other that (3.12) can
safely be applied. To this end, we note that, since |C (Kx)
K
| ≥ζKd′, at least
half of the sites in C (Kx)
K
will be farther from Zd \ ΛK(Kx) than ηK =
aζKd′−(d−1), where a is a constant depending only on the norm | · | and
the dimension. Moreover, an easy argument shows that if x,y ∈Zd are dis-
tinct and z ∈C (Kx)
K
and z′ ∈C (Ky)
K
are such that dist(z,Zd \ΛK(Kx)) ≥ηK
and dist(z′,Zd \ ΛK(Ky)) ≥ηK, then
ηK|x −y| ≤|z −z′| ≤2bK|x −y|.
(3.14)
Now our choice of d′ guarantees that, if K is suﬃciently large, we will have
ηK ≥Kd−d′ and the inequality on the left-hand side shows that (3.12) for pzz′
is in power. The bound on the right-hand side of (3.14) then allows us to
write
P(ΛK(Kx) ↮ΛK(Ky)) ≤exp

−α
(ζKd′/2)2
(|x −y|2bK)s

≤e−β|x−y|−s,
(3.15)
where we used (3.13) to derive the last inequality. This completes the proof.
□
Next we will focus on the proof of Lemma 3.7, which concerns the complete
graph. This lemma will be used to drive the induction argument in the next
section.
Proof of Lemma 3.7.
The proof starts by estimating the total number
of occupied sites. Once that number is known to be suﬃciently large, the
desired bound is a result of conditioning on occupied sites combined with
straightforward estimates concerning occupied bonds.
Fix r′ ∈(0,r) and p′ ∈[0,p), and let ρ = p′r′. Note that we can assume
that r′n > 1, because otherwise the right-hand side of (3.10) exceeds 1. Let
An denote the (random) number of occupied vertices of the complete graph.
Since An can be represented as a sum of independent random variables with
mean r, the event {An < r′n} forces An to deviate substantially from its
mean and the standard Chernoﬀbound implies
Pp,r
n (An < r′n) ≤e−ψ(r′,r)n,
(3.16)
where
ψ(q′,q) = sup
λ≥0
[−log(1 −q + qe−λ) −λq′].
(3.17)
As is easy to show, ψ(q′,q) > 0 for all q′ < q.
The bound (3.16) is responsible for the ﬁrst term on the right-hand side
of (3.10). It remains to show that the conditional probability given An ≥r′n

<!-- page 19 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
19
is bounded by the second term. Thus suppose that An ≥r′n and let Vn
denote the total number of unordered pairs of occupied sites not connected
by an occupied bond. The principal observation is that if |Cn| ≤ρn, then Vn
has to be rather large. More precisely, we claim that |Cn| ≤ρn implies Vn ≥
1
2(A2
n −ρnAn). Indeed, if we label the connected components of occupied
sites and bonds by index i and use ki to denote the number of sites in
the ith occupied component, then the number of vacant bonds certainly
exceeds
X
i<j
kikj = 1
2
 X
i
ki
!2
−1
2
X
i
k2
i .
(3.18)
On {|Cn| ≤ρn} we have ki ≤ρn for each i and since also P
i ki = An, the
second sum can be bounded by Anρn. The desired inequality Vn ≥1
2(A2
n −
ρnAn) follows.
In light of our previous reasoning, we are down to estimating the proba-
bility
Pp,r
n (Vn ≥1
2(A2
n −ρnAn)|An ≥r′n).
(3.19)
The estimate will be performed by conditioning on the set of occupied sites.
Once the set of occupied sites has been ﬁxed, Vn can be represented as
a sum of N = 1
2(A2
n −An) independent random variables, each of which
has mean 1 −p. Now, assuming An ≥r′n and recalling that ρ = p′r′ ≤r′
and r′n > 1, we can estimate
1
2(A2
n −ρnAn) = N

1 −ρn −1
An −1

≥N

1 −ρn −1
r′n −1

≥N(1 −p′).
(3.20)
(In the cases when ρn < 1, we just skip the intermediate inequality.) Now,
since p′ < p, the event {Vn ≥1
2(A2
n −ρnAn)} constitutes a large deviation
for the random variable Vn. Invoking again the Chernoﬀbound and a little
algebra, we thus get
Pp,r
n (Vn ≥1
2(A2
n −ρnAn)|An ≥r′n, 1
2(A2
n −An) = N) ≤e−Nψ(p′,p).
(3.21)
From here (3.10) follows by noting that on {An ≥r′n} we have N ≥1
2(n2(r′)2−n).
To verify the bound (3.11), we need to ﬁnd the minimizing λ in (3.17)
and use it to ﬁnd an explicit expression for ψ(q′,q). A computation gives
ψ(q′,q) = (1 −q′)log
1 −q′
1 −q

−q′ log
 q
q′

.
(3.22)
The second term—including the minus sign—can be split into two parts: the
term q′ log q′, which is bounded below by −1/e, and the term −q′ log q, which
is always positive. Moreover, if q′ →1, then (q′ log q′)/(1−q′) →1. From here
we infer that the second term is bounded below by (1 −α)(1 −q′)-times a
(negative) α-dependent constant. Using the bound 1 −q′ ≥(1 −q)α in the
ﬁrst term, (3.11) is proved.
□

<!-- page 20 -->

20
M. BISKUP
3.4. Site-bond percolation lemma.
Now we are ready to start proving
Lemma 3.6. The essential part of the proof relies on induction along a series
of scales deﬁned as follows. Fix an s′ ∈(s,2d), let ℓ0 be a positive integer
and consider an increasing sequence (ℓn) of integers such that ℓn ≥ℓ0 for
all n ≥1. Let N0 be another positive integer and let
Nn = N0
n
Y
k=1
ℓk.
(3.23)
Suppose that ℓn tend to inﬁnity so fast that
X
n≥1
ℓd−s
n
< ∞
and
X
n≥1
N s−s′
n
< ∞
(3.24)
but moderate enough that also
c0 = inf
n≥1
1
ℓs′
n+1
n
Y
k=1
ℓ2d−s′
k
> 0.
(3.25)
Next, let us deﬁne the sequences (rn) and (pn) by putting rn = 1 −6ℓd−s
n
and pn = 1 −N s−s′
n
, let us ﬁx a constant ρ0 ∈(0,1) and let
ρn = ρ0
n
Y
k=1
(rkpk).
(3.26)
Clearly, by (3.24) we have ρn →ρ∞> 0 as n →∞.
Remark 3.4.
An extreme example of a sequence (ℓn) satisfying these
constraints is
ℓn = e(1+a)n−1,
n ≥1,
(3.27)
where a = (2d −s′)/s′. Here a has been tuned in such a way that the term
in the inﬁmum in (3.25) is independent of n. Since a > 0, the bounds (3.24)
immediately follow.
The proof of Lemma 3.6 is based on the fact that a bound (3.9) for N = Nn
can be used to prove the same bound for (essentially) any N between Nn
and Nn+1. The proof works as soon as c0 and ρ∞are bounded away from zero
and ℓ0 and N0 are suﬃciently large; the precise form of (ℓn) is not important.
[In particular, we do not make explicit use of (3.24).] The induction step is
isolated into the following claim:
Lemma 3.8.
Suppose the assumptions of Lemma 3.6 hold and let s′ ∈
(s,2d), c′
0 > 0 and τ0 > 0. Then there exist two numbers ℓ′
0 < ∞and N ′
0 < ∞
and a constant c1 ∈(0,∞) such that for any N0 ≥N ′
0 and any sequence (ℓn)

<!-- page 21 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
21
with ℓ0 ≥ℓ′
0 and c0 ≥c′
0, the following holds: If τ ∈[τ0,c1ρ2
∞β] and if k is a
nonnegative integer such that
Pβ,r(|CNk| ≤ρkN d
k) ≤e−τN2d−s′
k
(3.28)
is true, then
Pβ,r(|CℓNk| ≤ρk+1(ℓNk)d) ≤e−τ(ℓNk)2d−s′
(3.29)
holds for all ℓ∈{ℓ0,...,ℓk+1}.
Proof.
The proof is based on the “complete-graph” Lemma 3.7. Our
preliminary task is to set up all the constants so that the bounds emerging
from this lemma are later easily converted to that on the right-hand side of
(3.29). Let s′ ∈(s,2d), c′
0 > 0 and τ0 > 0 be ﬁxed. First we will address the
choice of the constants ℓ′
0, N ′
0 and c1. We will assume that ℓ′
0 is so large that
(1 −6ℓd−s)2ℓ2d −ℓd ≥1
2ℓ2d
(3.30)
holds for all ℓ≥ℓ′
0. Then we choose N ′
0 so large that, for all N0 ≥N ′
0 and
all τ ≥τ0,
6ℓd−s > e−(1/2)τ(ℓN0)2d−s′
(3.31)
holds for all ℓ≥ℓ′
0 and that
τN 2d−s′
0
≥max{3C,log 2},
(3.32)
where C is the constant from Lemma 3.7 for α = 1/2. Moreover, we let b be
the constant such that Nb bounds the diameter of ΛN in the metric | · | for
all N ≥1 and choose c1 = b−s/16. Then we also require that N ′
0 be so large
that for any N0 ≥N ′
0,
βb−sρ2
nN 2d−s
n
ℓ−s −C ≥8τN 2d−s
n
ℓ−s
(3.33)
and
(ℓNn)s−s′ > e−(1/2)βb−sρ2
nN2d−s
n
ℓ−s
(3.34)
hold for all τ ≤c1ρ2
∞β, all n ≥0 and all ℓ∈{ℓ0,...,ℓn+1}. The ﬁrst bound
is veriﬁed for suﬃciently large N0 by noting the inequalities βb−sρ2
n ≥2τ0
and
N 2d−s
n
ℓ−s ≥N 2d−s
n
ℓ−s
n+1 ≥N 2d−s
0
 
ℓ−s′
n+1
Y
k≤n
ℓ2d−s′
k
!s/s′
≥N 2d−s
0
cs/s′
0
.
(3.35)
As to (3.34), we note that N 2d−s
n
ℓ−s/(ℓNn)s′−s ≥N 2d−s′
n
ℓ−s′
n+1, which again
can be made arbitrarily large by boosting up N0. The factors
1
2 in the
exponents of (3.31) and (3.34) have been inserted for later convenience.

<!-- page 22 -->

22
M. BISKUP
Now we are ready to prove (3.29) from (3.28). Let ℓ0 ≥ℓ′
0 and N0 ≥N ′
0.
Suppose that (3.28) holds for some τ ∈[τ0,c1ρ2
∞β] and some k ≥0. Pick an
ℓ∈{ℓ0,...,ℓk+1}. Viewing ΛℓNk as the disjoint union of ℓd translates Λ(i)
of ΛNk,
ΛℓNk =
ℓd
[
i=1
Λ(i),
(3.36)
let us call Λ(i) occupied if it contains a connected component of size at
least ρk|Nk|d. Choosing one such connected component in each occupied Λ(i)
(if necessary, employing lexicographic order on Zd to resolve potential ties),
we will call Λ(i) and Λ(j) with i̸ = j connected if their respective connected
components are directly connected by an occupied bond. Let
r = 1 −e−τN2d−s′
k
and
p = 1 −e−βb−sρ2
kN2d−s
k
ℓ−s
(3.37)
and note that (3.28) implies that r is a lower bound on the probability
that Λ(i) is occupied. Moreover, a simple calculation shows that p is a uni-
form lower bound on the probability that two distinct Λ(i) and Λ(j) are
connected.
Let us introduce the quantities
r′ = 1 −6ℓd−s
and
p′ = 1 −(ℓNk)s−s′
(3.38)
and let Ak,ℓbe the event that the largest connected component of the oc-
cupied Λ(i)’s comprises more than p′r′ℓd individual boxes. The assumption
(3.31) shows that r′ < r, while (3.34) in turn guarantees that p′ < p. Invoking
the fact that Ak,ℓis an increasing event, the probability of Ac
k,ℓis bounded by
the probability that, for site-bond percolation on a complete graph with ℓd
vertices and parameters r and p, the largest connected component involves
at most p′r′ℓd vertices:
Pβ,r(Ac
k,ℓ) ≤Pp,r
ℓd (|Cℓd| ≤p′r′ℓd).
(3.39)
Since the factors 1
2 in the exponents in (3.31) and (3.34) ensure that 1−r′ ≥
(1−r)1/2 and 1−p′ ≥(1−p)1/2, the right-hand side can be estimated using
Lemma 3.6 with α = 1
2 and C = C(1
2). To evaluate the ﬁrst term on the
right-hand side of (3.10), we estimate
ℓdψ(r′,r) ≥1
26ℓ2d−s(τN 2d−s′
k
−C) ≥2τ(ℓNk)2d−s′.
(3.40)
Here we used that, by our choice of α, we have ℓd(1 −α)(1 −r′) = 1
26ℓ2d−s
and then invoked (3.32) to show that τN 2d−s′
k
≥τN 2d−s′
0
≥3C. Similarly,
we get
1
2(ℓ2d(r′)2 −ℓd)ψ(p′,p) ≥1
4ℓ2d(ℓNk)s−s′(βb−sρ2
kN 2d−s
k
ℓ−s −C)
(3.41)
≥2τ(ℓNk)2d−s′

<!-- page 23 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
23
for the exponent in the second term in (3.10). Here we ﬁrst used (3.30)
to reduce the complicated ℓ-dependence on the extreme left and then we
inserted (3.11) and the deﬁnitions of p′ and p to produce the intermediate
inequality. Finally, we invoked (3.33).
By putting the bounds (3.10), (3.11) and (3.40), (3.41) together and recall-
ing that log 2 ≤τN 2d−s′
0
, the probability Pβ,r(Ac
k,ℓ) does not exceed the term
on the right-hand side of (3.29). But on Ak,ℓ, the box ΛℓNk contains a con-
nected component comprising (strictly) more than p′r′ℓd disjoint connected
components, each of which involves at least ρkN d
k sites. Using that r′ ≥rk+1
and p′ ≥pk+1, we have
|CℓNk| > p′r′ℓdρk|Nk|d ≥ρk+1(ℓNk)d
on Ak,ℓ,
(3.42)
and thus {|CℓNk| ≤ρk+1|ΛℓNk|} ⊂Ac
k,ℓ. From here (3.28) follows.
□
Lemma 3.8 encapsulates the induction step. However, we will also need
an estimate that allows us to start the induction. This is provided in the
following lemma.
Lemma 3.9.
Under the conditions of Lemma 3.6, for each c2 ∈(0,∞),
there exist numbers N0 < ∞and ϑ0 ∈(0,1) and, for each ϑ < ϑ0, there exists
a number β0 < ∞such that
Pβ,r(|CN| < 1
4|ΛN|) ≤e−c2ϑβN2d−s
(3.43)
holds once N ≥N0, β ≥β0 and r ≥1 −e−ϑβ.
Proof.
We will again apply the “complete-graph” Lemma 3.7. Let p′ =
r′ = 1
2, let b be a constant such that bN exceeds the diameter of ΛN for
any N and pick ϑ0 > 0 such that 128c2ϑ0 < b−s. Fix a number ϑ ∈(0,ϑ0).
Then the left-hand side of (3.43) is bounded by the left-hand side of (3.10)
with n = N d, p = 1 −exp{−β(bN)−s} and r = 1 −e−ϑβ. We will estimate
the right-hand side of (3.10) under the conditions when N0 is so large that
N s−d ≥16c2
and
N 2d(r′)2 −N ≥N 2d/8
(3.44)
are true for all N ≥N0 and, given such an N, the constant β0 is so large
that for all β ≥β0, we have 1 −r′ ≥(1 −r)1/2 and 1 −p′ ≥(1 −p)1/2, and
ϑβ ≥2C
and
β(bN)−s −C ≥128ϑβN −s.
(3.45)
Here, as before, C is the constant from Lemma 3.7 for α = 1
2.
In conjunction with these bounds, (3.11) with α = 1
2 shows that
ψ(r′,r) ≥1
4(ϑβ −C) ≥1
8ϑβ
(3.46)

<!-- page 24 -->

24
M. BISKUP
and
ψ(p′,p) ≥1
4(β(bN)−s −C) ≥32c2ϑβN −s.
(3.47)
Using the bounds in (3.44), we ﬁnd that both exponents on the right-
hand side of (3.10) exceed 2c2ϑβN 2d−s. This implies Pβ,r(|CN| ≤ϑ|ΛN|) ≤
2exp{−2c2ϑ×βN 2d−s}. Increasing β0 if necessary, the latter is no more than
exp{−c2ϑβN 2d−s}.
□
Equipped with the induction machinery from Lemmas 3.8 and 3.9, the
proof of the main site-bond percolation lemma is now easily concluded.
Proof of Lemma 3.6.
First we will adjust the parameters so that
Lemmas 3.8 and 3.9 can directly be applied. Let s′ ∈(s,2d) and let c′
0 > 0 and
τ0 > 0 be ﬁxed. Let ℓ′
0, N ′
0 and c1 be the constants from Lemma 3.8 and pick a
sequence (ℓ0) such that ℓ0 ≥ℓ′
0 and c0 ≥c′
0 are satisﬁed. Pick a number N0 ≥
N ′
0 so large that Lemma 3.9 holds for c2 = (2ℓ0)2d−s′ and N ≥N0, and
let ϑ0 be the corresponding constant from this lemma. Let ρ0 ∈(0, 1
4] and
deﬁne ρn and ρ∞as in (3.26). Let ϑ > 0 be such that ϑ ≤ϑ0, (2ℓ0)dϑ ≤ρ∞
and c2ϑ ≤c1ρ2
∞. Choose β0 so large that Lemma 3.9 holds for all β ≥β0 and
such that c2ϑβ0 > τ0. Note that τ = c2ϑβ necessarily satisﬁes τ ∈[τ0,c1ρ2
∞β]
as long as β ≥β0, which is needed in Lemma 3.8.
Now we are ready to run the induction argument: Since ρ0 ≤1
4 and ϑ ≤ϑ0,
Lemma 3.9 ensures that (3.28) holds for τ = c2ϑβ and k = 0. Applying the
induction step from Lemma 3.8, we recursively show that (3.29) is true for
all k ≥0 and all ℓ∈{ℓ0,...,ℓk+1}. Let N be a general integer and let k
be a nonnegative integer such that Nk+1 > N ≥Nk. Let ℓbe the maximal
number in {1,ℓ0,...,ℓk+1} such that ℓNk ≤N. A simple calculation now
shows that |ΛN| ≤max(2ℓ0)d|ΛℓNk| and, if we position ΛN and ΛℓNk so that
ΛℓNk ⊂ΛN,
{|CN| ≤ϑ|ΛN|} ⊂{|CℓNk| ≤(2ℓ0)dϑ|ΛℓNk|}.
(3.48)
By our previous conclusions and the fact that (2ℓ0)dϑ ≤ρk+1, the probability
of the event on the right-hand side is bounded by exp{−τ(ℓNk)2d−s′}. From
here, (3.9) for a general N follows by noting that, by our choice of c2, we
have τ(ℓNk)2d−s′ ≥ϑβN 2d−s′.
□
3.5. Proofs of Theorem 3.2 and Corollaries 3.3 and 3.4.
Now we are
ﬁnally ready to prove Theorem 3.2. After some preliminary arguments,
the proof follows a line of reasoning similar to the one just used to prove
Lemma 3.6.

<!-- page 25 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
25
Proof of Theorem 3.2.
Let s′ ∈(s,2d) and let ϑ, β0 and N0 be as
in Lemma 3.6. Pick numbers β ≥β0 and r ≥1 −e−ϑβ and let K and δ be
the corresponding constants from Lemma 3.5. First we will prove the claim
for L of the form L = KN, where N is a positive integer. To that end, let us
view ΛL as the disjoint union of N d translates ΛK(Kx) of ΛK, where x ∈
Zd. We will again call ΛK(Kx) occupied if |C (Kx)
K
| ≥δ|ΛK| and, similarly,
we will call two distinct ΛK(Kx) and ΛK(Ky) connected if the connected
components C (Kx)
K
and C (Ky)
K
—chosen with the help of lexicographic order
in case of a tie—are directly connected by a bond from G .
By Lemma 3.5 and our choice of δ and K, the box ΛK(x) is occupied
with probability exceeding r, while ΛK(Kx) and ΛK(Ky) are connected
with probability exceeding pxy in (3.8). Let AN,K be the event that the
box ΛL contains a connected component CN of boxes ΛK(Kx) such that at
least ϑ|ΛN| of the connected components in these boxes get joined in CN.
By Lemma 3.6, we know that
Pβ,r(Ac
N,K) ≤e−ϑβN2d−s′
.
(3.49)
On the other hand, on AN,K we have
|CL| ≥(ϑ|ΛN|)(δ|ΛK|) = ϑδ|ΛL|,
(3.50)
and thus {|CL| ≤ρ|ΛL|} ⊂Ac
N,K once ρ < ϑδ. If ρ is also less than ϑβKs′−2d,
this ﬁnishes the proof for L of the form NK. The general values of L are
handled by noting that if NK ≤L < (N + 1)K, then |ΛL| ≤2d|ΛNK| and,
if ΛNK ⊂ΛL, then also |CNK| ≤|CL|.
□
Proof of Corollary 3.3.
Fix s′ ∈(s,2d), let N0, ϑ and β0 be the
constants from Lemma 3.6, and let β ≥β0 and r ≥1 −e−ϑβ. Let ℓ> 3N0 be
an odd integer and let K be the constant from Lemma 3.5 for our choice
of β and r. Clearly, it suﬃces to show that (3.4) holds for L of the form
Ln = Kℓn and ρ proportional to the product of constants δ and ϑ from
Lemmas 3.5 and 3.6. All of the volumes ΛLn below are centered at x so we
omit that fact from the notation.
Our strategy is as follows: We pick an ε > 0 and show that, with probabil-
ity at least 1 −ε and some integer n, the largest connected component CLn
in ΛLn is connected to at least ρ|ΛLn′| sites in ΛLn′ , for every n′ ≥n. (Note
that this guarantees that CLn ⊂C∞.) Once this has been established, we
observe that |CLn′(x)| < ρ|ΛLn′| implies that x cannot be connected to CLn
within ΛLn′ . Assuming that x ∈C∞, the box ΛLn then contains at least
two distinct sites x,y ∈C∞which are not connected within ΛLn′ . By the
uniqueness of the inﬁnite cluster, the probability of the latter event can be
made smaller than ε by making n′ suﬃciently large. But then the limit in
(3.4) must be less than 2ε and, since ε was arbitrary, it must equal to zero.

<!-- page 26 -->

26
M. BISKUP
To make the proof complete, it remains to establish the ﬁrst claim in the
previous paragraph. Namely, we must show that the probability that CLn
is not connected to at least ρ|ΛLn′| sites in ΛLn′ for some n′ ≥n is less
than ε, provided n is suﬃciently large. To that end, let Bk, with k ≥0, be
a sequence of boxes (generally not centered at x) such that B0 = ΛLn and
that Bk is the maximal box in ΛLn+k that is centered on the (1,1,...,1)
half-axis, disjoint from all the previous Bk’s and with side a multiple of K.
Since ℓ> 3, it is easy to see that |Bk| grows proportionally to |ΛLn+k| (in
fact, |Bk|/|ΛLn+k| ≥3−d). Our goal is to show that the largest connected
components in all Bk’s are with overwhelming probability connected.
Invoking Lemmas 3.5 and 3.6 and choosing n suﬃciently large, the prob-
ability that each box Bk—viewed as the disjoint union of translates of ΛK—
contains a component comprising at least ϑ|Bk|/Kd maximal connected
components C (Kx)
K
of size at least δKd is bounded by ε/2. On the other
hand, the probability that the corresponding components in Bk and Bk+1
are not connected is bounded by
exp

−β ϑ2
K2d
|Bk||Bk+1|
(bLn+k+1/K)s

≤exp{−β′ℓ(2d−s)(n+k)},
(3.51)
where b is a constant such that bLn/K bounds the distance |x −y| for any
translates ΛK(Kx) and ΛK(Ky) contained in ΛLn, and where β′/β is a
constant that depends only on ℓand K. The right-hand side is summable
on k and the sum can be made smaller than ε/2 by increasing n. Thus with
probability as least 1 −ε, for each k ≥0 the component CLn is connected to
at least δϑ|Bk| sites in ΛLn+k. Choosing ρ < δϑ3−d, the above claim follows.
□
Proof of Corollary 3.4.
Let ℓ> 1 be an odd integer. Clearly, it
suﬃces to prove the result for D(ρ,2ℓ)
L
instead of D(ρ,ℓ)
L
and L a multiple of 3ℓ.
Viewing ΛL as a disjoint union of boxes Λ3ℓ(3ℓx) with x ∈Zd, let C (3ℓx)
ℓ
be
a maximal connected component in Λℓ(3ℓx). For ρ′ > 0, let Aℓ(x) be the
event that |C (3ℓx)
ℓ
| ≥ρ′|Λℓ| is true. Let NL,ℓdenote the number of x ∈Zd
with Λ3ℓ(3ℓx) ⊂ΛL such that Aℓ(x) occurs.
The events Aℓ(x) are independent and, if ℓis large enough and ρ′ > 0 is
suﬃciently small, Theorem 3.2 shows that Aℓ(x) occurs with probability at
least r = 1 −exp{−ρ′ℓ2d−s′}. This allows us to invoke the Chernoﬀbound
once again with the result
P(NL,ℓ< r′Ld/(3ℓ)d) ≤e−ψ(r′,r)Ld/(3ℓ)d,
(3.52)
where ψ(r′,r) is as in (3.17). Choosing r′ = 1
2, α = 1
2 and C = C(1
2) and
taking ℓso large that ρ′ℓ2d−s′ ≥2C, (3.11) gives us
ψ(r′,r) ≥1
4(ρ′ℓ2d−s′ −C) ≥1
8ρ′ℓ2d−s′.
(3.53)

<!-- page 27 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
27
On the other hand, on {NL,ℓ≥r′Ld/(3ℓ)d} we have
|D(ρ,2ℓ)
L
| ≥(ρ′ℓd)(r′Ld/(3ℓ)d) = ρ′r′3−dLd,
(3.54)
and so {|D(ρ,2ℓ)
L
| < ρLd} ⊂{NL,ℓ< r′Ld/(3ℓ)d} once ρ ≤ρ′r′3−d. Invoking
(3.52) and (3.53) and choosing ρ such that also ρ ≤1
83−dρ′, the desired
estimate follows.
□
4. Proof of main result.
In this section we will provide the proof of
Theorem 1.1. The arguments closely follow the outline presented in Sections
2.1 and 2.2. The reader may consider skimming through these sections once
again before plunging into the ultimate details of the proofs.
4.1. Upper bound.
The principal goal of this section is to establish the
upper bound in (1.6). By the continuity of s 7→∆(s,d) it suﬃces to prove the
upper bound for any number s exceeding the limit (1.3), so we will instead
assume that s obeys (3.1). The desired claim is then formulated as follows:
Proposition 4.1.
Let s ∈(d,2d) be such that (3.1) holds and let ∆=
∆(s,d) be as in (1.2). For each ∆′ > ∆and each ε > 0, there exists an
N0 < ∞such that
P(D(x,y) ≥(log |x −y|)∆′, x,y ∈C∞) ≤ε
(4.1)
holds for all x,y ∈Zd with |x −y| ≥N0.
As discussed in Section 2.1, the proof is conceptually rather simple: For
each pair of sites x and y we will construct a hierarchy of an appropriate
depth k connecting x and y, such that pairs (zσ01,zσ10) with σ ∈{0,1}k−2
are connected by paths of length (log |x −y|)o(1). The main diﬃculty stems
from the requirement that the bonds constituting the hierarchy be connected
in a prescribed (linear) order. This will be ensured by the condition that
all sites constituting the hierarchy are surrounded by a suﬃciently dense
connected component.
Recall our notation that ΛL(x) is a box of side L centered at x and CL(x)
is the set of sites in ΛL(x) connected to x by a path in ΛL(x). We will require
that the sites zσ are dense points according to the following deﬁnition:
Definition 4.1.
Given a number ρ ∈(0,1) and an odd integer ℓ> 1, we
will call x ∈Zd a (ρ,ℓ)-dense (or, simply, dense) point if |Cℓ(x)| ≥ρ|Λℓ(x)|.
For any real L > 0 suﬃciently large, let L+ be the minimal odd integer
larger than L and let L−be the minimal odd integer larger than L/2. Let
BL(x) = ΛL+(x) \ ΛL−(x).
(4.2)

<!-- page 28 -->

28
M. BISKUP
Given a number ρ ∈(0,1) and an odd integer ℓ> 1, let D(ρ,ℓ)
L
(x) denote the
set of all (ρ,ℓ)-dense points in BL(x). The input needed from Section 3 then
comes directly from Corollaries 3.3 and 3.4. By Corollary 3.4 and the fact
that BL(x) contains a box of side at least L/3, we know that
P(|D(ρ,ℓ)
L
(x)| ≤ρLd) ≤e−ρL2d−s′
(4.3)
once ρ is suﬃciently small, s′ ∈(s,2d) and ℓ0 ≤ℓ≤L/ℓ0. Corollary 3.3 in
turn shows that if x,y ∈C∞, then both x and y are dense points in the sense
that for each ε > 0 there exists an ℓ1 = ℓ1(ε) < ∞such that
P(|Cℓ(x)| ≤ρℓd,x ∈C∞) ≤ε
(4.4)
is true whenever ℓ≥ℓ1. A similar statement holds for y.
Now we can deﬁne the principal events: Let γ ∈(s/(2d),1) and let x and y
be two sites in Zd. Let N = |x −y| and deﬁne Nn = N γn. For each n ≥1,
let Bn = B(ρ,ℓ)
n,γ (x,y) be the event that there exists a hierarchy Hn(x,y) of
depth n connecting x and y subject to the following constraints: For all k =
0,1,... ,n −2 and all σ ∈{0,1}k,
zσ01 ∈D(ρ,ℓ)
Nk+1(zσ0)
and
zσ10 ∈D(ρ,ℓ)
Nk+1(zσ1).
(4.5)
The event Bn ensures that all sites of the hierarchy—except perhaps x and
y—are (ρ,ℓ)-dense points. To cover these exceptions, we also introduce the
event T = T (ρ,ℓ)(x,y) that both x and y are (ρ,ℓ)-dense in the above sense. In
the following, we will regard the number ρ as ﬁxed—such that (4.3) and (4.4)
hold—but ℓand γ will be adjustable.
The requirements (4.5) become appreciated in the proof of the following
bound:
Lemma 4.2.
For each ε ∈(0,1), each γ ∈(s/(2d),1) and each ∆′ satis-
fying
∆′ >
log 2
log(1/γ),
(4.6)
there exists a constant N ′ = N ′(ε,γ,∆′) < ∞such that the following is true
for all x,y ∈Zd with N = |x−y| ≥N ′: Let n be the maximal positive integer
such that
nlog(1/γ) ≤log log N −εlog log log N.
(4.7)
If ℓin the deﬁnition of the events Bn = B(ρ,ℓ)
n,γ (x,y) and T = T (ρ,ℓ)(x,y) is
an odd integer between Nn and 2Nn, then
P({D(x,y) ≥(log N)∆′} ∩Bn ∩T ) ≤ε.
(4.8)

<!-- page 29 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
29
The reason why we choose n as in (4.7) can be seen from the following
bounds:
2n ≤(log N)log 2/log(1/γ)
and
e(1/γ)(log log N)ε ≥N γn ≥e(log log N)ε.
(4.9)
These bounds will be important in the upcoming proof.
Proof of Lemma 4.2.
The main reason why N has to exceed a certain
constant is because we need the scales corresponding to successive levels of
the hierarchy to be clearly separated. To that end we observe that for all
k ≤n and Nk = N γk, we have
log Nk
Nk+1
≥(1 −γ)(log log N)ε,
(4.10)
which tends to inﬁnity as N →∞.
Introduce the abbreviation Bn = T ∩Bn. If Bn occurs, then there exists at
least one hierarchy Hn(x,y) of depth n connecting x and y such that (4.5)
is satisﬁed. Then (4.10) guarantees the existence of numbers N ′ < ∞and
b ∈(0,1) such that the following is true for any such hierarchy: If σ ∈{0,1}k
with k = 0,1,...,n −2, then
|zσ01 −zσ10| ≥bNk,
(4.11)
while if z ∈BNk+2(zσ0) and z′ ∈BNk+2(zσ01), then
b−1Nk+1 ≥|z −z′| ≥bNk+1,
(4.12)
whenever N ≥N ′. Similar statements hold for pairs z ∈BNk+2(zσ1) and
z′ ∈BNk+2(zσ10). Moreover, N ′ can be chosen so large that also the bounds
b−1Nn−1 < bNn−2
and
bNn−1 > diamΛℓ
(4.13)
hold true for all ℓbetween ℓ0 and 2Nn and all n satisfying (4.7).
Let An be the event that, for any hierarchy that would make Bn satisﬁed,
at least one of the “gaps” of the hierarchy, say (zσ0,zσ01) where σ ∈{0,1}n−2,
fails to have the components Cℓ(zσ0) and Cℓ(zσ01) connected by an occupied
bond. (Note that these components are quite large because both zσ0 and zσ01
are dense points.) We claim that
{D(x,y) ≥(log N)∆′} ∩Bn ⊂An ∩Bn.
(4.14)
Indeed, if all “gaps” do have the corresponding components connected, then
each zσ0 is connected to zσ01 by a path of no more than 1 + 2ℓd bonds [note
that ℓd bonds should be enough to get out of Cℓ(zσ0)], and similarly for the
pairs zσ1 and zσ10. Noting that a hierarchy of depth n involves only 2n−1
“gaps” and 2n−1 −1 bonds, we can use ℓ≤2Nn and (4.9) to write
D(x,y) ≤2n−1(1 + 2d+1N d
n) + 2n−2
(4.15)
≤2d+2(log N)log 2/log(1/γ)e(d/γ)(log log N)ε.

<!-- page 30 -->

30
M. BISKUP
In light of (4.6), this is not compatible with D(x,y) ≥(log N)∆′ if N is
suﬃciently large.
To ﬁnish the proof, we thus need to estimate the probability of An ∩Bn.
The above estimates show that the occurrence of Bn is determined by look-
ing only at the bonds longer than bNn−2 (to ensure the existence of a hi-
erarchy) or shorter than diamΛℓ[to ensure that all sites in the hierarchy
are (ρ,ℓ)-dense]. Explicitly, let F denote the σ-algebra generated by the
random variables (ωzz′) with |z −z′| ≥bNn−2 or |z −z′| ≤diamΛℓ. Then
(4.11) and (4.13) show that Bn ∈F. This allows us to prove (4.8) by condi-
tioning: Let β be a number such that pzz′ ≥1 −exp{−β|z −z′|−s} for any
pair (z,z′) of sites with |z −z′| ≥Nn. Then we have
P(An|F) ≤2n−1 exp

−βρ2
N 2d
n
(Nn−1/b)s

on Bn.
(4.16)
Here we used that, on Bn, the components Cℓ(zσ0) and Cℓ(zσ01) are both
larger than ρℓd ≥ρN d
n, while (4.12) dictates that the longest bond that
can connect them is not longer than Nn−1/b. The prefactor represents the
number of “gaps” in the hierarchy, which is the number of places where An
can fail. Inserting the upper bound on n from (4.7), the estimate (4.8) follows
once N is suﬃciently large.
□
Our next goal is to show that the event Bc
n is quite unlikely to occur:
Lemma 4.3.
Let γ ∈(s/(2d),1) and let s′ ∈(s,2dγ). Let N = |x −y| ≥
N ′ where N ′ is as in Lemma 4.2 and deﬁne Nk = N γk. Then there is a con-
stant c3 > 0 such that if ℓin the deﬁnition of Bk is an odd integer between Nn
and 2Nn, then
P(Bc
k+1 ∩Bk) ≤2k+1 exp{−c3N 2dγ−s′
k
},
(4.17)
for all k < n, where n is as in Lemma 4.2. In particular,
P(Bc
n) ≤2n+1 exp{−c3N 2dγ−s′
n
}.
(4.18)
Proof.
Clearly, (4.18) is a result of summing (4.17), so we just need
to prove (4.17) for all k = 0,1,...,n. By the fact that N ≥N ′, we can as-
sume that the scales Nk and Nk+1 are clearly separated in the sense of the
inequalities (4.11) and (4.13). Let B′
m be the event that there exists a hier-
archy Hm(x,y) of depth m connecting x and y such that for each k ≤m −2
and each σ ∈{0,1}k,
zσ01 ∈BNk+1(zσ0)
and
zσ10 ∈BNk+1(zσ1).
(4.19)
A comparison with (4.5) shows that Bk ⊂B′
k. Consider also the following
events:

<!-- page 31 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
31
1. The event A1 that, for any hierarchy Hk(x,y) that would make B′
k satis-
ﬁed, we have |D(ρ,ℓ)
Nk (zσ)| ≤ρN d
k for some σ ∈{0,1}k.
2. The event A2 that, for any hierarchy that would make B′
k satisﬁed, there
exists a pair of sites (z,z′) of the type (zσ0,zσ01) or (zσ1,zσ10) with σ ∈
{0,1}k−2 such that there is no occupied bond between the sets D(ρ,ℓ)
Nk (z)
and D(ρ,ℓ)
Nk (z′).
Now on Bk ∩Bc
k+1 there exists a hierarchy that would make B′
k satisﬁed,
but such that for some pair of sites as in the deﬁnition of A2, the sets
D(ρ,ℓ)
Nk (z) and D(ρ,ℓ)
Nk (z′) are not connected by an occupied bond. It follows
that Bk ∩Bc
k+1 ⊂B′
k ∩A2. The event A1 will be used to write A2 as the union
of A1 and Ac
1 ∩A2, whose probabilities are more convenient to estimate.
The proof now mimics the argument from the proof of Lemma 4.2. By the
fact that N ≥N ′, the event B′
k is determined by looking only at the bonds
that are longer than bNk−2. Let F ′ denote the σ-algebra generated by the
random variables (ωzz′) with |z −z′| ≥bNk−2. Then B′
k ∈F ′. On the other
hand, conditional on F ′, the event A1 is only determined by looking at the
bonds that are shorter than diamΛℓ. By (4.3), we have
P(A1|F ′) ≤2k exp{−ρN 2d−s′
k
}
on B′
k.
(4.20)
Here 2k counts the number of pairs where A1 can go wrong.
Concerning the event A2, we note that conditional on Ac
1 ∩Bk, the event
A2 is determined by the bonds of length between bNk−1 and b−1Nk−1, which
by (4.13) must be either longer than diamΛℓor shorter than bNk−2. Let
F be the σ-algebra generated by (ωzz′) with |z −z′| ≥bNk−2 or |z −z′| ≤
diamΛℓ. Then Ac
1 ∩B′
k ∈F and A2 is determined by bonds independent
of F. Let β be the same constant as in the proof of Lemma 4.2. Then we
have
P(A2|F) ≤2k exp

−βρ2
N 2d
k
(Nk−1/b)s

on Ac
1 ∩Bk.
(4.21)
Putting these bounds together and choosing c3 appropriately, (4.17) directly
follows.
□
Lemmas 4.2 and 4.3 ﬁnally allow us to prove Proposition 4.1.
Proof of Proposition 4.1.
Let ∆′ > ∆and let ε ∈(0,1). Choose
γ ∈(s/(2d),1) such that (4.6) holds true and pick an s′ ∈(s,2dγ). Suppose
N ≥N ′, where N ′ is the constant from Lemmas 4.2 and 4.3, and let n be
as in Lemma 4.2. Fix an odd integer ℓbetween Nn and 2Nn and let c3 be
the constant from Lemma 4.3.

<!-- page 32 -->

32
M. BISKUP
Invoking the inclusion
{D(x,y) ≥(log N)∆′}
(4.22)
⊂({D(x,y) ≥(log N)∆′} ∩Bn ∩T ) ∪Bc
n ∪T c,
we just need to estimate the probability of the three events on the right-hand
side. Lemma 4.2 shows that the probability of the ﬁrst event is less than ε.
Lemma 4.3 in conjunction with the bounds (4.9) shows that
P(Bn) ≤2exp{∆′(log log N) −c3e(2dγ−s′)(log log N)ε},
(4.23)
which can also be made less than ε by choosing N suﬃciently large. Finally,
the probability P(T c) is estimated from (4.4) where we assume that N is so
large that also ℓ≥Nn ≥ℓ1. From here the desired claim follows.
□
4.2. Lower bound.
The goal of this section is to prove the lower bound
in (1.6). As in Section 4.1, we formulate the relevant claim as a separate
proposition:
Proposition 4.4.
Suppose that (1.3) holds with an s ∈(d,2d) and let ∆=
∆(s,d) be as in (1.2). For each ∆′ < ∆and each ε > 0, there exists an
N0 < ∞such that
P(D(x,y) ≤(log |x −y|)∆′) ≤ε
(4.24)
holds for all x,y ∈Zd with |x −y| ≥N0.
In conjunction with Proposition 4.1, this result immediately implies The-
orem 1.1.
Proof of Theorem 1.1.
Let ε > 0 and let Dε(x) be the event in (1.6).
Choosing ∆′ such that |∆−∆′| ≤ε, Propositions 4.1 and 4.4 ensure that
lim|x|→∞P(Dε(x)c ∩{0,x ∈C∞}) = 0. Then (1.6) follows by noting that, by
FKG inequality and translation invariance, we have P(0,x ∈C∞) ≥P(0 ∈
C∞)2, which is positive by our assumption that there is percolation.
□
The remainder of this section will be spent on the proof of Proposition 4.4.
As discussed in Section 2.2, our strategy will be to show that each path
connecting x and y in less than (log |x−y|)∆steps contains a hierarchy whose
“gaps” obey the conditions (2.7). (As far as this claim is concerned, the
speciﬁc choice of the exponent plays no essential role; any positive number
will do.) This will be used to control the combined length of the paths needed
to span the “gaps” and show that it will eventually exceed (log |x−y|)∆′ for
any ∆′ < ∆.

<!-- page 33 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
33
We begin by deﬁning the relevant events. Let x,y ∈Zd be distinct (and
distant) sites and let N = |x −y|. Fix a number γ ∈(0,s/(2d)) and, for each
integer n ≥2, let En = En,γ(x,y) be the event that every hierarchy Hn(x,y)
of depth n connecting x and y such that
|zσ01 −zσ10| ≥|zσ0 −zσ1|(log N)−∆
(4.25)
holds for all k = 0,1,...,n−2 and all σ ∈{0,1}k will also satisfy the bounds
Y
σ∈{0,1}k
|zσ0 −zσ1| ∨1 ≥N (2γ)k
(4.26)
for all k = 1,2,... ,n −1. Here “∨” is a shorthand for maximum.
Remark 4.1.
Since we allow the possibility of “site collapse” in our
deﬁnition of a hierarchy (e.g., we do not forbid that zσ00 = zσ01), we must
use a “∨” on the left-hand side of (4.26). Note that (4.26) is a precise from
of (2.7) while (4.25) is a precise from of (2.6).
Our ﬁrst goal is to estimate the probability of Ec
n:
Lemma 4.5.
Let γ ∈(0,s/(2d)) and let s′ ∈(2dγ,s) be such that s′ > d.
Let En = En,γ(x,y) be as above. Then there exists a constant c4 ∈(0,∞) such
that for all x,y ∈Zd with N = |x −y| satisfying γn log N ≥2(s′ −d),
P(Ec
n+1 ∩En) ≤(log N)c42nN −(s′−2dγ)(2γ)n.
(4.27)
The proof of Lemma 4.5 requires certain combinatorial estimates whose
precise statements and proofs have been deferred to Lemmas A.1 and A.2
in the Appendix. We encourage the reader to skim through the statements
of these lemmas before plunging into the forthcoming proof.
Proof of Lemma 4.5.
On Ec
n+1 ∩En, there exists a hierarchy Hn(x,y)
such that the bound (4.26) holds for all k = 1,...,n−1 but does not hold for
k = n. In order to estimate the probability of such an event, let Θ(n) be the
collection of all 2n-tuples (zσ) of sites such that (4.25) holds for all σ ∈{0,1}k
with k = 0,1,... ,n −1, while (4.26) is true only for k = 1,...,n −1 but not
for k = n. Then we can write
P(Ec
n+1 ∩En) ≤
X
(zσ)∈Θ(n)
n−1
Y
k=0
Y
σ∈{0,1}k
p(zσ01,zσ10),
(4.28)
where p(z,z′) = 1 −e−q(z−z′) for z̸ = z′—see (1.4)—while p(z,z′) = 1 for z =
z′. As speciﬁed in the deﬁnition of the hierarchy, none of the bonds (zσ01,zσ10)
may appear more than once, whence (4.28) follows by invoking inclusion–
exclusion and independence.

<!-- page 34 -->

34
M. BISKUP
In order to estimate the right-hand side of (4.28), we will introduce a
convenient change of variables: For each k = 0,1,... ,n and each σ ∈{0,1}k,
let
tσ = zσ0 −zσ1.
(4.29)
(Thus, t∅is just x −y, while t0 represents the “gap” z00 −z01 and t1 rep-
resents the “gap” z10 −z11. Note that the Nσ’s from Section 2.2 are related
to the tσ’s via Nσ = |tσ|.) Clearly, once x and y are ﬁxed and tσ are deﬁned
for all σ ∈{0,1}k and all k = 1,...,n, all of zσ with σ ∈{0,1}n+1 can be
reconstructed from (4.29). In terms of the tσ’s, the conditions (4.26) can be
written as
Y
σ∈{0,1}k
|tσ| ∨1 ≥N (2γ)k,
(4.30)
which on Ec
n+1 ∩En is required to hold for all k = 1,2,...,n −1 and to fail
for k = n, while (4.25) can be rewritten as
|zσ01 −zσ10| ≥|tσ|(log N)−∆,
(4.31)
which is required to hold for all k = 0,1,...,n −1 and all σ ∈{0,1}k.
The latter condition—(4.31)—allows us to recast (4.28) entirely in terms
of the tσ’s. Indeed, let Ω(k) be the set of all collections (tσ), σ ∈{0,1}k,
of elements from Zd such that (4.30) holds true and let C ∈(0,∞) be a
constant so large that
p(z,z′) ≤
C
|z −z′|s′
(4.32)
is true for all distinct z,z′ ∈Zd. Then (4.32) and (4.30) allow us to write
P(Ec
n+1 ∩En)
(4.33)
≤
X
(tσ)∈Ω(1)
...
X
(tσ)∈Ω(n−1)
X
(tσ)/∈Ω(n)
n−1
Y
k=0
Y
σ∈{0,1}k
C(log N)s′∆
(|tσ| ∨1)s′ ,
where t∅= x −y and where we assumed that C is so large that the last
fraction exceeds 1 whenever |tσ| ≤1.
The right-hand side of (4.33) is now estimated as follows: First we will
extract the terms C(log N)s′∆and write the sequence of sums as a product
by grouping the corresponding tσ’s with their sum (and noting that |t∅| =
N). This gives
P(Ec
n+1 ∩En)
(4.34)
≤[C(log N)s′∆]2n
N s′
 
X
(tσ)/∈Ω(n)
1
! n−1
Y
k=1
 
X
(tσ)∈Ω(k)
Y
σ∈{0,1}k
1
(|tσ| ∨1)s′
!
.

<!-- page 35 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
35
Now s′ > d, which implies that the sum in the second parentheses can be
estimated using Lemma A.1. Explicitly, noting that |{z ∈Zd :n ≤|z| ∨1 <
n + 1}| ≤c(d)nd−1 for all n ≥1 and some ﬁxed constant c(d) < ∞, we intro-
duce a collection of positive integers (nσ) and ﬁrst sum over all (tσ) subject
to the constraint nσ ≤|tσ| ∨1 < nσ + 1. Then we are in a position to apply
Lemma A.1 with α = s′ −d, b = N γk and κ = 2k, which yields
X
(tσ)∈Ω(k)
Y
σ∈{0,1}k
1
(|tσ| ∨1)s′ ≤(C′ log N)2k
N (s′−d)(2γ)k ,
(4.35)
for some C′ < ∞independent of N and n. (Here we used that γk ≤1 to
bound log N γk by log N.) The sum in the ﬁrst parentheses can be estimated
in a similar fashion; the result of application of Lemma A.2 with α = d,
b = N γn and κ = 2n is
X
(tσ)/∈Ω(n)
1 ≤(C′′ log N)2nN d(2γ)n
(4.36)
for some C′′ < ∞. Combining these estimates with (4.34) and invoking the
identity (2.9), the desired bound (4.27) is proved.
□
Lemma 4.5 will be used to convert the event {D(x,y) ≤(log N)∆′} into
a statement about the total number of bonds needed to span the “gaps”
of a hierarchy identiﬁed within (one of ) the shortest paths connecting x
and y. Let Fn = Fn(x,y) be the event that, for every hierarchy of depth n
connecting x and y and satisfying (4.25), every collection of (bond) self-
avoiding and mutually (bond) avoiding paths πσ with σ ∈{0,1}n−1, such
that πσ connects zσ0 with zσ1 without using any bond from the hierarchy,
will obey the bound
X
σ∈{0,1}n−1
|πσ| ≥2n.
(4.37)
Then we have the following claim:
Lemma 4.6.
Let ∆′ < ∆. If N = |x −y| is suﬃciently large and
n > ∆′
log 2 log log N,
(4.38)
then
{D(x,y) ≤(log N)∆′} ∩Fn = ∅.
(4.39)
Proof.
We will show that on {D(x,y) ≤(log N)∆′} there exists a hier-
archy of a depth satisfying (4.38) such that (4.25) is true, and a collection of

<!-- page 36 -->

36
M. BISKUP
paths πσ “spanning” the “gaps” of this hierarchy such that (4.37) is violated.
Let π be a path saturating the distance D(x,y) between x and y. The path π
is necessarily (bond) self-avoiding. Since |π| ≤(log N)∆—by our restriction
to {D(x,y) ≤(log N)∆′}—and since x,y 7→|x −y| satisﬁes the triangle in-
equality, the path π must contain a bond whose length exceeds N/(log N)∆.
Let z01 be the endpoint of this bond on the x-side and let z10 denote the
endpoint on the y-side of the bond. Denoting z0 = x and z1 = y, we thus
have
|z01 −z10| ≥|z0 −z1|(log N)−∆,
(4.40)
that is, (4.25) for σ = ∅.
Similarly, we will identify the next level of the hierarchy. Let π0 be the
portion of π between z00 = x and z01, and let π1 be the portion of π
between z11 = y and z10. (Note that this agrees with our notation used
in the deﬁnition of the event Fn.) Again, we have |π0|,|π1| ≤(log N)∆
and thus both π0 and π1 contain at least one bond of length exceeding
|z00 −z01|/(log N)∆and |z10 −z11|/(log N)∆, respectively. The endpoints of
this bond in π0 identify the sites z001 and z010, and similarly for the end-
points of the bond in π1. (If π0 is empty, which can only happen if z01 = z0,
we let z001 = z010 = z0, and similarly for π1.) The very construction of these
bonds implies (4.25) for σ = 0,1. Proceeding in a similar way altogether n-
times, we arrive at a hierarchy of depth n connecting x and y and satisfying
(4.25).
The construction implicitly deﬁnes a collection of paths πσ with σ ∈
{0,1}n−1 such that πσ is the portion of π connecting the endpoints of the
“gap” (zσ0,zσ1). Now we are ready to prove (4.39). Indeed, if {D(x,y) ≤
(log N)∆′} occurs, the combined length of all πσ’s must be less than (log N)∆′,
which by (4.37) is strictly less than 2n. But then there exists a hierar-
chy of depth n and self-avoiding and mutually avoiding paths πσ “span-
ning” its “gaps” such that (4.37) is violated. Consequently, we must have
{D(x,y) ≤(log N)∆′} ⊂Fc
n.
□
In light of Lemma 4.6, to prove Proposition 4.4, we will need a bound
on the probability of Fc
n for some n obeying (4.38). However, invoking also
Lemma 4.5, we can as well focus just on the event Fc
n ∩En.
Lemma 4.7.
Let γ ∈(0,s/(2d)) and let s′ ∈(2dγ,s) be such that s′ > d.
Let En = En,γ(x,y) and Fn = Fn(x,y) be as above. Then there exists a con-
stant c5 ∈(0,∞) such that for all distinct x,y ∈Zd with N = |x−y| satisfying
γn log N ≥2(s′ −d),
P(Fc
n ∩En) ≤(log N)c52nN −s′(2γ)n−1.
(4.41)

<!-- page 37 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
37
Proof.
The proof will closely follow that of Lemma 4.5 so we will stay
rather brief. In fact, the only essential diﬀerence is that, instead of (4.36)—
which we cannot use because we are no longer in the complement of En+1—
the necessary decay for the last sum will have to be provided on the basis
of the containment in Fc
n.
We begin by noting that on Fc
n ∩En, the following events must occur:
1. There exists a hierarchy Hn(x,y) such that (4.25) and (4.26) hold.
2. There exists a collection of self-avoiding and mutually avoiding paths πσ,
with σ ∈{0,1}n−1, such that πσ connects zσ0 with zσ1 without using any
bonds from Hn(x,y).
3. The bound (4.37) fails.
As in Lemma 4.5, we will use brute force: First we will ﬁx a hierarchy
satisfying the desired condition and try to estimate the probability that,
for some collection of nonnegative integers (mσ), the length of the path πσ
is mσ for each σ ∈{0,1}n−1. The fact that the paths and the hierarchy are
all disjoint then allows us to write
P(Fc
n occurs and (|πσ|) = (mσ)|Hn(x,y)) ≤
Y
σ∈{0,1}n−1
Qmσ(zσ0,zσ1),
(4.42)
where
Qm(z,z′) =
X
π=(z0,...,zm)
z0=z, zm=z′
m−1
Y
i=0
C′
(|zi −zi+1| ∨1)s′ .
(4.43)
Here the sum runs over self-avoiding paths π of length m and C′ is so large
that the last quotient is an upper bound on the probability that zi and zi+1
are connected by an occupied bond.
To estimate (4.42), we ﬁrst need a bound on Qm(z,z′). To that end we
note that, in light of the inequality s′ > d, there exists a constant a ∈(1,∞)
such that for all x,y ∈Zd,
X
z∈Zd
1
(|x −z| ∨1)s′
1
(|y −z| ∨1)s′ ≤
a
(|x −y| ∨1)s′ .
(4.44)
From here we conclude that
Qm(z,z′) ≤
(C′a)m
(|z −z′| ∨1)s′ ;
(4.45)
that is, up to a multiplicative factor, Qm(z,z′) acts similarly as p(z,z′).
The paths still carry some entropy in the choice of the integers (mσ) which
amounts to counting the number #(m,n) of ordered partitions of a non-
negative integer m into 2n−1 nonnegative integers. A simple estimate shows

<!-- page 38 -->

38
M. BISKUP
that #(m,n) ≤2m+2n−1 and, noting that on Fc
n we only need to consider
m < 2n,
P(Fc
n|Hn(x,y)) ≤(4aC′)2n
Y
σ∈{0,1}n−1
1
(|zσ0 −zσ1| ∨1)s′ ,
(4.46)
because P
m<2n #(m,n) ≤42n.
Having dispensed with the paths πσ, we now start estimating the prob-
ability of Fc
n ∩En. Let Θ⋆(n) be the set of all collections (zσ), σ ∈{0,1}n,
obeying (4.25) for k = 0,1,... ,n −2 and (4.26) for k = 1,...,n −1. The
bounds (4.32), (4.25) and (4.46) then give
P(Fc
n ∩En) ≤(4aC′)2n
X
(zσ)∈Θ⋆(n)
n−1
Y
k=0
Y
σ∈{0,1}k
(C log N)s′2n
(|zσ0 −zσ1| ∨1)s′ .
(4.47)
Here C is the same constant as in (4.33) and the product still goes only up
to (n −1)—despite the insertion of the terms from (4.46)—because we are
now looking only at a hierarchy of depth n (and not n + 1 as in the proof of
Lemma 4.5). Passing again to the variables tσ = zσ0 −zσ1 and recalling the
deﬁnition of Ω(k) from the proof of Lemma 4.5, we now get
P(Fc
n ∩En) ≤[C′′(log N)s′∆]2n
N s′
n−1
Y
k=1
 
X
(tσ)∈Ω(k)
Y
σ∈{0,1}k
1
(|tσ| ∨1)s′
!
,
(4.48)
where C′′ < ∞. Each term in the product can now be estimated by (4.35).
Using
s′ + (s′ −d)
n−1
X
k=1
(2γ)k ≥s′(2γ)n−1
(4.49)
instead of (2.9), the estimate (4.41) directly follows.
□
Having assembled all necessary ingredients, we can now ﬁnish the proof
of Proposition 4.4.
Proof of Proposition 4.4.
Let ∆′ < ∆and, recalling that 2−1/∆=
s/(2d), choose an s′ ∈(d,s) such that 2−1/∆′ < s′/(2d). Pick a number γ
such that
2−1/∆′ < γ < s′
2d
(4.50)
and let δ = 1
2(s′ −2dγ).
By Lemma 4.6, we have {D(x,y) ≤(log N)∆′} ⊂Fc
n once n satisﬁes the
bound (4.38). On the other hand, if n also obeys the bound
nlog(1/γ) ≤log log N −2log log log N,
(4.51)

<!-- page 39 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
39
which is possible for large N by (4.50), then we have γn log N ≥(log log N)2.
This shows that, for N large enough, the right-hand side of the bound from
Lemma 4.5 is less than N −δ(2γ)n and similarly for the bound in Lemma 4.7.
Consequently, both bounds are summable on n and, increasing N if neces-
sary, the result can be made smaller than any number initially prescribed.
Hence, for any ε > 0 and N suﬃciently large, we will have
P(Fc
n) ≤P(Ec
n) + P(Fc
n ∩En) ≤2ε
(4.52)
once n satisﬁes both (4.38) and (4.51). By the inclusion {D(x,y) ≤(log N)∆′} ⊂
Fc
n, this ﬁnishes the proof.
□
APPENDIX
Here we establish the bounds needed in the proof of Lemmas 4.5 and 4.7.
To that end, let κ be a positive integer and, for b > 0 real, let
Ξκ(b) =
(
(ni) ∈Nκ :ni ≥1,
κ
Y
i=1
ni ≥bκ
)
.
(A.1)
We will also use Ξ⋆
κ(b) to denote a (formal) complement of this set, that is,
the set of all collections (ni) ∈Nκ of positive integers such that Q
i ni < bκ.
Lemma A.1.
For each ε > 0 there exists a constant g1 = g1(ε) < ∞such
that
X
(ni)∈Ξκ(b)
κ
Y
i=1
1
n1+α
i
≤(g1 b−α log b)κ
(A.2)
is true for all α > 0, all b > 1 and all positive integers κ satisfying
α −κ −1
κlog b ≥ε.
(A.3)
Proof.
As is common for this kind of estimates, we will turn the sum
into an integral. With each (ni) ∈Ξκ(b), we will associate a unique hypercube
h(ni) = (ni) + [−1
2, 1
2)κ in Rκ and note that whenever (xi) ∈h(ni), we have
xi ≥ni −1/2 ≥ni/2 and xi ≤ni + 1
2 ≤2ni for all i = 1,...,κ. This implies
that the product on the left-hand side of (A.2) can be bounded by the
product of (xi/2)−(1+α) and
[
(ni)∈Ξκ(b)
h(ni) ⊂
(
(xi) ∈Rκ :2xi ≥1,
κ
Y
i=1
(2xi) ≥bκ
)
.
(A.4)
Noting that the h(ni) are disjoint, we can now bound the sum in (A.2) by
the integral over the set on the right-hand side. Relabeling 2xi by xi, we

<!-- page 40 -->

40
M. BISKUP
thus get
left-hand side of (A.2) ≤2κ(1+2α)
Z
···
Z
xi≥1,i=1,...,κ
Qκ
i=1 xi≥bκ
dx1 ··· dxκ
κ
Y
i=1
1
x1+α
i
.
(A.5)
To evaluate the integral, we introduce the substitutions xi = eyi followed by
zj = y1 + ··· + yj for j = 1,...,κ. Since yi ≥0, the zj’s are ordered and since
the integrand depends only on zκ, the integrals over z1,...,zκ−1 can readily
be performed. The result is
right-hand side of (A.5) = 2κ(1+2α)
Z ∞
κ log b
dz
zκ−1
(κ −1)!e−αz,
(A.6)
where we have now dropped the subscript “κ” from z. Now the assump-
tion (A.3) ensures that for z ≥κlog b, the function z 7→zκ−1e−αz is strictly
decreasing and, in fact, its logarithm is concave. Applying (A.3), we easily
derive that for any z ≥κlog b,
zκ−1e−αz ≤(κlog b)s−1b−ακe−ε(z−κ logb).
(A.7)
Substituting this into (A.6), the integral is now easily performed. The cal-
culation is concluded by using Stirling’s formula to cancel the factor κκ−1
coming from the previous estimation against the leading order of (κ −1)! in
the denominator.
□
Our next claim concerns a similar sum over the indices in Ξ⋆
κ(a):
Lemma A.2.
There exists a constant g2 < ∞such that for each α ≥1,
each b ≥e/4 and any positive integer κ,
X
(ni)∈Ξ⋆κ(b)
κ
Y
i=1
nα−1
i
≤(g2 bα log b)κ.
(A.8)
Proof.
A moment’s thought reveals that we only have to address the
case α = 1. We will call upon the argument from Lemma A.1. Indeed, re-
placing (A.4) by
[
(ni)∈Ξ⋆κ(a)
h(ni) ⊂
(
(xi) ∈Rκ :2xi ≥1,
κ
Y
i=1
xi ≤(2b)κ
)
,
(A.9)
we easily ﬁnd out that
X
(ni)∈Ξ⋆κ(a)
1 ≤2−κ
Z
···
Z
xi≥1,i=1,...,κ
Qκ
i=1 xi≤(4b)κ
dx1 ··· dxκ.
(A.10)

<!-- page 41 -->

CHEMICAL DISTANCE IN LONG-RANGE PERCOLATION MODELS
41
Invoking the same substitutions as before, we then get that the right-hand
side of (A.10) equals
2−κ
Z κ log(4b)
0
dz
zκ−1
(κ −1)! ez ≤2−κ(4b)κ 1
κ!(κlog(4b))κ.
(A.11)
Here we used the bound ez ≤(4b)κ to get rid of the exponential in the
integral and then integrated out. Invoking Stirling’s formula, the desired
bound directly follows.
□
Acknowledgments.
This work was started when the author was a mem-
ber of Theory Group in Microsoft Research, Redmond, WA. Discussions
with Noam Berger and Itai Benjamini are gratefully acknowledged. The
author wishes to thank an anonymous referee for suggestions that led to
improvements in the presentation.
REFERENCES
[1] Aizenman, M. and Barsky, D. J. (1987). Sharpness of the phase transition in
percolation models. Comm. Math. Phys. 108 489–526. MR874906
[2] Aizenman, M. and Newman, C. (1986). Discontinuity of the percolation density in
one-dimensional 1/|x−y|2 percolation models. Comm. Math. Phys. 107 611–647.
MR868738
[3] Alexander, K., Chayes, J. T. and Chayes, L. (1990). The Wulﬀconstruction
and asymptotics of the ﬁnite cluster distribution for two-dimensional Bernoulli
percolation. Comm. Math. Phys. 131 1–51. MR1062747
[4] Antal, P. (1994). Trapping problems for the simple random walk. Dissertation ETH
No. 10759.
[5] Antal, P. and Pisztora, A. (1996). On the chemical distance for supercritical
Bernoulli percolation. Ann. Probab. 24 1036–1048. MR1404543
[6] Benjamini, I. and Berger, N. (2001). The diameter of long-range percolation clus-
ters on ﬁnite cycles. Random Structures Algorithms 19 102–111. MR1848786
[7] Benjamini, I., Kesten, H., Peres, Y. and Schramm, O. (2004). The geometry
of the uniform spanning forests: Transitions in dimensions 4, 8, 12, ... . Ann.
Math. To appear.
[8] Berger, N. (2002). Transience, recurrence and critical behavior for long-range per-
colation. Comm. Math. Phys. 226 531–558. MR1896880
[9] Bollob´as, B. (2001). Random Graphs. Cambridge Univ. Press. MR1864966
[10] Burton, R. M. and Keane, M. (1989). Density and uniqueness in percolation.
Comm. Math. Phys. 121 501–505. MR990777
[11] Campanino, M., Chayes, J. T. and Chayes, L. (1991). Gaussian ﬂuctuations in
the subcritical regime of percolation. Probab. Theory Related Fields 88 269–341.
MR1100895
[12] Campanino, M. and Ioffe, D. (2002). Ornstein–Zernike theory for the Bernoulli
bond percolation on Zd. Ann. Probab. 30 652–682. MR1905854
[13] Cerf, R. (2000). Large deviations for three dimensional supercritical percolation.
Ast´erisque 267 vi+177. MR1774341
[14] Coppersmith, D., Gamarnik, D. and Sviridenko, M. (2002). The diameter
of a long-range percolation graph. Random Structures Algorithms 21 1–13.
MR1913075

<!-- page 42 -->

42
M. BISKUP
[15] Gandolfi, A., Keane, M. S. and Newman, C. M. (1992). Uniqueness of the inﬁnite
component in a random graph with applications to percolation and spin glasses.
Probab. Theory Related Fields 92 511–527. MR1169017
[16] Grimmett, G. R. and Marstrand, J. M. (1990). The supercritical phase of perco-
lation is well behaved. Proc. Roy. Soc. London Ser. A 430 439–457. MR1068308
[17] Hara, T. and Slade, G. (2000). The scaling limit of the incipient inﬁnite cluster in
high-dimensional percolation. I. Critical exponents. J. Statist. Phys. 99 1075–
1168. MR1773141
[18] Hara, T. and Slade, G. (2000). The scaling limit of the incipient inﬁnite cluster in
high-dimensional percolation. II. Integrated super-Brownian excursion. J. Math.
Phys. 41 1244–1293. MR1757958
[19] Imbrie, J. and Newman, C. (1988). An intermediate phase with slow decay of
correlations in one-dimensional 1/|x −y|2 percolation, Ising and Potts models.
Comm. Math. Phys. 118 303–336. MR956170
[20] Menshikov, M. V. (1986). Coincidence of critical points in percolation problems.
Dokl. Akad. Nauk SSSR 288 1308–1311. (In Russian.) MR852458
[21] Milgram, S. (1967). The small-world problem. Psychology Today 1 61–67.
[22] Newman, C. M. and Schulman, L. S. (1986). One-dimensional 1/|j −i|s percolation
models: The existence of a transition for s ≤2. Comm. Math. Phys. 104 547–571.
MR841669
[23] Schulman, L. S. (1983). Long-range percolation in one dimension. J. Phys. A 16
L639–L641. MR723249
[24] Smirnov, S. (2001). Critical percolation in the plane. I. Conformal invariance and
Cardy’s formula. II. Continuum scaling limit. C. R. Acad. Sci. Paris S´er. I Math.
333 239–244. MR1851632
[25] Smirnov, S. and Werner, W. (2001). Critical exponents for two-dimensional per-
colation. Math. Res. Lett. 8 729–744. MR1879816
[26] Watts, D. J. and Strogatz, S. H. (1998). Collective dynamics of “small-world”
networks. Nature 303 440–442.
Department of Mathematics
University of California
Los Angeles, California 90095-1555
USA
e-mail: biskup@math.ucla.edu
url: http://www.math.ucla.edu/˜biskup/index.html
