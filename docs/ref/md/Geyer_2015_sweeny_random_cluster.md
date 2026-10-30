# Efficient simulation of the random-cluster model

> J. Stat. Mech. (2015) P08010 · [arXiv](https://arxiv.org/abs/1307.6647)
> Local copy: `docs/ref/Geyer_2015_sweeny_random_cluster.pdf` (14 pages)

<!-- page 1 -->

Eﬃcient simulation of the random-cluster model
Eren Metin El¸ci and Martin Weigel
Applied Mathematics Research Centre, Coventry University, Coventry, CV1 5FB, United Kingdom and
Institut f¨ur Physik, Johannes Gutenberg-Universit¨at Mainz, Staudinger Weg 7, D-55099 Mainz, Germany
(Dated: September 25, 2018)
The simulation of spin models close to critical points of continuous phase transitions is heavily
impeded by the occurrence of critical slowing down. A number of cluster algorithms, usually based on
the Fortuin-Kasteleyn representation of the Potts model, and suitable generalizations for continuous-
spin models have been used to increase simulation eﬃciency. The ﬁrst algorithm making use of this
representation, suggested by Sweeny in 1983, has not found widespread adoption due to problems in
its implementation. However, it has been recently shown that it is indeed more eﬃcient in reducing
critical slowing down than the more well-known algorithm due to Swendsen and Wang. Here, we
present an eﬃcient implementation of Sweeny’s approach for the random-cluster model using recent
algorithmic advances in dynamic connectivity algorithms.
I.
INTRODUCTION
In the vicinity of a continuous phase transition, par-
ticle or spin systems of statistical mechanics develop ex-
tended spatial correlations signaling the onset of long-
range translational order through spontaneous symme-
try breaking.
It has been realized early on that these
phenomena suggest a description of the ordering in geo-
metrical terms, using analogies to the percolation tran-
sition [1]. While Fisher’s droplet model initially consid-
ered simple clusters of like spins (geometrical clusters)
as the relevant quantities, it was only gradually real-
ized in the 1980s that the relevant collective degrees of
freedom, or “physical clusters”, are indeed of a diﬀer-
ent nature, and they can be constructed by breaking the
geometric clusters up following a suitable stochastic pre-
scription activating bonds within clusters with probabil-
ity p = 1 −exp(−J/kBT) < 1 [2]. Initially, this breakup
was only applied to one spin species (of the Ising model,
say), leading to inconsistencies in the symmetric high-
temperature phase.
It was in an independent line of
thought by Fortuin and Kasteleyn [3] that an equiva-
lence of the partition function of the q-state Potts model
[4] to a correlated bond-percolation problem known as
the random-cluster model with partition function
ZRC =
X
G′⊆G
vb(G′)qn(G′)
(1)
was established. Here, G′ denotes the set of b(G′) acti-
vated edges on a lattice graph G, resulting in n(G′) con-
nected components, and the bond weight v = p/(1 −p).
These results led Hu [5] to generalize the above bond ac-
tivation probability to all geometric clusters irrespective
of their orientation. As a consequence, the right choice of
physical droplets is now well understood, and the equiv-
alence of their percolation properties and thermal quan-
tities has been explicitly checked [6].
The growth without bounds of static correlations in
critical systems is accompanied in the time domain by
a divergence of relaxation times known as critical slow-
ing down [7]. While this is a physical phenomenon con-
nected, for instance, to the eﬀect of critical opalescence,
it is also of direct relevance for the pseudo-dynamics in
Monte Carlo simulations of near-critical systems.
As
this leads to an asymptotic ineﬃciency of Markov chain
Monte Carlo in producing independent samples, an im-
proved understanding of spatial correlations was hoped
to translate into suitable non-local updating procedures
allowing to precisely study near-critical systems. Initial
attempts in this direction, such as variants of the multi-
grid approach [8, 9], were based on renormalization group
ideas, and turned out to be only moderately successful.
The ﬁrst Monte Carlo algorithm based on the concept of
physical clusters was suggested by Sweeny in 1983 [10].
He considered a direct simulation of the bond variables of
Eq. (1), randomly suggesting state switches from active
to inactive and vice versa. As the relevant Boltzmann
weight depends on the number n(G′) of connected com-
ponents resulting from a given bond conﬁguration, cal-
culating the acceptance probability of the bond moves
needs up-to-date information about cluster connectiv-
ity.
Hence, a single update might require the expen-
sive traversal of large (possibly spanning) clusters, po-
tentially destroying the computational advantage of an
accelerated decorrelation of conﬁgurations through com-
putational critical slowing down [11].
An alternative suggestion made by Swendsen and
Wang [12] works directly on the spin conﬁguration, freez-
ing bonds between like spins with probability p = 1 −
exp(−J/kBT) and independently ﬂipping the resulting
spin clusters. Instead of working in the graph language
alone, this series of alternating updates of spin and bond
variables corresponds to a Markov chain in an augmented
state space [13]. The resulting algorithm (with its many
variants including, for instance, the single-cluster version
[14]) is rather straightforward to implement and turns
out to be very eﬃcient in beating critical slowing down,
reducing the dynamical critical exponent, e.g., of the 2D
Ising model from z ≈2 for local spin ﬂips to z ≈0.2
[15].
Owing to this success of the Swendsen-Wang al-
gorithm and related techniques as well as the delica-
cies of maintaining up-to-date connectivity information,
Sweeny’s approach was not used by many researchers.
Also, its reduction of critical slowing down was not pre-
arXiv:1307.6647v4  [cond-mat.stat-mech]  3 Oct 2013

<!-- page 2 -->

2
cisely investigated until, about 20 years after the original
work, it was claimed that a variant of the single-bond al-
gorithm was completely free of critical slowing down [16].
Although this was later shown to be incorrect [17], it was
not until recently that its dynamical critical behavior was
investigated in more detail [18, 19], revealing the surpris-
ing feature of critical speeding up, i.e. z < 0, for certain
ranges of q alongside generally smaller dynamical criti-
cal exponents than those found for the Swendsen-Wang
dynamics.
Besides it being a very elegant and direct sampling
procedure for the weights of Eq. (1), another favorable
feature of Sweeny’s approach is its general applicability
to arbitrary values of q: while the Potts model is only
deﬁned for integer q = 1, 2, 3, . . ., the random cluster
model of Eq. (1) is meaningful for any real value q ≥0,
serving as an analytic continuation of the Potts model
to real q [20]. The Swendsen-Wang algorithm, originally
working with a joint spin and bond representation mean-
ingful only for integer q, can be generalized to non-integer
q ≥1 [21]. The bond algorithm, however, is the only
approach for 0 < q < 1. This fact has prompted a num-
ber of researchers to use Sweeny’s approach to probe the
q < 1 regime, for instance to study fractal properties
of the cluster structure [22–24]. The main obstacle to a
more widespread adoption, however, has been the prob-
lem of expensive connectivity checks: inserting an edge
might join two previously unconnected clusters, delet-
ing a bond can lead to cluster fragmentation. A naive
approach without additional data structures appears to
require the tracing out of one (or two) randomly chosen
cluster(s) to check for connectivity. As the average clus-
ter size scales proportional to Lγ/ν [25] and γ/ν ≥1.75
for the random-cluster model, the cost of a full lattice
sweep is almost squared as compared to single spin ﬂips
or Swendsen-Wang. In his paper, Sweeny had suggested
a speciﬁc solution for the case of two-dimensional lat-
tices, replacing the traversal of clusters with a tracking
of boundary loops on the medial lattice [10, 11]. Irrespec-
tive of space dimension, a pair of interleaved breadth-ﬁrst
searches starting from both ends of the bond currently
examined can also dramatically improve the situation
[11, 26].
While these connectivity algorithms still ex-
hibit power-law scaling with the size of the system, fully
dynamic connectivity algorithms, where edge insertions
and removals can be performed in amortized times at
most (poly)logarithmic in the system size, are known in
computer science [27, 28]. Here, we compare a number
of diﬀerent implementations of Sweeny’s algorithm for
simulations of the random-cluster model to each other as
well as to the Chayes-Machta-Swendsen-Wang dynam-
ics [12, 21]. The combination of a polylogarithmic dy-
namic connectivity algorithm and Sweeny’s single-bond
approach is shown to be the more eﬃcient way, asymptot-
ically, to simulate the random-cluster model at criticality.
The rest of the paper is organized as follows. In Sec. II,
we introduce Sweeny’s algorithm in more detail and de-
scribe the three diﬀerent variants of connectivity checks
implemented here: breadth-ﬁrst search, union-and-ﬁnd,
and dynamic connectivities. Section III contains an in-
depth comparison of the scaling of properties of these ap-
proaches as compared to the Chayes-Machta-Swendsen-
Wang dynamics in terms of simulation as well as com-
puter time for the case of simulations on the square lat-
tice. Finally, Sec. IV contains our conclusions.
II.
MODEL AND ALGORITHMS
The random-cluster model (RCM) assigns weights to
(spanning) sub-graph conﬁgurations G′, i.e., subsets of
activated edges and the complete set of vertices, of the
underlying graph G according to [20]
wRC(G′) = qn(G′)vb(G′),
(2)
leading to the partition sum of Eq. (1). For integer val-
ues of the cluster weight q, the partition function (1)
is identical [3] to that of the q-state Potts model with
Hamiltonian
H = −J
X
b∈G
δσi,σj,
(3)
where b = (i, j) is an edge in the graph G, and σi ∈
{1, . . . , q}.
For the purposes of this study, we will re-
strict ourselves to graphs in two dimensions (2D), namely
compact L × L regions of the square lattice, applying
periodic boundary conditions. For this case, the order-
ing transition of the Potts model occurs at the coupling
J/kBT = ln(1 + √q), corresponding to the critical bond
weight vc = √q in (1). This transition is continuous for
q ≤4 and ﬁrst-order for q > 4 [4].
A.
Sweeny’s algorithm
Starting from the results of Fortuin and Kasteleyn [3],
Sweeny suggested to directly sample bond conﬁgurations
of the RCM according to the weight (2). For any sub-
graph G′, the basic update operation is then given by
the deletion of an occupied edge or the insertion of an
unoccupied edge. According to Eq. (2), the correspond-
ing transition probabilities depend on the changes ∆b of
the number of active edges and ∆n of the number of
connected components or clusters. While ∆b is trivially
determined to equal +1 for edge insertion and −1 for
edge removal, respectively, the change in cluster num-
ber depends on whether a chosen inactive edge is inter-
nal to one cluster (∆n = 0) or, instead, it is external
and hence amalgamates two existing clusters if activated
(∆n = −1). Likewise, removing an edge might lead to
∆n = 0 or ∆n = +1, depending on whether an alter-
native path exists connecting the end points of the re-
moved edge.
The construction and implementation of
data structures supporting the eﬃcient calculation of ∆n

<!-- page 3 -->

3
constitutes the intricacy of Sweeny’s algorithm and the
focus of the present work.
Importance sampling for the weight (2) can be con-
structed along well-known lines,
the most common
choices being the heat-bath and Metropolis schemes [29].
In both cases, a bond is randomly and uniformly selected
from the graph and a “ﬂip” of its occupation state from
inactive to active or vice versa is proposed. The heat-
bath acceptance ratio, used in the original approach of
Sweeny [10] is then given by
pHB
acc(∆b, ∆n) =
q∆nv∆b
1 + q∆nv∆b .
(4)
For Metropolis-Hastings, on the other hand, we have
pMH
acc (∆b, ∆n) = min(1, q∆nv∆b).
(5)
It is easily seen that
1
2 ≤pHB
acc
pMH
acc
= max(1, q∆nv∆b)
1 + q∆nv∆b
< 1.
Depending on q and v, we hence expect up to twice larger
acceptance rates for the Metropolis variant. At critical-
ity, v = √q, the minimal ratio pHB
acc/pMH
acc = 1/2 is reached
in the percolation limit q = 1. In contrast to Ref. [10],
our numerical experiments concentrate on Metropolis ac-
ceptance.
As, depending on the data structures used, the deter-
mination of the change ∆n in cluster number is the most
expensive operation, it is economic to only determine ∆n
if it is actually required for the update. In an update at-
tempt, one draws a random number uniformly in [0, 1[;
if r ≤pacc the move is accepted, otherwise it is rejected.
Given that
r ≤min
∆n pacc(∆b, ∆n),
where ∆n ∈{0, −1} for insertion and ∆n ∈{0, +1} for
deletion, respectively, the move can be unconditionally
accepted [16]. Conversely, for
r > max
∆n pacc(∆b, ∆n)
unconditional rejection occurs. At criticality, v = √q,
this results in a fraction
min(√q, 1/√q)
of moves which can be unconditionally accepted or re-
jected under the Metropolis dynamics [30]. Likewise, for
the heat-bath rate (4), a fraction
2
1 + √q min(1, √q)
of move attempts can be decided without actually work-
ing out ∆n.
Note that, in both cases, these fractions
tend to unity as q →1 which is a result of the cluster
weight (2) becoming independent of cluster number in
the uncorrelated percolation limit. Connectivity checks
are hence never required there.
TABLE I. Asymptotic run-time scaling at criticality of the
elementary operations of insertion or deletion of internal or
external edges, respectively, using sequential breadth-ﬁrst
search (SBFS), interleaved BFS (IBFS), union-and-ﬁnd (UF)
or the fully dynamic connectivity algorithm (DC) as a func-
tion of the linear system size L.
move
SBFS
IBFS
UF
DC
internal insertion
LdF −x2
LdF −x2
const.
log2 L
external insertion
Lγ/ν
LdF −x2
const.
log2 L
internal deletion
LdF −x2
LdF −x2
LdF −x2
log2 L
external deletion
Lγ/ν
LdF −x2
Lγ/ν
log2 L
dominant
Lγ/ν
LdF −x2
Lγ/ν
log2 L
B.
Connectivity algorithms
The main complication for an eﬃcient implementation
of the bond algorithm is to maintain the full connectiv-
ity information of the current sub-graph. Consider a ﬂip
attempt on a random edge; it can be currently in the
active or inactive (inserted or deleted) state. For each of
these cases, one needs to distinguish internal from exter-
nal edges, such that independent paths connecting the
two end points either exist (internal edge) or are absent
(external edge).
This leads to the four cases of inter-
nal/external insertions/deletions, each of which can ex-
hibit rather diﬀerent runtime scaling behavior depending
on the chosen implementation.
1.
Breadth-ﬁrst search
The simplest approach to the connectivity problem is
to not maintain any state information about clusters and
determine the value of ∆n for each individual bond move
from direct searches in the graph structure around the
current edge e = (i, j).
Such traversals are most nat-
urally implemented as breadth-ﬁrst search (BFS) [31]
starting from one of the end-points, say i, while not being
allowed to cross the edge e. In case of an external edge,
the cluster attached to i needs to be fully traversed. For
an internal edge, on the other hand, the search starting
at i terminates once it arrives at j, having found an al-
ternative path connecting i and j. Instead of the BFS
one could also use a depth-ﬁrst search (DFS) to achieve
the same result [31]. We found essentially no diﬀerences
in the run-time behavior of both variants, however, and
hence did not consider this possibility in more detail. To
determine the asymptotic run-time at criticality of these
operations, we note that the average number of clusters
with mass s per lattice site is [25]
ns ∼s−τe−cs,
(6)
where τ is the cluster-size or Fisher exponent. A ran-
domly picked site will therefore, on average, belong to a

<!-- page 4 -->

4
cluster of size
M2 =
X
s
s2ns ∼
Z
s2−τe−cs ∝cτ−3 ∼L(τ−3)/σν,
where the last identity follows from c ∼|p −pc|1/σ and
the standard ﬁnite-size scaling ansatz |p −pc| ∼L1/ν.
Since (3 −τ)/σ = γ, we arrive at a typical cluster size
M2 ∼Lγ/ν [25]. For operations on external edges, we
therefore expect an asymptotic scaling of run-times ∼
Lγ/ν. For internal edges, on the other hand, the relevant
eﬀort corresponds to the total number of visited sites of
a breadth-ﬁrst search starting from site i until it reaches
j.
In this case, the number of shells ℓin the BFS at
termination is just the shortest path between i and j.
As shown by Grassberger [32, 33], for bond percolation,
the probability of two nearby points on the lattice to
be connected by a shortest path of length ℓis p(ℓ) ∼
l−ψℓL(ℓ/Ldmin), where
ψℓ= 1 +
2β
νdmin
+
g1
dmin
.
Here, dmin is the shortest-path fractal dimension [34] and
g1 is the scaling exponent related to the density of growth
sites [33]. Ziﬀ[35] demonstrated that g1 = x2 −2β/ν,
where x2 is the two-arm scaling exponent [36, 37]. Hence
ψℓ= 1 + x2/dmin.
As a result, the average length of
shortest path between nearby points exhibits system-size
scaling according to
⟨ℓ⟩∼Ldmin−x2.
(7)
The number of sites touched by a BFS from i to j sep-
arated by a shortest path of length ℓis expected to be
ℓˆd, where ˆd = dF /dmin is known as spreading dimension
[34]. Here, dF = d −β/ν denotes the fractal dimension
of the percolating cluster. Hence, the average number of
sites touched by the BFS for an internal edge is
⟨ℓ
ˆd⟩∼LdF −x2.
(8)
Note that, while dF and x2 are exactly known [25, 38],
this is not the case for dmin [11, 39].
The asymptotic run-time scaling of the Sweeny update
using sequential BFS (SBFS) hence depends on the frac-
tions of internal and external edges encountered. These
are found to be asymptotically L independent numbers
0 < rint, rext < 1 which, however, vary with q [40, 41].
Comparing the scaling exponents for the operations on
internal and external edges, it is found that γ/ν > dF −x2
for the whole range 0 ≤q ≤4, cf. the data compiled in
Table III. As a consequence, the stronger scaling ∼Lγ/ν
of the operations on external edges will always dominate
the running time in the limit of large system sizes.
An improvement suggested in Refs. [11, 26, 42] con-
cerns the quasi simultaneous execution of both BFSs.
In practice, no hardware-level parallelism is needed here
and, instead, sites are removed from the BFS queues
FIG. 1.
(Color online) Equilibrium bond conﬁgurations for
a L = 8 system at criticality, vc = √q.
The left (right)
column shows a conﬁguration for q = 0.0005 (q = 2). The
ﬁrst (second) row corresponds to the level 0 (1) graph. Here,
tree edges are indicated in red, solid lines, while non-tree edges
are drawn with blue, dotted lines.
of the two searches in an alternating fashion, eﬀectively
leading to an interleaved structure of the cluster traver-
sals. To understand the beneﬁt of this modiﬁcation, con-
sider the insertion of an edge e. If it is external, i and
j belong to separate clusters C1 and C2 after the dele-
tion of e. In this case, the searches terminate as soon
as the smaller of the two clusters has been exhausted,
i.e., after C2,min ≡min(|C1|, |C2|) steps. Deng et al. [11]
have shown that, at criticality, this minimum scales as
LdF −x2, i.e., with the exponent already found above for
operations on internal edges in SBFS. For the case of
an internal edge, the interleaved searches terminate as
soon as they meet each other. As argued in Ref. [11] this
time again exhibits the same run-time scaling ∼LdF −x2
which is hence the relevant asymptotic behavior of the
critical bond algorithm using interleaved BFS (IBFS). In
Table I we compare the run-time scaling of the elemen-
tary operations between the diﬀerent implementations of
connectivity checks considered here.
2.
Union-and-ﬁnd
Sequential amalgamations of clusters through the addi-
tion of bonds can be handled eﬃciently using tree-based
data structures under a paradigm known as union-and-
ﬁnd. This is traditionally applied to set partitioning [31],
but has also been used in lattice models for highly eﬃ-
cient simulations of the bond percolation problem [43].
Each cluster is represented as a directed tree of nodes

<!-- page 5 -->

5
with pointers to their parent nodes; the root corresponds
to a designated site representing the cluster as a whole.
In this data structure, connectivity queries are answered
by path traversal to the root sites, such that two sites are
connected if and only if they have the same root. Using
path compression [31], where (most of) the node point-
ers directly link to the cluster root, as well as a balancing
heuristic that attaches the smaller cluster to the root of
the bigger in case of cluster fusion, allows to perform the
connectivity check with a worst-time scaling practically
indistinguishable from a constant [44].
Using this data structure for an implementation of the
bond algorithm for the RCM [46], edge insertion requires
a connectivity check. If the edge is identiﬁed as internal,
the cluster structure remains unchanged on its insertion
which hence can be performed in constant time. For an
external edge, insertion is realized through the attach-
ment of the cluster root of the smaller cluster to the
bigger which is, again, a constant-time operation. For
the deletion of an edge e = (i, j), the information about
alternate paths between i and j is not directly contained
in the data structure. We hence use interleaved BFS to
detect such paths, with a computational eﬀort asymp-
totically proportional to LdF −x2. For an internal edge,
this completes the deletion.
For the case of an exter-
nal edge, leading to fragmentation of the original cluster,
a complete re-labeling of both new clusters is required,
however, resulting in a total scaling of ∼Lγ/ν for this
step, cf. Table I.
The total eﬀective runtime of a bond simulation with
union-and-ﬁnd data structure depends on the frequency
of the individual operation types. The average number
of active bonds at criticality can be worked out from re-
sults for the square-lattice Potts model, where the criti-
cal internal energy density is found to be uc = 1 + 1/√q
[4].
Since ⟨b/N⟩= pu (see, e.g., Ref. [47]), one ﬁnds
⟨b/2N⟩= 1/2 for critical 0 ≤q ≤4. Here, N denotes
the total number of vertices. Thus, the bond occupation
probability of the critical RCM corresponds to the pure
bond percolation threshold, irrespective of q. For ran-
dom bond selection, this results in constant and equal
fractions of insertions and deletions. As mentioned above
in the context of the BFS technique, the fractions of in-
ternal versus external edges are diﬀerent from zero for
all values of q. Hence, it is the most expensive operation
which dominates the asymptotic scaling behavior, and we
hence expect ∼Lγ/ν scaling for the union-and-ﬁnd im-
plementation, although all insertion moves are performed
in constant time.
3.
Dynamic connectivity algorithm
While union-and-ﬁnd uses data structures that allow
for insertions and connectivity checks in constant time,
there exist modiﬁed data structures for which also edge
deletions are supported without the need of an expensive
∼Lγ/ν rebuild operation in case of an external edge.
A number of such “fully dynamic” graph algorithms has
been discussed in the computer-science literature [27, 28].
In the following, we refer to such approaches as dynamic
connectivity (DC) algorithms.
The advantage in run-
time for the deletion of edges is paid for in terms of in-
creased eﬀorts for edge insertion. We use the approach
suggested in Ref. [28] which is deterministic and features
amortized runtimes of O(log N) for connectivity queries
and O(log2 N) for deletions and insertions on graphs of
N nodes. The time complexity for connectivity queries
depends on the underlying binary search tree used to
encode the graphs. In our case we used splay trees [48]
which result in the amortized bound. An identical worst-
case bound holds for balanced binary search trees [28].
The algorithm of Holm et al. [28] is based on a re-
duction of the set of edges to be considered by focusing
on a spanning forest F(G) of a given graph G, which
encodes the same connectivity equivalence relation but
has less edges as no cycles occur [49]. This separates the
set of edges into “tree” edges e ∈F(G) and “non-tree”
edges e /∈F(G). Given this separation, one then stores
the spanning trees of all components and augments them
with information about incident non-tree edges.
This
edge separation already reduces the number of expensive
operations, as any edge e /∈F(G) must be an internal
edge, such that ∆n = 0. The only potentially expen-
sive cases remaining are the deletion of e ∈F(G) and
the insertion of an external edge into F(G). For the ﬁrst
case, one notes that the fact that the spanning tree is
fragmented by the removal of e ∈F(G) does not imply
that also the cluster on the original graph is split by this
operation as there might be non-tree edges still connect-
ing the parts. The concept of tree- and non-tree edges is
visualized in Fig. 1.
In order to eﬃciently search for a replacement edge in
the set of all non-tree edges, another edge separation is
introduced. This is achieved by associating a level in the
range 0, 1, · · · , lmax ≡⌊log N⌋to every edge (both tree
and non-tree).
Given the levels of all edges, one then
maintains a spanning forest Fi for i = 0, 1, · · · , lmax of
all edges with l(e) ≥i, thus F = F0 ⊇F1 · · · ⊇Flmax,
cf. the illustration in Fig. 1. An important part of the
algorithm is that the levels are not ﬁxed but they are
partially changed after every tree-edge deletion to ensure
the following two properties [28]:
1. The maximum cluster size at level i is ⌊N/2i⌋. This
implies lmax ≡⌊log N⌋.
2. If an edge at level i is deleted, then a possible re-
placement edge can only be found in levels l ≤i.
To eﬃciently save and manipulate spanning trees we rep-
resent each tree edge (i, j) by two directed edges (≡arcs)
(i →j) and (j →i) and every vertex i by a loop (i →i).
Given all directed edges we construct one of possibly
many Euler tours of this directed graph, i.e., a cycle that
traverses every edge exactly once. We split the tour at
one arbitrary point to save it as a sequence of traversed

<!-- page 6 -->

6
TABLE II.
Estimated dynamical critical exponents zint,O for the two-dimensional RCM at criticality and O = S2 as well
as O = N for a range of q values as compared to results reported in Refs. [15, 19, 45]. The values shown for zSW
int,S2 for the
Swendsen-Wang-Chayes-Machta algorithm are actually related to another observable E, but it was reported in Refs. [15, 45]
that the observables E′, N and S2 share the same dynamical critical exponent for this algorithm.
q
zint,S2
zint,S2 [19]
zSW
int,S2 [15, 45]
zint,N
α/ν
zint,N [19]
zSW
int,S2 −zint,S2
0.0005
−1.12(1)
-1.23
—
0.01(1)
-1.958
0
—
0.005
−1.09(1)
-1.21
—
0.01(1)
-1.868
0
—
0.05
−1.04(1)
-1.12
—
0.01(1)
-1.601
0
—
0.2
−0.86(1)
-1.01
—
−0.01(1)
-1.247
0
—
0.5
−0.63(1)
-0.71
—
0.00(1)
-0.878
0
—
1.0
−0.33(1)
-0.32
—
0.00(1)
-0.500
0
—
1.5
−0.11(2)
-0.16
0
0.06(2)
-0.227
0
2.0
0.03(3)
-0.08
0.143(3)
0.13(2)
0
0 (log)
0.11
3.0
0.44(4)
0.41
0.497(3)
0.45(4)
0.400
0.45(1)
0.06
4.0
0.75(6)
—
0.910(5)
0.73(6)
1
—
0.16
arcs and loops [50]. The fact that we now linearized the
tree by mapping it to a list or sequence allows to eﬃ-
ciently store it in a balanced binary search tree. Edge in-
sertions and deletions then translate into manipulations
using cuts and links on Euler tours.
By using a form
of self-adjusting binary search trees named “splay trees”
[48] we were able to do all the operations on trees in amor-
tized runtime of O(log N). Other types of trees might be
used alternatively, see e.g., [31, 51, 52]. Intuitively speak-
ing, we hence have O(log N) levels of edges with runtime
of O(log N) per level, resulting in the quoted O(log2 N)
amortized runtime for deletions and insertions.
A de-
tailed comparison and run-time analysis of diﬀerent DC
algorithms can be found in Ref. [53]. For simulations of
the RCM, we can therefore perform each operation in
a run-time asymptotically proportional to log2 L, hence
clearly outperforming the other approaches at criticality,
cf. Tab. I.
4.
Behavior oﬀcriticality
The run-time bounds summarized in Table I apply
to simulations of the RCM at criticality.
In the high-
temperature, non-percolating regime p < pc all clusters
are ﬁnite. Hence, the implementations based on BFS and
UF allow to perform insertions and deletions in asymp-
totically constant time there. For temperatures below the
transition, equivalent to p > pc, the behavior is a bit more
complicated. In this case clusters become compact, so
cluster masses scale proportional to L2. For SBFS, oper-
ations on external edges require complete cluster traver-
sal, leading to L2 run-time scaling there. Through the
dense nature of clusters, a re-connecting path replacing
an internal edge will, with probability one, be only of
ﬁnite length independent of L.
Similarly, if an exter-
nal edge connects two distinct clusters, at least one of
them will be small, i.e., of L independent size. Hence,
a constant run-time per operation is expected for IBFS
and p > pc. In the UF implementation again complete
cluster traversal is necessary in some cases, giving an L2
bound. The asymptotic scaling of DC equals log2 L, in-
dependent of p, although the absolute run-times will of
course depend on temperature.
III.
SIMULATIONS AND RESULTS
To gauge the eﬃciency of the connectivity implemen-
tations and conﬁrm the validity of the asymptotic analy-
sis presented above, we implemented Sweeny simulation
codes based on these diﬀerent approaches and subjected
them to a careful analysis of run-times as a function of
system size and q.
A.
Autocorrelation times and eﬃciency
Besides the run-time for individual bond operations
discussed above in Sec. II, the eﬃciency of the bond al-
gorithm is ultimately determined by the speed of decor-
relation of the Markov chain, i.e., the autocorrelation
times. Consider the autocorrelation function of the mea-
surements Os at time lag t,
ΓO(t) = ⟨OsOs+t⟩−⟨O⟩2
(9)
which, in equilibrium, is expected to be independent of
the initial time s due to stationarity. From the theory of
Markov chains, ΓO(t) is expected to decay exponentially
for large time lags, ΓO(t) ∼e−|t|/τexp,O, and one deﬁnes
the exponential autocorrelation time [54]
τexp,O = lim sup
t→±∞
|t|
−log |ρO(t)|,
(10)

<!-- page 7 -->

7
0
1
2
3
4
q
−2.0
−1.5
−1.0
−0.5
0.0
0.5
1.0
zint
α/ν
N
S2
C1
FIG. 2. (Color online) Estimates of the 2D dynamical critical
exponents zint,O for O = N, S2 and C1, respectively, as a
function of q.
The line corresponds to the Li-Sokal bound
zint > α/ν [56].
where ρO(t) = ΓO(t)/ΓO(0) is the normalized autocor-
relation function.
Unless the considered observable is
“orthogonal” to the slowest mode, one expects a result
independent of O, i.e., τexp ≈τexp,O.
The eﬃciency of sampling, on the other hand, is de-
termined by the integrated autocorrelation time,
τint,O = 1
2 +
M
X
t=1
ρO(t),
(11)
where M is the length of the time series. This is seen by
considering the variance of the average ¯O = (P
i Oi)/M
used as an estimator for ⟨O⟩,
σ2
¯
O ≈
σ2
O
M/2τint,O
.
(12)
Comparing this to the case of uncorrelated measure-
ments, where σ2
¯
O = σ2
O/M, it is seen that the eﬀective
number of independent measurements is reduced by a
factor of 1/2τint,O by the presence of autocorrelations.
In contrast to τexp, the integrated autocorrelation time
generically depends on the observable considered.
Close to criticality, dynamical scaling implies autocor-
relation times diverging according to τexp ∼ξzexp and
τint,O ∼ξzint,O. For ﬁnite systems, ξ eventually becomes
limited by L, resulting in the scaling form
τ ∼Lz.
(13)
For a purely exponential decay, exponential and inte-
grated autocorrelation times coincide.
More generally,
however, one has τint,O ≤τexp and hence zint,O ≤zexp.
Cases of a true inequality have been observed [55].
A number of diﬀerent techniques for the practical es-
timation of autocorrelation times have been discussed in
the literature. We use a direct summation of the auto-
correlation function according to Eq. (11). Since, for any
time lag, the estimated autocorrelation function adds a
constant amount of noise per term, a cutoﬀΛ needs to
be introduced to ensure a ﬁnite variance of the estima-
tor for the autocorrelation time [57]. A technique with
an adaptive summation window has originally been sug-
gested in Ref. [55].
We use a modiﬁcation of this ap-
proach discussed more recently in Ref. [58]. Determining
Λ involves a tradeoﬀbetween the bias for small cutoﬀs
and the exploding variance as Λ →M. This balance is
struck here by numerically minimizing the quantity
exp (−Λ/τ) + 2
p
Λ/M,
(14)
which is proportional to the sum of the relative statistical
and systematic errors. We ﬁnd this procedure to yield
stable results throughout. In particular, the estimates of
τint are consistent with those found from an alternative
jackkniﬁng technique [59, 60].
To achieve an appropriate judgment of Sweeny’s al-
gorithm in the diﬀerent implementations against other
approaches of simulating the RCM, we need to combine
the information contained in the autocorrelation times
with those of the run-times for the diﬀerent operations
involved.
We hence compare the eﬀective run-time to
create a statistically independent sample of a given ob-
servable, where independence is understood in the sense
of Eq. (12). To this end, we consider the eﬀective runtime
per edge,
TO ≡τint,O¯t,
(15)
where τint,O is measured in sweeps and ¯t is the runtime
per edge operation, averaged over a suﬃciently long sim-
ulation. Obviously, this quantity is hardware speciﬁc, yet
by looking at the ratio of two diﬀerent implementations
we expect the speciﬁc hardware-dependence to be small
[61].
B.
Observables and implementation
As τint,O depends on O, any conclusions about the
relative eﬃciency of the Sweeny and Swendsen-Wang-
Chayes-Machta dynamics for the RCM might depend on
the observable under consideration. It is therefore im-
portant to study a number of diﬀerent quantities cover-
ing the energetic and magnetic sectors of the model. For
the energetic sector, we studied the number N ≡b of
active bonds. In the magnetic sector, we considered the
sum of squares of cluster sizes S2 ≡Pn
i=1 |Ci|2, and the
size of the largest component C1 ≡maxn
i=1 |Ci|, where
Ci denotes the ith cluster resulting from the set of active
bonds. These observables of the RCM are related to stan-
dard observables of the Potts model [15, 47], namely the
internal energy per spin u = ⟨N⟩/(Ldp), the susceptibil-
ity χ = ⟨S2⟩/Ld, and the order parameter m = ⟨C1⟩/Ld.
We implemented a code for the Sweeny update us-
ing the four variants of connectivity algorithms discussed
above in Sec. II, i.e., sequential breadth-ﬁrst search

<!-- page 8 -->

8
101
102
103
L
10−4
10−3
10−2
10−1
100
101
102
103
104
105
τint,S2
q = 0.0005
q = 0.5
q = 1
q = 2
q = 3
q = 4
FIG. 3. (Color online) Integrated autocorrelation times as a
function of system size for the susceptibility S2 for q ranging
from 0.0005 up to q = 4. The lines show ﬁts of the functional
form (13) to the data.
(SBFS), interleaved BFS (IBFS), union-and-ﬁnd (UF),
as well as dynamic connectivities (DC). Our program was
implemented in C and compiled with the GNU Compiler
Collection (GCC) 4.5.0 at -O2 optimization level on an
Intel Xeon E4530 2.66 GHz. The available memory was
16 GB. These speciﬁcations resulted in speed ups through
caching eﬀects – mainly for the DC implementation, see
the discussion below – for system sizes L ≤48, which
approximately corresponds to a memory footprint of the
size of the L2 cache of 6 MB for the DC implementation.
The implementations based on BFS and union-and-
ﬁnd both require memory scaling as O(L2).
The DC
code, on the other hand, requires to maintain O(log L)
separate forests, leading to a total memory requirement
of O(L2 log L) [53].
In practise for a system of size
L = 1024 and q = 2 at criticality, the BFS and UF
implementations required around 10–15 MByte, whereas
the DC code used approximately 3.5 GByte.
The random number stream was generated by a GSL
implementation [62] of the Mersenne twister or MT19937
generator [63] which has a very large period of 219937 ≈
106000 and, more importantly, was shown to be equi-
distributed in less then 623 dimensions.
C.
Simulation results
In order to test the predicted asymptotic run-time be-
havior and determine the dynamic critical behavior of
Sweeny’s algorithm, we performed a series of simulations
of the RCM on the square lattice. To cover the complete
range of systems with continuous phase transitions, we
chose q = 0.0005, 0.005, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9,
1, 1.5, 2, 3 and 4. Simulations were performed for a series
of 14 system sizes ranging from L = 4 up to L = 1024.
All simulations were done at the exact critical coupling
vc = √q [4]. While, typically, observable measurements
0
20
40
60
80
100
t
10−3
10−2
10−1
100
ρC1(t)
q = 0.005
q = 0.05
q = 0.5
q = 2.0
FIG. 4.
(Color online) Normalized critical autocorrelation
function for the order parameter C1 and several values of q.
The linear system size is L = 64. Time is counted in units
of L single bond moves here, such that L = 64 time steps
correspond to a lattice sweep.
in Markov chain Monte Carlo are taken after every lat-
tice sweep of updates [29], for the case of the bond al-
gorithm and relatively small values of q, it turned out
that the fast decorrelation leads to autocorrelation times
way below a single sweep. To resolve these eﬀects, we
hence changed the measurement interval to multiples of
L, which turned out to be a good compromise between
the visibility of correlations and the resulting lengths of
time series. This setup resulted in at least 220 and up
to 223 measurements for some system sizes. As a result,
our simulations were at least 103, and at best 105 times
the relevant time scale τexp. Consistent with Sweeny’s
original observation, we were able to equilibrate all runs
in less then 200 sweeps.
1.
Dynamical critical behavior
We ﬁrst considered the behavior of the energy-like
observable N, determining (integrated) autocorrelation
times according to Ref. [58] with resulting cutoﬀparam-
eter ΛN (see the discussion in Sec. III A above). In or-
der to extract the corresponding dynamical critical ex-
ponents we ﬁtted a power-law of the form τint ∼Lzint
to the data, omitting some of the smallest system sizes
to account for scaling corrections. The quoted errors on
ﬁt parameters correspond to an interval of one standard
deviation. The resulting estimates of zint,N are shown in
Fig. 2 and the numerical values are summarized in Table
II. Note that zint,N ≈0 for q ≤2.
This is in agree-
ment with a key result due to Li and Sokal, providing a
lower bound for the autocorrelation times of N and its
corresponding dynamical exponents [15, 56]
τexp,N >∼τint,N ≥Cv ⇒zexp,N ≥zint,N ≥α/ν,
(16)
where Cv is the speciﬁc heat and α/ν the associated
ﬁnite-size scaling exponent. While this result was orig-

<!-- page 9 -->

9
101
102
L
101
102
103
Cmin,2(L)
q = 0.0005 : y(e)
IBFS = 1.24(1)
q = 1 : y(e)
IBFS = 0.67(1)
q = 2 : y(e)
IBFS = 0.47(2)
FIG. 5.
(Color online) Number of vertices Cmin,2 of the
smaller of two adjacent clusters as visited in a 2D simulation
using IBFS following (preceding) the deletion (insertion) of
an external edge. The lines show ﬁts of the power-law form
Cmin,2 = ALy(e)
IBFS to the data.
inally derived for the Swendsen-Wang dynamics, it was
also shown to hold for the Sweeny algorithm [19, 56].
As α/ν ≤0 for q ≤2 [4], our data are consistent with
this bound and indicate that it is close to being tight
for the Sweeny dynamics on the square lattice for q ≥2,
cf. Fig. 2. The values for q = 4 appear to violate this
bound, but we attribute these deviations to the logarith-
mic corrections expected for this particular value of q,
preventing us from seeing the truly asymptotic behavior
in the regime of system sizes considered here.
We next turned to the magnetic observables S2 and
C1. Deng et al. [19] ﬁrst showed that under the Sweeny
dynamics the susceptibility S2 exhibits a surprisingly fast
decorrelation on short timescales and, in particular, the
corresponding integrated autocorrelation time τint,S2 de-
creases with increasing system size, indicating a negative
value of the corresponding dynamical critical exponent
zint,S2. This phenomenon of critical speeding up is also
clearly seen in our data as is illustrated in the plot of
the L dependence of the autocorrelation times of S2 in
Fig. 3. A rather similar behavior is found for the order
parameter C1. The initial fast decay of correlations is
illustrated in Fig. 4, where it is clearly seen that C1 is
completely decorrelated in less than a single sweep for
q <∼2. Note that the measurements along the Markov
chain are still correlated, but with a system size scaling
weaker than L2 such that the impression of a complete
decorrelation appears on the scale of sweeps.
The values of zint,S2 resulting from power-law ﬁts to
the autocorrelation times are compiled in Table II. As is
clearly seen from the plot of the data in Fig. 3, we ﬁnd
zint,S2 ≤0 for q ≤2. Regarding the dynamical critical
exponents in the magnetic sector, we ﬁnd zint,C1 = zint,S2
within our error bars. The authors of Ref. [19] suggested
to determine zint,S2 from a data collapsing procedure us-
ing a two-time scaling ansatz combining the fast initial
0
500
1000
1500
2000
L
0
10
20
30
40
50
60
70
80
¯t in µs
q = 0.0005
q = 2
q = 4
FIG. 6. (Color online) Average run-time ¯t for several values of
q and L for the Sweeny update using a dynamic connectivity
(DC) algorithm based on splay trees. The lines correspond to
least-squares ﬁts of the model (18) to the data.
decay with a slower exponential mode for longer times.
In our studies, however, we found this approach to yield
rather unstable results and thus chose, instead, to per-
form more conventional ﬁts to a power law.
We note that, in line with Ref. [15] we used the summa-
tion cutoﬀΛN of the observable N for all observables, be-
cause the magnetic observables S2 and C1 with their fast
initial decay would lead to very small, sub-asymptotic
cutoﬀs if the rule of Ref. [58] would be directly applied.
We also checked that the estimators for τint,S2 and τint,C1
were on a plateau so that a change in the summation
window mainly inﬂuences the variance of the estimator,
which is monotonically increasing with Λ.
Comparing our results for the Sweeny dynamics to the
Swendsen-Wang algorithm, we note that, apart from the
slightly smaller dynamical critical exponent for the for-
mer, we also ﬁnd somewhat smaller amplitudes in the
τint = ALzint scaling for the bond algorithm.
Hence
for L = 256 we ﬁnd, e.g., τint,S2 ≈0.1 (q = 1) and
≈10 (q = 3) for the Sweeny update, while values of
τint,S2 ≈0.5 (q = 1) and ≈36 (q = 3) are found for the
Swendsen-Wang update.
2.
Run-time scaling
As discussed above, the relevant time scales for a com-
parison of the bond algorithm against other approaches
depend on both, the statistical decorrelation as well as
the run-time scaling of the elementary operations. We
therefore analyzed average run-times for bond updates
in the Sweeny algorithm with the diﬀerent implementa-
tions of connectivity updates discussed in Sec. II B.
For the techniques based on BFS, we studied the num-
ber of steps required to complete a connectivity check
for the case of operations on internal and external edges,
respectively, for the SBFS and SBFS implementations.

<!-- page 10 -->

10
TABLE III. Run-time scaling exponents in two dimensions according to Eq. (17) for the SBFS, IBFS and the UF implementa-
tion. The scaling exponents γ/ν and dF −x2 are shown for reference and comparison. The y exponents correspond to scaling
∼Ly of the number of vertices touched in a sequential (SBFS) and interleaved (IBFS) breadth-ﬁrst cluster traversal for internal
(i) and external (e) edges, respectively.
q
κSBFS
κIBFS
κUF
y(i)
IBFS
y(e)
IBFS
y(i)
SBFS
y(e)
SBFS
dF −x2
γ/ν
0.0005
1.81(1)
1.18(3)
1.54(6)
1.25(1)
1.24(1)
1.25(1)
1.99(1)
1.23407
1.99296
0.005
1.80(2)
1.14(2)
1.70(6)
1.22(1)
1.21(1)
1.22(1)
1.98(1)
1.20021
1.97823
0.05
1.74(1)
1.02(3)
1.77(4)
1.10(1)
1.11(1)
1.11(1)
1.93(1)
1.09783
1.93580
0.1
1.72(1)
0.97(2)
1.78(3)
1.05(1)
1.05(1)
1.06(1)
1.91(1)
1.03881
1.91284
0.5
1.65(2)
0.71(3)
1.77(2)
0.82(1)
0.82(1)
0.82(1)
1.83(1)
0.80768
1.83449
0.7
1.63(2)
0.69(4)
1.77(4)
0.75(2)
0.76(1)
0.76(2)
1.81(1)
0.73541
1.81407
1
0
0
1.74(2)
0.66(1)
0.67(1)
0.67(1)
1.79(1)
0.64583
1.79167
1.5
1.56(2)
0.43(2)
1.71(3)
0.55(2)
0.56(1)
0.57(2)
1.75(1)
0.52298
1.76644
2
1.57(2)
0.35(3)
1.68(3)
0.46(2)
0.47(2)
0.48(3)
1.73(1)
0.41667
1.75000
3
1.52(4)
0.19(1)
1.67(2)
0.32(3)
0.30(4)
0.35(3)
1.69(2)
0.21667
1.73333
4
1.42(4)
0.13(2)
1.64(7)
0.22(11)
0.23(1)
0.26(1)
1.68(1)
-0.12500
1.75000
For internal edges, this corresponds to the number of
vertices touched by the BFSs until a re-connecting path
is found. For external edges such a path is not found and
the search hence terminates after a number of steps corre-
sponding to the mass of either the ﬁrst cluster (SBFS) or
the smaller cluster (IBFS). Checking the number of steps
for operations on internal and external edges for SBFS
and IBFS, respectively, we used power-law ﬁts according
to ∼Ly to extract estimates of the four exponents y(i)
SBFS,
y(e)
SBFS, y(i)
IBFS, and y(e)
IBFS. The ﬁt results are collected in
Table III. The data and corresponding ﬁts for the case of
the number of steps C2,min relevant for the operation on
an external edge with IBFS are shown in Fig. 5. The ex-
ponents y follow the asymptotic values γ/ν and dF −x2,
respectively, derived above in Sec. II B 1 and also listed
in Table III for comparison.
For the total average run-time per bond operation, we
asymptotically expect power-law behavior as well,
¯t ∼Lκ.
(17)
This assumption in general describes well our data —
with only minor deviations for smaller system sizes due to
caching eﬀects. For SBFS, we expect the diﬀerent asymp-
totic scaling behavior for operations on internal and ex-
ternal edges, respectively, to result in an eﬀective run-
time exponent κ somewhere in between the exponents
dF −x2 and γ/ν relevant to operations on internal and
external edges, respectively (recall that internal and ex-
ternal edges occur in constant fractions). Our estimates
of κSBFS listed in Table III are in line with these expec-
tations. We have no doubt, however, that the asymp-
totically expected κSBFS = γ/ν ultimately holds for suf-
ﬁciently large systems. For the interleaved case, on the
other hand, all four operation types exhibit y = dF −x2,
and we hence ﬁnd κIBFS consistent with dF −x2 already
for the system sizes considered here, cf. Table III.
The analysis of the run-time behavior for the union-
and-ﬁnd approach is more subtle as the insertion of edges
is performed in constant time, whereas the deletion of
edges incurs an eﬀort proportional to LdF −x2 and Lγ/ν
for internal and external edges, respectively, cf. Table I.
As a consequence of the diﬀerent scaling of individual
operations, the eﬀective run-time scaling exponent κUF
according to Eq. (17) is again found to be smaller than
the expected limiting value γ/ν, see the values compiled
in Table III.
The scaling of run-times per step for our implemen-
tation of the dynamic connectivity algorithm and a rep-
resentative selection of q-values is shown in Fig. 6. We
ﬁnd a sub-algebraic growth and, according to the asymp-
totic run-time bounds derived in Ref. [28], we ﬁtted the
functional form
¯t(L) = a log2 L + b log L + c
(18)
to the data. The ﬁts resulted in c ≈0 such that we ﬁxed
c = 0 in the following. Somewhat surprisingly, our ﬁts
yield b < 0; we interpret this as a result of the presence
of correction terms and the amortized nature of the run-
time bounds leading to the asymptotic scaling only being
visible for very large system sizes. Similar observations
have been reported for general sets of inputs in Ref. [53].
Considering the ratio a/b, we ﬁnd that its modulus in-
creases with q, yielding a value of ≈0.3 for q = 0.0005
and ≈0.71 for q = 2. This corresponds to the increasing
fraction of non-tree edges for increasing q, resulting in an
increase of traversals of the edge level hierarchy with the
associated O(log2 L) complexity. Irrespective of that, as
a consequence of the larger number of cluster-splitting
operations the total run-time is found to be largest for
small q, cf. Fig. 6.
We also investigated the eﬀect of the unconditional ac-
ceptance of proposed updates for the Metropolis rule as
discussed in Sec. II A above.
This adds another q de-

<!-- page 11 -->

11
101
102
103
L
10−1
100
101
102
103
104
T/TSW
q = 1.0
DC
IBFS
UF
101
102
103
L
100
101
102
103
T/TSW
q = 2.0
DC
IBFS
UF
101
102
103
L
100
101
102
103
T/TSW
q = 3.0
DC
IBFS
UF
101
102
103
L
100
101
102
103
T/TSW
q = 4.0
DC
IBFS
UF
FIG. 7. (Color online) Eﬀective run-times T according to Eq. (15) for the diﬀerent implementations relative to the time TSW
of the Swendsen-Wang algorithm [12]. Dashed lines correspond to the run-time to generate an independent sample of the
observable N and dotted lines to samples of S2.
pendent, but system-size independent, element to the
run-time scaling. Such unconditional moves can save sig-
niﬁcant computational eﬀort in case no data structures
have to be updated after move acceptance. This is the
case for the algorithms based on BFS which are “state-
less” in the sense that no explicit record of connectivity
is kept. Unconditional insertion or removal of edges in
this case does not entail any further computational ef-
fort. On the contrary, unconditional insertion or removal
lead to data-structure updates for the union-and-ﬁnd and
DC implementations. As a consequence, we ﬁnd a con-
stant speed-up for the BFS based implementations pro-
portional to 1/[1 −min(√q, 1/√q)]. In the singular case
q = 1, BFS performs all edge updates in constant time
as all insertions and deletions can be performed uncon-
ditionally such that no cluster traversals are necessary.
On the contrary, no performance improvement from un-
conditional moves is observed for the more elaborate UF
and DC implementations.
We note that, for all implementations, the average run-
time per bond operation depends quite strongly on q.
This is, on the one hand, due to the q dependence of the
fraction rext of external edges reaching from rext = 1 for
q →0 down to rext = 0.33 for q = 4. For the case of
the BFS implementations, an additional q dependence is
introduced through the unconditional moves as discussed
above.
3.
Overall eﬃciency
As discussed above in Sec. III A, the relevant measure
for the overall eﬃciency of various implementations of
cluster algorithms is the total run-time for the genera-
tion of a statistically independent sample according to
Eq. (15).
We compared the eﬀective run-times of all
three implementations of the bond algorithm with a ref-
erence code for the Swendsen-Wang dynamics. As the
dynamical critical exponents zint for the Sweeny update
are found to be smaller than those of the Swendsen-
Wang-Chayes-Machta dynamics, see Refs. [15, 19, 45]
and Table II, we can expect asymptotically more eﬃ-
cient simulations for cases where the run-time exponent
κ < zSW
int,O−zint,O. Since any poly-logarithm is dominated
by Lϵ with ϵ > 0, this is clearly the case, asymptotically,
for the DC algorithm. From the data for zSW
int,O −zint,O
in Table II and those for κ in Table III, it is seen that
for 1 < q ≤4, this condition is not met for the imple-

<!-- page 12 -->

12
mentation based on union-and-ﬁnd. For the technique
based on (interleaved) breadth-ﬁrst-search, on the other
hand, such a case arises (for integer values q) only for
q = 4, where zSW
int,O −zint,O ≈0.16 and κIBFS = 0.13(2).
These observations are corroborated by the plots of the
relative eﬃciencies shown in Fig. 7. For comparison, we
here show the results for the two observables N and S2
with signiﬁcantly diﬀerent behavior for q ≤2. From the
plots for the integer values q = 2, 3 and 4, it is clear that
in absolute run-times IBFS is most eﬃcient for the range
of system sizes L ≤1024 considered here. Hence, the
asymptotic advantage of the DC algorithm only shows
for system sizes beyond this range.
The downturn of
the ratio T/TSW for the largest system sizes and q = 4
observed for the BFS and DC codes might be an indi-
cation of the asymptotic run-time advantage of Sweeny’s
algorithm over Swendsen-Wang with these connectivity
algorithms as discussed above. The results for the perco-
lation case q = 1 using IBFS, on the other hand, are of
exceptional nature as there the cost of bond operations
is completely independent of system size due to the ef-
fect of unconditional acceptance. For the case of S2, one
even ﬁnds a decrease of the relative cost of generating a
statistically independent sample as this observable prof-
its from the initial fast decorrelation or critical speeding
up.
For the most relevant case q ≤1, where Sweeny’s algo-
rithm provides the only means of simulation, we cannot
compare to another algorithm.
Instead, we present in
Fig. 8 a comparison of run times for the SBFS, IBFS, UF
and DC implementations for q = 0.005. As here κIBFS
is relatively unfavorable, we observe a clear advantage
for the DC algorithm which is signiﬁcantly more eﬃcient
than the other options in the full range of studied system
sizes 4 ≤L ≤1024. UF is found to be even less eﬃcient
than SBFS here which might be considered surprising
in view of the fact that all insertions are performed at
constant cost and deletions have the same asymptotic
run-time bounds as SBFS, see Table I. This is easily un-
derstood, however, noting that the factor of two gained
for UF from the 50% of operations performed in constant
time is spent again in having to traverse both clusters
fully in case of external edge deletions. Taking into ac-
count overheads for data-structure updates for UF, this
explains the slight disadvantage of UF over SBFS seen in
Fig. 8.
IV.
CONCLUSIONS
We have shown that it is possible to implement
Sweeny’s algorithm eﬃciently and in a lattice and dimen-
sionality independent way, using a dynamic connectivity
(DC) algorithm, in the sense that the runtime depen-
dence on the system size is poly-logarithmic and only con-
tributes a correction to the statistical Lzint contribution
of the runtime to create an eﬀectively uncorrelated sam-
ple.
Compared to the Swendsen-Wang-Chayes-Machta
101
102
103
L
100
101
102
103
¯t in µs
q = 0.005
UF
SBFS
IBFS
DC
FIG. 8.
(Color online) Run-time per edge operation of sim-
ulations of the q = 0.005 square-lattice RCM and the bond
algorithm employing the SBFS, IBFS, UF and DC connectiv-
ity implementations, respectively.
algorithm, we also ﬁnd somewhat smaller dynamical crit-
ical exponents, leading to an overall asymptotically more
eﬃcient simulation of the random-cluster or Potts model
with Sweeny’s approach. In addition, the bond algorithm
is the only known approach for simulations in the regime
q < 1 including, for instance, interesting q →0 limits
such as the maximally connected spanning sub-graphs of
Ref. [64].
We analyze in detail four implementations based on
(sequential and interleaved) breadth-ﬁrst searches, on
union-and-ﬁnd data structures, and on the fully dy-
namic connectivity algorithm suggested in Ref. [28], re-
spectively. For each implementation, we derive average
run-time bounds for insertions and deletions of internal
and external edges, respectively, and deduce the over-
all asymptotic run-time behavior. It is found that inter-
leaved breadth-ﬁrst searches, although relatively unfavor-
able as compared to union-and-ﬁnd and dynamic connec-
tivities at ﬁrst sight, perform rather well due to the lack
of an underlying data structure encoding the connectiv-
ity of the clusters, in particular if connectivity queries are
omitted whenever possible due to accepting moves for
which the drawn random number indicates acceptance
irrespective of the result of the connectivity query. The
union-and-ﬁnd based implementation, on the other hand,
although superior in asymptotic run-time in three out of
four cases of insertions or deletions of internal or exter-
nal edges, shows ultimately inferior performance due to
the run-time scaling for deletions of external edges that
require full traversals of the involved clusters. The dy-
namic connectivity algorithm, while asymptotically most
eﬃcient with a poly-logarithmic scaling of run-times per
operation, has rather large constants leading to some-
what weaker performance than breadth-ﬁrst search for
the considered lattice sizes L ≤1024 and q ≥1. For
q ≪1, on the other hand, where run-times are dominated
by operations on external edges, it outperforms the other

<!-- page 13 -->

13
implementations already for small systems. We see sig-
niﬁcant room for further run-time improvements for the
dynamic connectivity algorithm, however, for instance
by optimizations of the underlying tree data structure
or the implementation of additional heuristics as indi-
cated in the comparative study [53]. We note that due
to the lack of explicit connectivity information for the
breadth-ﬁrst search approach it becomes more expensive
than for the other techniques to perform measurements
of quantities such as cluster numbers or correlation func-
tions depending on the connectivity. As measurements
are typically taken at most once per sweep, however, any
cost of at most L2 operations for measurements results
in only O(1) amortized eﬀort per bond operation.
The observed fast initial decorrelation for q ≤2 and
quantities such as S2 and C1 depending on cluster con-
nectivity as illustrated in Figs. 2, 3 and 4 suggests that
there is an additional dynamical mechanism at play for
such observables. As argued in Ref. [19], this is due to a
larger number of operations on external edges for smaller
values of q which can lead to a large-scale change in the
connectivity structure through a single bond operation.
The concentration of external edges, bridges or fragment-
ing bonds drastically increases as q is decreased from 4
down to the tree limit q →0 [41] which is also clearly
expressed in a corresponding increase in the fractal di-
mension of “red” bonds [65]. It is currently not clear,
however, why this eﬀect only leads to a change in dy-
namical critical behavior for q ≤2.
While we have restricted our attention to simulations
on the square lattice, all implementations discussed here
are essentially independent of the underlying graph or
lattice, requiring only minimal adaptations for diﬀerent
situations. This makes our approach signiﬁcantly more
general than the implementation originally suggested by
Sweeny [10] which is based on tracing loops on the me-
dial lattice in two dimensions. Additionally, the latter
technique in two dimensions still suﬀers from polyno-
mial scaling of the run-time per edge operation [11], such
that a poly-logarithmic implementation is asymptotically
faster.
Until now we have only used the DC implementation
for canonical simulations but as proposed in Ref. [26]
an interesting application are generalized-ensemble sim-
ulations of the random cluster model where one directly
estimates the number of possible graphs g(k, b) with k
clusters and b edges which then allows for the calculation
of canonical ensemble averages as continuous functions
of temperature and the parameter q.
The source code of the implementations discussed here,
in particular including the dynamic connectivity algo-
rithm, is available on GitHub under a permissive license
[66].
ACKNOWLEDGMENTS
E.M.E. would like to thank T. Platini and N. Fytas for
fruitful discussions and U. Wolﬀfor providing an imple-
mentation of his automatic windowing method.
[1] M. E. Fisher, Physics 3, 255 (1967).
[2] A. Coniglio and W. Klein, J. Phys. A 13, 2775 (1980).
[3] C. M. Fortuin and P. W. Kasteleyn, Physica 57, 536
(1972); C. M. Fortuin, Physica 58, 393 (1972); Phys-
ica 59, 545 (1972).
[4] F. Y. Wu, Rev. Mod. Phys. 54, 235 (1982).
[5] C. K. Hu, Phys. Rev. B 29, 5103 (1984).
[6] M. D. De Meo, D. W. Heermann, and K. Binder, J. Stat.
Phys. 60, 585 (1990).
[7] P. C. Hohenberg and B. I. Halperin, Rev. Mod. Phys. 49,
435 (1977).
[8] K. E. Schmidt, Phys. Rev. Lett. 51, 2175 (1983).
[9] J. Goodman and A. D. Sokal, Phys. Rev. Lett. 56, 1015
(1986).
[10] M. Sweeny, Phys. Rev. B 27, 4445 (1983).
[11] Y. Deng, W. Zhang, T. M. Garoni, A. D. Sokal,
and
A. Sportiello, Phys. Rev. E 81, 020102 (2010).
[12] R. H. Swendsen and J. S. Wang, Phys. Rev. Lett. 58, 86
(1987).
[13] R. G. Edwards and A. D. Sokal, Phys. Rev. D 38, 2009
(1988).
[14] U. Wolﬀ, Phys. Rev. Lett. 62, 361 (1989).
[15] T. M. Garoni, G. Ossola, M. Polin, and A. Sokal, J. Stat.
Phys. 144, 459 (2011).
[16] F. Gliozzi, Phys. Rev. E 66, 016115 (2002).
[17] J. S. Wang, O. Kozan, and R. H. Swendsen, Phys. Rev.
E 66, 057101 (2002).
[18] X. Qian, Y. Deng, and H. W. J. Bl¨ote, Phys. Rev. E 71,
016709 (2005).
[19] Y. Deng, T. M. Garoni,
and A. D. Sokal, Phys. Rev.
Lett. 98, 230602 (2007).
[20] G. Grimmett, The random-cluster model
(Springer,
Berlin, 2006).
[21] L. Chayes and J. Machta, Physica A 254, 477 (1998).
[22] X. Qian, Y. Deng, and H. W. J. Bl¨ote, Phys. Rev. E 72,
056132 (2005).
[23] A. Zatelepin and L. Shchur, “Duality of critical in-
terfaces in Potts model:
Numerical check,”
(2010),
arXiv:1008.3573.
[24] F. Gliozzi and M. A. Rajabpour, JSTAT 2010, L05004
(2010).
[25] D. Stauﬀer and A. Aharony, Introduction to Percolation
Theory, 2nd ed. (Taylor & Francis, London, 1994).
[26] M. Weigel, Physics Procedia 3, 1499 (2010).
[27] M. R. Henzinger and V. King, J. ACM 46, 502 (1999).
[28] J. Holm, K. de Lichtenberg,
and M. Thorup, J. ACM
48, 723 (2001).
[29] K. Binder and D. P. Landau, A Guide to Monte Carlo
Simulations in Statistical Physics, 3rd ed. (Cambridge
University Press, Cambridge, 2009).
[30] Note that in the Metropolis dynamics there are no uncon-
ditional rejections, while unconditional acceptance and
rejection occur at equal rates for heat-bath rates.
[31] T. H. Cormen, C. E. Leiserson, R. L. Rivest,
and

<!-- page 14 -->

14
C. Stein, Introduction to Algorithms, 3rd ed. (MIT Press,
Cambridge, MA, 2009).
[32] P. Grassberger, J. Phys. A 25, 5867 (1992).
[33] P. Grassberger, J. Phys. A 32, 6233 (1999), arXiv:cond-
mat/9906309.
[34] P. Grassberger, J. Phys. A 25, 5475 (1992).
[35] R. M. Ziﬀ, J. Phys. A 32, L457 (1999), arXiv:cond-
mat/9907305.
[36] J. Cardy, J. Phys. A 31, L105 (1998).
[37] S. Smirnov and W. Werner, Math. Res. Lett. 8, 729
(2001).
[38] B. Nienhuis, in Phase Transitions and Critical Phenom-
ena, Vol. 11, edited by C. Domb and J. L. Lebowitz (Aca-
demic Press, London, 1987) p. 1.
[39] Z. Zhou, J. Yang, Y. Deng, and R. M. Ziﬀ, Phys. Rev.
E 86, 061101 (2012).
[40] M. F. Gyure and B. F. Edwards, Phys. Rev. Lett. 68,
2692 (1992).
[41] E. M. El¸ci and M. Weigel, in preparation.
[42] M. Weigel, Vertex Models on Random Graphs, Ph.d. the-
sis, University of Leipzig (2002).
[43] M. E. J. Newman and R. M. Ziﬀ, Phys. Rev. E 64, 016706
(2001).
[44] R. E. Tarjan, J. ACM 22, 215 (1975).
[45] Y. Deng, T. M. Garoni, J. Machta, G. Ossola, M. Polin,
and A. D. Sokal, Phys. Rev. Lett. 99, 055701 (2007),
arXiv:0705.2751.
[46] E. M. El¸ci, Sweeny-Algorithmus f¨ur Simulationen des
Potts-Modells, Bachelor’s thesis, Johannes Gutenberg-
Universit¨at Mainz (2011).
[47] M. Weigel, W. Janke, and C. K. Hu, Phys. Rev. E 65,
036109 (2002).
[48] D. D. Sleator and R. E. Tarjan, J. ACM 32, 652 (1985).
[49] A. Gibbons, Algorithmic Graph Theory (Cambridge Uni-
versity Press, Cambridge, 1985).
[50] R. E. Tarjan, Math. Program. 78, 169 (1997).
[51] D. E. Knuth, The Art of Computer Programming, Vol-
ume 3: Sorting and Searching, 2nd ed. (Addison Wesley,
1998).
[52] C. Mart´ınez and S. Roura, J. ACM 45, 288 (1998).
[53] R. Iyer, D. Karger, H. Rahul, and M. Thorup, J. Exp.
Algorithmics 6, 4 (2001).
[54] A. D. Sokal, in Functional Integration: Basics and Ap-
plications, Proceedings of the 1996 NATO Advanced
Study Institute in Carg`ese, edited by C. DeWitt-Morette,
P. Cartier,
and A. Folacci (Plenum Press, New York,
1997) pp. 131–192.
[55] N. Madras and A. D. Sokal, J. Stat. Phys. 50, 109 (1988).
[56] X. J. Li and A. D. Sokal, Phys. Rev. Lett. 63, 827 (1989).
[57] M. B. Priestley, Spectral Analysis and Time Series (Aca-
demic Press, London, 1996).
[58] U. Wolﬀ, Comput. Phys. Commun. 156, 143 (2006),
arXiv:hep-lat/0306017.
[59] M. Weigel and W. Janke, Phys. Rev. E 81, 066701 (2010).
[60] B. Efron and R. J. Tibshirani, An Introduction to the
Bootstrap (Chapman and Hall, Boca Raton, 1994).
[61] Caching eﬀects might lead to unexpected behavior in
some range of system sizes, but asymptotically these ef-
fects should be irrelevant.
[62] GNU
Scientiﬁc
Library,
http://www.gnu.org/software/gsl.
[63] M. Matsumoto and T. Nishimura, ACM Trans. Model.
Comput. Simul. 8, 3 (1998).
[64] J. Jacobsen, J. Salas, and A. Sokal, J. Stat. Phys. 119,
1153 (2005).
[65] H. E. Stanley, J. Phys. A 10, L211 (1977).
[66] “Eﬃcient
implementation
of
Sweeny’s
algorithm,”
https://github.com/ernmeel/sweeny.
