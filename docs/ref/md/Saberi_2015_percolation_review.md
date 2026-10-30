# Recent advances and open challenges in percolation

> Eur. Phys. J. B 89, 253 (2016) · [arXiv](https://arxiv.org/abs/1404.5325)
> Local copy: `docs/ref/Saberi_2015_percolation_review.pdf` (11 pages)

<!-- page 1 -->

arXiv:1404.5325v1  [cond-mat.stat-mech]  21 Apr 2014
Recent advances and open challenges in percolation
N. A. M. Ara´ujo,1, ∗P. Grassberger,2, † B. Kahng,3, ‡ K. J. Schrenk,4, § and R. M. Ziﬀ5, ¶
1Departamento de F´ısica, Faculdade de Ciˆencias, Universidade de Lisboa,
1749-016 Lisboa, Portugal, and Centro de F´ısica Te´orica e Computacional,
Universidade de Lisboa, 1749-016 Lisboa, Portugal
2JSC, FZ J¨ulich, D-52425 J¨ulich, Germany
3Department of Physics and Astronomy, Seoul National University 151-747, Korea
4Computational Physics for Engineering Materials, IfB,
ETH Zurich, Wolfgang-Pauli-Strasse 27, CH-8093 Zurich,
Switzerland and Department of Chemistry, University of Cambridge,
Lensﬁeld Road, Cambridge, CB2 1EW, United Kingdom
5Center for the Study of Complex Systems and Department of Chemical Engineering,
University of Michigan, Ann Arbor, Michigan 48109-2136 USA
Percolation is the paradigm for random connectivity and has been one of the most applied sta-
tistical models. With simple geometrical rules a transition is obtained which is related to magnetic
models. This transition is, in all dimensions, one of the most robust continuous transitions known.
We present a very brief overview of more than 60 years of work in this area and discuss several open
questions for a variety of models, including classical, explosive, invasion, bootstrap, and correlated
percolation.
I.
INTRODUCTION
Percolation is a classic problem in statistical physics
and, like a cat with nine lives, never seems to die. The
study by its proper name began with the work of engineer
Simon Broadbent and mathematician John Hammersley
in the 1950’s [1], inspired by the workings of activated
charcoal gas masks, but it was eﬀectively already consid-
ered by chemist Paul Flory in the early 1940’s in his study
of gelation in polymers [2–4]. The King’s College Lon-
don group of Cyril Domb, which included Michael Fisher,
John Essam, and M. F. Skyes, did a great deal to popu-
larize the percolation problem in the physics community
starting in the 1960’s [5]. In an early paper, Fisher and
Essam showed that the polymerization model of Flory
corresponds to percolation on the Bethe lattice [6, 7].
There followed a great deal of work by physicists, in-
cluding Rushbrooke, Stanley, Coniglio, Halperin, Her-
rmann, Stauﬀer, Aharony, Havlin, Duplantier, etc [8–10].
A major breakthrough was the demonstration in 1969 by
Fortuin and Kasteleyn [11] that percolation is a limiting
case of the general Potts model, which includes the Ising
model, and can be solved exactly. This paved the way to
many exact results in percolation, and allowed powerful
renormalization group ideas to be used [12]. Numerous
variations of percolation, such as invasion percolation,
ﬁrst-passage percolation, directed percolation, and boot-
strap percolation were introduced. Exact “static” expo-
nents in 2d were ﬁrst proposed in the 1980’s by den Nijs
[13–15], Pearson and others, but attempts to ﬁnd exact
∗nmaraujo@fc.ul.pt
† p.grassberger@fz-juelich.de
‡ bkahng@snu.ac.kr
§ kjs73@cam.ac.uk
¶ rziﬀ@umich.edu
exponents for conductivity (the Alexander-Orbach con-
jecture) [16], the backbone [17, 18], and the shortest path
(chemical distance) [19] have not succeeded.
Finding percolation thresholds both exactly and by
simulation has been an enduring subject of research in
this ﬁeld (see [20] and references therein), as well as
the development of algorithms such as those by Hoshen
and Kopelman [21], by Leath [22], and by Newman and
Ziﬀ[23].
Finding rigorous proofs of exact thresholds
and bounds has also been an enduring area of research
for mathematicians (Kesten [24], Wierman [25], Grim-
mett [26], Bollob´as and Riordan [27], etc.).
The derivation of an exact formula for the crossing
probability on a rectangle by Cardy in 1992 [28] (moti-
vated by numerical work of Langlands et al.
[29]) led
to a great deal of work on universal crossing properties
in two-dimensional system, such as those by Pinson [30],
Watts [31], Simmons [32], etc. Around 2000, Schramm
developed the Stochastic Loewner Evolution (SLE) the-
ory, which was soon applied to percolation hulls [33–35].
This caused once again percolation to be an active area of
research in mathematical probability theory, and led to
two Fields’ medals, to Stanislav Smirnov and Wendelin
Werner [36–38]. The results include rigorous derivation
of the static percolation exponents. Another infusion of
interest in percolation came from the surging ﬁeld of net-
work theory, which goes back to the study of random
and complete graphs by Erd˝os and R´enyi (1959), where
the formation of a “giant component” is exactly analo-
gous to percolation [39], but was revitalized by interest
in small-world and scale-free networks. In 2000, Newman
and Moore found the critical point for a random graph
in the limit of large size, in which case the system is ef-
fectively a Bethe lattice, and this result connects to the
early work of Flory, Fisher and Essam, but here with a
general degree distribution [40]. In the ﬁeld of random
networks, the model of explosive percolation was ﬁrst in-

<!-- page 2 -->

2
troduced, by Achlioptas, D’Souza and Spencer [41], and
this has been another problem that has led to a wave of
new interest in percolation.
In this article, we present some perspective on open
problems and challenges in percolation. The ﬁeld is so
vast that we can only touch on a subset of them, more
aligned with our own interests.
II.
CLASSICAL UNCORRELATED
PERCOLATION
The most basic question in percolation is the value
of the threshold for a given system.
Exact thresholds
in two dimensions for the square, triangular, honeycomb
and related lattices (kagom´e and (3-122) site thresholds)
were found by Sykes and Essam in 1963, using the star-
triangle transformation [42]. In 1984 Wierman [25] gen-
eralized this transformation to ﬁnd the threshold of the
bow-tie lattice, and in 2006 Ziﬀand Scullard [43] showed
that thresholds can be found for any lattice that can be
represented as a self-dual 3-hypergraph (that is, decom-
posed into triangles that form a self-dual arrangement).
Just recently, Grimmett and Manolescu [44] showed that
thresholds can be found for any lattice that can be rep-
resented geometrically as an isoradial graph, yielding a
broad new class of exact thresholds and providing a proof
[45] of Wu’s 1979 conjecture [46] for the threshold of the
checkerboard lattice. So the widely held belief that there
are just a handful of lattices where the threshold is known
exactly is no longer true. However, many systems of long
interest (such as site percolation on the square and hon-
eycomb lattices, and bond percolation on the kagom´e
lattice) do not ﬁt into the self-dual hypergraph or iso-
radial forms and cannot presently be solved.
Is there
another approach that can give exact thresholds to these
and other lattices? Or can a proof be made that thresh-
olds for some lattices can never be found in a closed form?
New methods to obtain very high precision estimates of
the threshold with unprecedented precision (over 12 dig-
its for some models) have been developed by Jacobsen
and Scullard [47, 48], but none of these seem to suggest
any obvious closed-form values. For higher dimensions,
exact thresholds would be very beneﬁcial but seems un-
likely. In four dimensions, a claim for an exact threshold
of plaquettes has been made [49].
Likewise, there are no exactly known exponents for
more dynamical properties such as the shortest-path ex-
ponent dmin, dynamic exponent z, conductivity exponent
t, spectral dimension, etc., even in two dimensions (only
the related exponent g is known [50, 51]). For directed
percolation, which can be interpreted as a dynamical
form of percolation also related to the contact process, no
exact results are known for any dimension. After years of
trying, it seems unlikely that these other exponents will
ever be found exactly, but perhaps a new breakthrough
will be found.
Many questions relating to scaling functions and am-
plitude ratios were looked at extensively in the 1980’s,
but several questions remain unanswered, and with new
computers and algorithms these problems are worth
revisiting [52].
Amplitude ratios are useful for con-
ﬁrming universality of diﬀerent percolation models [53,
54].
Near the percolation threshold the size distribu-
tion of percolation clusters satisﬁes the scaling form
ns(p) ∼As−τf(B(p−pc)sσ) where A and B are system-
dependent scale factors, while τ and σ are universal ex-
ponents and f(z) is a universal function for a given di-
mensionality. The universality of f(z) implies the uni-
versality of amplitude ratios (such as the ratio of the
mean cluster size an ǫ above and below pc), which were
studied extensively in the past [55]. With a precise de-
termination of the scaling function, many of these uni-
versal ratios should be able to be determined, but as far
as we know the scaling function has never been deter-
mined precisely and used for this. For the correlation-
length ξ, the amplitude ratio has been predicted to be
exactly 2 [56], but because of questions of how the corre-
lation length should be measured, it seems that this value
has never been veriﬁed numerically.
In 2001, Seaton
[57] predicted an exact amplitude ratio for percolation,
R+
ξ = [α(2−α)(1−α)A]1/dξ+
0 = [40/27
√
3]
1
2 = 0.9248 . . .
where α = −2/3, A+ is the coeﬃcient of the singular part
|p −pc|2−α of the free energy (number of clusters) and
ξ+
0 is the coeﬃcient to the divergence of the correlation
length ξ ∼ξ0(pc −p)ν for p < pc. While this ratio agreed
with previous results by Delﬁno and Cardy [58], questions
about the correct deﬁnition of the correlation function to
use remains.
A way to put the size distribution into a universal form
is by considering the enclosing area or volume; in that
case the area of the n-th ranked cluster satisﬁes Zipf’s
law An ∼C/n, with C known exactly in the case of
two-dimensional enclosed area [59, 60]. How this quan-
tity behaves away from pc and how that relates to reg-
ular scaling functions has not been detailed. This same
approach should apply to any critical system, including
directed percolation.
An early exact result in percolation is the density n of
clusters of wet sites in bond percolation at the threshold,
found by Temperley and Lieb [61] for the square lattice
and evaluated simply as (3
√
3 −5)/2 [62], and is also
known for bond percolation on the triangular lattice [63].
Surprisingly, while many things have been proven exactly
for site percolation on the triangular lattice, this quantity
is known only numerically there: n = 0.0176255(5) [62].
Can this quantity also be found exactly?
III.
CLASSICAL PERCOLATION WITH
LONG-RANGE CORRELATED DISORDER
Percolation theory and related models have been ap-
plied to study transport and geometrical properties of
disordered systems [64, 65]. Frequently, such disorder is
characterized by a power-law long-ranged spatial correla-

<!-- page 3 -->

3
tion. This fact has motivated studies of percolation mod-
els where the sites of the lattice are not occupied indepen-
dently, but instead with long-range spatial correlation,
in a process named correlated percolation [64–78]. The
qualitative picture that emerged is that, in the presence
of long-range correlations, percolation clusters become
more compact and their transport properties change ac-
cordingly. These ﬁndings have also been conﬁrmed by
experimental studies of the transport properties of clus-
ters in correlated invasion percolation [79, 80].
The typical strategy to investigate correlated percola-
tion is to work with a landscape of random heights h,
where h(x) is the height of the landscape at the lattice
site at position x [64, 66–68, 81, 82]. Recently, ranked
surfaces have been introduced, providing the adequate
framework to tackle this problem [78, 83]. The ranked
surface of a discrete landscape is constructed as follows:
One ﬁrst ranks all sites in the landscapes according to
their height, from the smallest to the largest value. Then,
a ranked surface is constructed where each site has a
number corresponding to its position in the rank. The
following percolation model can then be deﬁned: Ini-
tially, all sites of the ranked surface are unoccupied. The
sites are occupied one-by-one, following the ranking. At
each step, the fraction of occupied sites p increases by the
inverse of the total number of sites in the surface. By this
procedure, a conﬁguration of occupied sites is obtained,
the properties of which depend on the landscape. For
example, if the heights are distributed uniformly at ran-
dom, classical percolation with fraction of occupied sites
p is obtained [23, 84, 85]. To generate ranked surfaces
with power-law correlated disorder the Fourier ﬁltering
method is used [67, 86–93].
The focus is usually on the dependence of percolation
properties on how the correlations decay with the spa-
tial distance, typically characterized by the Hurst ex-
ponent H.
Schmittbuhl et al.
have shown that for
H > 0, there is no percolation transition [68].
In-
stead, compact clusters grow and merge in such a way
that the size of the largest cluster grows linearly with
the fraction of occupied sites. By contrast, for H < 0
the percolation transition is critical and critical expo-
nents can be found.
According to the extended Har-
ris criterion, the critical exponents should change with
H [66, 68, 77, 78, 94, 95]. It was found that the frac-
tal dimension of the largest cluster, its perimeter, as well
as the dimension of its shortest path, backbone, and red
sites depend on H [67, 68, 78, 81, 96]. Strong dependence
on H is also found for the electrical conductivity expo-
nent of the largest cluster and the growth of bridge sites
in the correlated percolation model. Schrenk et al. have
proposed a functional dependence of several exponents on
the Hurst exponent H [78]. While these relations were
found as the simplest rational expressions giving the con-
sidered critical exponents as function of H, they do not
have theoretical support. It is still an open question if
these expressions can be found by scaling arguments or
ﬁeld theory. Concerning the perimeter fractal dimension
of the largest cluster, it was found numerically that the
duality relation of Duplantier [97] is fulﬁlled in the full
range of Hurst exponents [78]. A more theoretical argu-
ment supporting this ﬁnding is also still missing.
IV.
VARIANTS OF PERCOLATION
“Classical” (ordinary) percolation deals [9, 10] with the
situation where small bridges are either established or
removed randomly, irreversibly, and independently, until
large scale connectivity is either established or lost. Here,
a “bridge” might be a link or a node in a network, with
regular lattices being very regular class of special net-
works. Apart from this no other processes of relevance go
on, which means that each bridge is dealt with only once,
whence the randomness or the phenomenon can be either
considered as frozen (a view mostly taken in condensed
matter) or as characterizing a stochastic contagion pro-
cess, in which case it is describing the spreading of some
epidemic [98], computer malware, or public opinion. In
the latter context, it is also called the SIR (“susceptible-
infected-removed”) epidemic model [54, 99, 100], to be
contrasted to the SIS (“susceptible-infected-susceptible”)
model where bridges are removed after a ﬁnite life time
(of infectivity), and which corresponds to directed per-
colation in space-time (the directness being that of the
time arrow). Both SIR and SIS (i.e., undirected and di-
rected percolation) show transitions between phases with
and without long range connectivity, which are continu-
ous (“second order”) in the sense that a suitable order
parameter (which usually is the density of the largest
connected cluster) is zero at the transition point and in-
creases continuously (but of course with singular deriva-
tives) as one enters the supercritical (connected) phase.
Deviations from this classical model can involve any of
the basic ingredients:
1. “Bridges” (links, infections) themselves can already
be long ranged. Models with power behaved long
range interactions have been studied in statistical
physics since the late 1960’s , when it was real-
ized that the 1-dimensional Ising model can have a
non-trivial phase transition if interactions are suf-
ﬁciently long ranged, with a decay ∼1/r2 being
the turning point [101–103]. Very early, Thouless
[103] argued that the transition in this model is
of Kosterlitz-Thouless [104] type, i.e.
it is of in-
ﬁnite order, the correlation length does not scale
with a power of the distance from the critical point,
and there are generic power laws with continuously
varying exponents in the subcritical phase.
The
same should be true for percolation. But his ar-
guments were soon shown to be incomplete [102].
Indeed, when it was proven that in this case the
transition is ﬁrst order (discontinuous) [105], the
predictions of Kosterlitz seemed to be obsolete. But
recent simulations [106] show that there are indeed

<!-- page 4 -->

4
continuously varying critical exponents (in the su-
percritical phase), and that the correlation length
diverges exponentially. Related to this is a model
by Boettcher et al. [107], where long range links
are not established in a completely random fash-
ion. Rather, in the Boettcher et al. model all long
range bonds must have length 2k, k = 1, 2, 3, . . .,
bonds with length 2k connect only to sites with
coordinates i ≡0
mod
2k, and they are all es-
tablished with the same probability p. This leads
to the same average density of long range links as
in the random model with P(r) ∼1/r2, and both
models show very similar behavior.
How similar
they are in details remains to be seen.
2. In two dimensions, there is no similarity with the
Kosterlitz-Thouless transition, but there is a long-
standing controversy. When increasing the power
in the contact probability ∼1/rσ, there comes a
point where the critical behavior changes from a
regime with σ-dependent critical exponents to or-
dinary (short-range) percolation. The controversy
deals with the question when this happens. Diﬀer-
ent ﬁeld theoretical arguments suggest that it either
happens when also the scaling of the supercritical
process changes, or already for smaller σ. Exten-
sive simulations (both for percolation and for the
Ising model, see [108, 109]) seemed to have settled
this in favor of the latter option, but very recent
simulations suggest that neither is right and things
might be more complicated [110–112].
3. Bridges can be established in a not entirely ran-
dom way. In particular, when establishing a new
bridge, one may take two random (and not yet
established) “candidates,” and establish that can-
didate that leads to less long range connectivity
(i.e., to a smaller cluster). This was proposed by
Achlioptas et al.
[41], who termed it “explosive
percolation” and claimed it to have a discontinuous
transition, where the transition is delayed so long
that it ﬁnally takes place in a sudden, “explosive”
way. This paper has led to several similar models,
some of which actually do have ﬁrst-order transi-
tions (see e.g. [113, 114]) — although the original
Achlioptas model does not [115–117]. Notice that
the choice between the two candidates introduces a
severe non-locality which goes far beyond anything
that can be described by a local ﬁeld theory.
It
is presumably for this reason that the Achlioptas
model seems to show (there is not yet any analytic
proof) ﬁnite size scaling that is qualitatively diﬀer-
ent from that expected from the renormalization
group [117]. This family of models is discussed in
detail in Sec. V.
4. Bridges can be established in a partly synchro-
nized way. In ordinary bond percolation, links in
a graph are randomly selected one after the other.
In “agglomerative percolation” [118–120] this is re-
placed by the rule that a cluster is selected ran-
domly, and then all links connected to it are es-
tablished immediately.
On non-bipartite lattices
in two and three dimensions this is in the uni-
versality class of ordinary percolation, but not in
one-dimensional lattices [119], on a class of trees
[118], and on Erd˝os-R´enyi graphs [121]. Notice that
this model also introduces some non-locality (al-
beit much weaker than in the Achlioptas model),
if very large clusters are chosen during late stages
of the process. A special feature happens on bipar-
tite lattices like the square and simple cubic lattice
[120, 122]. There, at the critical point also a Z2
symmetry is spontaneously broken, which leads at
least in two dimensions to a new universality class.
5. In the standard scenario, bridges are newly estab-
lished in a pre-existing “virtual” network e.g., in
bond percolation all nodes and all possible links
are already deﬁned, while in site percolation the
new sites are chosen from a pre-existing set and all
links are already established. In contrast to this,
Callaway et al.
[123] considered the case where
links are established in a growing network. More
precisely, each time a node is added to the exist-
ing network, also a link is added with probability
δ < 1. In contrast to the Barab´asi-Albert model
[124], where the new links are added so that nodes
with higher degree get attached preferentially, new
links are added completely at random. Again this
leads to a Kosterlitz-Thouless type transition when
δ passes through a critical value, but now the tran-
sition is (as the original KT transition) of inﬁnite
order instead of being discontinuous.
6. The establishment of bridges can be a cooperative
eﬀect. This seems to be the most interesting case.
Models of this type have been studied as bootstrap
percolation (or as the closely related k-core percola-
tion) since the 1970’s [125, 126], while co-epidemics
(see below) are a major health hazard.
• In bootstrap percolation with index m one
considers a connected cluster on some graph,
and removes recursively all nodes which have
less than m neighbors in the cluster.
For
m = 1 nothing is done, as a connected graph
has always at least one neighbor. For m = 2
all leaves are removed in a ﬁrst step.
This
might turn some nodes, which originally had
≥2 neighbors, into leaves. They are removed
in a second step. If this again turns nodes into
leaves, they are removed again, etc. In this
way one ends up either with an empty clus-
ter, or with a cluster where all nodes have ≥2
neighbors. In k−core percolation one starts
with the cluster which contains all nodes and
keeps at the end only the largest connected

<!-- page 5 -->

5
component, i.e. one ends up ﬁnally with the
largest connected subgraph in which all nodes
have degree at least k. Bootstrap percolation
can show either continuous or discontinuous
transitions, depending on the type of networks
(lattice or random), and on the dimensional-
ity when the network is a regular lattice. For
details see [126].
• Instead of deleting “poor” nodes, let us now
consider the opposite process of complex
contagion [127, 128].
Here nodes are at-
tached to a growing cluster and nodes with
more neighbors in the cluster are more likely
to be attached. A possible application is po-
litical opinion spreading. We assume that no
person is so eloquent that (s)he alone is able to
convince his or her neighbor. But if there is al-
ready a small group of early adapters, then the
combined argumentation of several of them
is able to convince others with probability q.
When people get convinced if and only if they
have ≥m convinced neighbors, this leads ex-
actly to bootstrap percolation.
More inter-
esting, however, is the case where the chance
qm to be ﬁnally convinced (or infected, or at-
tached to the cluster, depending on the appli-
cation) depends non-trivially on m, the num-
ber of neighbors in the attacking cluster. Take
e.g. site percolation on a square lattice, where
any site is incorporated with probability p into
a growing cluster, if it has at least one infected
neighbor. In this case qm = p for all m ≥1. In
bond percolation, on the other hand, each site
has the same chance p to convince (or infect)
any not yet convinced (or healthy) neighbor.
A site with m infected neighbors succumbs
then with probability qm = p + (1 −p)p +
. . . (1 −p)m−1p = 1 −(1 −p)m.
In general
we expect qm (the probability that m neigh-
bors together will lead to infection) to grow
with m. If this growth is moderate (as e.g. in
bond percolation) we expect the transition to
be continuous and in the universality class of
ordinary percolation. This is no longer true if
qm grows too fast. In that case one ﬁnds ﬁrst-
order transitions, where the epidemic (politi-
cal opinion, fashion, computer virus, ...) ei-
ther does not take oﬀat all or takes oﬀex-
plosively. In between the ﬁrst and second or-
der regimes there is a tricritical point whose
properties on ﬁnite-dimensional lattices where
studied by Jansen et al.
by ﬁeld theoretic
methods [129].
• For random sparse (e.g.
Erd˝os-R´enyi) net-
works the tricritical point is particularly easy
to ﬁnd [127, 130]: If q1 > 0, the transition
from second to ﬁrst order happens exactly
when q2 = 2q1, independent of its degree dis-
tribution.
• One important physical application of this
model is to critically pinned driven interfaces
in isotropic random media [131].
It is well
known that such interfaces are often self-aﬃne
(i.e., show nontrivial scaling), while the bulk
behind them is not fractal (as it would be
for critical percolation). Thus in this problem
the phase transition is actually hybrid, i.e. it
combines a jump in the order parameter (the
density of the infected phase) with non-trivial
scaling related to the interface. The fact that
any ﬁrst order transition in k−core percola-
tion is actually hybrid was stressed in partic-
ular in [132].
• Notice that interfaces in the present model
have in general overhangs (as in real isotropic
media, but in contrast to the standard model
[131, 133] used to describe them). This has
important consequences. First of all, such in-
terfaces in 1+1 dimension are always in the
ordinary percolation universality class [130],
as conjectured already in [134] (see, however,
[135, 136]).
This is related to the proof of
[137] that there cannot be ﬁrst-order transi-
tions in 2-dimensional random media.
Sec-
ondly, it means that there will always be some
holes left in a growing cluster, and these holes
in general will show non-trivial scaling at the
pinning point (as required by its hybrid na-
ture). Finally, it might imply that the point
where the zero temperature random ﬁeld Ising
model changes from being hysteresis-free to
having hysteresis is not, as claimed in [138],
a critical point with upper critical dimension
dc = 6, but a tricritical point with dc = 5.
The ﬁeld theoretical predictions of [129] for
this tricritical point were not well conﬁrmed
by the simulations of [130] for d = 3, but they
were very well conﬁrmed for d = 4, where they
are supposed to be much more precise.
• Instead of cooperativity between attacking
neighbors, we can also have cooperativity be-
tween diﬀerent modes of attack. This refers
in particular to cooperative co-epidemics
or syndemics [139]. Take e.g. the examples
of the 1918 Spanish ﬂu and TB [140, 141], or of
HIV and a host of other diseases like diabetes,
hepatitis, TB and others [142, 143]. In these
cases one disease weakens the victim so much
that the victim falls prey to the other one.
A simple SIR type model for this with per-
manently increased susceptibility for the other
disease was developed in [144]. Simulations of
this model on random graphs and on lattices
[145] show a host of diﬀerent behaviors, de-
pending on the topology of the graph but also
on the detailed of the implementation (and, in

<!-- page 6 -->

6
one case, even on the initial conditions). But
whenever a ﬁrst-order transition is found, it is
in general hybrid.
7. The last model we shall discuss here are interde-
pendent networks. Consider e.g. a country like
Italy with towns connected by a network of elec-
tric power lines and another network of computer
connections. The latter are needed for controlling
the power stations, while the former are needed
to provide power to the computers. If some node
fails, then this can lead to a huge cascade of fail-
ures as happened indeed in Italy on September 28,
2003. As shown by Buldyrev et al. in a seminal pa-
per [146], this failure would have been a ﬁrst-order
transition, if it had happened on a locally loopless
network like e.g.
a sparse Erd˝os-R´enyi network.
But, as shown in [119], it still would be a contin-
uous transition on a 2-dimensional regular lattice.
Since Italy is 2-dimensional but not quite regular,
it is unclear what should apply in this case. But in
any case, in a series of papers (see e.g. [147–149])
several other types of interdependencies were stud-
ied. In all cases with ﬁrst-order transitions, these
are indeed again hybrid.
Finally, let us point out a possible connection between
interdependent networks and cooperative co-epidemics.
People have livers and lungs, and both can become sick.
Someone with a liver disease will be more likely to get
also a lung infection and vice versa. Let us now assume
that liver and lung infections are both lethal.
In this
extreme case we are dealing exactly with two interde-
pendent networks as described in [146]. Details of this
correspondence have not yet been worked out.
V.
EXPLOSIVE PERCOLATION
Conventionally percolation transitions are known to
be continuous. Recently, however, interest in discontin-
uous percolation transitions (DPT) has been triggered
and boosted by i) the explosive percolation model [41]
and ii) the cascading failure model in interdependent net-
works [146]. Actually the subject of DPT was initiated by
the k-core percolation model [125, 126, 150, 151] a long
time ago. However, the mechanism of the DPT arising in
the k-core model was unusual, so that further studies of
the DPT had not proceeded very much. The two models
i) and ii) were designed for DPTs; however, the governing
processes to the PT and their mechanisms are diﬀerent:
cluster coagulation and fragmentation processes for the
model i) and ii), respectively. Moreover, the order pa-
rameters are also diﬀerent. Thus, we need to discuss the
two problems separately, and ﬁnd common features in a
uniﬁed framework.
We ﬁrst discuss the DPT driven by cluster merging
processes. The explosive percolation (EP) model was in-
troduced motivated by a mathematical invention, the so-
called Achlioptas process: at each step, we take the one
(among multiple options) which is the optimal one for
the formation of a given target pattern. To generate a
DPT in random graph, at each step, two pairs of nodes
that are not connected yet are selected at random, and
one pair of nodes among them is taken and connected by
a new link, which is the one that produces a smaller size
of connected cluster than the size of the cluster produced
by the other option. Then, the growth of large clusters is
suppressed and clusters with medium sizes are abundant.
As the number of added links is increased, such medium-
sized clusters merge and create a macroscopic-sized clus-
ter, which emerges drastically at a delayed percolation
threshold.
Following this pioneering paper [41], many related
models have been introduced [113, 152–158]. One of no-
ticeable issues addressed in following works was that the
EP transition in random graph is not indeed discontin-
uous but continuous in the thermodynamic limit.
Ac-
tually, due to extremely slow convergence to asymptotic
behavior as the system size is increased, whether the EP
is indeed continuous or discontinuous was hard to be de-
termined numerically and became a controversial issue.
This claim was addressed ﬁrstly based on numerical re-
sults for a speciﬁc model [115]; however, since the claim
was not ﬁrmly supported by analytic solution, the claim
had not been ﬁrmly accepted ﬁrst time.
Later, based
on the numerically observed fact that the average size
of medium-sized cluster is not suﬃciently large [154], a
mathematical proof [116] was presented, which is the fol-
lowing: the number of clusters participated to generate a
macroscopic-sized cluster is not sub-extensive to system
size in order of magnitude, and thus, it cannot bring out
a DPT, but a continuous PT.
The EP transition in Euclidean space has been also ex-
tensively studied [152, 159]. However, whether the perco-
lation transition is discontinuous or continuous had not
been clariﬁed, either. In spite of such extensive studies,
understanding of the EP transition in Euclidean space
and on random graph in a uniﬁed scheme had been ab-
sent. For the continuous PT case, the Erd˝os and R´enyi
model on random graph was regarded as the mean-ﬁeld
model of the ordinary percolation model.
To resolve
this goal for the DPT case, a model called the spanning
cluster-avoiding (SCA) model was introduced [114]. In
this model, the target pattern in the Achlioptas process
is taken as a spanning cluster, which is actually stan-
dard in Euclidean space. This model used previous re-
sults, the scaling behavior of the number of bridge bonds
above a percolation threshold.
Here the bridge bonds
are those that would form a spanning cluster if occu-
pied [83]. Using the scaling formula for the bridge bonds,
the percolation thresholds of the SCA model were calcu-
lated analytically for various numbers of potential links
in the Achlioptas process. This analytic result leads to
the following conclusion: the EP transition can be either
continuous or discontinuous, depending on the number of
multiple options, if the spatial dimension is less than the

<!-- page 7 -->

7
upper critical dimension, and is always continuous other-
wise. Subsequently, it was concluded that the transition
of the ordinary EP model is continuous as a mean-ﬁeld
solution of the SCA model.
Further intriguing and fundamental questions remain
open. For example, the DPT occurring in the SCA model
is rather trivial in the sense that the order parameter in-
creases drastically all the way to unity at a trivial per-
colation threshold, which is equal to the system size in
the thermodynamic limit. This behavior is similar to the
one of the DPT governed by the Smoluchowski coagula-
tion equation with the reaction kernel kikj ∼(ij)ω with
0 ≤ω < 0.5 [160].
Therefore, it would be intriguing
to study a non-trivial DPT model in which the order pa-
rameter is increased all the way to a ﬁnite value less than
unity at a ﬁnite threshold, and gradually beyond that
point [107, 161–167]. Investigation of the origin of such
a non-trivial DPT would be also interesting. Moreover,
understanding of other discontinuous phase transitions in
non-equilibrium systems such as synchronization transi-
tions [168] and jamming transitions [169] in the perspec-
tive of the DPT under cluster merging processes would
be interesting.
Next we discuss the DPT occurring in inter-dependent
networks. In this DPT, the order parameter is taken as
the fraction of nodes in the mutually connected largest
component, which is not the same as the standard ones.
The model associated with this DPT was introduced to
study cascading-failure dynamics in the inter-dependent
systems. Even though the original model was understood
analytically, a simpler analytic solution was derived from
the viewpoint of epidemic spreading [149]. In this per-
spective, the DPT can be understood as the emergence
of a mutually connected giant cluster in cluster merging
processes. Thus, it would be interesting to compare the
mechanism of such a DPT with that of the EP model.
This investigation would provide a clue to enhance the
robustness of inter-dependent networks.
Percolation transitions have served as a platform for
understanding phase transitions in non-equilibrium sys-
tems. Likewise, the theoretical framework developed to
understand the origin of the DPTs occurring in the mod-
els i) and ii) is anticipated to serve as a basis for further
research on drastic phase transitions in non-equilibrium
complex systems.
VI.
RETENTION CAPACITY AND
WATERSHEDS OF LANDSCAPES
Percolation has also been applied to study the prop-
erties of real and artiﬁcial landscapes. A landscape is
typically represented as a digital elevation map (DEM),
which consists of a two-dimensional array of regular cells
(sites) to which average heights can be associated. By
mapping the DEM to a ranked surface, it is possible
to identify the sequence of ﬂooded sites when water is
dripped over the landscape, ﬁlling it from the valleys to
the mountains, and letting the water ﬂow out through the
open boundaries [83]. In this spirit, one can ask what the
maximum volume of water is that can be retained by the
landscape. This problem can indeed be mapped to stan-
dard percolation [170, 171]. For example, consider the
simplest case of a two-level random landscape of L × L
sites, with a fraction p sites of unity level and 1 −p of
level zero. The water retention of such a landscape is
R(L)
2
= L2(p −P∞), where P∞(p) is the fraction of sites
belonging to the percolation cluster [170]. Since this clus-
ter touches the borders, the water that falls on its sites
ﬂows out of the landscape. In the general case of a land-
scape with equal number of n levels, it was argued [170]
that its retention R(L)
n
can be expressed as,
R(L)
n
=
n−1
X
i=1
R(L)
2
(i/n) .
(1)
This expression was also shown to hold for correlated
landscapes [172].
Some exact expressions have been
found for the retention capacity of ﬁnite (and very small)
lattice sizes [170, 172] but, a general analytic expression
is still missing.
In hydrology it is also important to identify the water-
shed lines dividing the landscape into diﬀerent drainage
basins. When the landscape is ﬂooded from the valleys
lakes are formed. As the level of water rises, those lakes
start to merge and form even larger lakes. If one merges
the lakes under the constraint that no lake percolates,
i.e., no lake connects two predeﬁned opposite boundaries
of the landscape, one ends up with only two lakes sep-
arated by the main watershed line [173].
For random
uncorrelated landscapes this line is a fractal of dimen-
sion 1.2168±0.0005 [174] and its statistical properties are
consistent with SLEκ, with κ = 1.734 ± 0.005 [175]. For
correlated landscapes, the fractal dimension of the water-
shed line decreases with the Hurst exponent [174, 176].
The determination of the watershed line can be mapped
to a percolation problem where sites are sequentially oc-
cupied (according to their rank in the ranked surface) un-
der the constraint that percolation is suppressed. Such
description reveals a tricritical point at a critical fraction
of occupied (ﬂooded) sites and unveils how several seem-
ingly unrelated physical models (e.g., optimal path, op-
timal path cracks, and polymers in strongly disordered
media) fall into the same universality class [83].
The
generalization of ranked percolation in three dimensions
provides the framework to determine the eﬀective shares
when diﬀerent companies or nations extract either oil,
gas, or water, from the same porous formation [177].
VII.
DIRECTED PERCOLATION
Let us add ﬁnally two comments on models related to
directed percolation.
One deals with the (in-)famous “pair contact process
with diﬀusion” (PCPD) in one dimension of space. For

<!-- page 8 -->

8
a long time this was believed to have continuously vary-
ing critical exponents, and yet two recent simulations,
Ref. [178] and Ref. [179], still come to opposite conclu-
sions. While the former claim that it is in the directed
percolation universality class, the second rules this out.
The
second
deals
with
the
“parity
conserving
branching-annihilating random walk” (pcBARW) [180],
which mainly diﬀers from directed percolation and the
contact process by preserving the “parity” P(t)
=
(−1)N(t), where N(t) is the number of particles at time
t. When starting with two initial particles (i. e. in the
“even” sector), simulations (the most recent and careful
ones being in [181]) suggest that the exponent η deﬁned
by N(t) ∼tη is exactly equal to zero. This cries out for
a theoretical explanation, but so far none seems in sight.
VIII.
FINAL REMARKS
Percolation is a vast subject with more than 80 thou-
sand papers published over the last 60 years according to
the ISI Web of Knowledge, and about one paper posted
per day on the arXiv related to percolation. Thus, in-
stead of a comprehensive review (which would be an epic
journey and impossible in the space available), we de-
cided to rather give a ﬂavor of this fascinating and active
ﬁeld and oﬀer a glimpse at the extensive list of still open
questions, following our own interests. We foresee many
more years of interesting ﬁndings, constantly raising new
questions.
ACKNOWLEDGMENTS
We acknowledge ﬁnancial support from: (NA) the Por-
tuguese Foundation for Science and Technology (FCT),
Contracts No.
IF/00255/2013 and No.
EXCL/FIS-
NAN/0083/2012; (BK) the Korean NRF, Grant No.
2010-0015066; and (KJS) the Swiss National Science
Foundation, Grant No.
P2EZP2-152188.
PG thanks
Golnoosh Bizhani for interesting discussion and collab-
oration on cooperative percolation.
[1] S.
R.
Broadbent
and
J.
M.
Hammersley,
Proc. Cambridge Philos. Soc. 53, 629 (1957).
[2] P. J. Flory, J. Am. Chem. Soc. 63, 3083 (1941).
[3] P. J. Flory, J. Am. Chem. Soc. 63, 3091 (1941).
[4] P. J. Flory, J. Am. Chem. Soc. 63, 3096 (1941).
[5] J. W. Essam, Rep. Prog. Phys. 43, 833 (1980).
[6] M. E. Fisher and J. W. Essam, J. Math. Phys. 2, 609
(1961).
[7] M. E. Fisher, J. Math. Phys. 2, 620 (1961).
[8] D. Stauﬀer, Phys. Rep. 54, 1 (1979).
[9] D. Stauﬀer and A. Aharony, Introduction to Percolation
Theory, 2nd ed. (Taylor and Francis, London, 1994).
[10] M. Sahimi, Applications of Percolation Theory (Taylor
and Francis, London, 1994).
[11] P. W. Kasteleyn and C. M. Fortuin, J. Phys. Soc. Japan
(Suppl.) 26, 11 (1969).
[12] J. Cardy, Scaling and Renormalization in Statisti-
cal Physics (Cambridge University Press, Cambridge,
1996).
[13] M. P. M. den Nijs, J. Phys. A 12, 1857 (1979).
[14] B.
Nienhuis,
E.
K.
Riedel,
and
M.
Schick,
J. Phys. A 13, L189 (1980).
[15] R. B. Pearson, Phys. Rev. B 22, 2579 (1980).
[16] S. Alexander and R. Orbach, J. Physique (Paris) Lett.
43, 625 (1982).
[17] P. Grassberger, J. Phys. A 25, 5867 (1992).
[18] P. Grassberger, Physica A 262, 251 (1999).
[19] H. J. Herrmann, D. C. Hong,
and H. E. Stanley,
J. Phys. A 17, L261 (1984).
[20] K. J. Schrenk, N. A. M. Ara´ujo, and H. J. Herrmann,
Phys. Rev. E 87, 032123 (2013).
[21] J.
Hoshen
and
R.
Kopelman,
Phys. Rev. B 14, 3438 (1976).
[22] P. L. Leath, Phys. Rev. B 14, 5046 (1976).
[23] M.
E.
J.
Newman
and
R.
M.
Ziﬀ,
Phys. Rev. Lett. 85, 4104 (2000).
[24] H. Kesten, Commun. Math. Phys. 74, 41 (1980).
[25] J. C. Wierman, J. Phys. A 17, 1525 (1984).
[26] G. R. Grimmett, Percolation, Grundlehren der mathe-
matischen Wissenschaften, Vol. 321 (Springer, 1999).
[27] B. Bollob´as and O. Riordan, in An Irregular Mind,
Bolyai Society Mathematical Studies, Vol. 21 (Springer,
Berlin, 2010).
[28] J. L. Cardy, J. Phys. A 25, L201 (1992).
[29] R. P. Langlands, C. Pichet, Ph. Pouliot, and Y. Saint-
Aubin, J. Stat. Phys. 67, 553 (1992).
[30] H. T. Pinson, J. Stat. Phys. 75, 1167 (1994).
[31] G. M. T. Watts, J. Phys. A 29, L363 (1996).
[32] J. J. H. Simmons, P. Kleban, and R. M. Ziﬀ, J. Phys.
A 40, F771 (2007).
[33] O. Schramm, Israel J. Math. 118, 221 (2000).
[34] O. Schramm, Electron. Commun. Prob. 6, 115 (2001).
[35] O.
Schramm
and
S.
Sheﬃeld,
Acta Math. 202, 21 (2009).
[36] S. Smirnov, C. R. Acad. Sci. Paris I 333, 239 (2001).
[37] S. Smirnov and W. Werner, Math. Res. Lett. 8, 729
(2001).
[38] S. Smirnov, in Proceedings of the International Congress
of Mathematicians, Madrid, Spain, 2006, edited by
M. Sanz-Sol´e, J. Soria, J. L. Varona,
and J. Verdera
(European Mathematical Society, Z¨urich, 2006) p. 1421.
[39] P. Erd˝os and A. R´enyi, Publ. Math. Inst. Hung. Acad.
Sci. 5, 17 (1960).
[40] C. Moore and M. E. J. Newman, Phys. Rev. E 62, 7059
(2000).
[41] D. Achlioptas,
R. M. D’Souza,
and J. Spencer,
Science 323, 1453 (2009).
[42] M.
F.
Sykes
and
J.
W.
Essam,
Phys. Rev. Lett. 10, 3 (1963).
[43] R.
M.
Ziﬀ
and
C.
R.
Scullard,
J. Phys. A 39, 15083 (2006).

<!-- page 9 -->

9
[44] G.
R.
Grimmett
and
I.
Manolescu,
Prob. Thoery Relat. Fields (2013), 0.1007/s00440-013-0507-y,
publ. online.
[45] R. M. Ziﬀ, C. R. Scullard, J. C. Wierman, and M. R. A.
Sedlock, J. Phys. A 45, 494005 (2012).
[46] F. Y. Wu, J. Phys. C 12, L645 (1979).
[47] J.
L.
Jacobsen
and
C.
R.
Scullard,
J. Phys. A 46, 075001 (2013).
[48] J. L. Jacobsen, J. Phys. A 47, 135001 (2014).
[49] S. Miyazima, Bulletin, Institute of Science and Technol-
ogy Research 17, 119 (2005).
[50] P. Grassberger, J. Phys. A 32, 6233 (1999).
[51] R. M. Ziﬀ, J. Phys. A 32, L457 (1999).
[52] R. M. Ziﬀ, arXiv:1103.3243.
[53] A.
Aharony
and
D.
Stauﬀer,
J. Phys. A 30, L301 (1997).
[54] D.
R.
de
Souza,
T.
Tom´e,
and
R.
M.
Ziﬀ,
J. Stat. Mech. , P03006 (2011).
[55] V. Privman, P. C. Hohenberg,
and A. Aharony,
in Phase transitions and critical phenomena, Vol. 14,
edited by C. Domb and J. L. Lebowitz (Academic, NY,
1991).
[56] J. T. Chayes, L. Chayes, G. R. Grimmett, H. Kesten,
and R. H. Schonmann, Ann. Probab. 17, 1277 (1989).
[57] K. A. Seaton, J. Phys. A 34, L759 (2001).
[58] G.
Delﬁno
and
J.
L.
Cardy,
Nucl. Phys. B 519, 551 (1998).
[59] J. Cardy and R. M. Ziﬀ, J. Stat. Phys. 110, 1 (2003).
[60] R.
M.
Ziﬀ,
C.
D.
Lorenz,
and
P.
Kleban,
Physica A 266, 17 (1999).
[61] H.
N.
V.
Temperley
and
E.
H.
Lieb,
Proc. R. Soc. Lond. A 322, 251 (1971).
[62] R. M. Ziﬀ,
S. R. Finch,
and V. S. Adamchik,
Phys. Rev. Lett. 79, 3447 (1997).
[63] R. J. Baxter, H. N. V. Temperley,
and S. E. Ashley,
Proc. R. Soc. Lond. A 358, 535 (1978).
[64] M. B. Isichenko, Rev. Mod. Phys. 64, 961 (1992).
[65] M. Sahimi, Rev. Mod. Phys. 65, 1393 (1993).
[66] A. Weinrib, Phys. Rev. B 29, 387 (1984).
[67] S. Prakash, S. Havlin, M. Schwartz, and H. E. Stanley,
Phys. Rev. A 46, R1724 (1992).
[68] J.
Schmittbuhl,
J.-P.
Vilotte,
and
S.
Roux,
J. Phys. A 26, 6115 (1993).
[69] M. Sahimi, J. Phys. I France 4, 1263 (1994).
[70] C.
Du,
C.
Satik,
and
Y.
C.
Yortsos,
AIChE Journal 42, 2392 (1996).
[71] H.
A.
Makse,
S.
Havlin,
and
H.
E.
Stanley,
Nature 377, 608 (1995).
[72] H. A. Makse, J. S. Andrade Jr., M. Batty, S. Havlin,
and H. E. Stanley, Phys. Rev. E 58, 7054 (1998).
[73] M. Sahimi, Phys. Rep. 306, 213 (1998).
[74] H. A. Makse, J. S. Andrade Jr.,
and H. E. Stanley,
Phys. Rev. E 61, 583 (2000).
[75] A.
D.
Ara´ujo,
A.
A.
Moreira,
H.
A.
Makse,
H.
E.
Stanley,
and
J.
S.
Andrade
Jr.,
Phys. Rev. E 66, 046304 (2002).
[76] A. D. Ara´ujo, A. A. Moreira, R. N. Costa Filho,
and
J. S. Andrade Jr., Phys. Rev. E 67, 027102 (2003).
[77] N.
Sandler,
H.
R.
Maei,
and
J.
Kondev,
Phys. Rev. B 70, 045309 (2004).
[78] K. J. Schrenk, N. Pos´e, J. J. Kranz, L. V. M. van
Kessenich, N. A. M. Ara´ujo,
and H. J. Herrmann,
Phys. Rev. E 88, 052102 (2013).
[79] H. Auradou, K. J. M˚aløy, J. Schmittbuhl, A. Hansen,
and D. Bideau, Phys. Rev. E 60, 7224 (1999).
[80] J. Schmittbuhl, A. Hansen, H. Auradou,
and K. J.
M˚aløy, Phys. Rev. E 61, 3985 (2000).
[81] J. Kalda, EPL 84, 46003 (2008).
[82] J.
Kondev,
C.
L.
Henley,
and
D.
G.
Salinas,
Phys. Rev. E 61, 104 (2000).
[83] K. J. Schrenk, N. A. M. Ara´ujo, J. S. Andrade Jr., and
H. J. Herrmann, Sci. Rep. 2, 348 (2012).
[84] M.
E.
J.
Newman
and
R.
M.
Ziﬀ,
Phys. Rev. E 64, 016706 (2001).
[85] H.
Hu,
H.
W.
J.
Bl¨ote,
and
Y.
Deng,
J. Phys. A 45, 494006 (2012).
[86] H.-O. Peitgen and D. Saupe, eds., The Science of Fractal
Images (Springer, New York, 1988).
[87] K. B. Lauritsen, M. Sahimi,
and H. J. Herrmann,
Phys. Rev. E 48, 1272 (1993).
[88] H. A. Makse, S. Havlin, M. Schwartz, and H. E. Stanley,
Phys. Rev. E 53, 5445 (1996).
[89] H.
G.
Ballesteros
and
G.
Parisi,
Phys. Rev. B 60, 12912 (1999).
[90] B.
D.
Malamud
and
D.
L.
Turcotte,
J. Stat. Plan. Infer. 80, 173 (1999).
[91] E. A. Oliveira, K. J. Schrenk, N. A. M. Ara´ujo,
H.
J.
Herrmann,
and
J.
S.
Andrade
Jr.,
Phys. Rev. E 83, 046113 (2011).
[92] B.
Ahrens
and
A.
K.
Hartmann,
Phys. Rev. B 84, 144202 (2011).
[93] P. A. Morais,
E. A. Oliveira,
N. A. M. Ara´ujo,
H.
J.
Herrmann,
and
J.
S.
Andrade
Jr.,
Phys. Rev. E 84, 016102 (2011).
[94] A.
Weinrib
and
B.
I.
Halperin,
Phys. Rev. B 27, 413 (1983).
[95] W.
Janke
and
M.
Weigel,
Phys. Rev. B 69, 144208 (2004).
[96] I.
Mandre
and
J.
Kalda,
Eur. Phys. J. B 83, 107 (2011).
[97] B. Duplantier, Phys. Rev. Lett. 84, 1363 (2000).
[98] P. Grassberger, Math. Boisci. 63, 157 (1983).
[99] D. Mollison, J. R. Statist. Soc. 39, 283 (1977).
[100] J. D. Murray, Mathematical Biology, 3rd ed. (Springer,
Berlin, 2005).
[101] F. J. Dyson, Commun. Math. Phys. 91, 212 (1969).
[102] P. W. Anderson, G. Yuval, and D. R. Hamann, Phys.
Rev. B 1, 4464 (1970).
[103] D. J. Thouless, Phys. Rev. 187, 732 (1969).
[104] J. M. Kosterlitz and D. J. Thouless, J. Phys. C 6, 1181
(1973).
[105] M. Aizenman and C. M. Newman, Commun. Math.
Phys. 107, 611 (1986).
[106] P. Grassberger, J. Stat. Mech. , P04004 (2013).
[107] S.
Boettcher,
V.
Singh,
and
R.
M.
Ziﬀ,
Nat. Commun. 3, 787 (2012).
[108] F. Linder, J. Tran-Gia, S. R. Dahmen, and H. Hinrich-
sen, J. Phys. A 41, 185005 (2008).
[109] E.
Luijten
and
H.
W.
J.
Bl¨ote,
Phys. Rev. Lett. 89, 025703 (2002).
[110] P. Grassberger, J. Stat. Phys. 153, 289 (2013).
[111] M. Picco, arXiv:1207.1018.
[112] T. Blanchard,
M. Picco,
and M. A. Rajabpour,
EPL 101, 56003 (2013).
[113] N.
A.
M.
Ara´ujo
and
H.
J.
Herrmann,
Phys. Rev. Lett. 105, 035701 (2010).

<!-- page 10 -->

10
[114] Y. S. Cho, S. Hwang, H. J. Herrmann, and B. Kahng,
Science 339, 1185 (2013).
[115] R. A. da Costa, S. N. Dorogovtsev, A. V. Goltsev, and
J. F. F. Mendes, Phys. Rev. Lett. 105, 255701 (2010).
[116] O. Riordan and L. Warnke, Science 333, 322 (2011).
[117] P. Grassberger, C. Christensen, G. Bizhani, S.-W. Son,
and M. Paczuski, Phys. Rev. Lett. 106, 225701 (2011).
[118] G. Bizhani, V. Sood, M. Paczuski, and P. Grassberger,
Phys. Rev. E 83, 036110 (2011).
[119] S.-W. Son, G. Bizhani, C. Christensen, P. Grassberger,
and M. Paczuski, Europhys. Lett. 95, 58007 (2011).
[120] C. Christensen, G. Bizhani, S.-W. Son, M. Paczuski,
and P. Grassberger, EPL 97, 16004 (2012).
[121] G. Bizhani, P. Grassberger,
and M. Paczuski, Phys.
Rev. E 84, 066111 (2011).
[122] H.
W.
Lau,
M.
Paczuski,
and
P.
Grassberger,
Phys. Rev. E 86, 011118 (2012).
[123] D.
S.
Callaway,
J.
E.
Hopcroft,
J.
M.
Klein-
berg,
M.
E.
J.
Newman,
and
S.
H.
Strogatz,
Phys. Rev. E 64, 041902 (2001).
[124] A.-L. Barab´asi and R. Albert, Science 286, 509 (1999).
[125] J. Chalupa, P. L. Leath, and G. R. Reich, J. Phys. C
12, L31 (1979).
[126] J. Adler, Physica A 171, 453 (1991).
[127] P. S. Dodds and D. J. Watts, Phys. Rev. Lett. 92,
218701 (2004).
[128] P. S. Dodds and D. J. Watts, J. Theor. Biology 42, 232
(2005).
[129] H.-K. Janssen, M. M¨uller, and O. Stenull, Phys. Rev.
E 70, 026114 (2004).
[130] G.
Bizhani,
M.
Paczuski,
and
P.
Grassberger,
Phys. Rev. E 86, 011128 (2012).
[131] A.-L. Barab´asi and H. E. Stanley, Fractal Concepts
in Surface Growth (Cambridge University Press, New
York, 1995).
[132] A. V. Goltsev, S. N. Dorogovtsev, and J. F. F. Mendes,
Phys. Rev. E 73, 056101 (2006).
[133] P. Le Doussal, K. J. Wiese, and P. Chauve, Phys. Rev.
B 66, 174201 (2002).
[134] B. Drossel and K. Dahmen, Eur. Phys. J. B 3, 485
(1998).
[135] N. J. Zhou and B. Zheng, Phys. Rev. E 82, 031139
(2010).
[136] X. P. Qin, B. Zheng,
and N. J. Zhou, J. Phys. A 45,
345005 (2012).
[137] M. Aizenman and J. Wehr, Phys. Rev. Lett. 62, 2503
(1989).
[138] J. P. Sethna, K. A. Dahmen,
and O. Perkovi´c,
arXiv:cond-mat/0406320.
[139] M. Singer, Introduction to Syndemics: A Critical Sys-
tems Approach to Public and Community Health (John
Wiley & Sons, Hoboken, 2009).
[140] J. F. Brundage and J. D. Shanks, Emerging Infectious
Diseases 14, 1193 (2008).
[141] W. Oei and H. Nishiura, Comput. Math. Methods in
Med. 2012, 124861 (2012).
[142] M. S. Sulkowski, J. Hepatology 48, 353 (2008).
[143] S. K. Sharma, A. Mohan,
and S. Kadhivaran, Indian
J. Med. Res. 121, 550 (2005).
[144] L. Chen, W. Cai, F. Ghanbarnejad, and P. Grassberger,
Europhys. Lett. 104, 50001 (2013).
[145] W. Cai, L. Chen, F. Ghanbarnejad, and P. Grassberger,
(2014), to be published.
[146] S. V. Buldyrev, R. Parshani, G. Paul, H. E. Stanley,
and S. Havlin, Nature 464, 1025 (2010).
[147] R. Parshani, S. V. Buldyrev, and S. Havlin, Phys. Rev.
Lett. 105, 048701 (2010).
[148] J. Gao, S. V. Buldyrev, H. E. Stanley, and S. Havlin,
Nat. Phys. 8, 40 (2011).
[149] S.-W. Son, G. Bizhani, C. Christensen, P. Grassberger,
and M. Paczuski, Europhys. Lett. 97, 16006 (2012).
[150] B. Bollob´as, Graph Theory and Combinatorics: Proc.
Cambridge Combinatorial Conference in Honour of
Paul Erd˝os (Academic Press, New York, 1984) p. 35.
[151] S. N. Dorogovtsev, A. V. Goltsev, and J. F. F. Mendes,
Phys. Rev. Lett. 96, 040601 (2006).
[152] R. M. Ziﬀ, Phys. Rev. Lett. 103, 045701 (2009).
[153] Y. S. Cho, J. S. Kim, J. Park, B. Kahng, and D. Kim,
Phys. Rev. Lett. 103, 135702 (2009).
[154] E. J. Friedman and A. S. Landsberg, Phys. Rev. Lett.
103, 255701 (2009).
[155] R. M. D’Souza and M. Mitzenmacher, Phys. Rev. Lett.
104, 195702 (2010).
[156] J.
Nagler,
A.
Levina,
and
M.
Timme,
Nat. Phys. 7, 265 (2011).
[157] R. M. Ziﬀ, Science 339, 1159 (2013).
[158] Y. S. Cho and B. Kahng, arXiv:1404.4470.
[159] R. M. Ziﬀ, Phys. Rev. E 82, 051105 (2010).
[160] Y. S. Cho, B. Kahng,
and D. Kim, Phys. Rev. E 81,
030103(R) (2010).
[161] T. Bohman, A. Frieze,
and N. C. Wormald, Random
Struct. Algorithms 25, 432 (2004).
[162] W.
Chen
and
R.
M.
D’Souza,
Phys. Rev. Lett. 106, 115701 (2011).
[163] W.
Chen,
Z.
Zheng,
and
R.
M.
D’Souza,
EPL 100, 66006 (2012).
[164] K. J. Schrenk,
A. Felder,
S. Deﬂorin,
N. A. M.
Ara´ujo,
R.
M.
D’Souza,
and
H.
J.
Herrmann,
Phys. Rev. E 85, 031103 (2012).
[165] K. Panagiotou, R. Sp¨ohel, A. Steger, and H. Thomas,
Electron. Notes Discrete Math. 38, 699 (2011).
[166] R. Zhang, W. Wei, B. Guo, Y. Zhang, and Z. Zheng,
Physica A 392, 1232 (2013).
[167] W.
Chen,
J.
Nagler,
X.
Cheng,
X.
Jin,
H.
Shen,
Z.
Zheng,
and
R.
M.
D’Souza,
Phys. Rev. E 87, 052130 (2013).
[168] J. G´omez-Garde˜nes,
S. G´omez,
A. Arenas,
and
Y. Moreno, Phys. Rev. Lett. 106, 128701 (2011).
[169] P. Echenique, J. G´omez-Garde˜nes, and Y. Moreno, Eu-
rophys. Lett. 71, 325 (2005).
[170] C. L. Knecht, W. Trump, D. ben-Avraham, and R. M.
Ziﬀ, Phys. Rev. Lett. 108, 045703 (2012).
[171] S.
K.
Baek
and
B.
J.
Kim,
Phys. Rev. E 85, 032103 (2012).
[172] K. J. Schrenk, N. A. M. Ara´ujo, R. M. Ziﬀ, and H. J.
Herrmann, arXiv:1403.2082v1.
[173] E. Fehr, J. S. Andrade Jr., S. D. da Cunha, L. R. da
Silva, H. J. Herrmann, D. Kadau, C. F. Moukarzel, and
E. A. Oliveira, J. Stat. Mech. , P09007 (2009).
[174] E. Fehr, K. J. Schrenk, N. A. M. Ara´ujo, D. Kadau,
P. Grassberger, J. S. Andrade Jr., and H. J. Herrmann,
Phys. Rev. E 86, 011117 (2012).
[175] E.
Daryaei,
N.
A.
M.
Ara´ujo,
K.
J.
Schrenk,
S.
Rouhani,
and
H.
J.
Herrmann,
Phys. Rev. Lett. 109, 218701 (2012).
[176] E. Fehr, D. Kadau, J. S. Andrade Jr., and H. J. Her-
rmann, Phys. Rev. Lett. 106, 048501 (2011).

<!-- page 11 -->

11
[177] K. J. Schrenk, N. A. M. Ara´ujo, and H. J. Herrmann,
Sci. Rep. 2, 751 (2012).
[178] R.
D.
Schram
and
G.
T.
Barkema,
J. Stat. Mech. , P03009 (2012).
[179] S.-C. Park, arXiv:1202.6024v2.
[180] J. Marro and R. Dickman, Nonequilibrium Phase Tran-
sitions in Lattice Models (Cambridge University Press,
Cambridge, 1998).
[181] S.-C. Park, J. Korean Phys. Soc. 62, 469 (2013).
