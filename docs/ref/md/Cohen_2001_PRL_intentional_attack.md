# Breakdown of the Internet under Intentional Attack

> Phys. Rev. Lett. 86, 3682 (2001) · [arXiv](https://arxiv.org/abs/cond-mat/0010251)
> Local copy: `docs/ref/Cohen_2001_PRL_intentional_attack.pdf` (4 pages)

<!-- page 1 -->

arXiv:cond-mat/0010251v2  [cond-mat.dis-nn]  29 Mar 2001
Breakdown of the Internet under intentional attack
Reuven Cohen1 ∗, Keren Erez1, Daniel ben-Avraham2, and Shlomo Havlin1
1Minerva Center and Department of Physics, Bar-Ilan university, Ramat-Gan, Israel
2Department of Physics,Clarkson University, Potsdam NY 13699-5820, USA
We study the tolerance of random networks to intentional attack, whereby a fraction p of the
most connected sites is removed. We focus on scale-free networks, having connectivity distribution
P(k) ∼k−α (where k is the site connectivity), and use percolation theory to study analytically and
numerically the critical fraction pc needed for the disintegration of the network, as well as the size
of the largest connected cluster. We ﬁnd that even networks with α ≤3, known to be resilient to
random removal of sites, are sensitive to intentional attack. We also argue that, near criticality, the
average distance between sites in the spanning (largest) cluster scales with its mass, M, as
√
M,
rather than as logk M, as expected for random networks away from criticality. Thus, the disruptive
eﬀects of intentional attack become relevant even before the critical threshold is reached.
02.50.Cw, 05.40.a, 05.50.+q, 64.60.Ak
The question of stability of scale-free random networks
to removal of a fraction of their sites, especially in the
context of the Internet, has been recently of interest [1–3].
The Internet can be viewed as a special case of a ran-
dom, scale-free network, where the probability of a site
to be connected to k other sites follows a power-law:
P(k) ∼k−α (α ≈2.5, for the Internet). It is now well
established that if a fraction p of the sites is removed ran-
domly, then for α > 3 there exists a critical threshold,
pc, such that for p > pc the network disintegrates; net-
works with α ≤3 are more resilient and do not undergo
this transition, although ﬁnite networks (such as the In-
ternet) may be eventually disrupted when nearly all of
their sites are removed, as shown numerically in [1,2],
and analytically in [2].
Albert et al., [1] have introduced a model for inten-
tional attack, or sabotage of random networks: the re-
moval of sites is not random, but rather sites with the
highest connectivity are targeted ﬁrst. Their numerical
simulations suggest that scale-free networks are highly
sensitive to this kind of attack. In this Letter we study
the problem of intentional attack in scale-free networks.
Our study focuses on the exact value of the critical frac-
tion needed for disruption and the size of the remaining
largest connected cluster.
We also study the distance
between sites on this cluster near the transition.
We
ﬁnd, both analytically and numerically, that scale-free
networks are highly sensitive to sabotage of a small frac-
tion of the sites, for all values of α, lending support to
the view of Albert et al., [1].
In a recent paper [2] we have studied the properties of
the percolation phase transition in scale-free random net-
works, and applied a general criterion for the existence
of a spanning cluster (a cluster whose size is proportional
to the size of the network) [4,2]:
κ ≡⟨k2⟩
⟨k⟩= 2 .
(1)
Here k is the site connectivity, and averages, indicated by
angular brackets, are taken over all sites of the network.
When a fraction p of the sites are randomly removed (or
a fraction p of the links are removed, or lead to deleted
sites), the distribution of site connectivity is changed
from the original P(k) to a new distribution ˜P(k):
˜P(k) =
K
X
k0≥k
P(k0)
k0
k

(1 −p)kpk0−k .
(2)
Using this criterion together with Eq. (1), the critical
threshold p = pc is found to be:
1 −pc =
1
κ0 −1 ,
(3)
where κ0 = ⟨k2
0⟩/⟨k0⟩is calculated from the original con-
nectivity distribution, before the removal of any sites [2].
A wide range of networks, including the Internet, have
site connectivities which follow a power-law distribu-
tion [1,5,6]:
P(k) = ck−α,
k = m, m + 1, ..., K,
(4)
where k = m is the minimal connectivity and k = K is
an eﬀective connectivity cutoﬀpresent in ﬁnite networks.
For the distribution (4), κ0 can be approximated by [6]:
κ0 =
2 −α
3 −α
K3−α −m3−α
K2−α −m2−α .
(5)
This, together with Eq. (2), was used to show that net-
works with α ≤3, which have a divergent second mo-
ment, are resilient to random deletion of sites [2]. Indeed,
∗e-mail: cohenr@shoshi.ph.biu.ac.il
1

<!-- page 2 -->

when the number of sites in such networks N →∞, then
the upper cutoﬀK →∞, and there exists a spanning
cluster for all values of p < 1. Another approach, based
on generating functions, was introduced in [8] and was
used to study a similar problem in [3].
Consider now intentional attack, or sabotage [1],
whereby a fraction p of the sites with the highest con-
nectivity is removed.
(The links emanating from the
sites are removed as well.)
This has the following ef-
fect: (a) the cutoﬀconnectivity K reduces to some new
value, ˜K < K, and (b) the connectivity distribution of
the remaining sites is no longer scale-free, but is changed,
because of the removal of many of their links. The upper
cutoﬀK before the attack may be estimated from
∞
X
k=K
P(k) = 1
N ,
(6)
where N is the total number of sites in the network. Simi-
larly, the new cutoﬀ˜K, after the attack, can be estimated
from
K
X
k= ˜
K
P(k) =
∞
X
k= ˜
K
P(k) −1
N = p .
(7)
If the size of the system is large, N ≫1/p, the original
cutoﬀK may be safely ignored. We can then obtain ˜K
approximately by replacing the sum with an integral [6]:
˜K = mp1/(1−α) .
(8)
We estimate the impact of the attack on the distribu-
tion of the remaining sites as follows. The removal of a
fraction p of the sites with the highest connectivity re-
sults in a random removal of links from the remaining
sites — links that had connected the removed sites with
the remaining sites. The probability ˜p of a link leading
to a deleted site equals the ratio of the number of links
belonging to deleted sites to the total number of links:
˜p =
K
X
k= ˜
K
kP(k)
⟨k0⟩,
(9)
where ⟨k0⟩is the initial average connectivity. With the
usual continuous approximation, and neglecting K, this
yields
˜p =
 ˜K
m
2−α
= p(2−α)/(1−α) ,
(10)
for α > 2. For α = 2, ˜p →1, since just a few nodes of
very high connectivity control the entire connectedness
of the system. Indeed, consider a ﬁnite system of N sites
and α = 2. The upper cutoﬀK ≈N must then be taken
into account, and approximating Eq. (9) by an integral
yields ˜p = ln(Np/m). That is, for α = 2, very small val-
ues of p are needed to destroy an arbitrarily large fraction
of the links as N →∞.
With these results known we can compute the eﬀect
of intentional attack, using the theory previously devel-
oped for random removal of sites [2].
Essentially, the
network after attack is equivalent to a scale-free network
with cutoﬀ˜K, that has undergone random removal of a
fraction ˜p of its sites. This can be seen as the result of
two processes: (a) Removal of the highest connectivity
sites reduces the upper cutoﬀ. Since this eﬀect changes
the connectivity distribution, κ0 needs to be recalculated
accordingly. (b) Removal of the links leading to the re-
moved sites. The probability of removing a link is ˜p —
the probability of a randomly chosen link to lead to one
of the removed sites — and all links have the same proba-
bility of being deleted. Since this eﬀect has the inﬂuence
on the probability distribution described in Eq. (2), the
result in Eq. (3) can be used, with ˜p replacing p. (Notice
that for random site deletion the probability of a link
leading to a deleted site is identical to the fraction of
deleted sites.)
Although the number of nodes removed in intentional
attack is diﬀerent than in the random breakdown model,
this aﬀects the size of the spanning cluster (see below)
but not the critical point. This is because the transition
point is deﬁned as the point where the spanning cluster
becomes a ﬁnite fraction of the whole network. A ﬁnite
fraction of the remaining nodes is also a ﬁnite fraction
of the original network, so the diﬀerence has no eﬀect on
pc.
We therefore use Eqs. (3) and (5), but with ˜p =
( ˜K/m)2−α and ˜K replacing pc and K. This yields the
equation:
( ˜K/m)2−α −2 = 2−α
3−αm[( ˜K/m)3−α −1] ,
(11)
which can be solved numerically to obtain ˜K(m, α), and
then pc(m, α) can be retrieved from Eq. (8). In Fig. 1 we
plot pc — the critical fraction of sites needed to be re-
moved in the sabotage strategy to disrupt the network —
computed in this fashion, and compared to results from
numerical simulations. A phase transition exists (at a
ﬁnite pc) for all α > 2. The decline in pc for large α is
explained from the fact that as α increases the spanning
cluster becomes smaller in size, even before attack. (Fur-
thermore, for m < 2 the original network is disconnected
for some large enough α.) The decline in pc as α →2
results from the critically high connectivity of just a few
sites: their removal disrupts the whole network.
This
was already argued in [1]. We note that for inﬁnite sys-
tems pc →0 as α →2. The critical fraction pc is rather
sensitive to the lower connectivity cutoﬀm. For larger
m (the case of m = 1 is shown in Fig. 1) the networks
are more robust, though they still undergo a transition
at a ﬁnite pc.
The size of the spanning cluster as a fraction of the
number of undeleted sites, P∞(p), can be calculated us-
ing the methods introduced in [7] and developed in [3,8].
2

<!-- page 3 -->

Following closely the derivation in [8], a generating func-
tion is built for the connectivity distribution:
G0(x) =
∞
X
k=0
P(k)xk.
(12)
The probability of reaching a site with connectivity k
by following a speciﬁc link is kP(k)/⟨k⟩[2–4,8], and the
corresponding generating function is
G1(x) =
P kP(k)xk−1
P kP(k)
= d
dxG0(x)/⟨k⟩.
(13)
Hence, the generating function for the probability of
reaching a branch of a given size by following a link is
H1(x) = xG1(H1(x)) ,
(14)
while the generating function for the size of a component
is
H0(x) = xG0(H1(x)) .
(15)
Then, P∞= 1−H0(1), since H0 contains only the ﬁnite-
size clusters. It follows that
P∞(p) = 1 −
∞
X
k=0
˜P(k)uk ,
(16)
where u ≡H1(1) is the smallest positive root (found nu-
merically) of
⟨k⟩u =
∞
X
k=0
kP(k)uk−1 .
(17)
We argue that the same holds true after attack, but the
sums in Eqs. (16) and (17) should run only up to k = ˜K,
and the original distribution P(k) should be replaced by
the new connectivity distribution [2]:
˜P(k) =
˜
K
X
k0=m
P(k0)
k0
k

(1 −˜p)k ˜pk0−k .
(18)
The actual fraction of removed sites, p, is inconsequential
since we calculate the size of the inﬁnite cluster relative
to the number of undeleted sites. P∞evaluated in this
fashion, and compared to numerical simulations, is shown
in Fig. 2.
Since intentional attack leads always to a ﬁnite cutoﬀ
which does not scale with the system size, all the mo-
ments are ﬁnite and Eqs. (16) and (17) are well behaved.
Therefore, there always exists a linear term in the se-
ries expansion of P∞(p). Hence, near the critical point
P∞∼|pc −p|β, where β = 1 [10].
Finally, we consider the distance between sites in the
spanning cluster. The behavior of this quantity in diluted
networks is diﬀerent from that at the highly connected
regime. The average distance between two sites in the
spanning cluster of a highly connected network is propor-
tional to logk N, where k is the average connectivity and
N is the number of sites [9]. This has also been shown
to hold for scale-free and general networks [8].
How-
ever, the diluted case is essentially the same as inﬁnite-
dimensional percolation. In this case, there is no notion
of geometrical distance (since the graph is not embedded
in an Euclidean space), but only of a distance along the
graph (which is the shortest distance along bonds). It is
known from inﬁnite-dimensional percolation theory that
the fractal dimension at criticality is df = 2 [11]. There-
fore the average (chemical) distance d between pairs of
sites on the spanning cluster at criticality behaves as
d ∼
√
M ,
(19)
where M is the number of sites in the spanning clus-
ter. This is analogous to percolation in ﬁnite dimensions,
where in lengthscales smaller than the correlation length
the cluster is a fractal with dimension df and above the
correlation length the cluster is homogeneous and has
the dimension of the embedding space. In our inﬁnite-
dimensional case, the crossover between these two behav-
iors occurs around the correlation length ξ ≈|pc −p|, as
can be seen in Fig. 3.
In summary, we have shown that scale-free networks
are highly sensitive to intentional attack. This is true
even for networks with α < 3, which are known to be
resilient to random removal of their sites. The high sen-
sitivity near α = 2 results from the presence of just a
few sites with connectivity comparable to the size of the
system: their removal disrupts the whole network. We
note that while the cutoﬀ˜K must reach a typically small
number before the network is disrupted, this is achieved
with a modest removal fraction p (Eq. 8). The eﬀect of
sabotage on the connectivity distribution of the remain-
ing sites after the attack and thus the relation between
˜p, ˜K, and p (Eqs. (7) and (9) ) is found explicitly in our
approach. This eﬀect is particularly important near the
borderline case of α = 2. We have also shown that the
average distance between pairs of sites in the spanning
cluster grows dramatically near criticality. This makes
communication very ineﬃcient, even before the spanning
cluster is completely disrupted.
ACKNOWLEDGMENTS
We thank the National Science Foundation for support,
under grant PHY-9820569 (D.b.-A.).
3

<!-- page 4 -->

2
2.5
3
3.5
4
α
0
0.02
0.04
0.06
0.08
pc
FIG. 1. Critical probability, pc, as a function of α, for
networks of size N = 500, 000 (circles) and N = 64, 000
(squares). Lines represent the analytical solution, obtained
from Eqs. (7) and (9).
0
0.05
0.1
p
0
0.2
0.4
0.6
0.8
1
P∞(p)
FIG. 2. Fraction of sites belonging to the spanning clus-
ter, P∞, as a function of the fraction of removed sites, p,
for networks with α = 2.5 (circles), α = 2.8 (squares), and
α = 3.3 (diamonds).
Lines represent the analytical result
from Eqs. (16) and (17). Both the simulation and analysis
are for system size N = 500, 000.
10
−2
10
−1
10
0
10
1
d/ξ
10
−2
10
0
10
2
10
4
10
6
M/ξ
2
slope=2
FIG. 3. Mass (number of sites), M, as a function of dis-
tance, d, on the spanning cluster. The correlation length is
ξ = |p −pc|−1. Note that for d/ξ < 1, the slope is 2, cor-
responding to the behavior in the critical regime, while for
d/ξ > 1, M grows exponentially with d, corresponding to the
well connected regime.
[1] R. Albert, H. Jeong, and A. L. Barab´asi, Nature, 406,
6794, 378 (2000).
[2] R. Cohen, K. Erez, D. ben-Avraham, and S. Havlin,
Phys. Rev. Lett. 85, 4626 (2000).
[3] D. S. Callaway, M. E. J. Newman, S. H. Strogatz, and
D. J. Watts, Phys. Rev. Lett. 85, 5468 (2000)..
[4] M. Molloy and B. Reed, Random Structures and Algo-
rithms 6, 161 (1995).
[5] M. Faloutsos, P. Faloutsos, and C. Faloutsos, ACM SIG-
COMM ’99 Comput. Commun. Rev. 29, 251 (1999).
[6] In practice, in our simulations we generate a continu-
ous variable k, with distribution P(k) = ck−α (k ≥m),
and pick up the closest integer as the actual connectiv-
ity. Thus, we really have P(k) =
R k+1/2
k−1/2 ck′−αdk′. This
has the eﬀect of rendering the continuum approximation,
used in our analysis, more accurate.
[7] M. Molloy and B. Reed, Combinatorics, Probability and
Computing 7, 295 (1998).
[8] M. E. J. Newman, S. H. Strogatz, and D. J. Watts cond-
mat/0007235.
[9] B. Bollob´as, Random Graphs. pp. 123-136 (Academic
Press, London, 1985).
[10] However, for random breakdown when the upper cutoﬀ
K →∞, there are diverging moments. Using Abelian
theories for calculating the sums in Eqs. (16) and (17),
it can be shown that for 3 < α < 4, β = 1/(α −3)
while for α > 4, β = 1 as in regular inﬁnite-dimensional
percolation.
[11] A. Bunde, and S. Havlin (editors), Fractals and Disor-
dered System (Springer, New York, 1996).
4
