# Sharp Metastability Threshold for Two-Dimensional Bootstrap Percolation

> Probab. Theory Relat. Fields 125, 195-224 (2003) · [arXiv](https://arxiv.org/abs/math/0206132)
> Local copy: `docs/ref/Holroyd_2003_PTRF_bootstrap_sharp_threshold.pdf` (40 pages)

<!-- page 1 -->

arXiv:math/0206132v1  [math.PR]  12 Jun 2002
Sharp Metastability Threshold for
Two-Dimensional Bootstrap Percolation
Alexander E. Holroyd
May 8, 2002 (updated May 29, 2018)
Abstract
In the bootstrap percolation model, sites in an L by L square are
initially independently declared active with probability p.
At each
time step, an inactive site becomes active if at least two of its four
neighbours are active. We study the behaviour as p →0 and L →
∞simultaneously of the probability I(L, p) that the entire square is
eventually active. We prove that I(L, p) →1 if lim inf p log L > λ,
and I(L, p) →0 if lim sup p log L < λ, where λ = π2/18. We prove
the same behaviour, with the same threshold λ, for the probability
J(L, p) that a site is active by time L in the process on the inﬁnite
lattice.
The same results hold for the so-called modiﬁed bootstrap
percolation model, but with threshold λ′ = π2/6. The existence of
the thresholds λ, λ′ settles a conjecture of Aizenman and Lebowitz [3],
while the determination of their values corrects numerical predictions
of Adler, Stauﬀer and Aharony [2].
1
Introduction
We consider the bootstrap percolation process in two dimensions. Let Z2 =
{x = (x1, x2) : x1, x2 ∈Z} be the set of all 2-vectors of integers. Elements of
Key words: bootstrap percolation, cellular automaton, metastability, ﬁnite-size scaling
2000 Mathematics Subject Classiﬁcations: Primary 60K35; Secondary 82B43
Address: UCLA Department of Mathematics, CA 90095-1555, USA
Research funded in part by NSF Grant DMS–0072398
1

<!-- page 2 -->

Z2 are called sites. The neighbourhood N(x) of a site is
N(x) = {y ∈Z2 : ∥x −y∥= 1},
where ∥· ∥is Euclidean distance. So |N(x)| = 4 for all x. Let K be a subset
of Z2. Deﬁne B(K) by
B(K) = K ∪{x ∈Z2 : |N(x) ∩K| ≥2},
and ⟨K⟩by
⟨K⟩= lim
t→∞Bt(K),
where Bt denotes the t-th iterate. (If K is the active set at time 0, then Bt
is the active set at time t and ⟨K⟩is the active set at time ∞).
Now ﬁx p ∈[0, 1] and let X be a random subset of Z2 in which each site is
independently included with probability p. More formally, denote by Pp the
product probability measure with parameter p on the product σ-algebra of
{0, 1}Z2, and deﬁne the random variable X by X(ω) = {x ∈Z2 : ω(x) = 1}
for ω ∈{0, 1}Z2. A site x ∈Z2 is said to be occupied if x ∈X.
We say that a set K ⊆Z2 is internally spanned if ⟨X ∩K⟩= K. A
rectangle is a set of sites of the form
R(a, b; c, d) := {a, . . . , c} × {b, . . . , d},
where a ≤c, b ≤d are integers. We also write R(c, d) = R(1, 1; c, d). Deﬁne
the function
I(L, p) = Pp

R(L, L) is internally spanned

.
Our main result is the following.
Theorem 1 Let Ln, pn be sequences such that Ln →∞and pn →0. Then:
(i) if lim inf
n→∞
pn log Ln > λ then lim
n→∞I(Ln, pn) = 1;
(ii) if lim sup
n→∞
pn log Ln < λ then lim
n→∞I(Ln, pn) = 0;
where
λ = π2
18.
The main step in the proof of Theorem 1 will be Theorem 2 below, concerning
the probability of much smaller squares being internally spanned.
2

<!-- page 3 -->

Theorem 2
(i) lim sup
p→0
sup
m≥1
−p log I(m, p) ≤2λ.
(ii) lim inf
B→∞lim inf
p→0
−p log I(⌊B/p⌋, p) ≥2λ.
(Here ⌊·⌋denotes the integer part).
We also establish the following result about the time-evolution of the
process on Z2. Deﬁne
J(t, p) = Pp

(0, 0) ∈Bt(X)

,
(the probability the origin is active by time t).
Theorem 3 Let tn, pn be sequences such that tn →∞and pn →0. Then:
(i) if lim inf
n→∞
pn log tn > λ then lim
n→∞J(tn, pn) = 1;
(ii) if lim sup
n→∞
pn log tn < λ then lim
n→∞J(tn, pn) = 0.
The modiﬁed bootstrap model is a variant of the above model, in which
the deﬁnition of B is replaced with
B′(K) = K ∪
n
x ∈Z2 : |{x + ei, x −ei} ∩K| ≥1 for each of i = 1, 2
o
,
where e1 = (1, 0) and e2 = (0, 1). We deﬁne I′(·, ·) and J′(·, ·) accordingly.
The arguments in this article can be used to prove the following.
Theorem 4 For the modiﬁed bootstrap model, the analogues of Theorems
1,2,3 hold, with λ replaced by
λ′ = π2
6 .
Our results answer the main question posed in [3] in the case of two-
dimensional bootstrap percolation. In that paper is was proved (for a wider
class of models) that for p →0 and L →∞, we have I(L, p) →1 if
lim inf p log L > c1 and I(L, p) →0 if lim sup p log L < c2, for diﬀerent
constants c1, c2, and similarly for J(L, p). Other aspects of the model have
been subsequently studied in detail (see [4],[10] and the references therein, for
3

<!-- page 4 -->

example), but the natural question of whether c1, c2 could be replaced with a
single sharp threshold remained open. Our results answer this aﬃrmatively,
as well as establishing the precise value of the threshold.
Predictions for the thresholds λ, λ′ based on simulation are not in good
agreement with our rigorous result. In [1],[2], the estimates 0.245 ± 0.015
for λ and 0.47 ± 0.02 for λ′ are reported, whereas λ = π2/18 = 0.548311 · · ·
and λ′ = π2/6 = 1.644934 · · ·. The likely reason for this discrepancy is that
it is necessary to take L extremely large in order to “see” the true limiting
behaviour. (These simulations used values of L up to 28, 800). Similar phe-
nomena have been noted for several other variants of bootstrap percolation;
for details see [8],[13],[14],[17],[18].
Bootstrap percolation in three and higher dimensions presents a new
set of challenges; for details see [5],[6],[14].
In particular for the three-
dimensional model in which a site becomes active if at least three of its six
neighbours are active, it was established in [5] that the form of the threshold
regime is diﬀerent: I(L, p) →1 if lim inf p log log L > c1 and I(L, p) →0 if
lim sup p log log L < c2; the extension to other dimensions is treated in [6].
The arguments required for these results are more sophisticated than those
in [3], and it seems that the problem of establishing a sharp threshold here
is likely to be correspondingly harder.
Aside from their intrinsic mathematical interest, bootstrap percolation
models ﬁnd numerous applications, both directly and as tools in the analysis
of more complicated systems. See for example [1],[7],[8],[12],[19].
In [3] it was also conjectured that stochastic Ising models show similar
metastability behaviour. The analogue of the result in [3] was proved (in
arbitrary dimension) in [15], and the analogue of our result was proved in
[16].
We omit the proof of Theorem 4 regarding the modiﬁed bootstrap model.
The proof is almost identical to those given for the earlier theorems, and in
fact a few simpliﬁcations are possible. The main diﬀerence is that the deﬁni-
tion of a rectangle being “horizontally (respectively vertically) traversable”
in Section 3 should be replaced with the statement that all the columns (re-
spectively rows) are occupied, and the function g (see Sections 2, 3, 4) should
be replaced with f.
The article is organized as follows. In Section 2 we introduce and solve a
deﬁnite integral which gives rise to the constant λ = π2/18. In Section 3 we
give notation and basic results. In Section 4 we prove Theorem 2 (i), and in
Section 5 we deduce Theorems 1 and 3 from Theorem 2. These proofs are
4

<!-- page 5 -->

relatively standard; the rest of the article is devoted to the much harder task
of proving Theorem 2 (ii). In Sections 6, 7, 8 and 9 we assemble the tools
needed for this; the reader may prefer to skip the proofs in these sections on
a ﬁrst reading. Finally in Section 10 we complete the proof of Theorem 2
(ii).
We now present a sketch of the ideas behind the proofs in this article.
Theorem 1 may be deduced from Theorem 2 as follows. Roughly speaking,
Theorem 2 states that if we take B large and then let p →0, then I(B/p, p)
behaves roughly as e−2λB/p, where λB →λ as B →∞. Now suppose L ≈ec/p
for some c. It may be shown that the event that R(L, L) is internally spanned
is roughly equivalent to the event that it contains some internally spanned
square of side length B/p. This is because for B large, most intervals of
length B/p contain an occupied site, so when a growing active cluster reaches
approximately this size it is almost certain to grow to ﬁll the whole square
(such a growing cluster is sometimes referred to as a “critical droplet” in the
literature). The expected number of internally spanned squares of side B/p
in R(L, L) is approximately L2e−2λB/p ≈e2(c−λB)/p, which is ≫1 if c > λB
and ≪1 if c < λB. However, we may choose B so that λB is as close to λ as
desired. A more careful application of these ideas proves Theorem 1, and a
similar approach gives Theorem 3.
We now turn to the proof of Theorem 2. The lower bound on I in (i) is
relatively straightforward; the upper bound in (ii) is much harder. For (i),
assume m = B/p and consider one way in which R(m, m) can be internally
spanned. Suppose that (1, 1) is occupied, and that for each k = 2, . . . , m,
the column with horizontal coordinate k and height k −1 contains at least
one occupied site, and similarly the row with vertical coordinate k and width
k −1. Then it is easily seen that each square with bottom-left corner (1,1)
is internally spanned, so the active set grows to ﬁll whole of R(m, m). This
picture gives the correct lower bound for the modiﬁed bootstrap model, but
for the bootstrap model one can do better. Since a site becomes active if it
has active neighbours on two opposite sides, we can allow some of the rows
and columns described above to be vacant. Roughly speaking, the growth
will still take place provided there are never two adjacent vacant rows or
columns (actually, a little more care is needed, see Section 4). A calculation
shows that the probability of this event behaves roughly as e−(2/p)
R B
0 g where
g is a certain function which satisﬁes
R ∞
0 g(z) dz = λ = π2/18. This proves
Theorem 2 (i), and also shows that in some sense the “natural” length unit
5

<!-- page 6 -->

for the problem is 1/p.
For the upper bound, Theorem 2 (ii), the basic idea is to consider all
other possible ways in which R = R(B/p, B/p) could be internally spanned,
and ﬁnd upper bounds on the number of ways and the probability of each
one.
One way in which R could be internally spanned is for the growth
described above to start from the center of the square R (say), and for an
active square to grow in all four directions rather than just two. It turns out
that the probability of this is essentially the same, e−(2/p)
R B
0 g, roughly because
growth by one unit each to the left and to the right is equivalent to growth
by two units to the right. Hence there is in some sense no loss of generality
in assuming that growth starts from the corner. However, there are many
other ways in which this growth could occur. For example two adjacent rows
as described above might be vacant, in which case the active region could
still continue to grow horizontally for a while, until it encounters an occupied
site in the vertical direction, at which point vertical growth can resume. It
can be shown that such a possibility has much smaller probability than the
one considered earlier. Indeed, each such growth history corresponds to an
oriented path γ in [0, ∞)2 (where the length is rescaled by a factor 1/p), and
it may be shown that the probability of such a history is roughly e−w(γ)/p
where w is a certain functional on paths deﬁned as a path integral involving
the function g. By the convexity of g, it can be shown (see Section 6) that
this functional is minimized along the main diagonal x1 = x2, with minimum
2
R B
0 g, which corresponds to growth as a square as described above.
At this point it is natural to try to get an upper bound on the probability
of R being internally spanned by summing the above probability over all
possible paths. However, the number of such paths is large, roughly 2B/p,
so this does not give the correct bound.
The solution is to introduce a
new “coarse-graining” length scale T/p, which is small compared with B/p
but large compared with the lattice spacing 1. It may be shown that the
probability of seeing some path which coincides with γ at a set of points
of spacing T is also approximately e−w(γ)/p. This is because all such paths
require roughly the same occupied sites. Actually, for this argument to work,
we must assume that the growth starts not from the point (1, 1) but from
a larger rectangle, of size A/p, where T ≪A ≪1 ≪B; the details are
in Section 7. The number of possible choices of a set of points at spacing
T/p is much less than before (in fact, bounded in p), so now summing over
all possibilities gives an upper bound of approximately e−(2/p)
R B
A g (and the
6

<!-- page 7 -->

integral here is close to
R ∞
0 g for A small and B large, since the latter integral
converges).
However, there is another diﬃculty. There are other ways for R to be
internally spanned which do not involve only growth starting from a single
“seed”. For example, the growing square described above might encounter
another internally spanned rectangle S within R, and the two would combine
to give a larger internally spanned rectangle without any of the intervening
growth taking place. It seems unlikely that such an event could have higher
probability, since the probability of a small internally spanned rectangle S
should be much less than the probability of ﬁnding at least one occupied
site in each of the corresponding rows and columns required for growth as
described previously. In fact, most of the work in the proof is to rule out
possibilities such as this.
The idea is as follows. Suppose R is internally spanned. Then we expect
to be able to ﬁnd a sequence of successively slightly smaller internally spanned
rectangles inside R (corresponding to the growth picture above). However,
at some point there may be a “split” into two separate internally spanned
rectangles. Now we apply the same reasoning to each of these, and so on.
In this way we obtain a “hierarchy” - a tree structure of nested rectangles
each of which is internally spanned.
We can arrange that when there is
no splitting, the sizes of consecutive rectangles diﬀer by approximately the
coarse-graining scale T/p; when there is a split, the two oﬀspring rectangles
might be much smaller than the parent, but they should have the property
that they “span” the parent, in the sense that bootstrap percolation starting
with the two oﬀspring rectangles completely ﬁlled results in the parent being
completely ﬁlled. Also, for the coarse-graining argument to work, we should
stop whenever the rectangles in a line of descent become too small (perhaps
smaller than A/p), and declare such a rectangle to be a seed.
Now, another diﬃculty arises.
Since there may be many splits, there
may also be many seeds; indeed R could be “almost all seeds” in which case
there is no room left for the growth estimates used previously. The solution
is to use the fact that seeds are small, together with an a priori bound on
the probability a small rectangle is internally spanned, to show that there
are “not too many” seeds. More precisely, we introduce yet another length
Z/p ≪A/p, and declare a rectangle a seed if its size is less than Z/p (the
exact meaning of “size” turns out to be a delicate issue here).
We show
that the probability that the total perimeter of all the seeds is more than
A/p is very small. Then, provided the total perimeter is less than A/p, a
7

<!-- page 8 -->

geometrical argument together with a variational principle can be used to
show that the probability of the whole hierarchy is approximately at most
e−(2/p) R B
A g. The order of choosing the constants is important here: they must
be chosen in the order B, A, Z, T, then ﬁnally p →0. The total number of
hierarchies is large, but it can be shown to be at most p−K where K depends
on B but not p. Hence summing gives a bound of the correct form.
At the heart of any such proof must be a (deterministic) result giving
rigorous necessary conditions for a rectangle to be internally spanned. In
[3], this role was played by a result which we shall also make use of, Lemma
10. For the construction of hierarchies, the required condition is provided
by an apparently new result, Proposition 31, which roughly speaking states
that any internally spanned rectangle must be spanned by a pair of smaller
rectangles which are themselves internally spanned disjointly (in the sense of
the Van den Berg-Kesten inequality).
2
Integrals
We deﬁne the functions f : (0, ∞) →(0, ∞), β : (0, 1) →(0, 1), and g :
(0, ∞) →(0, ∞) by
f(z) = −log(1 −e−z),
β(u) = u +
p
u(4 −3u)
2
,
g(z) = −log β(1 −e−z).
Note that f and g are continuously diﬀerentiable, positive, decreasing and
convex.
To see the latter for g note that the function −z 7→g(z) is the
composition of the increasing, convex functions z 7→ez, z 7→−β(1 −z) and
z 7→−log(−z). Observe also that β(u) > u so f(z) > g(z) for z ∈(0, ∞).
Finally note that f(z), g(z) →∞as z →0, and f(z), g(z), zf(z), zg(z) →0
as z →∞.
Proposition 5
(i)
Z ∞
0
f(z) dz = π2
6 ;
8

<!-- page 9 -->

(ii)
Z ∞
0
g(z) dz = π2
18.
We deﬁne λ to equal the second integral.
Proposition 5 (ii) is somewhat remarkable; the integral does not appear in
standard tables such as [9], and seems not to be directly solvable by computer
programs such as Mathematica. The proof which we shall give is mysterious,
and uses many special features of the function g.
Proof of Proposition 5.
In (i) we substitute x = e−z to obtain
Z ∞
0
−log(1 −e−z) dz =
Z 1
0
−log(1 −x)dx
x = π2
6 ,
by [9], number 4.291.2.
The integral in (ii) converges because g ≤f. Substituting y = e−z and
then taking out a factor 1 −y, we have
Z ∞
0
g(z) dz
=
Z ∞
0
−log 1 −e−z +
p
(1 −e−z)(1 + 3e−z)
2
dz
=
Z 1
0
−log 1 −y +
p
(1 −y)(1 + 3y)
2
dy
y
=
Z 1
0
−log(1 −y) dy
y +
Z 1
0
−log
1 +
q
1+3y
1−y
2
dy
y .
Note that the ﬁrst of these two integrals converges by the above, hence the
second must also.
In the ﬁrst integral we make the substitution x = 1 −y. In the second
we make the substitution
x =
2
1 +
q
1+3y
1−y
which yields
y =
1 −x
1 −x + x2
9

<!-- page 10 -->

and
dy
y =
(2 −x)x
(x −1)(1 −x + x2)dx,
and also interchanges the limits 0 and 1. Thus we obtain
Z ∞
0
g(z) dz
=
Z 1
0
−log x
dx
1 −x +
Z 1
0
−log x
(2 −x)x
(x −1)(1 −x + x2) dx
=
Z 1
0
log x

−
1
1 −x −
(2 −x)x
(x −1)(1 −x + x2)

dx
=
Z 1
0
log x
2x −1
1 −x + x2 dx
=
h
log x log(1 −x + x2)
i1
0 −
Z 1
0
log(1 −x + x2) dx
x
(by parts)
=
0 −
Z 1
0
log
1 + x3
1 + x
 dx
x
=
Z 1
0
log(1 + x) dx
x −
Z 1
0
log(1 + x3) dx
x
=
Z 1
0
log(1 + x) dx
x −
Z 1
0
log(1 + u) du
3u
(substituting u = x3)
=
2
3
Z 1
0
log(1 + x) dx
x
=
2
3
π2
12
(by [9], 4.291.1)
=
π2
18.
✷
3
Basic Notation and Results
In this section we introduce basic tools which will be used throughout.
The following elementary properties of ⟨·⟩will be fundamental: ⟨K⟩⊇K;
⟨⟨K⟩⟩= ⟨K⟩; if K ⊆K′ then ⟨K⟩⊆⟨K′⟩.
These have the important
consequence that if K ⊆K′ ⊆⟨K⟩then ⟨K′⟩= ⟨K⟩.
10

<!-- page 11 -->

Let R = R(a, b; c, d) be a rectangle. By the dimensions of R we mean the
2-vector
dim(R) = (c −a + 1, d −b + 1).
If dim(R) = (m, n) we deﬁne the short side short(R) = min{m, n}, the long
side long(R) = max{m, n} and the semi-perimeter φ(R) = m + n of R.
A site x ∈Z2 is occupied if x ∈X. A set of sites K ⊆Z2 is full if every
site in K is occupied, and occupied if at least one site in K is occupied. It
will be convenient to write
q = −log(1 −p).
Note that q ≥p and q ∼p (that is, q/p →1) as p →0. The advantage of
this notation is that
Pp(K is occupied) = 1 −e−|K|q = e−f(|K|q)
(6)
(where f was deﬁned in Section 2).
We say that a sequence of events (Aj) has a double gap if there is an
adjacent pair Ai, Ai+1 neither of which occurs.
Lemma 7 In a sequence of k independent events each with probability u ∈
(0, 1), the probability ak(u) that there are no double gaps satisﬁes
β(u)k ≤ak(u) ≤β(u)k−1
where
β(u) = u +
p
u(4 −3u)
2
.
Proof.
By induction on k, on noting that a0 = a1 = 1, ak+2 = uak+1 +
(1 −u)uak, 0 < β < 1, and β2 = uβ + (1 −u)u.
✷
Let R = R(a, b; c, d). For a ≤i ≤c (respectively b ≤j ≤d), column i (re-
spectively row j) of R is the rectangle R(i, b; i, d) (respectively R(a, j; c, j)).
We say that R is horizontally (respectively vertically) traversable if the se-
quence (column i is occupied)c
i=a (respectively (row j is occupied)d
j=b) has no
double gaps. The following deﬁnitions will also be convenient. The rectangle
R is East-traversable (respectively North-traversable) if it is horizontally (re-
spectively vertically) traversable and in addition column c (respectively row
d) is occupied.
Recall the deﬁnitions of f and g from Section 2.
11

<!-- page 12 -->

Lemma 8 If R is a rectangle with dimensions (m, n) then
(i)
e−mg(nq) ≤Pp(R is horizontally traversable) ≤e−(m−1)g(nq);
(ii)
e−mf(nq) ≤e−(m−1)g(nq)−f(nq) ≤Pp(R is East-traversable) ≤e−mg(nq).
And similar inequalities hold for vertical and North- traversability, with m
and n exchanged.
Proof.
Part (i) follows immediately from (6), Lemma 7 and the deﬁnition
of g. For the upper bound in (ii), note that if R(a, b; c, d) is East-traversable
then R(a, b; c + 1, d) is horizontally traversable, and use (i).
The second
inequality in (ii) is straightforward, and the ﬁrst follows because f ≥g.
✷
To see the usefulness of the above concepts note the following.
Lemma 9
(i) If R is internally spanned then R is East- and North- traversable.
(ii) If R1 = (a, b; c, d) is internally spanned and R2 = (c + 1, b; e, d)
is East-traversable then R1 ∪R2 is internally spanned. And a similar
statement holds for North-traversability.
Proof.
For part (i), note that if two adjacent columns or the East-most
column of R contains no occupied sites, then no site in these columns can be
in ⟨R ∩X⟩. A similar remark applies to rows. For (ii), it is easy to see that
each successive column in R2 (moving in the East direction) is in ⟨R ∩X⟩.
✷
An event A of {0, 1}Z2 is called increasing if whenever ω ∈A and ω′ ≥ω
we have ω′ ∈A. The FKG inequality states that for increasing events A, B
we have Pp(A ∩B) ≥Pp(A)Pp(B) (see [11] p34 for example).
12

<!-- page 13 -->

r
r
r
r
r
r
m
Figure 1: An illustration of the event A. The arrows indicate East- and
North- traversability.
4
Lower Bound
Proof of Theorem 2 (i).
Let r = ⌊p−1/2⌋, and let A = A(m, p) be the
event that all the following occur:
R(1, 1; 1, r) and R(1, 1; r, 1) are full,
the sites (m, 1) and (1, m) are occupied,
and for all integers k ≥1:
R(kr + 1, 1; kr + r, kr) is East-traversable,
R(1, kr + 1; kr, kr + r) is North-traversable.
See Figure 1.
Using Lemma 9 (ii), it is easily seen that if A occurs then R(m, m) is
internally spanned, hence using Lemma 8 (ii) and the FKG inequality,
I(m, p) ≥Pp(A) ≥p2r+1
∞
Y
k=1

e−(r−1)g(krq)−f(krq)2
,
13

<!-- page 14 -->

for all m. Hence,
sup
m≥1
−p log I(m, p)
≤
−(2r + 1)p log p + 2(r −1)p
rq
∞
X
k=1
g(krq)rq + 2p
rq
∞
X
k=1
f(krq)rq
≤
−3rp log p + 2p
q
Z ∞
0
g(z) dz + 2p
rq
Z ∞
0
f(z) dz
where in the last step we have use the fact that f and g are decreasing. Now,
recalling that q ∼p, r = ⌊p−1/2⌋, and
R ∞
0 f(z) dz < ∞, we see that as p →0
the above expression converges to 0 + 2
R ∞
0 g(z) dz + 0 = 2λ, as required. ✷
5
Metastability
In this section we deduce the (i) parts of Theorems 1 and 3 from the (i) part
of Theorem 2, and similarly for the (ii) parts.
Proof of Theorem 1 (i).
It is clearly suﬃcient to prove that for
any ǫ > 0, if pn →0 and Ln →∞are such that pn log Ln ≥λ + ǫ then
I(Ln, pn) →1.
Suppose p log L ≥λ + ǫ and p < 1/2. Let m = ⌊p−3⌋, and let S =
S(L, p) be the event that R(L, L) contains at least one internally spanned
rectangle of dimensions (m, m). By dividing R(L, L) into disjoint rectangles
of dimensions (m, m) we see that
Pp(S) ≥1 −

1 −I(m, p)
⌊L/m⌋2
,
so
−log(1 −Pp(S)) ≥L2p6
2
I(m, p).
Therefore
p log
h
−log(1 −Pp(S))
i
≥2p log L + 6p log p −p log 2 + p log I(m, p).
Hence, by Theorem 2 (i) we have
lim inf
n→∞pn log
h
−log(1 −Ppn(Sn))
i
≥2(λ + ǫ) + 0 −0 −2λ = 2ǫ,
14

<!-- page 15 -->

where Sn = S(Ln, pn), so in particular Ppn(Sn) →1.
The following is proved in [3]. There exists a sequence of increasing events
Hn deﬁned in terms of the states of sites in R(Ln, Ln) such that Ppn(Hn) →1,
and if Sn and Hn occur then R(Ln, Ln) is internally spanned. This proves
the result, since by the FKG inequality,
I(Ln, pn) ≥Ppn(Sn)Ppn(Hn) →1.
Here is a sketch of the construction of Hn. If pn log Ln →a > λ then
we may take Hn to be the event that every rectangle of dimensions (m, 1) or
(1, m) in R(Ln, Ln) is occupied. If Ln grows faster than this then we must
also use a renormalization argument, dividing R(Ln, Ln) into disjoint squares
of size e(λ+δ)/pn, so that the probability each one is internally spanned exceeds
the critical probability for site percolation on Z2.
✷
Proof of Theorem 3 (i).
Note that, in contrast with I, J(t, p) is clearly
increasing in t, so we may assume that pn log tn →λ + ǫ, where ǫ > 0. Let
mn = ⌊p−3
n ⌋and Ln = ⌊tn/(3mn)⌋; then pn log Ln →λ + ǫ. It is easily seen
that the events Hn in the proof of Theorem 1 (i) above may be chosen in
such a way that if Sn and Hn occur then the whole of R(Ln, Ln) becomes
active in time at most
|R(mn, mn)| + 2mnLn ≤tn
(The inequality holds if pn is suﬃciently small). Hence,
Pn

(1, 1) ∈Btn(X ∩R(Ln, Ln))

→1,
therefore J(tn, pn) →1.
✷
To prove Theorem 1 (ii) we need the following result from [3].
Lemma 10 Let k be a positive integer and let R be a rectangle with long(R)
≥2k + 1. If R is internally spanned then there exists an internally spanned
rectangle T ⊆R with long(T) ∈[k, 2k + 1].
The idea of the proof of Lemma 10 is to construct the bootstrap percola-
tion process by an algorithm which sequentially replaces a pair of internally
spanned rectangles with a larger internally spanned rectangle; the result then
follows because at each step the long side of the longest rectangle increases at
most as n 7→2n + 1. For the details see [3]. We shall use a similar argument
in Section 8 to prove Proposition 31.
15

<!-- page 16 -->

a
b
K
K
T
S1
S2
Figure 2: The rectangles T, S1, S2. The arrows indicate North- and East-
traversability.
Proof of Theorem 1 (ii).
It is clearly suﬃcient to show that for any
ǫ > 0, if pn →0 and Ln →∞with pn log Ln ≤λ −ǫ then I(Ln, pn) →0.
We write λB =
1
2 infB′≥B lim infp→0 −p log I(⌊B′/p⌋, p), so that Theorem 2
(ii) states limB→∞λB ≥λ.
Fix ǫ > 0 and B > 0, and suppose that p log L ≤λ −ǫ. We write R =
R(L, L) and let K = ⌊B/p⌋and k = ⌊B/(2p)⌋−1, so that k ≤2k + 1 ≤K.
Assume that p is suﬃciently small that 1 < k < K < L. We claim ﬁrst that
if T is any rectangle satisfying long(T) ∈[k, 2k + 1] then
Pp(T is internally spanned) ≤e2Kf(kq)I(K, p).
(11)
To prove this, suppose without loss of generality that T is of the form T =
R(a, b) where a ∈[k, 2k + 1] and b ≤K. Let S1 = R(1, b + 1; a, K) and S2 =
R(a + 1, 1; K, K); see Figure 2. By Lemma 9 (ii), if T is internally spanned,
and S1 is North-traversable, and S2 is East-traversable, then R(K, K) is
internally spanned. Hence by Lemma 8 (ii) we have
I(K, p) ≥Pp(T is internally spanned)e−(K−b)f(aq)−(K−a)f(Kq)
which yields (11) since f is decreasing.
Hence using Lemma 10 and (11) we have
I(L, p) ≤L2K2e2Kf(kq)I(K, p).
Here L2K2 is an upper bound on the number of possible choices for the
rectangle T. Hence
p log I(L, p) ≤2p log L + 2p log K + 2Kpf(kq) + p log I(K, p).
16

<!-- page 17 -->

Therefore, recalling that pn log Ln ≤λ −ǫ, and the deﬁnitions of k, K, λB,
we have
lim sup
n→∞
pn log I(Ln, pn)
≤
2(λ −ǫ) + 0 + 2Bf(B/2) −2λB.
Since limB→∞Bf(B/2) = 0 and limB→∞λB ≥λ, we may choose B suﬃ-
ciently large (depending on ǫ) that the above expression is negative, so in
particular I(Ln, pn) →0.
✷
Proof of Theorem 3 (ii).
In consequence of the results (5.4),(5.5) of
[3], if I(L, p) →0 then J(L, p) →0. Hence the result follows from Theorem
1 (ii).
✷
6
Variational Principles
We write R2
+ = (0, ∞)2 = {a = (a1, a2) : a1, a2 ∈(0, ∞)} for the set of
all 2-vectors of positive reals. Boldface letters will denote elements of R2
+
unless stated otherwise. We write a ≤b if a1 ≤b1 and a2 ≤b2. Let g
be a continuously diﬀerentiable, positive, decreasing, convex function from
(0, ∞) to (0, ∞) (g is otherwise arbitrary for the purposes of this section). A
(piecewise-linear, oriented) path γ from a to b is a subset of R2
+ consisting of
the union of a ﬁnite sequence of points a = u0 ≤u1 ≤· · · ≤uk = b (called
the vertices of the path) together with the line segments {αui + (1 −α)ui+1 :
α ∈(0, 1)} joining the vertices in order. A path γ may be parameterized as
γ = {(x(t), y(t)) : t ∈[a, b]} where x(t) + y(t) = t. We deﬁne the functional
w(γ) as the path integral
w(γ) =
Z
γ

g(y)dx + g(x)dy

=
Z b
a

g(y)dx
dt + g(x)dy
dt

dt,
and for a ≤b we deﬁne
W(a, b) = inf
γ:a→b w(γ)
where infγ:a→b denotes the inﬁmum over all paths γ from a to b.
The purpose of this self-contained section is to prove the following four
properties of W. Propositions 12, 13 and 14 are fairly natural; Proposition
15 is tailored to a speciﬁc application in Section 10.
17

<!-- page 18 -->

Proposition 12 If a ≤b ≤c then W(a, b) + W(b, c) ≥W(a, c).
Proposition 13 W(a, b) ≤(b1 −a1)g(a2) + (b2 −a2)g(a1).
Proposition 14 If a1 + a2 = A and b = (B, B) where A ≤B then
W(a, b) ≥2
Z B
A
g(z) dz.
Proposition 15 Suppose that a ≤b; c ≤d; r ≥b; r ≥d; r ≤b+d+(q, q)
and r ≥(2Z, 2Z) where q < Z. Then there exists s satisfying s ≤r and
s ≤a + c such that
W(a, b) + W(c, d) ≥W(s, r) −2qg(Z).
Indeed we may take s = a ∨[(a + c) ∧(a + c + r −b −d)].
In the above, ∨and ∧denote coordinate-wise maximum and minimum re-
spectively.
Proof of Proposition 12.
We have
W(a, c)
=
inf
γ:a→c w(γ)
≤
inf
γ:a→b→c w(γ)
=
inf
γ1:a→b w(γ1) +
inf
γ2:b→c w(γ2)
=
W(a, b) + W(b, c),
where infγ:a→b→c denotes the inﬁmum of all paths from a to c containing b.
✷
Proof of Proposition 13.
For any path γ from a to c, since g is
decreasing and u ≥a for all u ∈γ, we have
w(γ)
=
Z
γ

g(y)dx + g(x)dy

≤
Z
γ

g(a2)dx + g(a1)dy

=
(b1 −a1)g(a2) + (b2 −a2)g(a1).
✷
18

<!-- page 19 -->

b
a
∆
V
U
Figure 3: The regions U and V .
In order to prove Proposition 14 we need Lemma 16 below.
For sets
A, B ⊆R2
+ we say A lies Northwest of B and write A ⪰B if for any a ∈A
and b ∈B satisfying a1 + a2 = b1 + b2, we have a2 ≥b2. Let ∆be the
“main diagonal” of R2
+:
∆:= {u ∈R2
+ : u1 = u2}.
Lemma 16 If γ1, γ2 are paths from a to b, and either γ1 ⪰γ2 ⪰∆or
∆⪰γ2 ⪰γ1, then w(γ1) ≥w(γ2).
Proof.
Without loss of generality (since the deﬁnition of w is symmetric
in the two coordinates), we may assume γ1 ⪰γ2 ⪰∆.
First note that for any path γ from a to b, we may express w(γ) as
w(γ) =
ZZ
U
g′(y) dx dy + (b1 −a1)g(a2) +
ZZ
V
g′(x) dx dy + (b2 −a2)g(a1)
(17)
where the regions U, V are deﬁned by
U
=
{u : a ≤u ≤b and γ ⪰{u}},
V
=
{u : a ≤u ≤b and {u} ⪰γ};
see Figure 3. To check (17) note that by performing the y integral we have
ZZ
U
g′(y) dx dy =
Z
γ

g(y) −g(a2)

dx =
Z
γ
g(y) dx −(b1 −a1)g(a2),
and similarly for the second integral.
19

<!-- page 20 -->

a
u
∆
b
γ
γ0
Figure 4: A path γ from a to b, and the optimizing path γ0.
Applying (17) to γ1 and γ2 and subtracting, we have
w(γ1) −w(γ2) =
ZZ
H

g′(y) −g′(x)

dx dy,
where H is the region between γ1 and γ2:
H = {u : a ≤u ≤b and γ1 ⪰{u} ⪰γ2}.
Since γ1 ⪰γ2 ⪰∆, we have H ⪰∆, and hence y ≥x for (x, y) ∈H, and
since g is convex this implies that g′(y)−g′(x) ≥0 on H, so w(γ1)−w(γ2) ≥0.
✷
Proof of Proposition 14.
Let γ be a path from a to b. We claim that
w(γ) ≥w(γ0)
where γ0 is the path with vertices a, u, b, where u1 = u2 = max{a1, a2}
(see Figure 4). Thus W(a, b) = w(γ0). To check the above claim, split γ
into sections separated by the intersections of γ and γ0, and observe that
by Lemma 16, each section of γ has a value of w at least as large as the
corresponding section of γ0.
Now since g is positive, w(γ0) ≥w(γ′
0) where γ′
0 is the straight path with
vertices u and b. From the deﬁnition of w we have
w(γ′
0) = 2
Z B
u1
g(z) dz ≥
Z B
A
g(z) dz.
✷
20

<!-- page 21 -->

To prove Proposition 15 we need Lemmas 18, 19 below.
Lemma 18 If a ≤b and k ∈[0, ∞)2 then W(a, b) ≥W(a + k, b + k).
(Note that here k might not be an element of R2
+).
Proof.
We have
W(a + k, b + k)
=
inf
γ:a+k→b+k w(γ)
=
inf
γ:a→b w(γ + k)
≤
inf
γ:a→b w(γ)
since g is decreasing
=
W(a, b),
where γ + k denotes the shifted path obtained by adding k to each point in
γ.
✷
Lemma 19 If a ≤b ≤c then W(a, c) ≥W(b, c).
We prove Lemma 19 via the following.
Lemma 20 If a ≤c and either b = (a1, c2) or b = (c1, a2) then W(a, c) ≥
W(b, c).
Proof.
Without loss of generality suppose that b = (a1, c2). Let γ be any
path from a to b, and let δ be the unique (straight, horizontal) path from b
to c. Since g is decreasing and positive we have
W(b, c) = w(δ) =
Z
δ
g(y)dx ≤
Z
γ
g(y)dx ≤w(γ).
✷
Proof of Lemma 19.
Let γ be any path from a to b. Without loss of
generality suppose that {b} ⪰γ. Let u = (a1, b2) and v = (v1, b2) where
v1 = inf{t : (t, b2) ∈γ}.
21

<!-- page 22 -->

a
c
b
v
u
γ1
γ2
δ1
δ2
Figure 5: Illustration of the proof of Lemma 19.
Write γ = γ1 ∪γ2 where γ1 is a path a to v and γ2 is a path v to b (both
of which are thus uniquely deﬁned). Let δ = δ1 ∪δ2 where δ1 is the unique
path from u to b and δ2 is the unique path from b to v. See Figure 5.
Now
w(γ)
=
w(γ1) + w(γ2)
≥
w(δ) + w(γ2)
by Lemma 20
≥
w(δ2) + w(γ2)
=
w(δ2 ∪γ2)
≥
W(b, c)
since δ2 ∪γ2 is a path from b to c.
✷
We are now ready to prove Proposition 15; we start with a version without
the “error term”.
Lemma 21 Suppose that a ≤b; c ≤d; r ≥b; r ≥d and r ≤b + d. Then
there exists s satisfying s ≤r and s ≤a + c such that
W(a, b) + W(c, d) ≥W(s, r).
Indeed we may take s = a ∨(a −b + c −d + r).
Proof.
Deﬁne s as indicated, and note that s ≥a; s ≤a + c; and s ≤r.
Then we have
W(a, b) + W(c, d)
≥
W(s, s −a + b) + W(r + c −d, r)
by Lemma 18 twice
≥
W(s, s −a + b) + W(s −a + b, r)
by Lemma 19
≥
W(s, r)
by Lemma 12.
22

<!-- page 23 -->

R
R4
R5
R6
R8
R7
R′
R2
R1
R3
Figure 6: The rectangles R1, . . . , R8.
In the ﬁrst inequality we have used the facts that s −a ≥(0, 0) and r −d ≥
(0, 0). For the second (and third) we must check that r+c−d ≤s−a+b ≤r,
which is achieved as follows:
r −(s −a + b) = (r + a −b) −s = (r −b) ∧(d −c) ≥(0, 0),
and
(s−a+b)−(r+c−d) = s−(a−b+c−d+r) = (0, 0)∨(b−c+d−r) ≥(0, 0).
✷
Proof of Proposition 15.
Let r′ = r ∧(b + d), and note that r′ ≥b;
r′ ≥d; r′ ≤b + d; r′ ≤r and r′ ≥r −(q, q) ≥(2Z −q, 2Z −q) ≥(Z, Z).
Hence we may apply Lemma 21 to obtain
W(a, b) + W(c, d)
≥
W(s, r′)
≥
W(s, r) −W(r′, r)
by Proposition 12
≥
W(s, r) −2qg(Z)
by Proposition 13,
where s = a ∨(a −b + c −d + r′) = a ∨[(a + c) ∧(a + c + r −b −d)], and
s ≤r′ ≤r and s ≤a + c.
✷
7
Border Events
Let R, R′ be two rectangles satisfying R ⊆R′. Deﬁne rectangles R1, . . . , R8
(some of which may be empty) according to Figure 6, so that R′ is the disjoint
union of R, R1, . . . , R8. Deﬁne D(R, R′) to be the event that each of the two
23

<!-- page 24 -->

rectangles R1 ∪R8 ∪R7 and R3 ∪R4 ∪R5 is horizontally traversable, and each
of the two rectangles R1 ∪R2 ∪R3 and R7 ∪R6 ∪R5 is vertically traversable.
One may think of D(R, R′) as a necessary condition for the “border” between
R and R′ to be “traversable” from R to R′. More precisely, the event has
the following properties.
(i) If R′ is internally spanned then D(R, R′) occurs.
(ii) D(R, R′) is deﬁned in terms of the states of sites in R′ \ R.
Property (i) holds because if R′ is internally spanned then it must be hor-
izontally and vertically traversable. The purpose of this section is to prove
the following upper bound on the probability of D(R, R′).
Proposition 22 For any Z > 0 and c ∈(0, 1/2), there exist Q = Q(c, Z) <
∞and T = T(c, Z) ∈(0, Z/2) such that for any rectangles R ⊆R′ with
dimensions (m, n) and (m + s, n + t) respectively, and any q > 0, provided
m, n ≥Z/q and s, t ≤T/q we have
Pp(D(R, R′)) ≤Q exp

−(1 −2c) [g(nq)s + g(mq)t]

.
In the applications of Proposition 22, it will be essential that Q, T do not
depend on q.
Proof.
Let H be the “corner region”, H = R1 ∪R3 ∪R5 ∪R7, and let Y be
the number of occupied sites in H, Y = |X ∩H|. The idea of the proof is as
follows: if the four events in the deﬁnition of D were independent, the proof
would be easy. If Y is small, then the events are nearly independent, and if
s, t are suﬃciently small compared with m, n, the probability Y is large can
be made smaller than the bound we are trying to obtain for Pp(D).
Without loss of generality suppose that s ≤t. We have
P(D)
=
P(D | Y ≤cs)P(Y ≤cs) + P(D | cs < Y ≤ct)P(cs < Y ≤ct)
+P(D | Y > ct)P(Y > ct)
≤
P(D | Y ≤cs) + P(D | cs < Y ≤ct)P(Y > cs) + P(Y > ct)
(23)
24

<!-- page 25 -->

We claim that the terms appearing on the right side of (23) may be bounded
as follows
P(D | Y ≤cs)
≤
e4g(Z) exp −(1 −2c)[g(nq)s + g(mq)t]
(24)
P(D | cs < Y ≤ct)
≤
e2g(Z) exp −(1 −2c)g(mq)t
(25)
P(Y > cs)
≤
exp −cs(log c −log T −1)
(26)
P(Y > ct)
≤
exp −ct(log c −log T −1)
(27)
for a suitable choice of T.
The inequality (26) is an instance of the Chernoﬀbound, as follows. Y
is a binomial random variable with parameters st and p, so
P(Y > cs)
=
P(ea(Y −cs) > 1)
for any a > 0
≤
E(ea(Y −cs))
=
exp[−csa + st log(1 −p + pea)]
≤
exp[−csa + stpea]
=
exp[−cs log(c/(tp)) + cs]
taking ea = c/(tp)
≤
exp −cs(log c −log T −1),
provided tp ≤tq ≤T. (T will be a function of Z and c to be chosen later).
The bound (27) follows similarly.
To prove (25), we condition further on which sites in H are occupied. If
for each choice of this set of sites (satisfying cs < Y ≤ct) the conditional
probability of D is bounded above by the right side of (25), the desired bound
will follow. For simplicity, suppose ﬁrst that Y = ct, and also that the ct
occupied sites in H all lie in diﬀerent horizontal rows. These sites (together
with R) split R2 ∪R6 into ct + 2 “horizontal strips”, some of which may be
empty. (More precisely, if we remove from R2 ∪R6 every row which contains
an occupied site in H, we are left with ct + 2 rectangles of width m, some of
which may be empty). In order for D to occur, a necessary condition is that
each of these strips is vertically traversable. By Lemma 8 (i), the probability
of this is at most
exp −[t −ct −(ct + 2)]g(mq)
(since the sum of the vertical heights of the strips is t −ct, hence the sum of
their heights minus one is t −ct −(ct + 2)). Provided mq ≥Z, the above
expression is at most
e2g(Z) exp −(1 −2c)g(mq)t,
(28)
25

<!-- page 26 -->

since g is decreasing. Now, if Y < ct, or if some of the occupied sites in H
lie in the same horizontal rows, clearly the conditional probability of D will
be even smaller, hence we have proved (25).
The bound (24) is proved similarly. If Y ≤cs then Y ≤ct also, and by
conditioning on the occupied sites in H we obtain a collection of horizontal
strips in R2 ∪R6 together with a collection of vertical strips in R4 ∪R8.
For D to occur, each of horizontal strips must be vertically traversable, and
each of vertical strips must be horizontally traversable, but these two events
are independent (conditional on the set of occupied sites in H), so the two
bounds corresponding to (28) are multiplied, to obtain the right side of (24).
Now choose T > 0 suﬃciently small that
c(log c −log T −1) ≥2(1 −2c)g(Z)
(and also T < Z/2). Since g is decreasing this ensures that
exp −cs(log c −log T −1) ≤exp −(1 −2c)g(nq)s
(29)
and also
exp −ct(log c −log T −1) ≤exp −(1 −2c)[g(nq)s + g(mq)t].
(30)
In the latter we have used the fact that s ≤t so 2g(Z)t ≥g(Z)s + g(Z)t.
Now substituting (24)–(27) into (23) and using (29),(30) we obtain
P(D) ≤3e4g(Z) exp −(1 −2c)[g(nq)s + g(mq)t].
✷
8
Disjoint Spanning
For a collection of increasing events A1, . . . , Ak on {0, 1}Z2, the event A1 ◦
· · · ◦Ak that A1, . . . , Ak occur disjointly is deﬁned as the event that there
exist pairwise disjoint full sets of sites K1, . . . , Kk such that for each i, Ai
occurs whenever Ki is full. The BK inequality states that if A1, . . . , Ak are
deﬁned in terms the states of a ﬁnite set of sites then
Pp(A1 ◦· · · ◦Ak) ≤Pp(A1) · · ·Pp(Ak).
For more details see for example [11] p37.
26

<!-- page 27 -->

Figure 7: In this example the rectangles R′, R′′ must overlap.
Proposition 31 Let R be a rectangle with |R| ≥2.
If R is internally
spanned then there exist two distinct non-empty rectangles R′, R′′ such that
(i) the strict inclusions R′ ⊂R, R′′ ⊂R hold,
(ii) ⟨R′ ∪R′′⟩= R,
(iii) {R′ is internally spanned} ◦{R′′ is internally spanned} occurs.
To see the subtlety of Proposition 31, note that we cannot in general take the
two rectangles R′, R′′ to be disjoint; in Figure 7 for example, the whole square
R is internally spanned, but the only possible choice for R′, R′′ is the pair
of 6 by 6 squares indicated. It is for this reason that the concept of disjoint
occurrence is important. The idea of the proof of Proposition 31 is simple:
we run an algorithm which produces successively larger internally spanned
rectangles, combining two rectangles into one at each step and eventually
obtaining R; then we consider the last step. The details of the proof require
a little more care.
Proof of Proposition 31.
If K is any ﬁnite set of sites, we may con-
struct ⟨K⟩via the following algorithm. For each time step t = 0, 1, . . . , τ, we
shall construct a collection of mt rectangles Rt
1, . . . , Rt
mt, and corresponding
sets of sites Kt
1, . . . , Kt
mt, with the following properties:
(i) Kt
1, . . . , Kt
mt are pairwise disjoint;
(ii) Kt
i ⊆K;
(iii) Rt
i = ⟨Kt
i⟩;
(iv) if i̸ = j then Rt
i̸ ⊆Rt
j;
27

<!-- page 28 -->

(v) K ⊆Ht ⊆⟨K⟩, where
Ht :=
mt
[
i=1
Rt
i.
Initially, the rectangles and sets of sites are just the individual sites of
K. That is, let K be enumerated as K = {x1, . . . , xk}, and set m0 = k and
R0
i = K0
i = {xi}, so that in particular
H0 = K.
The ﬁnal set of rectangles Rτ
1, . . . Rτ
mτ will have the property that
Hτ = ⟨K⟩.
(32)
Before describing the algorithm we make the following observation. Let
R1, R2 be two distinct rectangles neither of which is a subset of the other,
and consider ⟨R1 ∪R2⟩. The following three possibilities exist.
(a) ⟨R1 ∪R2⟩= R1 ∪R2, and R1 ∪R2 is not connected.
(b) ⟨R1 ∪R2⟩= R1 ∪R2, and R1 ∪R2 is a rectangle.
(c) ⟨R1 ∪R2⟩⊃R1 ∪R2, and R1 ∪R2 is a rectangle.
(As usual, a set of sites is said to be connected if any two sites can be joined
via a sequence of sites at Euclidean distances 1).
The algorithm proceeds as follows. Suppose Rt
1, . . . , Rt
mt and Kt
1, . . . , Kt
mt
have already been constructed.
Step (I). If there do not exist a pair of rectangles Rt
i, Rt
j with i̸ = j such
that ⟨Rt
i ∪Rt
j⟩is a rectangle (that is, if case (a) above holds for all
pairs), then STOP, and set τ = t.
Step (II). If there do exist Rt
i, Rt
j with i̸ = j such that ⟨Rt
i∪Rt
j⟩is a rect-
angle (case (b) or (c) above), then choose one such pair of rectangles.
Denote the rectangle ⟨Rt
i ∪Rt
j⟩by R′. Also let K′ = Kt
i ∪Kt
j.
Step (III). Construct the state (Rt+1
1
, Kt+1
1
), . . . , (Rt+1
mt+1, Kt+1
mt+1) at time
t+1 as follows. From the list (Rt
1, Kt
1), . . . , (Rt
mt, Kt
mt) at time t, delete
every pair (Rt
l, Kt
l ) for which Rt
l ⊆R′.
This includes (Rt
i, Kt
i) and
(Rt
j, Kt
j), and may include others. Then add (R′, K′) to the list.
28

<!-- page 29 -->

Step (IV). Increase t by 1 and return to Step (I).
It is straightforward to see that properties (i)–(v) are preserved by this
procedure. Also mt is strictly decreasing with t, so the algorithm must stop
eventually. To check that (32) is satisﬁed, observe that if there exists a site
x ∈⟨K⟩\ Hτ, then there must exist a site y ∈⟨K⟩\ Hτ having at least
two neighbours in Hτ, but these neighbours must lie in distinct rectangles
Rτ
i , Rτ
j , and hence the algorithm should not have stopped at time τ.
Furthermore, observe that if ⟨K⟩is a single rectangle R, then we must
have mτ = 1 and Rτ
1 = R. If not, since Hτ = R, there must exist two distinct
rectangles Rτ
i , Rτ
j whose union is connected, and again this means that the
algorithm should not have stopped.
Finally, to prove the proposition, note that if R is internally spanned then
running the algorithm on the set of sites K = R ∩X results in mτ = 1 and
Rτ
1 = R. If |R| ≥2 then there must have been at least one step, τ ≥1. Now
considering the last time step of the algorithm (from time τ −1 to time τ)
we obtain the two rectangles R′ = Rτ−1
i
, R′′ = Rτ−1
j
with all the required
properties.
✷
9
Hierarchies
A directed graph is a set of vertices V together with a set of ordered pairs of
vertices E. If (u, v) ∈E then we say there is an edge from u to v and write
u ❀v. The children of a vertex u are all the vertices v such that u ❀v.
A hierarchy H is a ﬁnite directed graph in which every vertex v is labeled
with a non-empty rectangle Rv (where the rectangles corresponding to diﬀer-
ent vertices are not necessarily distinct), with the following properties. The
graph is a tree. There is a special vertex r called the root, and all edges are
directed away from the root. (So for any vertex v there is a unique directed
chain of vertices r = u0 ❀u1 ❀· · · ❀uk = v). If u ❀v then we have the
strict inclusion Ru ⊃Rv. Every vertex has 0, 1 or 2 children. A vertex with
no children is called a seed. If u has exactly one child v, we call u normal
and write u ⇒v. If u has two children v, w we call u a splitter and write
u ⇒(v, w). If u ⇒(v, w) then we have ⟨Rv ∪Rw⟩= Ru. See Figure 8 for an
example of a hierarchy.
29

<!-- page 30 -->

R4
R5
R3
R6
R1
R2
Rr
r
1
2
3
4
5
6
Figure 8: An example of a hierarchy.
We say that the hierarchy H occurs if all the following events occur dis-
jointly:
D(Rv, Ru)
for each pair u, v such that u ⇒v,
and
{Rw is internally spanned}
for each seed w.
(Note in particular that if for example Rv = Rw for two distinct seeds v, w
then we require that the event {Rv is internally spanned} occurs twice dis-
jointly; that is, there are two disjoint sets of occupied sites each of which
internally spans Rv).
Fix q, T, Z satisfying 0 < 4q ≤2T ≤Z ≤1/2 (in fact we shall be
concerned with the case 0 < q ≪T ≪Z ≪1). We call a hierarchy H good
if it satisﬁes all the following.
(i) If w is a seed then
short(Rw) < 2Z/q;
(ii) if u is normal or a splitter then
short(Ru) ≥2Z/q;
(iii) if u ⇒v and v is a seed or a normal vertex then
φ(Ru) −φ(Rv) ∈[T/(2q), T/q];
30

<!-- page 31 -->

(iv) if u ⇒v and v is a splitter then
φ(Ru) −φ(Rv) ≤T/q;
(v) if u ⇒(v, w) then
φ(Ru) −φ(Rv)
≥
T/(2q)
and
φ(Ru) −φ(Rw)
≥
T/(2q).
Proposition 33 Let q, T, Z satisfy 0 < 4q ≤2T ≤Z ≤1/2, and let R be
a non-empty rectangle. If R is internally spanned then some good hierarchy
with root-label Rr = R occurs.
Proof.
The proof is by induction on the size of the rectangle. Let R be
a rectangle and suppose the proposition holds for all rectangles with semi-
perimeter less than φ(R).
If short(R) < 2Z/q, then the good hierarchy having only one vertex r
(which is the root and a seed), and Rr = R, occurs.
If short(R) ≥2Z/q then we construct a sequence of non-empty internally
spanned rectangles R = S0 ⊃S1 ⊃· · · ⊃Sm by an algorithm as follows. The
idea is that S1, . . . , Sm are successive attempts to ﬁnd a rectangle S such that
φ(R) −φ(S) is in the range [T/(2q), T/q]; the attempt may succeed in which
case the root will be a normal vertex, or we may “overshoot” in which case we
need to introduce a splitter. Here are the details. Given Si, apply Proposition
31 to obtain two rectangles S′
i, S′′
i . Let Si+1 be the one of S′
i, S′′
i with the
larger semi-perimeter φ (choosing according to an arbitrary rule if they are
equal). Stop, after m ≥1 steps, the ﬁrst time φ(R) −φ(Sm) ≥T/(2q).
(This must occur eventually because the sequence of rectangles is strictly
decreasing, and φ(R) ≥4Z/q > T/(2q) + 2, so that a rectangle S containing
only one site does satisfy φ(R) −φ(S) ≥T/(2q); and we must have m ≥1
because T/(2q) > 0).
Now consider the following three possible cases.
(i) If φ(R) −φ(Sm) ≤T/q,
then we have that φ(R) −φ(Sm) ∈[T/(2q), T/q], Sm ⊂R, and R and
Sm are internally spanned. By the inductive hypothesis, there exists a
good hierarchy H′ with root r′ and root-label Rr′ = Sm. Furthermore,
31

<!-- page 32 -->

the events D(Sm, R) and {H′ occurs} occur disjointly, since they are
deﬁned in terms of disjoint sets of sites. We construct a hierarchy H as
follows: start with H′, and add a new vertex r, with Rr = R, and an
edge from r to r′ so that r ⇒r′. It follows from the above observations
that H is good and occurs.
(ii) If φ(R) −φ(Rm) > T/q and m = 1,
then there exist disjointly internally spanned rectangles S′
0, S′′
0 such
that ⟨S′
0 ∪S′′
0⟩= R, and φ(R) −φ(S′
0), φ(R) −φ(S′′
0) > T/q > T/(2q).
By Proposition 31 and the inductive hypothesis there exist disjointly
occuring good hierarchies H′, H′′ with roots r′, r′′ and root-labels Rr′ =
S′
0, Rr′′ = S′′
0, and we may take the vertex sets of H′, H′′ to be disjoint.
Now we construct H from H′ ∪H′′ (this last object being a labeled
directed graph deﬁned in the obvious way) by adding a new vertex r,
with Rr = R, and edges so that r ⇒(r′, r′′). It is easily seen that H is
a good hierarchy and occurs.
(iii) If φ(R) −φ(Rm) > T/q and m ≥2,
then we have internally spanned rectangles Sm−1, S′
m−1, S′′
m−1 satisfying
R ⊃Sm−1 and Sm−1 = ⟨S′
m−1 ∪S′′
m−1⟩, where φ(R) −φ(Sm−1) <
T/(2q) and φ(R) −φ(S′
m−1), φ(R) −φ(S′′
m−1) > T/q, and therefore we
have φ(Sm−1)−φ(S′
m−1), φ(Sm−1)−φ(S′′
m−1) ≥T/(2q). By Proposition
31 and the inductive hypothesis there exist disjointly occuring good
hierarchies H′, H′′ with roots r′, r′′ and root-labels Rr′ = S′
m−1, Rr′′ =
S′′
m−1, and we may take the vertex sets of H′, H′′ to be disjoint. Now
we construct H from H′ ∪H′′ by adding new vertices r, y with Rr =
R, Ry = Sm−1, and new edges so that r ⇒y and y ⇒(r′, r′′). Then it
is easily seen that H is a good hierarchy and occurs.
✷
10
Upper Bound
We are now ready to prove Theorem 2 (ii). Fix B > 2, and let A = c = 1/B.
We shall prove that
lim inf
p→0
−p log I(⌊B/p⌋, p) ≥2(1 −2c)
Z B
A
g(z) dz,
32

<!-- page 33 -->

from which Theorem 2 (ii) follows. The approach is to use Proposition 33,
and obtain upper bounds on the number of possible good hierarchies, and on
the probability that each one occurs.
Choose Z > 0 suﬃciently small that Z < A/2 and
g(2Z) ≥4λ
A
(34)
where λ = π2/18. Recall that g(z) →∞as z →0, so this is indeed possible;
the reason for this particular choice of Z will become clear later. Finally
choose T = T(c, Z), Q = Q(c, Z) according to Proposition 22. It will also be
convenient to assume that q < T/2. Thus we have
16q < 8T < 4Z < 2A < 1 < B/2.
We shall be concerned with “large” B and “very small” q, in which case we
have
q ≪T ≪Z ≪A ≪1 ≪B.
Later we will let q →0 while keeping B ﬁxed.
It will be important to
distinguish between quantities which depend only on B (such as A, c, Z, T, Q)
and those which also depend on q.
If R is a rectangle we deﬁne
V (R) = q long(R) g

q short(R)

,
and for rectangles R ⊆R′ we deﬁne
U(R, R′) = W

q dim(R), q dim(R′)

,
where W is deﬁned as in Section 6. We claim that for any R we have
Pp(R is internally spanned) ≤exp −q−1V (R),
(35)
while for any R ⊆R′ such that short(R′) ≥2Z/q and φ(R′) −φ(R) ≤T/q
we have
Pp(D(R, R′)) ≤Q exp −q−1(1 −2c)U(R, R′).
(36)
Inequality (35) follows from Lemma 8 (ii) and Lemma 9 (i) (using which-
ever of East- and North- traversability gives the better bound). For (36),
note that dim(R′) −dim(R) ≤(T/q, T/q), and short(R) ≥2Z/q −T/q ≥
33

<!-- page 34 -->

Z/q, so R, R′ satisfy the conditions of Proposition 22. Combining this with
Proposition 13 and the deﬁnition of U above yields (36).
Now if H is any good hierarchy, by BK inequality and the deﬁnition of
H occuring, we have
Pp(H occurs)
≤
Y
u⇒v
Pp

D(Rv, Ru)
 Y
w seed
Pp(Rw internally spanned)
≤
QNnorm exp −q−1
"
(1 −2c)
X
u⇒v
U(Rv, Ru) +
X
w seed
V (Rw)
#
,
(37)
where the ﬁrst product and sum are over all pairs of vertices u, v such that
u ⇒v, the second product and sum are over all seeds w, and Nnorm is the
number of normal vertices of H. The second inequality in (37) follows from
(35),(36) above and from properties (ii),(iii),(iv) of a good hierarchy.
Next we derive lower bounds for the two sums in (37).
Lemma 38 For any good hierarchy H with root-label Rr = R, there exists a
rectangle S = S(H) ⊆R satisfying
dim(S) ≤
X
w seed
dim(Rw)
such that
X
u⇒v
U(Rv, Ru) ≥U(S, R) −Nsplit2qg(Z),
where Nsplit is the number of splitters of H.
We call S as above the pod of H. The idea is that, if the total size of the
seeds is not too big, then the pod is not too big, and therefore the ﬁrst sum
in (37) is large enough to give a good bound. Note that the location of the
pod is actually immaterial - only its dimensions are ever used. It is deﬁned
to be a rectangle rather than just a 2-vector as a notational convenience only.
The pod is function of the hierarchy H only, not of B or q (even though the
deﬁnitions of a good hierarchy and of U(·, ·) do depend on B and q).
The following will be needed in the proof of Lemma 38.
Lemma 39 If ⟨R′ ∪R′′⟩= R then dim(R′) + dim(R′′) ≥dim(R) −(1, 1).
34

<!-- page 35 -->

Proof.
This is proved in a similar manner to Lemma 9 (i). If the sum
of the widths of R′, R′′ is less than the width of R minus one, then either R
has two adjacent columns which do not intersect R′ ∪R′′, or the East-most
or West-most column of R does not intersect R′ ∪R′′. In either case, no site
in such a column can be in ⟨R′ ∪R′′⟩.
✷
Proof of lemma 38.
The proof is by induction on the number of vertices
of H. Suppose the lemma holds for all hierarchies with fewer vertices than
H. We consider three cases according to whether the root r is a seed, normal,
or a splitter.
If r is a seed (so it is the only vertex), then we take S = R, and the result
holds trivially.
If r is normal, so that r ⇒y say, then we apply the inductive hypothesis
to the sub-hierarchy H′ rooted at y (that is, the hierarchy obtained by taking
all vertices and edges in directed chains y = v0 ❀v1 ❀v2 ❀· · · away from
y, together with the associated rectangles). Let S = S(H) = S(H′), and
note that H′ has the same number of splitters as H, to obtain
X
u⇒v
U(Rv, Ru)
≥
U(Ry, R) + U(S, Ry) −Nsplit2qg(Z)
≥
U(S, R) −Nsplit2qg(Z),
by Proposition 12 and the deﬁnition of U(·, ·).
If r is a splitter, so that r ⇒(y1, y2) say, we apply the inductive hypothesis
to the sub-hierarchies H1, H2 rooted at y1, y2, and denote their pods S1 =
S(H1), S2 = S(H2). We also write R1 = Ry1, R2 = Ry2. Since the total
number of splitters in H1 and H2 is one less than Nsplit, we obtain
X
u⇒v
U(Rv, Ru) ≥U(S1, R1) + U(S2, R2) −(Nsplit −1)2qg(Z).
(40)
Now we apply Proposition 15 with
a
=
q dim(S1),
b
=
q dim(R1),
c
=
q dim(S2),
d
=
q dim(R2),
r
=
q dim(R).
35

<!-- page 36 -->

We choose the pod S = S(H) to be a rectangle satisfying S ⊆R and
s = q dim(S),
where s is as in Proposition 15. The formula for s ensures that the dimensions
of S are indeed integers and that S depends only on H. Furthermore, since
s ≤a + c we have
dim(S) ≤dim(S1) + dim(S2) ≤
X
w seed
dim(Rw),
by the inductive hypothesis, since the set of seeds of H is the disjoint union
of the sets of seeds of H′ and H′′. It is easy to check that the conditions of
Proposition 15 are met, by Lemma 39 and property (ii) of a good hierarchy,
so we obtain
U(S1, R1) + U(S2, R2) ≥U(S, R) −2qg(Z).
Combining this with (40) gives the required bound.
✷
Now we derive a lower bound on the second sum in (37). If w is a seed
then
V (Rw)
qφ(Rw) = long(Rw) g(q short(Rw))
long(Rw) + short(Rw) ≥g(2Z)
2
,
by property (i) of a good hierarchy. Hence
X
w seed
V (Rw) ≥g(2Z)
2
X
w seed
qφ(Rw) ≥g(2Z)
2
qφ(S),
(41)
since the pod S satisﬁes dim(S) ≤P
w seed dim(Rw) by Lemma 38.
Substituting from Lemma 38 and (41) into (37), for any good hierarchy
H we have
Pp(H occurs) ≤
QNnormQ
Nsplit
1
exp −q−1

(1 −2c)U(S, R) + g(2Z)
2
qφ(S)

,
(42)
where Q1 = e2g(Z).
Now suppose that H is a good hierarchy with root-label Rr = R =
R(⌈B/q⌉, ⌈B/q⌉). Let us ﬁnd an upper bound on the number of vertices
36

<!-- page 37 -->

of H. By properties (iii),(v) of a good hierarchy, in any directed chain of
vertices r = v0 ❀v1 ❀· · · ❀vk away from the root, at least half the edges
have φ decreasing by at least T/(2q), so the number of vertices in such a
chain is at most
22⌈B/q⌉
T/(2q) + 1 ≤20B
T .
Hence, since the graph underlying H is a binary tree, the total number of
vertices in H is at most
M := 220B/T .
All that matters is that this number depends only on B, not on q.
We now divide hierarchies into two diﬀerent types according to the semi-
perimeter of the pod. If qφ(S) ≤A then by Lemma 14 we have
U(S, R) ≥2
Z q⌈B/q⌉
qφ(S)
g(z) dz ≥2
Z B
A
g(z) dz.
On the other hand, if qφ(S) > A then by the choice of Z, (34), we have that
g(2Z)
2
qφ(S) ≥g(2Z)A
2
≥2λ ≥2(1 −2c)
Z B
A
g(z) dz,
since λ =
R ∞
0 g(z) dz. In both cases we obtain from (42) that
Pp(H occurs) ≤QM
2 exp −q−12(1 −2c)
Z B
A
g(z) dz,
(43)
where Q2 = max{Q, Q1}.
We now bound the total number of possible good hierarchies with root-
label R. The number of abstract directed graphs with at most M vertices
is at most M2M2, and the number of diﬀerent rectangles in R is at most
(B/q + 1)4, so the number of possible good hierarchies is at most
M2M2(2B/q)4M.
(44)
From Proposition 33 and (43), (44) we deduce the existence of constants
K1, K2 ∈(0, ∞) depending only on B, such that when q is suﬃciently small,
Pp(R is internally spanned) ≤K1q−K2 exp −q−12λB,
37

<!-- page 38 -->

where
λB = (1 −2c)
Z B
A
g(z) dz.
(45)
Hence, recalling that q ∼p, we have
lim inf
p→0
−p log I(⌈B/q⌉, p) ≥2λB.
This implies that
lim inf
p→0
−p log I(⌊B/p⌋, p) ≥2λB.
(46)
To check this note that for any ǫ > 0 we may write ⌊B/p⌋= ⌈(B + ǫ)/q′⌉,
where q′ = −log(1 −p′) and p′ ∼p(B + ǫ)/B as p →0; since I(L, p) is
increasing in p this implies that the left side of (46) is at least 2λB+ǫB/(B+ǫ),
establishing (46). Finally, since A = c = 1/B we have from (45) that
λB →
Z ∞
0
g(z) dz = λ
as B →∞.
Hence we have proved Theorem 2 (ii).
Acknowledgments
I thank Marek Biskup, Lincoln Chayes, Nathaniel Grossman, Tom Liggett
and Roberto Schonmann for valuable conversations.
References
[1] J. Adler. Bootstrap percolation. Physica A, 171:453–470, 1991.
[2] J. Adler, D. Stauﬀer, and A. Aharony. Comparison of bootstrap perco-
lation models. Journal of Physics A, 22:L297–L301, 1989.
[3] M. Aizenman and J. L. Lebowitz.
Metastability eﬀects in bootstrap
percolation. J. Phys. A, 21(19):3801–3813, 1988.
[4] E. D. Andjel, T. S. Mountford, and R. H. Schonmann. Equivalence of
exponential decay rates for bootstrap-percolation-like cellular automata.
Ann. Inst. H. Poincar´e Probab. Statist., 31(1):13–25, 1995.
38

<!-- page 39 -->

[5] R. Cerf and E. N. M. Cirillo. Finite size scaling in three-dimensional
bootstrap percolation. Ann. Probab., 27(4):1837–1850, 1999.
[6] R. Cerf and F. Manzo. The threshold regime of ﬁnite volume bootstrap
percolation. To appear.
[7] L. R. Fontes, V. Sidoravicius, and R. H. Schonmann. Stretched expo-
nential ﬁxation in stochastic ising models at zero temperature. Commu-
nications in Mathematical Physics. To appear.
[8] K. Frob¨ose. Finite-size eﬀects in a cellular automaton for diﬀusion. J.
Statist. Phys., 55(5-6):1285–1292, 1989.
[9] I. S. Gradshteyn and I. M. Ryzhik. Table of integrals, series, and prod-
ucts. Academic Press, New York, 1965.
[10] J. Gravner and D. Griﬀeath. Scaling laws for a class of critical cellular
automaton growth rules. In Random walks (Budapest, 1998), pages 167–
186. J´anos Bolyai Math. Soc., Budapest, 1999.
[11] G. R. Grimmett. Percolation. Springer-Verlag, second edition, 1999.
[12] S. Kirkpatrick, W. Wilcke, R. Garner, and H. Huels.
Percolation in
dense storage arrays. Physica A. To appear.
[13] R. H. Schonmann. Finite size scaling behavior of a biased majority rule
cellular automaton. Phys. A, 167(3):619–627, 1990.
[14] R. H. Schonmann. On the behavior of some cellular automata related
to bootstrap percolation. Ann. Probab., 20(1):174–193, 1992.
[15] R. H. Schonmann.
Slow droplet-driven relaxation of stochastic Ising
models in the vicinity of the phase coexistence region. Comm. Math.
Phys., 161(1):1–49, 1994.
[16] R. H. Schonmann and S. B. Shlosman. Wulﬀdroplets and the metastable
relaxation of kinetic Ising models. Comm. Math. Phys., 194(2):389–462,
1998.
[17] A. C. D. van Enter, J. Adler, and J. A. M. S. Duarte. Finite-size eﬀects
for some bootstrap percolation models. J. Statist. Phys., 60(3-4):323–
332, 1990.
39

<!-- page 40 -->

[18] A. C. D. van Enter, J. Adler, and J. A. M. S. Duarte.
Addendum:
“Finite size eﬀects for some bootstrap percolation models”. J. Statist.
Phys., 62(1-2):505–506, 1991.
[19] S. Wolfram, editor. Theory and applications of cellular automata. World
Scientiﬁc Publishing Co., Singapore, 1986. Including selected papers
1983–1986.
Alexander E. Holroyd
UCLA Department of Mathematics
405 Hilgard Avenue
Los Angeles
CA 90095-1555
U. S. A.
holroyd@math.ucla.edu
40
