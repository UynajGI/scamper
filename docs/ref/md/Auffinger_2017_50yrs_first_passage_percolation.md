# 50 years of first passage percolation

> AMS Univ. Lecture Series 68 (2017) · [arXiv](https://arxiv.org/abs/1511.03262)
> Local copy: `docs/ref/Auffinger_2017_50yrs_first_passage_percolation.pdf` (160 pages)

<!-- page 1 -->

50 years of ﬁrst passage percolation
Antonio Auﬃnger ∗
Northwestern University
Michael Damron †
Georgia Tech
Jack Hanson ‡
CUNY
September 26, 2016
Abstract
We celebrate the 50th anniversary of one the most classical models in probability
theory. In this survey, we describe the main results of ﬁrst passage percolation,
paying special attention to the recent burst of advances of the past 5 years. The
purpose of these notes is twofold. In the ﬁrst chapters, we give self-contained proofs
of seminal results obtained in the ’80s and ’90s on limit shapes and geodesics, while
covering the state of the art of these questions. Second, aside from these classical
results, we discuss recent perspectives and directions including (1) the connection
between Busemann functions and geodesics, (2) the proof of sublinear variance
under 2+log moments of passage times and (3) the role of growth and competition
models. We also provide a collection of (old and new) open questions, hoping to
solve them before the 100th birthday.
Contents
1
Introduction
4
1.1
The model of ﬁrst passage percolation and its history . . . . . . . . . . .
4
1.2
Acknowledgments . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
7
2
The time constant and the limit shape
8
2.1
Subadditivity and the time constant . . . . . . . . . . . . . . . . . . . . .
8
2.2
The time constant through a homogenization problem . . . . . . . . . . .
14
2.3
The limiting ball: Cox-Durrett shape theorem . . . . . . . . . . . . . . .
15
∗The research of A. A. is supported by NSF grant DMS-1597864.
†The research of M. D. is supported by NSF grant DMS-1419230 and an NSF CAREER award.
‡The research of J. H. is supported by the AMS-Simons travel grant and NSF grant DMS-1612921.
0MSC2000: Primary 60K35, 82B43.
0Keywords: First-passage percolation, shape ﬂuctuations, Busemann functions, Richardson’s growth
model, graph of infection.
1
arXiv:1511.03262v2  [math.PR]  23 Sep 2016

<!-- page 2 -->

2.4
Other limit shapes
. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
19
2.4.1
Shell passage times . . . . . . . . . . . . . . . . . . . . . . . . . .
20
2.4.2
FPP in the super-critical percolation cluster . . . . . . . . . . . .
20
2.5
Properties of the limit shape . . . . . . . . . . . . . . . . . . . . . . . . .
21
2.5.1
Flat edges for limit shapes . . . . . . . . . . . . . . . . . . . . . .
21
2.6
The subadditive ergodic theorem revisited
. . . . . . . . . . . . . . . . .
24
2.7
Pointed Gromov-Hausdorﬀconvergence . . . . . . . . . . . . . . . . . . .
27
2.8
Strict convexity of the limit shape . . . . . . . . . . . . . . . . . . . . . .
30
2.9
Simulations
. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
32
3
Fluctuations and concentration bounds
34
3.1
Variance bounds
. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
34
3.2
Log improvement to upper bound for d ≥2 . . . . . . . . . . . . . . . . .
37
3.3
Log improvement to lower bound for d = 2 . . . . . . . . . . . . . . . . .
42
3.4
Concentration bounds
. . . . . . . . . . . . . . . . . . . . . . . . . . . .
47
3.4.1
Subdiﬀusive concentration . . . . . . . . . . . . . . . . . . . . . .
47
3.4.2
Talagrand’s theorem via the entropy method . . . . . . . . . . . .
49
3.5
Convergence of the mean for sub-additive ergodic processes . . . . . . . .
52
3.5.1
Non-random ﬂuctuations for subadditive sequences
. . . . . . . .
53
3.5.2
Asymptotic geodesicity and rate of convergence of the expected ball 57
3.5.3
Alexander’s method
. . . . . . . . . . . . . . . . . . . . . . . . .
59
3.6
Large deviations . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
65
3.7
Cases where Gaussian ﬂuctuations appear
. . . . . . . . . . . . . . . . .
68
3.7.1
Critical ﬁrst-passage percolation . . . . . . . . . . . . . . . . . . .
68
3.7.2
Passage time in small cylinders
. . . . . . . . . . . . . . . . . . .
72
4
Geodesics
73
4.1
Existence of ﬁnite geodesics, their sizes, and the geodesic tree
. . . . . .
73
4.2
The wandering exponent . . . . . . . . . . . . . . . . . . . . . . . . . . .
78
4.2.1
Upper bounds on ξ . . . . . . . . . . . . . . . . . . . . . . . . . .
79
4.2.2
Licea-Newman-Piza lower bound on ξ . . . . . . . . . . . . . . . .
83
4.3
The scaling relation χ = 2ξ −1
. . . . . . . . . . . . . . . . . . . . . . .
88
4.4
Inﬁnite geodesics
. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
90
4.4.1
Existence of geodesic rays
. . . . . . . . . . . . . . . . . . . . . .
90
4.4.2
Directions and coalescence . . . . . . . . . . . . . . . . . . . . . .
92
4.5
Absence of geodesic lines and connection to the two-dimensional Ising
ferromagnet . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
98
4.5.1
Heuristic argument . . . . . . . . . . . . . . . . . . . . . . . . . .
98
4.5.2
Rigorous results . . . . . . . . . . . . . . . . . . . . . . . . . . . .
101
2

<!-- page 3 -->

5
Busemann functions
107
5.1
Basics of Busemann functions . . . . . . . . . . . . . . . . . . . . . . . .
107
5.2
Hoﬀman’s argument for multiple geodesics . . . . . . . . . . . . . . . . .
108
5.3
Directions of geodesics via Busemann functions
. . . . . . . . . . . . . .
110
5.4
Busemann increment distributions and geodesic graphs . . . . . . . . . .
113
5.5
Busemann functions along boundaries in Z2
. . . . . . . . . . . . . . . .
119
5.6
Nonexistence of Bigeodesics in ﬁxed directions . . . . . . . . . . . . . . .
123
6
Growth and infection models, competition interface, and more
128
6.1
Eden Model and the limit shape in high dimensions . . . . . . . . . . . .
128
6.2
Growth and competition models . . . . . . . . . . . . . . . . . . . . . . .
130
6.3
Competition with same speed . . . . . . . . . . . . . . . . . . . . . . . .
131
6.4
Competition with diﬀerent growth speeds . . . . . . . . . . . . . . . . . .
134
6.5
The competition interface
. . . . . . . . . . . . . . . . . . . . . . . . . .
135
7
Related models and questions
138
7.1
The maximum ﬂow . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
138
7.2
Variants . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
140
7.2.1
Diﬀerent edge weights
. . . . . . . . . . . . . . . . . . . . . . . .
140
7.2.2
FPP on diﬀerent graphs . . . . . . . . . . . . . . . . . . . . . . .
142
7.2.3
Last Passage Percolation . . . . . . . . . . . . . . . . . . . . . . .
143
8
Summary of open questions
146
3

<!-- page 4 -->

1
Introduction
1.1
The model of ﬁrst passage percolation and its history
First passage percolation (FPP) was originally introduced by Hammersley and Welsh [98]
in 1965 as a model of ﬂuid ﬂow through a random medium.
It has been a stage of
research for probabilists since its origins but despite all eﬀorts through the past decades,
most of the predictions about its important statistics remain to be understood. Most of
the beauty of the model lies in its simple deﬁnition (as a random metric space) and the
property that several of its fascinating conjectures do not require much eﬀort to be stated.
During these 50 years, FPP brought attention of theoretical physicists, biologists, and
computer scientists and also gave birth to some of the most classical tools in mathematics,
the sub-additive ergodic theorem as one of the main examples. Here, we will focus on
the model deﬁned on the lattice Zd; some variants will be discussed in Section 7.
The model is deﬁned as follows. We place a non-negative random variable τe, called
the passage time of the edge e, at each nearest-neighbor edge in Zd. The collection (τe)
is assumed to be independent, identically distributed with common distribution F and
probability measure ν. The random variable τe is interpreted as the time or the cost
needed to traverse edge e.
A path Γ is a ﬁnite or inﬁnite sequence of edges e1, e2, . . . in Zd such that for each
n ≥1, en and en+1 share exactly one endpoint. For any ﬁnite path Γ we deﬁne the
passage time of Γ to be
T(Γ) =
X
e∈Γ
τe.
Given two points x, y ∈Rd one then sets
T(x, y) = inf
Γ T(Γ),
(1.1)
where the inﬁmum is over all ﬁnite paths Γ that contain both x′ and y′, and x′ is the
unique vertex in Zd such that x ∈x′ + [0, 1)d (similarly for y′). The random variable
T(x, y) will be called the passage time between points x and y. In the original interpre-
tation of the model, T(x, y) represents the time that a ﬂuid with source in x takes to
reach a location y.
For each t ≥0 let
B(t) = {y ∈Rd : T(0, y) ≤t}.
In the case that F(0) = 0, the pair (Zd, T(·, ·)) is a metric space and B(t) ∩Zd is
the (random) ball of radius t around the origin.
The ultimate goal of ﬁrst passage
percolation is to understand this metric as the observer moves away from Zd or as we
make the Euclidean length of the edges small. A variety of questions comes to mind
almost immediately. We write | · | for the ℓ1 norm in Rd.
4

<!-- page 5 -->

1. What is the typical distance between two points that are far from each other in the
lattice? Or in other words, what can we say about T(x, y) as |x −y| →∞? Does
it converge? What is the rate of convergence?
2. How does a ball of large radius look? Do we have a scaling limit and ﬂuctuation
theory for the set B(t)?
3. What is the geometry of geodesics (time-minimizing paths) between two distant
points? How diﬀerent they are from straight lines?
4. What role does the distribution of the passage times play in describing the metric?
In this manuscript, we will discuss progress on these and related questions.
The
purpose is twofold. First, we hope that this set of notes will serve as a quick guide for
readers who are not necessarily experts in the ﬁeld. We will try to provide not only
the main results but also the main techniques and a large collection of open problems.
Second, the ﬁeld had a burst of activity in the past 5 years and the most complete survey
is more than a decade old. We hope that these notes will ﬁll this gap. At least, we hope
to share some of the beautiful mathematical ideas and constructions that arise through
FPP and which have enchanted many throughout these years.
Let’s now go back to questions 1 to 4. The original paper of Hammersley and Welsh
[98] considered question 1 for a class of passage times in Z2.
If we write e1 for the
ﬁrst coordinate vector, they showed that T(0, ne1) grows linearly in n. Their result was
extended in the famous paper of Kingman [36,118,119]. It was also the building block
for the classical “shape theorem” of Richardson [141], improved by Cox and Durrett [53]
and Kesten [115] that gives the analogue of the law of large numbers for the random ball
B(t). It roughly says that B(t) grows linearly in t and, when normalized, it converges to
a deterministic subset Bν of Rd, called the limit shape. The set Bν is not universal and
depends on the distribution of the passage times. Section 2 is devoted to explaining the
shape theorem and certain properties of the limit shape Bν.
In Section 3 we discuss the variance and the order of ﬂuctuations of the passage
time T. In two dimensions, it is expected that under certain assumptions on Bν the
ﬂuctuations are governed by the predictions of physicists, including Kardar, Parisi,
Zhang [66, 110, 111, 122].
In higher dimensions, the picture is less clear and some of
the predictions disagree. After stating what is conjectured, we focus our attention on
presenting a simple proof of sub linear variance valid under minimal assumptions on the
passage time.
Section 4 is devoted to the study of geodesics. We brieﬂy discuss the existence of ﬁnite
geodesics between any two points, then move to the study of geodesic rays. We present
results on coalescence, directional properties and we sketch the proof of the absence
of geodesic lines (or bigeodesics) in the upper half plane. The important connection
between geodesic lines and ground states of the two-dimensional Ising ferromagnet is
also presented.
5

<!-- page 6 -->

(a) t = 0
(b) t = 70
(c) t = 135
(d) t = 250
Figure 1: Simulation of the ball B(t) at t = 0, t = 70, t = 135 and t = 250. The passage
times τe have exponential mean 1 distribution. Simulations by Si Tang.
In Section 5 we describe the modern role of Busemann functions.
We explain a
beautiful argument by Hoﬀman for the existence of 2 or more geodesic rays. We then
focus on Busemann function limits and the relation to limiting geodesic graphs. Section
6 tries to cover the vast relation between FPP, growth processes and infection models.
We focus on questions of coexistence of multiple species and the limiting interface.
Section 7 is our attempt to show the reader what this survey is not about. In the
literature, there are thousands of pages of related (equally fascinating) questions and
models, similar to or inspired by FPP. We collect a few of these directions and try to
point the right references. In particular, we brieﬂy discuss FPP on diﬀerent graphs, the
maximum ﬂow problem, the exactly solvable models for last passage percolation and the
positive temperature version of the model. Section 8 is just a recollection of the open
questions spread throughout this manuscript for easy reference.
6

<!-- page 7 -->

1.2
Acknowledgments
We thank the American Institute of Mathematics and its staﬀfor helping us organize a
workshop on this subject. A.A. thanks the hospitality of the Mathematics department
at Indiana University, and its visiting professor program, where part of this survey was
completed. He also thanks Elizabeth Housworth who, introduced him to the japanese
eraser and chalk dust cleaner. M.D. thanks the hospitality of the Mathematics Depart-
ment at Northwestern University. We thank Phil Sosoe for several fruitful discussions
on this topic, especially on concentration estimates. We thank Si Tang for making the
simulations used in Figure 1 available to us, Sven Erick Alm and Maria Deijfen for the
simulations and Figures 2, 4, and 5 and Xuan Wang for discussions on non-random ﬂuc-
tuations. We also thank Daniel Ahlberg, Wai-Kit Lam, and Phil Sosoe for spotting typos
and giving comments.
7

<!-- page 8 -->

2
The time constant and the limit shape
2.1
Subadditivity and the time constant
The ﬁrst-order of growth of the passage time T(0, ne1) is described by the following
theorem.
Theorem 2.1 (Theorem 2.18 in [113]). Assume that
E min[t1, . . . , t2d] < ∞
(2.1)
where t′
is are i.i.d copies of τe. Then there exists a constant µ(e1) ∈[0, ∞) (called the
time constant) such that
lim
n→∞
T(0, ne1)
n
= µ(e1)
a.s. and in L1.
The proof of Theorem 2.1 is a classic application of the subadditive ergodic theorem
that we now state. The version that we write here is due to Liggett [129] and suﬃces
for our purposes. Several versions, with diﬀerent hypotheses, including the one with
Kingman’s original assumptions can be found in [129].
Theorem 2.2 (Subadditive Ergodic Theorem [129]). Let (Xm,n)0≤m<n be a family of
random variables that satisﬁes:
(a) X0,n ≤X0,m + Xm,n, for all 0 < m < n.
(b) The distribution of the sequences (Xm,m+k)k≥1 and (Xm+1,m+k+1)k≥1 is the same
for all m ≥0.
(c) For each k ≥1, the sequence (Xnk,(n+1)k)n≥0 is stationary.
(d) EX0,1 < ∞and EX0,n > −cn for some ﬁnite constant c.
Then
lim
n→∞
X0,n
n
exists a.s. and in L1.
(2.2)
Furthermore, if the stationary sequence in (c) is also ergodic, then the limit in (2.2) is
constant almost surely and equal to
lim
n→∞
EX0,n
n
= inf
n
1
nEX0,n.
One can ﬁnd diﬀerent proofs of the subadditive ergodic theorem in the literature and
some of them are in standard books of probability theory, for instance [68, Section 6.6].
We will visit and discuss the original proof given by Kingman in Section 2.6. We now
show how Theorem 2.1 is an easy consequence of Theorem 2.2.
8

<!-- page 9 -->

Γ(x, y)
x
z
y
Γ(y, z)
Figure 2: The geodesics Γ(x, y) and Γ(y, z) and their concatenation. As a geodesic from
x to z does not need to pass through y, we have T(x, z) ≤T(x, y) + T(y, z).
Proof of Theorem 2.1. We
apply
the
sub-additive
ergodic
theorem
to
Xm,n
=
T(me1, ne1).
It is not diﬃcult to verify conditions (a) to (d).
To check (a) (which
is just the triangle inequality) note that a path from 0 to ne1 does not necessarily need
to go through me1 while a concatenation of paths from 0 to me1 and from me1 to ne1
gives us a path from 0 to ne1. Therefore, see Figure 2,
T(0, ne1) ≤T(0, me1) + T(me1, ne1).
Items (b), (c) and ergodicity follow directly from the fact that the environment is i.i.d.,
thus it is invariant under horizontal shifts of Zd. EX0,n > −cn holds as passage times are
non-negative. The rest of item (d) follows from assumption (2.1) as there are 2d disjoint
deterministic paths Γ1, . . . , Γ2d in Zd joining 0 to e1 and therefore
T(0, e1) ≤min{T(Γ1), . . . , T(Γ2d)}.
Order them in such a way that Γ1 is the path with the largest number of edges and
denote by B the number of edges of Γ1. Then
P(T(0, e1) > s) ≤
2d
Y
i=1
P(T(Γi) > s) ≤P(T(Γ1) > s)2d,
and
P(T(Γ1) > s) ≤BP(τe > s/B).
The second inequality comes from the fact that if T(Γ1) > s then at least one of the
edges in Γ1 must have a passage time larger than s/B. Combining the previous two
inequalities and setting Y = min[t1, . . . , t2d] we obtain
P(T(0, e1) > s) ≤B2dP(τe > s/B)2d = B2dP(Y > s/B)
(2.3)
which proves the desired result.
In fact, the following lemma comes directly from (2.3) and the observation that any
path from 0 to x must contain at least one edge incident to 0.
9

<!-- page 10 -->

Lemma 2.3. Let k ≥1. Then E min[t1, . . . , t2d]k < ∞if and only if ET(0, x)k < ∞for
all x ∈Zd.
The convergence in probability of the normalized passage time was ﬁrst proved in two
dimensions in the original paper of Hammersley and Welsh [98] under the assumption
of ﬁnite mean of the random variable τe. Still under the assumption of Eτe < ∞, this
result was strengthened to almost sure and L1 convergence by Kingman, using his sub-
additivity theorem.
The condition (2.1) is necessary to have almost sure or L1 convergence. Indeed, if
(2.1) does not hold, denoting by τ n
1 , . . . , τ n
2d the edge weights of edges incident to the
vertex 2ne1, then for any c > 0, the events
An := {min[τ n
1 , . . . , τ n
2d] > cn}
are independent and satisfy
X
n
P(An) = ∞.
Thus,
an application of Borel-Cantelli Lemma shows that with probability one
lim supn T(0, ne1)/n > c.
Without assuming (2.1), Kesten [113, Page 137], Cox-Durrett [53] and Wierman [160]
establish the existence of a constant ˜µ(e1) such that
lim
n→∞T(0, ne1)/n = ˜µ(e1) in probability.
Clearly, if (2.1) holds, then ˜µ(e1) = µ(e1).
We now gather more information on the time constant µ(e1). It is clear that µ(e1)
satisﬁes
0 ≤µ(e1) ≤Eτe,
by considering the direct path from 0 to ne1 and using the law of the large numbers.
That strict inequality does not always hold is seen by taking τe = 1 almost surely.
However, if the distribution F of the passage times has at least two points in its support,
we can prove:
Theorem 2.4 (Hammersley-Welsh [98]). If F is not a trivial distribution, we have
µ(e1) < Eτe.
Proof. Choose a < b such that 0 < F(a) ≤F(b) < 1.
Pick n > 2a/(b −a).
Let
e(1), . . . , e(n) be the n edges from 0 to ne1 and f(1), . . . , f(n + 2) the edges from 0 to
e2, e2 to e2 + ne1, and e2 + ne1 to ne1. On the positive probability event
{τe(i) ≥b, τf(j) ≤a,
∀i, j}
one has T(0, ne1) < τe(1) + . . . + τe(n). Thus, ET(0, ne1) < nEτe and
µ(e1) ≤ET(0, ne1)
n
< Eτe.
10

<!-- page 11 -->

Let’s now look at lower bounds for the time constant. If F(0) > 0; that is, if we have
edges that are cost-free to cross, one may wonder if µ is equal to 0 and the growth of the
time constant is in fact not linear. This issue is handled in the next theorem. Let pc(d)
be the critical probability for bond percolation in Zd.
Theorem 2.5 (Theorem 6.1 in [113]). For FPP on Zd, µ(e1) > 0 if and only if F(0) <
pc(d).
We now make a very simple but important remark. One can extend Theorem 2.1
with a similar proof to arbitrary directions with rational coordinates. Then, we deﬁne a
homogeneous function µ : Qd →R such that, for any x ∈Qd,
lim
n→∞
T(0, nx)
n
= µ(x)
a.s. and in L1.
(2.4)
The reader should see (2.4) as the analogue of a law of large numbers. In Section 3,
we will discuss the ﬂuctuations of T(0, nx) around nµ(x), and we will discover that, in
general, a central limit theorem with Gaussian ﬂuctuations does not hold. Before moving
to the study of ﬂuctuations, we continue to gather more information on µ(x). It is not
diﬃcult to establish the following properties of µ for x, y ∈Qd and c ∈Q:
1. µ(x + y) ≤µ(x) + µ(y).
2. µ(cx) = |c|µ(x)
3. µ is invariant under symmetries of Zd that ﬁx the origin.
4. µ is uniformly continuous and Lipschitz on bounded subsets of Qd, so it has a
unique continuous extension to Rd.
Furthermore, these properties imply the following:
Theorem 2.6. µ(e1) > 0 if and only if µ(x) > 0 for all x ∈Rd \ {0}.
Subadditivity is a nice argument to show the existence of the time constant. However
it gives no insight, nor closed expression for µ. The determination of µ(e1) as a function
of F is a fundamental and diﬃcult problem in ﬁrst passage percolation.
Question 1. Find a non-trivial explicit distribution for which we can actually de-
termine µ(e1).
Although fundamental and old, the question above is perhaps not as interesting as
one addressing geometric properties of the limit shape discussed later in this section.
Furthermore, it is not even clear what one should mean by “determine” µ(e1). In FPP,
one does not know any distribution that allows exact computations, and maybe there
11

<!-- page 12 -->

are none. In Section 7, we discuss similar solvable models where explicit computations
are possible. A more interesting related question at this point is:
Given distributions F and ˜F, how diﬀerent are their respective time constants µ(e1)
and ˜µ(e1)?
In general, we lack strong information about how the limit shape changes under small
perturbations of the edge-weight distribution. If one could derive strong results in this
direction, perhaps the establishment of various conjectures about the limit shape (e.g.,
curvature) could be made easier, or reduced to ﬁnding some special class of distributions
for which the properties are explicitly derivable. The best current results on stability,
dating back over thirty years, say simply that the time constant is a continuous function
of the edge-weight distribution.
Theorem 2.7 (Cox-Kesten [54], Kesten [113]). The time constant µ(e1) is continuous
under weak convergence of i.i.d. distribution. That is, if (Fn)n is a sequence of distri-
bution functions for the edge weight τe with Fn ⇒F, and if (µn(e1)), µ(e1) denote the
respective time constants, then
lim
n µn(e1) = µ(e1) .
While the preceding result is not strong enough to preserve curvature, it does guaran-
tee a certain semicontinuity property of the set of extreme points of B. In [58], this was
used to establish the existence of limit shapes with arbitrarily many extreme points for
some nonatomic edge weight distributions; improvements to Theorem 2.7 could be useful
for similar constructions. One existing improvement of Theorem 2.7 is the recent work of
Garet, Marchand, Procaccia and Th´eret [83], which establishes an analogous continuity
result in the case that the edge weights are allowed to assume the value τe = +∞.
On the other hand, when comparing distributions F, ˜F which obey certain stochastic
orderings, much more can be said [54, Theorem 3] (see also [98, Section 6.4] and [52]).
When ˜F(t) ≤F(t) for all t ∈R, it is possible to provide an easy answer to this question, as
we can straightfowardly couple the passage times to obtain µ(e1) ≤˜µ(e1). The question
above was considered by several authors [113, 132, 164]. A gorgeous answer came with
the work of van den Berg and Kesten [164] proving the strict inequality µ(e1) < ˜µ(e1) if
F is strictly more variable than ˜F. Here we follow an extension of the van den Berg -
Kesten comparison theorem provided by Marchand [132].
Deﬁnition 2.8. Let F and ˜F be two distributions on R. We say that F is more
variable than ˜F if for every concave increasing function φ : R →R,
Z
φ dF ≤
Z
φ d ˜F
when the two integrals exist. If in addition F̸ = ˜F then we say that F is strictly more
variable than ˜F.
12

<!-- page 13 -->

Example 2.9. If ˜F stochastically dominates F then F is more variable than ˜F.
Example 2.10. For any non-negative random variable X and any constant M the
distribution of X ∧M is more variable than the distribution of X.
Example 2.11. For 0 < p < 1 and ρ a probability measure, consider the probability
measures ˜ν = pδa + (1 −p)ρ and ν = (p/2)δa−ϵ + (p/2)δa+ϵ + (1 −p)ρ. Writing Fπ for the
probability distribution associated to a probability measure π, then Fν is strictly more
variable than F˜ν.
Theorem 2.12 (van den Berg-Kesten [164], Marchand [132]). Assume that d = 2 and
suppose that F(0) < pc. If F is strictly more variable than ˜F then µ(e1) < ˜µ(e1).
For x̸ = e1 a version of µ(x) < ˜µ(x) can also be found in [132,164]. Van den Berg-
Kesten-Marchand’s comparision theorem also comes with a nice remark that gets rid of
some na¨ıve intuition. Suppose that the distribution of the passage time has unbounded
support.
One may (wrongly) imagine that there exists a threshold M such that an
optimal path between 0 and ne1, n large, never takes a linear fraction of edges with
weights above M. However, the theorem above implies that, if one truncates the passage
time at level M, the time constant µM(e1) of the truncated model is strictly less than the
original µ(e1). Thus, it is more eﬃcient for the model to use a certain positive proportion
of edges with very large passage time than to try to always avoid them.
This remark leads also to the following question, which seems to be open.
Question 2. Suppose that the support of the distribution of τe equals R+. Let
Xn = max{τe : e is in a geodesic from 0 to ne1}.
How does Xn scale with n?
The assumption F(0) < pc in Theorem 2.12 cannot be removed because of Theorem
2.5. Now, let r be the inﬁmum of the support of F. In higher dimensions, a version of
the theorem above is known [164] if F(r) <⃗pc(d), where⃗pc(d) is the critical probability
for directed edge percolation. The only missing case at this point is:
Question 3. Extend the comparision theorem to the case d > 2, r > 0 and⃗pc(d) ≤
F(r).
Remark 2.13. Marchand’s original approach does not work in higher dimensions, as
she uses large deviations for supercritical oriented percolation available only in dimension
2 (see [132, Section 7]) at the time of her paper (see the estimate in [132, Page 1014]).
New estimates were obtained in d > 2 in [84, Proposition 2] recently. These combined
with Marchand’s arguments may provide a solution to Question 3.
13

<!-- page 14 -->

2.2
The time constant through a homogenization problem
Another way to interpret the time constant was recently explored by Krishnan [120] in
FPP and by Georgiou, Rassoul-Agha and Sepp¨al¨ainen [85] in last-passage percolation.
We brieﬂy describe it here.
The idea is to interpret the passage time as a solution of an optimal-control problem.
Deﬁne
A := {±e1, . . . , ±ed}.
We think of A as the collection of possible directions to exit a vertex. We now write
τ(v, α) to refer to the weight τe at v ∈Zd along the direction α ∈A. It is now possible
to check that
T(0, x) = inf
α∈A{T(0, x + α) + τ(x, α)}.
This suggests that one could think of the problem as a homogenization problem for
metric Hamilton-Jacobi equations in Rd. The advantage of such a perspective is to allow
us to use the work of Lions-Souganidis [130,131] on certain homogenization problems in
Rd to give a diﬀerent characterization of the time constant. To see this we need a few
deﬁnitions.
Deﬁnition 2.14. For a function f : Zd →R, let
Df(x, α) = f(x + α) −f(x)
be its discrete derivative at x in the direction α ∈A.
We write Ωfor the probability space and σv : Ω→Ωfor the shift that translates the
random variables τe by v. Deﬁne
F :=
(
f : Zd × Ω→R

Df(x + v, α)(ω) = Df(x, α)(σvω), ∀v, x ∈Zd, ∀ω ∈Ω
E[Df(x, α)] = 0 ∀x ∈Zd and α ∈A
)
.
For p ∈Rd, x ∈Zd and f ∈F deﬁne
H(f, p, x) = sup
α∈A
−Df(x, α) −(p, α)
τ(x, α)

,
(2.5)
where (·, ·) is the standard inner product in Rd. The following is the main result of [120].
Theorem 2.15. Assume that the passage times are bounded and bounded away from 0,
that is, there exist m, M such that 0 < m < τe < M almost surely. Then µ(x) solves the
following Hamilton-Jacobi equation
¯H(Dµ(x)) = 1,
µ(0) = 0
14

<!-- page 15 -->

where ¯H(p) is a convex, coercive, Lipschitz continuous function given by
¯H(p) = inf
f∈F ess sup
ω∈Ω
sup
x∈Zd H(f, p, x)(ω).
Furthermore, ¯H(p) is the dual norm of µ(x) on Rd, deﬁned by
¯H(p) = sup
µ(x)=1
(p, x).
Although the theorem above gives a diﬀerent characterization for the time constant,
it has not yet been used to get a better understanding of µ(e1) as a function of F. An
algorithm for ﬁnding a minimizer for the above variational formula is explained in part II
of the Ph.D. thesis of Krishnan [121]. It would be extremely nice to extend these ideas to
say more geometric information on the limit ball (see questions in the next subsections).
2.3
The limiting ball: Cox-Durrett shape theorem
For each unit vector x ∈Rd we deﬁne the time constant µ(x) in direction x through (2.4).
In this section, we will see how the function µ describes the ﬁrst order approximation of
the random ball B(t) as t goes to inﬁnity. The main result of the section is the world
famous shape theorem, Theorem 2.16.
Let M be the set of Borel probability measures on [0, ∞) satisfying
E min{td
1, . . . , td
2d} < ∞,
(2.6)
where ti, i = 1, . . . 2d, are independent copies of τe and with
F(0) < pc(d),
(2.7)
where pc(d) is the threshold for bond percolation in Zd. If S is a subset of Rd and r ∈R
we write rS = {rs : s ∈S}.
Theorem 2.16 (Cox and Durrett [53]). For each ν ∈M, there exists a deterministic,
convex, compact set Bν in Rd such that for each ε > 0,
P

(1 −ε)Bν ⊂B(t)
t
⊂(1 + ε)Bν for all large t

= 1 .
(2.8)
Furthermore, Bν has non-empty interior and is symmetric about the axes of Rd.
Remark 2.17. If (2.7) does not hold, edges with zero passage time percolate, cre-
ating several instantaneous ’highways’. In this case, Theorem 2.5 says that the time
constant µ(e1) = 0.
As for the limit shape, one has Bν = Rd [113, Theorem 1.10].
Precisely, under assumption (2.6), we have F(0) ≥pc(d) if and only if for every M > 0
{x : |x| ≤M} ⊂B(t)
t
eventually w.p. 1.
15

<!-- page 16 -->

Remark 2.18. For x ∈Zd, let m(x) be the minimum of τe over all edges e incident
to x. If (2.6) fails, then for any C > 0,
X
x∈Zd
P

m(x) > C|x|

= ∞.
Now note that T(0, x) ≥m(x) and the random variables {m(x) : x ∈(2Z)d} are inde-
pendent. Then we can apply the Borel-Cantelli Lemma to obtain
T(0, x)
|x|
> C for inﬁnitely many x ∈Zd, a.s.
Thus, if (2.6) fails, (2.8) also does not hold.
Remark 2.19. With Theorem 2.16, Cox and Durrett provided suﬃcient and neces-
sary conditions for the shape theorem to hold in FPP. The ﬁrst shape theorem, however,
was proven in the seminal work of Richardson [141] in 1973. We will come back to discuss
Richardson’s model in Section 6.
Figure 3: In red, simulation of the 2 dimensional ball B(t), intersected with the ﬁrst
quadrant, when t = 15000. The passage times τe are distributed according to an ex-
ponential random variable of mean 1 plus a constant C. Here, C = 0 (left), C = 0.5
(middle), and C = 4 (right). In blue, the shape of a circle of radius T(0, te1). Simulations
and ﬁgures are courtesy of S.-E. Alm and M. Deijfen.
The idea of the proof of Theorem 2.16 is to ﬁrst use subadditivity to demonstrate the
linear growth of B(t) in a ﬁxed rational direction. This implies that, with probability
one, we have the right growth rate in a countable dense set of directions simultaneously.
To obtain the full result from this, we need some bound which allows us to interpo-
late between these directions, insuring that the convergence occurs along all rays with
probability one.
16

<!-- page 17 -->

There are diﬀerent ways to implement this interpolation step.
Here we sketch a
method which ﬁrst appeared in [99] for establishing the interpolative bound just men-
tioned. This method has the advantage of applying (with modiﬁcations) to the stationary
ergodic case of FPP, and is inspired by [27]. This proof method was further used for the
Busemann shape theorem (see Lemma 5.11) in [59].
Lemma 2.20 (Diﬀerence estimate). Let ν ∈M. Then there exists a constant κ < ∞
such that, for any x ∈Zd,
P


sup
z∈Zd
z̸=x
T(x, z)
∥x −z∥1
< κ


> 0 .
(2.9)
Idea of proof. This proof follows the line of a similar estimate in [53]. If Eτ d
e < ∞, the
proof follows by noting that there exist 2d edge-disjoint paths {ri} between x and z,
each of length order ∥x −z∥1. For the event {T(x, z) > κ∥x −z∥1} to occur, each of
these paths must have T(ri) > κ∥x −z∥1. Using standard estimates for i.i.d. sums,
this probability is small for κ suﬃciently large–in fact, we can get an estimate which is
summable in z. Using the Borel-Cantelli Lemma (and adjusting the constant κ upwards
if necessary) allows us to complete the proof.
In the general case, we consider a sparse lattice CZd such that between “neighboring”
vertices of this lattice, there exist 2d disjoint paths lying inside cells of side length of order
C. In particular, paths corresponding to well-separated pairs of “neighboring” vertices
are disjoint, and their passage times are independent. The smallest passage time among
these 2d paths is clearly an upper bound for the passage time between “neighboring”
vertices; we treat this as the passage time of a “renormalized” edge of the sparse lattice.
Mimicking the argument of the preceding paragraph, then extending the bound to the
rest of Zd, completes the proof.
Proof of Theorem 2.16. We will call an x in Zd for which the event appearing in (2.9)
occurs a “good” vertex. We can immediately leverage the information in that lemma to
show
Claim 1. Let ζ ∈Zd \ {0}. For a given realization of edge-weights, denote by (nk)
the sequence of natural numbers such that nkζ is a good vertex. Then with probability
one, the sequence (nk) is inﬁnite and limk(nk+1/nk) = 1.
To see that the claim is correct, note that the ergodic theorem implies that the
sequence (nk) is inﬁnite almost surely.
Let Bm denote the event that mζ is a good
vertex. Then
k
nk
= 1
nk
nk
X
i=1
1Bi;
17

<!-- page 18 -->

the right side converges to the probability in (2.9) by the ergodic theorem. Thus,
nk+1
nk
=
 nk+1
k + 1
  k
nk
 k + 1
k

→1
almost surely. This proves the claim.
Let Ξ1 denote the event that limn T(0, nq)/n = µ(q) for all q having rational coor-
dinates; let Ξ2 denote the event that for every ζ ∈Zd, the sequence (nk) deﬁned in
Claim 1 is inﬁnite and that the ratio of successive terms tends to one. From here, the
proof of Theorem 2.16 proceeds by contradiction. Assume the Shape Theorem does not
hold. Then there exists δ > 0 and a collection of edge-weight conﬁgurations Dδ with
P(Dδ) > 0 such that, for every outcome in Dδ, there are inﬁnitely many vertices x ∈Zd
with
|T(0, x) −µ(x)| > δ∥x∥1.
(2.10)
Since P(Ξ1) = 1 = P(Ξ2), the event Dδ ∩Ξ1 ∩Ξ2 contains some outcome ω; we
claim that ω has contradictory properties. On outcome ω, there must exist a sequence
(xi) ⊂Zd satisfying the condition in (2.10). We can assume that xi/∥xi∥1 converges
to some y with ∥y∥1 = 1 by compactness of the unit sphere. Let δ′ > 0 be arbitrary;
we will ﬁx its value at the end of the proof. We ﬁrst choose some large N such that
∥xn/∥xn∥1 −y∥1 < δ′ and such that
|µ(xn) −∥xn∥1µ(y)| < δ∥xn∥1/2
for n > N. Then we have for n > N (using our assumption (2.10):
|T(0, xn) −∥xn∥1µ(y)| > δ∥xn∥1/2.
(2.11)
Next, we set up a sequence of approximating good vertices. We ﬁnd some z ∈Rd,
∥z∥1 = 1 such that ∥z −y∥1 < δ′, with the additional property that z = x/M for some
x ∈Zd and some positive integer M. This can be done because vectors with rational
coordinates are dense in the unit sphere. On ω, there must exist a sequence (nk) such
that nkMz is a good vertex and such that nk+1/nk tends to one. For any n, there exists
a value of k such that
nk+1M ≥∥xn∥1 ≥nkM;
denote this value by k(n). Finally, ﬁx K > 0 such that nk+1 < (1 + δ′)nk and
|T(0, nkMz)/(nkM) −µ(z)| < δ′
for all k > K. We now let n > N be large enough that k(n) > K.
Before completing the calculation here, it is worth considering where the contradiction
will arise. We have (essentially by assumption) that T(0, ny) −nµ(y) is of order n for
inﬁnitely many n.
Since µ is a norm, µ(y) and µ(z) are arbitrarily close and since
18

<!-- page 19 -->

inﬁnitely many of the {nz} are good vertices, T(0, ny) and T(0, nz) are arbitrarily close.
Thus |T(0, nz) −nµ(z)| is large – but this is counter to the properties assumed for ω.
To turn the above into a rigorous estimate, write k for k(n) and expand

T(0, xn)
∥xn∥1
−µ(y)
 ≤

T(0, xn) −T(0, nkMz)
∥xn∥1
 + T(0, nkMz)
nkM

1 −nkM
∥xn∥1

+

T(0, nkMz)
nkM
−µ(z)
 + |µ(z) −µ(y)|.
There are four terms on the right side of the above, which we number from left to right
and bound individually in terms of δ′.
Term 1.
Since n > N and k > K, one has nkM ≤∥xn∥1 ≤(1 + δ′)nkM, that
∥nkMy −nkMz∥1 ≤δ′nkM, and that ∥xn/∥xn∥1 −y∥1 < δ′. Therefore, ∥xn −nkMz∥1 ≤
2δ′∥xn∥1. Using the fact that nkMz is a good vertex yields
|T(0, xn) −T(0, nkMz)| ≤κ∥xn −nkMz∥1 ≤2κδ′∥xn∥1.
Term 2. The relationship between nkM and ∥xn∥1 given in the Term 1 estimates yield
an upper bound for the second factor of Term 2. By the fact that k > K, we can bound
the ﬁrst factor. The overall bound is
[µ(z) + δ′]
 1 −(1 + δ′)−1
.
Term 3. By the fact that k is chosen greater than K, this term is bounded above by δ′.
Term 4. If µ is identically zero, this term is trivially zero. If µ is not identically zero, it
is a norm on Rd and is thus bounded by the ∥· ∥1 norm:
cL∥· ∥1 ≤µ(·) ≤cU∥· ∥1.
Since ∥z −y∥1 < δ′, Term 4 is bounded above by δ′ times a constant depending only on
µ.
We have therefore bounded the left side of (2.11) by an expression of the form
f(δ′)∥xn∥1, where f tends to zero as δ′ →0. Since δ′ was arbitrary, we can choose
it such that f(δ′)∥xn∥1 is smaller than the right side of (2.11). This contradiction proves
the theorem.
2.4
Other limit shapes
In this section, we brieﬂy discuss a few extensions of Theorem 2.16.
19

<!-- page 20 -->

2.4.1
Shell passage times
As we saw in Remark 2.18, when (2.6) fails we have with probability one
lim sup
v→∞
T(0, v)
|v|
= ∞.
Nevertheless, without any moment condition, one can deﬁne a modiﬁed passage time
ˆT(u, v) such that the family of random variables ˆT(u, v) −T(u, v) is tight and one has
a limit shape for the modiﬁed ˆT. This was ﬁrst done by Cox-Durrett in dimension 2
and later extended by Kesten to all dimensions. Their construction goes as follows. Let
M ∈R+ be large enough so that F([0, M]) is very close to 1. The collection of edges e
such that τe ≤M is a super-critical percolation process, so if we denote by CM its inﬁnite
cluster, each point u ∈Zd is a.s. surrounded by a small contour (or shell) S(u) ⊂CM.
They deﬁne ˆT(u, v) = T(S(u), S(v)) for u, v ∈Zd. The times ˆT(0, u) have good moment
properties; thus their limit shape can be deﬁned using the a.s. and L1 limit given by the
previous arguments. The details of this construction in d > 2 requires certain topological
properties of the exterior boundary of a subset of Zd. These properties were derived by
Kesten and generalized in the work of Tim´ar [155].
2.4.2
FPP in the super-critical percolation cluster
Another direction where shape theorems have been proven is where we allow passage
times to be inﬁnite. This is equivalent to considering FPP on a super-critical Bernoulli
percolation performed independently. When the edge weights are either 1 or ∞, the
passage time is also known as the chemical distance.
In this setting, the benchmark is the work of Gar´et-Marchand [80], where the ana-
logues of Theorems 2.5, 2.16 were proven under a moment condtion Eτ α(d)
e
1{τe<∞} < ∞,
where α(d) = 2(d2 + d −1) + ϵ; see hypothesis Hα on page 4 in [80]. Their results are
also valid for stationary ergodic passage times, in the spirit of the work of Boivin [27].
In the i.i.d.
case, in two independent works, Cerf-Th´eret [39] and Mourrat [134]
recently removed all moment assumptions of [80], by proving a weak shape theorem. In
his paper, Mourrat considers a model of a random walk in a random potential, but he
discusses how his theorems easily extend to our setting (See Section 11 of [134]). We
describe their results below, as they are a nice compromise between the results of Gar´et-
Marchand and shell passage times of Cox-Durrett and Kesten. Let C∞be the inﬁnite
cluster for the Bernoulli percolation. For all x ∈Zd, let x∗∈Zd be the random point of
C∞such that ∥x −x∗∥1 is minimal, with a deterministic rule to break any possible ties.
Deﬁne
T ∗(x, y) = T(x∗, y∗)
and, for t ≥0,
B∗
t = {z + u|z ∈Zd, T ∗(0, z) ≤t, u ∈[−1/2, 1/2]d}.
20

<!-- page 21 -->

The time constant satisﬁes:
Theorem 2.21 (Theorem 4 [39], Theorem 1.2 [134]). Suppose that F([0, ∞)) > pc(d).
Then there exists µ∗(x) such that for all x ∈Zd,
lim
n
T ∗(0, nx)
n
= µ∗(x) in probability,
and
lim
n
T(0, nx)
n
= Z
in law,
where the distribution measure of Z is given by θ2δµ∗(x) + (1 −θ2)δ∞and θ = P(0 ∈C∞).
The limit shape obeys:
Theorem 2.22 (Theorem 5 [39], Theorem 1.2 [134]). Suppose that F([0, ∞)) > pc(d)
and F({0}) < pc(d). Then there exists a compact set B∗such that almost surely
lim
t λd
B∗
t
t △B∗

= 0,
where λd denotes the Lebesgue measure in Rd and A △B is the symmetric diﬀerence of
sets A and B.
2.5
Properties of the limit shape
2.5.1
Flat edges for limit shapes
In this section, we address one of the questions presented in the introduction.
Which compact convex sets C are realizable as limit shapes?
(2.12)
The question above is
completely open in the i.i.d. setting. A partial expected
answer is given by the following conjecture.
Question 4. Show that if F is a continuous distribution then the limit shape is
strictly convex.
Surprisingly, not even the following is known:
Question 5. Show that the d-dimensional cube is not a possible limit shape for a
FPP model with independent, identically distributed passage times.
Remark 2.23. Interestingly, question (2.12) is solved by H¨aggstr¨om and Meester [94]
in the case of stationary (not necessarily i.i.d.) passage times. They establish that any
non-empty compact, convex set C that is symmetric about the coordinate axes is a limit
shape for some FPP model with weights distributed according to a stationary (under
translations of Zd) and ergodic measure. This is in sharp contrast with the i.i.d. case
explained above.
21

<!-- page 22 -->

However, there is one class of weights where the limit shape is known in some direc-
tions. This collection was introduced by Durrett and Liggett [70] and further studied by
Marchand [132], Zhang [172, 173] and by Auﬃnger and Damron [16]. Its main feature
is the presence of a ﬂat edge for the limit shape, as we describe below. We will stick to
dimension 2 in what follows.
Write supp(ν) for the support of ν where ν[0, x] = F(x) is the probability distribution
of τe. Let Mp be the set of measures ν that satisfy the following:
1. supp(ν) ⊆[1, ∞) and
2. ν({1}) = p ≥⃗pc,
where⃗pc is the critical parameter for oriented percolation on Z2 (see, e.g., [69]). In [70],
it was shown that if ν ∈Mp then Bν has some ﬂat edges. The precise location of these
edges was found in [132]. To describe this, write B1 for the closed ℓ1 unit ball:
B1 = {(x, y) ∈R2 : |x| + |y| ≤1}
and write int B1 for its interior. For p >⃗pc let αp be the asymptotic speed of oriented
percolation [69], deﬁne the points
Mp =
1
2 −αp
√
2, 1
2 + αp
√
2

and Np =
1
2 + αp
√
2, 1
2 −αp
√
2

(2.13)
and let [Mp, Np] be the line segment in R2 with endpoints Mp and Np. For symmetry
reasons, the following theorem is stated only for the ﬁrst quadrant.
Theorem 2.24 (Durrett-Liggett [70], Marchand [132]). Let ν ∈Mp.
1. Bν ⊂B1.
2. If p <⃗pc then Bν ⊂int B1.
3. If p >⃗pc then Bν ∩[0, ∞)2 ∩∂B1 = [Mp, Np].
4. If p =⃗pc then Bν ∩[0, ∞)2 ∩∂B1 = (1/2, 1/2).
The angles corresponding to points in the line segment [Mp, Np] are said to be in the
percolation cone; see Figure 4 below.
Let βp := 1/2 + αp/
√
2, that is, deﬁne βp as the x coordinate of Np. Convexity and
symmetry of the limit shape imply that 1/µ(e1) ≥βp. A non-trivial statement about the
edge of the percolation cone came in 2002 when Marchand [132, Theorem 1.4] proved
that this inequality is in fact strict:
1/µ(e1) > 1/2 + αp/
√
2.
22

<!-- page 23 -->

Percolation cone
@B⌫
Mp
Np
Figure 4: Pictorial description of the limit shape Bν. If ν ∈Mp, the shape has a ﬂat
edge between the points Mp and Np. The limit shape is diﬀerentiable at Mp and Np.
Outside the percolation cone the limit shape is unknown.
In other words, Marchand’s result says that the line that goes through Np and is or-
thogonal to the x-axis is not a tangent line of ∂Bν. The following theorem builds on
Marchand’s result and technique and says that at the edge of the percolation cone, one
cannot have a corner.
Theorem 2.25 (Auﬃnger-Damron [16]). Let ν ∈Mp for p ∈[⃗pc, 1). The boundary ∂Bν
is diﬀerentiable at Np.
Remark 2.26. Theorem 2.25 is stated for the single point Np but, due to symmetry,
it is clearly valid for Mp and for the reﬂections of these two points about the coordinate
axes.
The theorem above shows that any measure in Mp has a non-polygonal limit shape.
The question of ﬁnding a single distribution where the limit shape is non-polygonal
was raised by H. Kesten [113]. The ﬁrst example of a non-polygonal limit shape was
discovered by Damron-Hochman [58].
23

<!-- page 24 -->

The ﬂat part of the percolation cone ends at Mp and Np; however, that does not
exclude the limit shape from having further ﬂat spots. This is not expected though.
Question 6. Show that for any measure ν ∈Mp the boundary of the limit shape is
not ﬂat outside the percolation cone.
A similar, but perhaps more ambitious question is to show:
Question 7. Show that if the limit shape of a measure ν has a ﬂat piece then ν ∈Mp,
for p > pc and the ﬂat piece is delimited by the percolation cone.
More approachable may be the following two open questions:
Question 8. Show that for any measure ν, in direction e1 the boundary of the
boundary of the limit shape does not contain any segment parallel to the e2 axis.
A solution to 8 immediately implies a solution to Question 5.
Question 9. Find a non-trivial example of a measure ν such that in direction e1 the
limit shape boundary is not parallel to the e2 axis.
One may reasonably guess that to solve Question 9, it suﬃces to consider small
perturbations of the trivial measure ν = δ1 (every edge has non-random passage time
equal to 1) as the limit shape in this case is just the ℓ1 ball, a diamond. This observation
was explored by Basdevant and co-authors by considering FPP in the highly super-critical
bond percolation cluster. Precisely, let τe be equal to 1 with probability p ∈(0, 1) and
inﬁnite with probability 1 −p.
Theorem 2.27 (Basdevant et al. [20]). For all 0 ≤λ ≤1, on the event that the origin
is connected to inﬁnity, almost surely we have
lim
n→∞
T(0, n(1, λ(1 −p))
n
= 1 + (1 −p)1 + λ2
2
+ O
 (1 −p)2
.
The theorem above roughly says that when p is close to 1, the four corners of the ℓ1
ball are replaced by curves that resemble parabolas.
2.6
The subadditive ergodic theorem revisited
Recall that the main tool to prove the existence of the time constant (Theorem 2.1)
was the subadditive ergodic theorem (Theorem 2.2). The fact that Theorem 2.2 requires
few assumptions makes it very powerful and widens its scope far beyond FPP. This
strength also comes with two main drawbacks.
First, it only provides the existence
of the limiting object.
The characterization as inf EX0,n/n is not always helpful to
obtain further information on the time constant. Second, if one tries to obtain diﬀerent
24

<!-- page 25 -->

results (ﬂuctuations, concentration, large deviations) for ergodic processes under the
same assumptions, one fails miserably. Nearly anything can happen.
The goal of this section is to dissect a major tool in Kingman’s original proof of the
subadditive ergodic theorem. Kingman’s proof provides a non-trivial decomposition of
the subadditive process. One term of this decomposition has the same properties as the
Busemann function in FPP, an object that we will study in detail in Section 5. Later
on, we will pursue this direction and provide extra assumptions that allow one to derive
further information about the subadditive process X0,n. In this chapter, we stick with
the ﬁrst task.
An array X = (Xm,n) satisfying the assumptions of Theorem 2.2 but also X0,n =
X0,m + Xm,n is called an additive ergodic sequence. The proof of Kingman’s theorem
depends on the following decomposition:
Theorem 2.28. If Xm,n is a subadditive ergodic sequence then there exist arrays Ym,n
and Zm,n such that
(a) Ym,n is an additive ergodic sequence with EY0,1 = µ = infn EX0,n/n.
(b) Zm,n is a non-negative subadditive ergodic sequence with time constant equal to 0.
(c) Xm,n = Zm,n + Ym,n.
The decomposition above is not necessarily unique. Here is an easy counter-example.
Let (yk)k≥1 be a sequence of independent standard Gaussians and set Sm,n = ym+1 +
. . . + yn. Then S is additive, with time constant 0. If we put
Xm,n = max(Sm,n, 0),
we see that X is subadditive, with EX0,n = (n/2π)1/2, so µ = 0. Two decompositions of
X are given by Y = 0, Z = X and Y = S, Z = max(−S, 0).
To see some implications of the above decomposition, suppose that the following limit
exists almost surely:
Bm,n := lim
N (Xm,N −Xn,N)
(2.14)
and satisﬁes E|B0,n| ≤Cn for some positive constant C.
We claim that Bm,n is an
additive process in the decomposition above. To see this, note that
B0,n = lim
N (X0,N −Xn,N) = lim
N (X0,N −Xm,N + Xm,N −Xn,N) = B0,m + Bm,n
and (b) and (c) from the deﬁnition in Theorem 2.2 follow from the stationarity of the
process (Xm,n). For instance,
Bm,n = lim
N (Xm,N −Xn,N)
=
lim
N (Xm,N−1 −Xn,N−1)
d=
lim
N (Xm+1,N −Xn+1,N) = Bm+1,n+1,
25

<!-- page 26 -->

where we used the fact that if XN
d= YN and XN, YN converge in distribution to X and
Y respectively then X
d= Y . The moment bounds follows from assumption.
Now set
Zm,n = Xm,n −Bm,n.
For any integer N > n, subadditivity gives Xm,n ≥Xm,N −Xn,N almost surely. Thus,
we also have Zm,n ≥0 almost surely. As B is additive, Zm,n is sub-additive and thus it
satisﬁes (b).
The importance of the limit in (2.14) will become clear in Section 5. At this point,
it would be interesting to determine when we can actually use (2.14).
Question 10. Find conditions that guarantee the existence of the limit (2.14).
The way that Kingman avoided the problem of existence of the limit in (2.14) was to
construct weak averaged limits. It will be beneﬁcial to explain his idea here. Let Ωbe
the space of all subadditive functions. A subadditive ergodic process is a random element
on Ω, inducing a probability measure P on this space. For any x ∈Ωwe deﬁne the shift
θx as the element in Ωgiven by
θx(m, n) = x(m + 1, n + 1).
Now let f be an element of L1(Ω, P) and deﬁne T : L1(Ω, P) →L1(Ω, P) as the bounded
linear operator given by
(Tf)(x) = f(θx).
Kingman’s magic was to construct a function f in L1(Ω, P) such that
f + Tf + . . . + T n−1f ≤X0,n
and Ef = g. Once this function is constructed, the reader can easily check that the
representation follows by taking
Ym,n = T mf + T m+1f + . . . + T n−1f
(2.15)
and Zm,n = Xm,n −Ym,n. Note that both sides in (2.15) are random variables in Ω. To
construct such f he considered the process:
fk = 1
k
k
X
i=1
(X0,i −X1,i)
(2.16)
and for each n ≥1 its iterates
fk + Tfk + . . . + T n−1fk = 1
k
k+n−1
X
i=1
(Xa,i −Xb,i)
(2.17)
26

<!-- page 27 -->

where a = max(i −k, 0) and b = min(i, n). The reader can now see the connection:
as k goes to inﬁnity, (2.16) and (2.17) play the role of B0,1 and B0,n, respectively. The
existence of the k limit turns out to follow from the Bourbaki-Alaoglu theorem. We will
come back to this point in Section 5. When the limit in (2.14) exists, we will call it a
generalized Busemann function for the subadditive process X.
2.7
Pointed Gromov-Hausdorﬀconvergence
FPP is a model of a random (pseudo)metric space. There is a classical way of deﬁning
convergence of a sequence of metric spaces, using the Gromov-Hausdorﬀdistance on the
space of metric spaces. We follow this route in this section and rephrase the limit shape
in this context. The reader is invited to check [32, 35] for a detailed explanation and
historic motivation of the basic topics we touch here. Although this approach brings a
diﬀerent perspective, these methods have not yet provided signiﬁcant new progress in
FPP on Zd. However, they were successfully used to extend the limit shape to FPP in
Cayley graphs with polynomial growth, where the subadditive ergodic theorem does not
immediately apply [24,154]. We will need a few deﬁnitions before we start, and most of
the following is taken from [32].
Deﬁnition 2.29. A subset S of a metric space X is said to be ϵ-dense if every point
of X lies in the ϵ-neighborhood of S.
Deﬁnition 2.30. An ϵ-relation between two (pseudo)metric spaces X1 and X2 is a
subset R ⊂X1 × X2 such that
1. For i = 1, 2, the projection of R to Xi is ϵ-dense.
2. If (x1, x2), (x′
1, x′
2) ∈R then |dX1(x1, x′
1) −dX2(x2, x′
2)| < ϵ.
If there exists an ϵ-relation between the metric spaces X1 and X2, we say that X1 and
X2 are ϵ-related and use the notation X1 ∼ϵ X2. When the projection of an ϵ-relation
is onto in both of its coordinates, we say that the relation is surjective and we write
X1 ≃ϵ X2. It is an exercise to show that if X1 ∼ϵ X2 then X1 ≃3ϵ X2.
The Gromov-Hausdorﬀdistance between X1 and X2 is deﬁned as:
DH(X1, X2) := inf{ϵ|X1 ≃ϵ X2}.
If there is no ϵ such that X1 ≃ϵ X2, then DH(X1, X2) is inﬁnite. It is not diﬃcult to see
that DH satisﬁes the triangle inequality on the space of metric spaces and thus it is a
pseudo-metric which may take the value inﬁnity.
Deﬁnition 2.31. We say that a sequence of (pseudo)metric spaces Xn converges to
X in the Gromov-Hausdorﬀmetric, and write Xn →X if and only if DH(Xn, X) →0
as n →∞.
27

<!-- page 28 -->

Gromov-Hausdorﬀconvergence works well in contexts where one deals with sequences
of compact metric spaces, but it is a less satisfactory concept when applied to non-
compact spaces. One disadvantage is that the distance between a compact space and an
unbounded set is always inﬁnity. Since the spaces that we care about are not compact,
we will need the following alternative deﬁnition of convergence. The magic here is that
the intuitive sense of convergence comes from observations from a ﬁxed point.
Deﬁnition 2.32. A pointed space (X, d, x) is a pair of a metric space (X, d) and a
point x ∈X. The point x is called the basepoint of the pointed space (X, d, x).
Deﬁnition 2.33. A sequence of pointed spaces (Xn, dn, xn) converges to a pointed
space (X, d, x) if for every r > 0 the sequence of closed balls B(xn, r) (with induced
metrics) converges to the closed ball B(x, r) in the Gromov-Hausdorﬀmetric.
One of the nice features of pointed Gromov-Hausdorﬀconvergence is that it preserves
several properties of the sequence of metric spaces in the limit. We will comment on this
at the end of the section.
Now, let’s go back to FPP. Let (Ω, F, P) be a probability space where all the τe’s are
deﬁned. For each ω ∈Ωwe deﬁne a sequence of pseudometric spaces:

Xn, dn(x, y)

=
1
nZd, T(nx, ny)
n

.
The pseudometric space Xn is just the original lattice rescaled by n with the normal-
ized pseudometric; that is, dn(0, x) = T(0, nx)/n. The origin 0 is a point of Xn for all
n ≥1. Now recall the construction from Section 2.3. Given an edge distribution on Zd,
there exists a norm µ on Rd where the unit ball in that norm is the limit shape of the
FPP model. The pair (Rd, µ) is a normed vector space with distance d(x, y) = µ(y −x)
for x, y ∈Rd. We will assume that the passage times have ﬁnite exponential moments.
This assumption is to make sure the metric satisﬁes a concentration bound given by
Lemma 2.35 below.
The limit shape theorem translates to the following statement.
Theorem 2.34. Assume that F(0) < pc(d) and
R
eαx dν < ∞for some α > 0. Al-
most surely, the sequence (Xn, dn, 0) converges in the pointed Gromov-Hausdorﬀsense to
(Rd, µ, 0).
Proof. Fix r > 0 rational. We ﬁrst show that almost surely the balls Bn := Bn(0, r) ⊂Xn
converge in the Gromov-Hausdorﬀsense to the ball
B := B(0, r) = {x ∈Rd|µ(x) ≤r}.
For this, it suﬃces to show that for any ϵ > 0 there exists n0 ∈N so that, for any n ≥n0
there is an ϵ-relation between Bn and B.
28

<!-- page 29 -->

We construct such a relation as follows. Fix ϵ > 0. Given 0 < ϵ′ < ϵ, use Theorem
2.16 to choose n0 = n0(ω) so that for n ≥n0
P

(1 −ϵ′)B ⊆¯Bn ⊆(1 + ϵ′)B for all n ≥n0

= 1,
where ¯Bn = Bn + [−1/(2n), 1/(2n)]d. The set R is deﬁned as the union of two sets,
R = R1 ∪R2 where
R1 :=

(x1, x2) ∈Bn × B : x2 = x1/(1 + ϵ′)

,
R2 :=

(x1, x2) ∈Bn × B : (1 −ϵ′)x2 ∈x1 + [−1/(2n), 1/(2n)]d

.
Note that R is surjective. Indeed, since Bn ⊆(1 + ϵ′)B, every element of Bn is related
to some element of B through R1 while (1 −ϵ′)B ⊆¯Bn implies that every element of B
is related to at least one element of Bn through R2. Now take (x1, x2) and (x′
1, x′
2) in R.
We have
|dn(x1, x′
1) −d(x2, x′
2)| = |T(nx1, nx′
1)/n −µ(x2 −x′
2)|
≤|T(nx1, nx′
1)/n −µ(x1 −x′
1)| + |µ(x1 −x′
1) −µ(x2 −x′
2)|.
(2.18)
Let’s ﬁrst look at the second term in the right side of (2.18). If (x1, x2) and (x′
1, x′
2) are
both in R1 we have
µ(x1 −x′
1) −µ(x2 −x′
2) = (1 + ϵ′)µ(x2 −x′
2) −µ(x2 −x′
2) = µ(x2 −x′
2)ϵ′ < 2rϵ′. (2.19)
If both pairs of points are in R2, then x2 = (1 −ϵ′)−1(x1 + zn), x′
2 = (1 −ϵ′)−1(x′
1 + z′
n)
for some zn, z′
n in [−1/(2n), 1/(2n)]d. We thus obtain
µ(x2 −x′
2) = (1 −ϵ′)−1µ(x1 −x′
1 + zn −z′
n).
(2.20)
Since µ is a norm, we can ﬁnd n1 ∈N so that for any n ≥max(n0, n1) and any |w| <
2
√
d/n, |µ(z + w) −µ(z)| ≤ϵ′. As |zn −z′
n| ≤2
√
d/n, we have by (2.20) and the triangle
inequality
|µ(x1 −x′
1) −µ(x2 −x′
2)| ≤((1 −ϵ′)−1 −1)µ(x1 −x′
1) + ϵ′(1 −ϵ′)−1 ≤2ϵ,
(2.21)
for ϵ′ small. If (x1, x2) ∈R1 and (x′
1, x′
2) ∈R2, then similarly to (2.19) and (2.21) we
obtain
|µ(x1 −x′
1) −µ(x2 −x′
2)| ≤2ϵ,
(2.22)
for suﬃciently small ϵ′. Thus a combination of (2.19), (2.21) and (2.22) tells us that if
we choose ϵ′ small enough
|µ(x1 −x′
1) −µ(x2 −x′
2)| < ϵ/2
for all n ≥max{n0, n1}.
(2.23)
The ﬁrst term of (2.18) is controlled by the following concentration bound.
29

<!-- page 30 -->

Lemma 2.35. Given ϵ > 0 and R > 0 there exists C1 > 0 such that for any n ≥1
P

∃x, y ∈[−Rn, Rn]d with T(x, y) > n(µ(y −x) + ϵ/2)

≤exp(−nC1).
Proof. See Theorem 3.11.
Combining (2.18), (2.23) and Lemma 2.35, we see that
P

R is an ϵ-relation between B and Bn

≥1 −exp(−nC1),
and thus by taking a countable sequence of ϵn →0 and using Borel-Cantelli we obtain
the desired result for each r ∈Q. However, if R is a ϵ-relation between Bn(0, r) and
B(0, r) then (the restriction of) R is also an ϵ-relation between Bn(0, r′) and B(0, r′) for
any r′ < r. This last observation suﬃces to end the proof of Theorem 2.34.
2.8
Strict convexity of the limit shape
In this section, we explore in more detail the conjecture that, under mild assumptions on
F, the limit shape (see Question 4) is strictly convex. We also introduce the deﬁnition
of uniform positive curvature, a concept related to strict convexity. In Sections 3 and 4,
we will discuss important results of Newman where this unproven property of uniform
positive curvature will play a major role. Strict convexity also plays an important part
in certain questions regarding the evolution of multi-type stochastic competition models
discussed in Section 6.
Recall that we call a subset B of Rd strictly convex if every line segment connecting
any two points of B is entirely contained, except for its endpoints, in the interior of B.
Let u be a unit vector of Rd and let H0 be a hyperplane such that u+H0 is supporting
hyperplane for µ(u)Bν at u (this means that u + H0 contains u and µ(u)Bν intersects
only one of the two halfspaces determined by u + H0). We introduce an exponent that
captures the nature of the boundary of Bν in direction u, called the curvature exponent,
as follows.
Deﬁnition 2.36 (Curvature Exponent). Assume that ∂Bν is diﬀerentiable.
The
curvature exponent κ(u) in the direction u is a real number such that there exist positive
constants c, C and ε such that for any z ∈H0 with |z| < ε, one has
c|z|κ(u) ≤µ(u + z) −µ(u) ≤C|z|κ(u).
(2.24)
Deﬁnition 2.37 (Uniformly curved). We say that Bν is uniformly curved if for every
unit vector u ∈Rd,
κ(u) ≥2,
with constants in (2.24) that are uniform in u.
30

<!-- page 31 -->

In the case that ∂Bν is not necessarily diﬀerentiable, Newman [135] gave a general
deﬁnition of uniform curvature: there exists C > 0 such that for all z1, z2 ∈∂Bν and
z = (1 −λ)z1 + λz2 with λ ∈[0, 1],
1 −µ(z) ≥C min{µ(z −z1), µ(z −z2)}2.
Either of these deﬁnitions is suitable for the results in this survey.
x
Figure 5: Representation of the limit shape (in red) in direction x with curvature expo-
nent κ(x) = 2. The limit shape curve ﬁts two tangent parabolas, one inside (orange),
other outside (black).
In two dimensions, the exponent κ(u) tells us that it is possible to trace two curves of
the form y = xκ(u) that are tangent to the limit shape, one inside and the other outside
of Bν (see Figure 5). For instance, a Euclidean ball is uniformly curved with κ(u) = 2
for every u. The ℓ1 ball is not uniformly curved as outside its corners one has κ(u) = 1.
Unfortunately, uniform curvature has not been proved for the limit shape of any FPP
model with i.i.d. passage times. Any advance in the direction of the following question
would be a major contribution.
Question 11. Show that for continuous distributions of passage times, the limit
shape is uniformly curved.
The importance of the notion of curvature will be revealed in the next two sections.
One characterization of curvature of the limit shape is to establish that B is a strict
convex set of Rd. A conditional proof of strict convexity was obtained by Lalley [126].
The two hypotheses of Lalley’s result, however, seem to be out of reach at this moment.
Hypothesis two may not be valid as for instance the Tracy-Widom distribution does not
have mean 0. We describe them now.
For u a ﬁxed nonzero vector in Rd let Lu be the ray through u emanating from the
origin.
(H1). For any convex cone A of Rd containing the vector u in its interior, and for
each δ > 0 there exists R = R(δ, A) < ∞such that the following is true: For each point
v ∈Zd ∩A at distance ≤2 from the line Lu, the probability that the time-minimizing
path from the origin to v is contained in A ∪{x : ∥x∥≤R} is at least 1 −δ.
31

<!-- page 32 -->

The second assumption requires a ﬂuctuation theorem for the normalized passage
times.
(H2). There exists a mean-zero probability distribution Gu on the real line and a
scalar sequence a(n) →∞such that as n →∞
T(0, nu) −nµ(u)
a(n)
d→Gu.
Theorem 2.38 (Theorem 1, [126]). Let u and v be linearly independent vectors in Rd
and assume hypothesis (H1) and (H2) for both u and v. Then for each λ ∈(0, 1),
µ(λu + (1 −λ)v) < λµ(u) + (1 −λ)µ(v).
2.9
Simulations
Although we are celebrating the ﬁftieth anniversary of the model, simulation studies on
ﬁrst passage percolation were somewhat limited until very recently. Initial work is due
to Richardson [141] in 1973, where the model with exponentially distributed weights was
analyzed. In [141] the limit shape Bν seemed to be curved, with a shape resembling
a circle. As one could imagine, these simulations were restricted due to limitations in
computer power. Further investigation (also in the Eden model) came in the work of
Zabolitzky and Stauﬀer [169], in 1986, and by Durrett and Liggett in 1981. In particular,
the numbers obtained in [169] indicate the predicted ﬂuctuation exponents for ξ = 2/3
and χ = 1/3 by theoretical physicists [110,111,122] (see next two sections for the study
of these exponents).
A major contribution was done recently in the beautiful and extensive work of Alm
and Deifjen [9] for FPP in two dimensions. Running 19 years of CPU time, in a cluster
of 28 Linux machines, they investigate the value of the time constant and the limit shape
for several continuous distributions. Their results are consistent with most of the famous
conjectures of the model. Their numerical simulations show strict convexity of the limit
shape, with a limit shape diﬀerent from a circle in all cases. The exponents for the
standard deviation of hitting times and for the ﬂuctuations of hitting points on lines also
matched the predicted values 1/3 and 2/3, respectively.
The paper of Alm and Deifjen also brought new ﬁndings to the table. It seems that
the time constant depends primarily on the mean of the minimum edge weight adjacent
to the origin, that is,
E4τe := E min{t1, t2, t3, t4},
where the t′
is are independent copies of τe, at least for continuous distributions which
are not too concentrated. They reported that the time constants along the axis and
along the diagonal for all simulated distributions have an almost perfect linear relation
to E4τe. The distributions are not scaled to have the same mean, but most of them have
inf supp ν = 0. They also suggest that if F(0) = 0 then µ ≥E4τe.
32

<!-- page 33 -->

Question 12 (Alm-Deifjen). Assume F(0) = 0. Show that µ ≥E4τe.
Figure 6: In red, simulation of the ball B(t) for passage times distributed according to
a Gamma random variable Γ(k, k) with parameters k = 1 (left) and k = 2 (right) with
t = 20000. In blue, the shape of a circle of radius T(0, te1). Simulations and ﬁgures are
courtesy of S.-E. Alm and M. Deijfen.
33

<!-- page 34 -->

Figure 7: In red, simulation of the ball B(t) for passage times distributed according to
a Gamma random variable Γ(k, k) with parameters k = 3 (left) and k = 4 (right) with
t = 20000. In blue, the shape of a circle of radius T(0, te1). Simulations and ﬁgures are
courtesy of S.-E. Alm and M. Deijfen.
3
Fluctuations and concentration bounds
The passage time between 0 and a vertex x ∈Zd can be approximated (almost surely)
as
T(0, x) = µ(x) + o(∥x∥1),
due to the shape theorem. Quantifying the error term is the main subject of this section.
It has been traditionally analyzed in two pieces:
o(∥x∥1) = T(0, x) −ET(0, x)
|
{z
}
random fluctuations
+
ET(0, x) −µ(x)
|
{z
}
non−random fluctuations
.
The reason this splitting is useful is that the ﬁrst term can typically be treated using
techniques from concentration of measure, whereas the second term is analyzed using, in
part, bounds on the ﬁrst.
3.1
Variance bounds
The most basic control on the random ﬂuctuation term is a variance bound.
It has
been predicted in the physics literature by simulation [162,169] and some scaling theory
[111, 122] (see also in Kesten [115]) that there is a dimension-dependent exponent χ =
χ(d) such that
Var T(0, x) ∼∥x∥2χ
1
.
(3.1)
34

<!-- page 35 -->

This exponent is expected to be universal; it should not depend on the underlying envi-
ronment, as long as F satisﬁes some mild moment conditions and the limit shape has no
ﬂat edges. The meaning of “∼” has not been made clear. For instance, it could be that as
x →∞, the ratio of both sides converges to a constant, is bounded away from 0 and ∞,
or even that the variance has the expression ∥x∥2χ+o(1)
1
, as in the current case of Bernoulli
percolation exponents. Regardless, the following dependence on d is predicted for χ: In
d
χ
1
1/2
2
1/3
3
?
·
·
dc
0 ?
d = 1, the passage time T(0, x) is just a sum of ∥x∥1 i.i.d. random variables, and so one
has χ(1) = 1/2 under any reasonable deﬁnition of “∼”. For d ≥2, the inﬁmum in the
deﬁnition of T(0, x) is predicted to produce sub-diﬀusive ﬂuctuations, giving χ < 1/2. It
is clear that χ should decrease with dimension, but there is not agreement on whether
χ(d) = 0 for all d at least equal to some dc, and some even debate whether χ →0.
The history of rigorous variance bounds begins with Kesten’s work [113, Theo-
rem 5.16], showing that Var T(0, ne1) ≤C(n/ log1/p n)2 for p = 9d + 3. Although this
bound is only logarithmically better than a trivial bound (say in the case of bounded
weights), the proof is far from trivial. In ’93, Kesten introduced the “method of bounded
diﬀerences” to FPP, and with this he was able to prove the best current bounds for χ:
0 ≤χ(d) ≤1/2 for all d ≥1 .
We will begin by giving a sketch of his argument using the Efron-Stein inequality.
Theorem 3.1 (Kesten [115]). Assume Eτ 2
e < ∞, P(τe = 0) < pc(d) and that the
distribution of τe is not concentrated at one point. There exist C1, C2 > 0 such that
for all non-zero x ∈Zd,
C1 ≤Var T(0, x) ≤C2∥x∥1 .
The proof will use the following inequality for functions of independent random vari-
ables. We use the notation x+ = max{0, x}.
Lemma 3.2 (Efron-Stein’s inequality). Let X1, X2, . . . be independent and let X′
i be an
independent copy of Xi, for i ≥1. If f is an L2 function of (X1, X2, . . .) then
Var f ≤
∞
X
i=1
E[(Zi −Z)+]2 ,
where Z = f(X1, X2, . . .) and
Zi = f(X1, . . . , Xi−1, X′
i, Xi+1, . . .) .
35

<!-- page 36 -->

Proof. Letting Σi = σ(X1, . . . , Xi), the proof consists of writing Var f = P∞
i=1 E∆2
i ,
where ∆i = E[f | Σi] −E[f | Σi−1] and then applying Jensen’s inequality to ∆2
i , along
with symmetry. See [72,151].
Proof of Theorem 3.1. We now apply Efron-Stein to the passage time, noting that the
condition Eτ 2
e < ∞implies that T has two moments. Then
Var T(0, x) ≤
∞
X
i=1
E[(Ti(0, x) −T(0, x))+]2 ,
where e1, e2, . . . is any enumeration of the edges and Ti is the passage time in the edge-
weight conﬁguration (τe) but with the weight τei replaced by an independent copy τ ′
ei.
Note that Ti(0, x) > T(0, x) only when both ei is in Geo(0, x), the intersection of all
geodesics from 0 to x in the original edge-weight conﬁguration (τe), and τ ′
ei > τei. Fur-
thermore, in this case, Ti(0, x) −T(0, x) ≤τ ′
ei. Therefore we obtain the bound
∞
X
i=1
E(τ ′
ei)21{e∈Geo(0,x)} .
By independence, this equals
Eτ 2
e E#Geo(0, x) .
Therefore we can conclude the upper bound for Var T(0, x) with the following lemma.
Lemma 3.3. Assume F(0) < pc(d). There exists C > 0 such that for all x ∈Zd,
E#Geo(0, x) ≤CET(0, x) .
Proof. Note that if τe ∈[a, b] almost surely, where 0 < a < b < ∞, then the statement is
easy to prove. Indeed, letting Γ be a geodesic from 0 to x, one has
a#Γ ≤T(Γ) = T(0, x) ≤b∥x∥1 ,
giving
#Geo(0, x) ≤(b/a)∥x∥1.
In the general case, one can modify the idea from [56, Cor. 1.4]. Setting G(0, x) to be
the maximal number of edges in any self-avoiding geodesic from 0 to x and
Yx = G(0, x)1{T(0,x)<aG(0,x)} ,
then [56, Prop. 1.3] gives existence of a > 0 such that for some C2 > 0, one has
P(Yx ≥n) ≤e−C2n
for all x ∈Zd and n ≥1. Then write
#Geo(0, x) ≤G(0, x) ≤a−1T(0, x) + Yx .
Taking expectations gives the result.
36

<!-- page 37 -->

The lower bound in Kesten’s theorem is easier. Setting Σ to be the sigma-algebra
generated by the 2d edge-weights for edges adjacent to 0, one has
Var T(0, x) ≥E (E[T(0, x) | Σ] −ET(0, x))2 = Var

E[T(0, x)|Σ]

.
Let t′
1, . . . , t′
2d be independent copies of the passage times of edges adjacent to 0 and
set Σ′ = σ(t′
1, . . . , t′
2d). Let T ′(0, x) be the passage time from 0 to x where we replaced
the edge weights of adjacent edges to the origin by t′
is. Then the right side of the display
above is equal to
1
2E

E[T(0, x)|Σ] −E[T ′(0, x)|Σ′]
2
.
Now pick a < b such that P(τe < a) > 0 and P(τe > b) > 0. By considering the
events A = {t{0,y} < a for all y ∈Zd with ∥y∥1 = 1} and B the same event but with < a
replaced by > b and with t{0,y}’s replaced by t′
is, one has the lower bound
E

E[T(0, x)|Σ] −E[T ′(0, x)|Σ′]

1A∩B
2
≥C(b −a)2 > 0
independently of x.
We end by restating three open questions discussed at the beginning of the section.
Question 13. Show that for any d ≥2, under suitable conditions on F, χ < 1/2. If
d = 2, show that χ = 1/3.
Question 14. Determine whether or not
lim
d→∞χ(d) = 0.
Question 15. For suitable d, show that χ > 0.
3.2
Log improvement to upper bound for d ≥2
The main tools used in the proof of Kesten’s bound were (a) E#Geo(0, x) ≤C∥x∥1 and
(b) Ti(0, x) −T(0, x) ≤τ ′
ei. To improve the variance upper bound we will need to use
one more piece of information: for most edges e, the probability that e is in a geodesic
from 0 to x is small (in x). Another way to say this is that each edge has small inﬂuence
on the variable T(0, x). This statement requires d ≥2 since in d = 1, each edge has high
inﬂuence (there is only one path from 0 to x in that case).
The ﬁrst proof of sublinear variance for T(0, x) was due to Benjamini-Kalai-Schramm
[23] in ’03 and applied only to τe that are Bernoulli: there exist 0 < a < b < ∞such that
τe takes values a or b with probability 1/2. This speciﬁc distribution was needed to take
advantage of Talagrand’s inﬂuence inequality on the hypercube [152]. The second proof
37

<!-- page 38 -->

was due to Bena¨ım-Rossignol [22] in ’08 and applied to distributions in the nearly-Gamma
class: those distributions that satisfy a log-Sobolev inequality similar to that for the
Gamma distribution. Their methods were based on entropy and used an inequality due
to Falik-Samorodnitsky [73] which replaced Talagrand’s inequality (a similar inequality
was derived by Rossignol [142]).
The most recent proof is due to Damron-Hanson-Sosoe [57] in ’14 and applies to all
distributions with 2+log moments. Their method follows that of Bena¨ım-Rossignol, but
replaces the representation of an edge-weight as a push-forward of a Gaussian variable
with a representation using a Bernoulli encoding. That is, each edge-weight is encoded as
an inﬁnite sequence of 0/1-valued random variables, and the Gross two-point inequality
[93] is used to bound the entropy.
Below we state the sublinear variance bound from [57], but we will sketch the proof
only in the simplest case (uniform weights), following [22], and indicating where compli-
cations arise in extending the argument.
Theorem 3.4. For d ≥2, suppose P(τe = 0) < pc and Eτ 2
e (log τe)+ < ∞. There exists
C > 0 such that for all x ∈Zd with ∥x∥1 > 1,
Var T(0, x) ≤C
∥x∥1
log ∥x∥1
.
To prove the theorem above, we will use Falik-Samorodnitsky’s inequality. Recall
that the entropy of a nonnegative random variable X is deﬁned as
Ent X = EX log X −EX log EX .
Now, let νn := ν ×. . .×ν be the uniform measure on [0, 1]n. In what follows, expectation
is with respect to νn. For 1 ≤k ≤n, deﬁne Σk to be the sigma-algebra generated by
the ﬁrst k coordinates in Rn, with Σ0 the trivial sigma-algebra. Let f : Rn →R be such
that Ef 2 < ∞. Last, deﬁne the martingale diﬀerence
∆k = E[f | Σk] −E[f | Σk−1] .
Theorem 3.5. (Falik-Samorodnitsky) Let f : Rn →R be nonconstant and such that
Ef 2 < ∞. Then
Var f log
"
Var f
Pn
k=1 (E|∆k|)2
#
≤
n
X
k=1
Ent ∆2
k .
(3.2)
Before we prove the above inequality, a few words of comment are needed. First, note
that by the martingale decomposition of the variance
Var f =
n
X
k=1
E∆2
k .
38

<!-- page 39 -->

Applying Jensen’s inequality,
Var f ≥
n
X
k=1
(E|∆k|)2 ,
(3.3)
so the term inside the logarithm in (3.2) is greater than equal to 1.
The inequality (3.2) is useful to obtain log-sublinear bounds if one can show that the
right side of (3.3) is lower order of Var f. We will see that this is the case in our setting
as the ratio can be shown to be at least of order nα, α > 0. Cutting the story short, the
log ∥x∥1 improvement in Theorem 3.4 comes from the appearance of the log term in the
left side of Falik-Samorodnitsky’s inequality.
Last, Equation (3.2) ﬁrst appeared in a paper of Falik-Samorodnitsky [73, Equation
(3)] as a functional version of an edge-isoperimetric inequality for Boolean functions. We
provide a diﬀerent proof than the one given in the Appendix of [73].
We start the proof of (3.2) with the following lemma.
Lemma 3.6. Let f be a nonnegative function on a probability space (Ω, F, P) such that
Entf 2 < ∞. Then,
Ef 2 log Ef 2
(Ef)2 ≤Ent(f 2).
(3.4)
If f is identically zero we interpret the left side of (3.4) to be 0.
Proof. As the inequality is preserved if we multiply f by any positive constant, we can
assume that Ef 2 = 1. In this case, the inequality reads
−log(Ef)2 ≤Ef 2 log f 2,
or, as f ≥0,
0 ≤Ef 2 log(fEf).
However, on the event f > 0, we can use the fact that 1 −x ≤log x−1 with x = (fEf)−1
to obtain
0 = Ef 2 −1 = E

f 2

1 −
1
fEf

; f > 0

≤Ef 2 log(fEf).
Therefore, (3.4) holds.
Proof of Theorem 3.5. We use Lemma 3.6 on each ∆k:
n
X
k=1
Ent(∆2
k) ≥
n
X
k=1
E(∆2
k) log E(∆2
k)
(E|∆k|)2 = −(Var f)
n
X
k=1
E(∆2
k)
Var f log (E|∆k|)2
E(∆2
k) .
39

<!-- page 40 -->

As Var f = Pn
k=1 E(∆2
k), we can apply Jensen’s inequality with the function x 7→−log x
to get the lower bound
−(Var f) log
n
X
k=1
E(∆2
k)
Var f
(E|∆k|)2
E(∆2
k) .
The above equation is exactly the left-side of (3.2).
With Falik-Samorodnitsky’s inequality in our hands, we turn back to the proof of
Theorem 3.4. As mentioned before, we will prove this theorem under the assumption
that the passage times are uniformly distributed in the interval [0, 1]. This assumption
allow us to use the fact that this probability measure satisﬁes a log-Sobolev inequality.
Precisely, for any n ≥1, under uniform measure on [0, 1]n, there exists C > 0 such that
for any f : [0, 1]n →R that is smooth, one has
Ent(f 2) ≤CE∥∇f∥2
2 .
The above inequality combined with Theorem 3.5 leads to
Var f log
"
Var f
Pn
k=1 (E|∆k|)2
#
≤C
n
X
k=1
E∥∇∆k∥2
2 = CE∥∇f∥2
2 ,
(3.5)
where in the last equality we used the fact that ∂i∆k = 0 if i̸ = k and Pn
k=1 E|∂k∆k|2 =
E∥∇f∥2
2.
Proof of Theorem 3.4. We now apply (3.5) to the passage time T(0, x), noting that since
we assumed that our edge-weights are uniform, T(0, x) is bounded, and so we can extend
the above inequality with n →∞as
Var T(0, x) log
"
Var T(0, x)
P∞
k=1 (E|∆k|)2
#
≤C
∞
X
i=1
E
 ∂
∂τei
T(0, x)
2
.
(3.6)
The derivative on the right is in the sense of distributions (since the passage time is not
a smooth function of the edge-weights) and it is relative to the i-th edge-weight, where
we have enumerated the edges in the lattice as e1, e2, . . .. When the edge-weights are not
bounded, one needs to argue that n can be taken to ∞more carefully, imposing that T
has at least 2 + ϵ moments (and this is guaranteed by existence of 2 + log moments for
τe), to exploit uniform integrability.
Next we use the fact that
∂
∂τei
T(0, x) = 1{ei∈Geo(0,x)} ,
40

<!-- page 41 -->

where we recall that Geo(0, x) is the intersection of all geodesics from 0 to x. This holds
Lebesgue almost surely, which is ok for us since the weights are uniform. So we obtain
an upper bound for the right side of (3.6) of
CE#Geo(0, x) ≤C∥x∥1 .
Here we have used Lemma 3.3. Note that this is the same bound we obtained from
Efron-Stein (in Kesten’s method in the last section) but now the advantage is that we
have an extra factor of log[· · · ] on the left of (3.6).
We are now left to show that P∞
k=1 (E|∆k|)2 is at most ∥x∥a
1 for some a < 1. If
we succeed in this, then (3.6) implies sub-linear variance. Indeed, in that case, either
Var T(0, x) is already ≤C∥x∥1−ϵ
1
for some ϵ > 0, or it is not, in which case, the term
log[· · · ] is at least order log ∥x∥1, and we divide it to the other side to complete the proof.
Unfortunately, it is not known how to show the bound on P∞
k=1 (E|∆k|)2. The reason
is that it is at most of order P∞
k=1(P(ek ∈Geo(0, x)))2, and the only information we
have on these probabilities is
∞
X
k=1
P(ek ∈Geo(0, x)) = E#Geo(0, x) ≍∥x∥1 .
If the geodesic prefers to take certain nearly deterministic edges (say in a small tube
centered on an ℓ2-geodesic from 0 to x), then the sum of squares can be of order ∥x∥1.
The work of Benjamini-Kalai-Schramm [23] introduced an averaging trick to get around
this. The main realization is that if the system is translation-invariant, we can give the
appropriate inequality. First, we can restrict attention to a box around 0 of size C∥x∥1
for a large constant C. If all the probabilities are equal, they must be of order ∥x∥1−d
1
.
Plugging this in gives the correct bound.
So we consider an averaged passage time (this form of averaging was used in [7,150])
Fm =
1
Zm
X
∥z∥1≤m
T(z, z + x) ,
where the sum is over integer sites only and m =
j
∥x∥1/4
1
k
. Zm is the number of terms
in the sum. By Jensen’s inequality, we can still obtain the same upper bound in (3.5),
using f = Fm:
∞
X
i=1
E

∂
∂ei
1
Zm
X
∥z∥1≤m
T(z, z + x)


2
≤
1
Zm
X
∥z∥1≤m
∞
X
i=1
P(ei ∈Geo(z, z + x)) ,
which is bounded by C∥x∥1. Furthermore, Fm is not too diﬀerent from T(0, x): by the
triangle inequality,
∥Fm −T(0, x)∥2 ≤
1
Zm
X
∥z∥1≤m
∥T(z, z + x) −T(0, x)∥2 ≤2 max
∥z∥1≤m ∥T(0, z)∥2 .
41

<!-- page 42 -->

By our bound on the weights,
∥Fm −T(0, x)∥2 ≤2∥x∥1/4
1
,
which is o(∥x∥1/2/ log ∥x∥1), so it suﬃces to bound Var Fm. By the arguments in the
beginning of the proof, we need only show that P∞
k=1 (E|∆k|)2 ≤∥x∥a
1 for some a < 1,
where ∆k is the martingale diﬀerence associated to Fm.
By applying Jensen’s inequality (as in the proof of Efron-Stein),
E|∆k| ≤
2
Zm
X
∥z∥1≤m
E|Ti(z, z + x) −T(z, z + x)|+ ≤
2
Zm
X
∥z∥1≤m
P(ek ∈Geo(z, z + x)) .
By translation invariance, we obtain the upper bound
2
Zm
X
∥z∥1≤m
P(ek −z ∈Geo(0, x)) ≤Cm/Zm .
Here we have used an extension of Lemma 3.3 found in [57]. So
∞
X
k=1
(E|∆k|)2 ≤2Cm
Z2
m
X
∥z∥1≤m
∞
X
k=1
P(ek ∈Geo(z, z + x)) ≤Cm
Zm
∥x∥1 .
This is bounded by Cm1−d∥x∥1 ≤C∥x∥3/4
1
and completes the proof.
In the general case (assuming only Eτ 2
e (log τe)+ < ∞), the diﬀerence in the proof is
in the entropy bound. By writing each τe as the push-forward of an inﬁnite sequence
ωe = (ωe,1, ωe,2, . . .) of Bernoulli random variables, one can apply the Gross two-point
entropy bound for
∞
X
k=1
Ent ∆2
k ≤C
∞
X
k=1
X
e
E (∆e,kFm) ,
where ∆e,k is the discrete derivative operator relative to ωe,k. To give the upper bound
C∥x∥1 for the right side, one needs a careful analysis of these discrete derivatives and
tools from the theory of greedy lattice animals. See [56] for details.
3.3
Log improvement to lower bound for d = 2
The following theorem by Newman-Piza in ’95 represents the state of the art for the
lower bound on the variance of passage times. It improves on Kesten’s lower bound by
a factor of log n.
Theorem 3.7 (Newman-Piza [137]). Let λ = inf{s : P(τe ≤s) > 0}. Assume Eτ 2
e < ∞
and Var τe > 0. Assume in addition that one of the following two conditions is satisﬁed:
λ = 0 and P(τe = 0) < pc
42

<!-- page 43 -->

λ > 0 and P(τe = λ) <⃗pc .
If d = 2 then there is a constant B > 0 such that
Var T(0, nx) ≥B log n
(3.7)
for n ≥1 and all unit vectors x ∈R2.
In the case when τe is exponential mean 1, the log n lower bound was also obtained,
using diﬀerent methods, by Pemantle and Peres [138]. Theorem 3.7 was extended to
distributions in Mp (whose limit shape has a ﬂat edge – see Section 2.5.1) ﬁrst in the
e1 direction by Zhang [173] and then for all directions outside the percolation cone by
Auﬃnger and Damron [16] (see also Kubota [124, Corollary 1.4], who reduced the moment
condition of [16]). It is important to note, however, that inside the percolation cone, the
variance of the passage time is of order constant [173], and this is a strong version of
χ = 0 in those directions.
Newman-Piza used a martingale method to obtain the desired lower bound in (3.7).
Their technique goes as follows. Enumerate the edges (ei) of Z2 in a spiral order starting
from the origin and deﬁne a ﬁltration (Σi), where Σi is the sigma-algebra generated by
the weights τe1, . . . , τei and Σ0 is the trivial sigma-algebra. Writing T = T(0, nx), we
may use L2-orthogonality of martingales to ﬁnd
Var T =
∞
X
i=1
E [E[T | Σi] −E[T | Σi−1]]2 .
The i-th term in the sum represents a part of the contribution to the variance given by
ﬂuctuations of the passage time of the i-th edge ei. The idea of [137] is that if ei is in a
geodesic from 0 to nx then if we lower its passage time (while keeping all other weights
ﬁxed), the variable T will decrease linearly, and therefore τei will have inﬂuence on the
ﬂuctuations of T. They consequently argue (see Theorem 3.8 below) that one has a lower
bound for Var T of
C
∞
X
i=1
P(Fi)2 ,
(3.8)
where Fi is the event that ei is in a geodesic from 0 to nx.
Precisely, one can summarize the Newman-Piza lower bound method in the following
statement.
Let I be a countable set. Consider the probability space (Ω, F, P) where Ω= RI =
{ω = (ωi
:
i ∈I)}, F is the Borel sigma-ﬁeld BI and P is a probability measure.
Suppose T is a random variable with E(T 2) < ∞. For W ⊆I, write F(W) for the
Borel sigma-ﬁeld BW. Let U1, U2, . . . be disjoint subsets of I and express ω for each k as
(ωk, ˆωk), where ωk (resp. ˆωk) is the restriction of ω to Uk (resp. to I \ Uk). For each k,
let D0
k and D1
k be disjoint events in BUk. Deﬁne
Hk(ω) = T 1
k (ˆωk) −T 0
k (ˆωk) ,
43

<!-- page 44 -->

where
T 0
k (ˆωk) = sup
ωk∈D0
k
T((ωk, ˆωk)) and T 1
k (ˆωk) =
inf
ωk∈D1
k
T((ωk, ˆωk)) .
Theorem 3.8 (Newman-Piza [137]). Assume the setting just described and the following
three hypotheses about P, the Uk’s, the Dδ
k’s and T:
1. Conditional on F(I \ ∪kUk), the F(Uk)’s are mutually independent.
2. There exist a, b > 0 such that, for any k,
P
 ωk ∈D0
k | F(U c
k)

≥a and P
 ωk ∈D1
k | F(U c
k)

≥b a.s.
3. For every k, Hk ≥0 a.s.
Suppose that, for some ε > 0 and each k, Fk ∈F is a subset of the event {Hk ≥ε}.
Then
Var T ≥abε2 X
k
P(Fk)2 .
(3.9)
We will show in detail how (3.9) is used. We will use the approach of [16]. Before
doing so, let us sketch how Newman-Piza proceeded in their original paper.
When Fi is the event that ei is in a geodesic from 0 to nx, the sum on the right-side
of (3.9) is exactly the expected overlap (the number of common edges) of two geodesics
from 0 to nx sampled independently. However, this overlap is diﬃcult to control and
Newman-Piza replace the right-hand side by an expression involving the partial sums:
k
X
j=1
P(Fj) = E(number of minimizing edges among the ﬁrst k edges).
(3.10)
Once (3.10) is established, the logarithmic lower bound follows after showing that
(3.10) is bounded below by B
√
k for k ≤Cn2. This step is where the 2 dimensional
assumption comes into play. They use the fact that within a box of ℓ1-distance L to the
origin, there are order L2 edges. Thus, the right side of (3.10) is bounded below by the
expected passage time to the boundary of a box with radius
√
k, which is, by the shape
theorem, of order B
√
k.
A solution of the following question, combined with the Newman-Piza method, would
provide an improvement on the lower bound of the variance.
Question 16. Show that in d = 2 there exists a > 0 such that the expected overlap
of two geodesics from 0 to nx sampled independently is at least na.
The above technique works well if we are allowed to lower the edge weights along a
geodesic. In the class Mp, when the passage time distribution has an atom at 1, the
bottom of its support, this turns out to be a problem. Therefore if we are to use the
44

<!-- page 45 -->

same technique, we must show that geodesics use many edges with weight above 1. It is
not enough only to know this though; if one repeats the computations above, one ﬁnds
only a lower bound of a constant (with no logarithm term). It is essential also to know
information about the location of these non-one edges on the lattice. In particular, if they
are heavily concentrated enough near the origin, we can extract a logarithmic bound. To
do this, we need to know something about the geometry of geodesics (for instance, that
they avoid certain regions of the plane). The proof that we will give will also allow us to
extend Theorem 3.8 as follows.
Recall that Np is the edge of the percolation cone, and let θp be the unique angle such
that the line segment connecting 0 and Np has angle θp with the e1-axis. Let wθ be the
vector
wθ = (cos θ, sin θ).
Theorem 3.9. Let µ ∈Mp for p ∈[⃗pc, 1) and θ ∈[0, θp). Suppose that Eτ 2
e < ∞. Then
there exists Cθ > 0 such that for all n,
Var T(0, nwθ) ≥Cθ log n .
(3.11)
Fix y > 1 such that
q = ν([y, ∞)) > 0.
Assume that there exists ρ > 0 such that given ε > 0 we may ﬁnd 0 < ζ < 1 and n0 ∈N
such that n > n0 and nζ ≤M ≤n/2 imply that with probability at least 1 −ε/2,
(A) every geodesic from 0 to nwθ contains at least ρM edges in the set MBν with
weights at least equal to y.
Also assume that we can ﬁnd c > 1 and L0 > 0 such that with probability at least
1 −ε/2,
(B) for all L > L0 every geodesic from 0 to nwθ contains at most cL edges in LBν.
These assumptions are not diﬃcult to verify given the information that we already
know about the limit shape. Indeed, (B) is a straight-forward consequence of the shape
theorem and the fact that 1 is the inﬁmum of the support. To verify (A), we use the
following Lemma, taken from [16, Lemma 9].
Lemma 3.10. Let C be a closed subset of ∂Bν that does not intersect the ﬂat edge [Mp, Np]
or any of its reﬂections about the axes. Given ε > 0, there exists M1 and ρ > 0 such
that with probability at least 1 −ε, the following holds. For all M > M1, x ∈MC and
for every geodesic Γ from 0 to x, at least ρM edges of Γ have passage times ≥y.
In the rest of this subsection, we explain how to prove Theorem 3.9 using Theorem
3.8 and assumptions (A) and (B).
45

<!-- page 46 -->

Proof of Theorem 3.9. Recall that p = ν({1}). We use Theorem 3.8 with T = T(0, nwθ),
Uk = {ek}, D0
k the event that τek = 1 and D1
k the event that τek ≥y. Last, set ε = y −1,
a = p and b = q. The reader may verify that all of the hypotheses of Theorem 3.7 are
satisﬁed. Therefore, setting
Fk =

τek ≥y and ek is in a geodesic from 0 to nwθ
	
,
we ﬁnd that
Var T(0, nwθ) ≥pq(1 −y)2 X
k
P(Fk)2 .
(3.12)
We now deﬁne a sequence of numbers (Mi) by M1 = nζ and, setting J = 2c/ρ, for
i = 1, . . . , I (for I = ⌊(1 −ζ)/ log J⌋log(n/2)),
Mi+1 = JMi .
From (3.12) it follows that if we write Ei, i = 0, . . . , I −1 for the set of edges in
Mi+1Bν but having no endpoints in MiBν (for i = 0 we take all edges with an endpoint
in MiBν), then
Var T(0, nwθ) ≥pq(1 −y)2
I−1
X
i=0
X
k:ek∈Ei
P(Fk)2 .
The idea of the above decomposition is that if the aspect ratio J of the annuli is large
enough then the inner sum will be always at least of order constant. Since there is a
logarithmic number of such annuli, we will obtain the desired lower bound.
Indeed, using Jensen’s inequality, we get a lower bound of
pq(1 −y)2
I−1
X
i=0
1
|Ei|
 X
k:ek∈Ei
P(Fk)
!2
.
(3.13)
We will now give a lower bound for the inner sum. Call Xi+1 the event that (A) and
(B) hold for M = Mi+1. On this event the number of edges on any geodesic Γ from 0 to
nwθ in the set MiBν is at most (ρ/2)Mi+1 and the number of edges on Γ in Mi+1Bν with
weight at least y is at least ρMi+1. From this it follows that on the event Xi+1, we have
the lower bound
X
k:ek∈Ei
I(Fk) ≥(ρ/2)Mi+1 .
Since Bν is bounded, there exists d > 0 such that for all i, Mi ≥d
p
|Ei|. So (3.13) is
bounded below by
pq(1 −y)2 [d(ρ/2)(1 −ε)]2 I ≥C log n
for C independent of n. This completes the proof.
46

<!-- page 47 -->

3.4
Concentration bounds
3.4.1
Subdiﬀusive concentration
In addition to sub-linear variance bounds, there has been work to establish concentra-
tion inequalities for T(0, x) on the scale (∥x∥/ log ∥x∥)1/2. These have so far only been
exponential inequalities, not Gaussian ones.
As in the case of sub-linear variance, the ﬁrst such exponential inequality was not for
general distributions, only those in the “nearly Gamma” class. The result we present
below is from [56], and only moment conditions are needed. Note that the condition for
the lower-tail inequality is weaker than that for the upper-tail inequality.
Theorem 3.11. Let d ≥2 and suppose that P(τe = 0) < pc. If Eeατe < ∞for some
α > 0 then there exist c1, c2 > 0 such that for all x ∈Zd with ∥x∥1 > 1,
P
 
T(0, x) −ET(0, x) ≥λ
s
∥x∥1
log ∥x∥1
!
≤c1e−c2λ for λ ≥0.
If Eτ 2
e (log τe)+ < ∞, then for all x ∈Zd with ∥x∥1 > 1,
P
 
T(0, x) −ET(0, x) ≤−λ
s
∥x∥1
log ∥x∥1
!
≤c1e−c2λ for λ ≥0.
Problem. Is the bound on the right side optimal?
The main strategy is again due to Bena¨ım-Rossignol and follows the same lines as
their proof of sub-linear variance. The proof of the general case in [56] again involves
a Bernoulli encoding and estimating discrete derivatives after applying the two-point
entropy estimate. Either way, one deﬁnes the averaged passage time
Fm =
1
Zm
X
∥z∥1≤m
T(z, z + x),
where the sum is over integer sites only, m =
j
∥x∥1/8
1
k
, and Zm is the number of terms
in the sum. One can show [56, Section 2.1] that it suﬃces to derive the concentration
inequality for Fm.
The main idea is to obtain a variance estimate for an exponential function of Fm
analogous to the one obtained by T. By following the sub-linear variance strategy, with
more technical diﬃculty, one obtains the following inequality: for some C > 0,
Var eλFm/2 ≤Kλ2EeλFm < ∞for |λ| <
1
2
√
K
and ∥x∥1 > 1,
(3.14)
where K =
C∥x∥1
log ∥x∥1. Note the similarity to the entropy bound obtained in (3.17) in the
proof of Talagrand’s theorem. The above inequality can be thought of as weaker, due
47

<!-- page 48 -->

to the presence of the variance instead of entropy, but stronger due to the logarithmic
factor.
The above variance estimate is turned into an exponential concentration bound using
the “iteration method.” The following comes from [31, p. 70-71].
Proposition 3.12. If EZ = 0 and for some constants B, C satisfying 0 < C ≤B,
Var etZ/2 ≤Ct2EetZ < ∞for t ∈(0, B−1/2) ,
then putting ψZ(t) = log EetZ, one has
ψZ(t) ≤−2 log(1 −Ct2) for t ∈(0, B−1/2) .
Proof. Beginning with
EetZ −
 EetZ/22 = Var etZ/2 ≤Ct2EetZ ,
we obtain
EetZ ≤
1
1 −Ct2
 EetZ/22 ,
or
ψZ(t) ≤−log(1 −Ct2) + 2ψZ(t/2) .
By induction, for n ≥1,
ψZ(t) ≤−
n
X
k=1
2k−1 log(1 −C(t/2k−1)2) + 2nψZ(t/2n) .
Because EZ = 0,
2nψZ(t/2n) = tψZ(t/2n) −ψZ(0)
t/2n
→tψ′
Z(0) = 0 ,
so
ψZ(t) ≤−
∞
X
k=1
2k−1 log(1 −C(t/2k−1)2) .
(3.15)
To bound these terms, we use the fact that −u−1 log(1−u) is non-decreasing in u ∈(0, 1).
Therefore
−C−1(t/2k−1)−2 log(1 −C(t/2k−1)2) ≤−C−1t−2 log(1 −Ct2) for k ≥1 .
Rewritten,
−2k−1 log(1 −C(t/2k−1)2) ≤−2−(k−1) log(1 −Ct2) .
Placing this in (3.15),
ψZ(t) ≤−log(1 −Ct2)
∞
X
k=1
2−(k−1) = −2 log(1 −Ct2) .
48

<!-- page 49 -->

To make the last bound more clear, we can use the inequality
−log(1 −u) ≤
u
1 −u for u ∈[0, 1) ,
which follows from the mean value theorem.
So we obtain
ψZ(t) ≤
2Ct2
1 −Ct2 ≤
4Ct2
2(1 −
√
Ct)
for t ∈(0, B−1/2).
In our case, B = C = K in (3.14), so we get
log Eet|Fm−EFm| ≤
4Kt2
2(1 −
√
Kt)
for t ∈(0, K−1/2).
By using the bound
P(X ≥x) = P(eλX ≥eλx) ≤e−λxEeλX,
which is valid for λ ≥0, and optimizing over λ, one can then complete the ﬁrst inequality
of Theorem 3.11.
3.4.2
Talagrand’s theorem via the entropy method
Here we will give the concentration argument from [56, Corollary A.5]. The goal will
be to give an exponential concentration inequality for the passage time about its mean
assuming certain moment conditions for the edge weights. This is not the main result of
that paper, but an auxiliary one used to obtain the main one. The lower tail inequality
comes from [61]. Our aim will be to prove the following result. It was initially established
by Talagrand in [153], using diﬀerent methods.
Theorem 3.13. Let d ≥2. Assuming P(τe = 0) < pc and Eeατe < ∞for some α > 0,
there exist C1, C2 > 0 such that
P

T(0, x) −ET(0, x) ≥t
p
∥x∥1

≤e−C1t2 for t ∈

0, C2
p
∥x∥1

.
If P(τe = 0) < pc and EY 2 < ∞, where Y is the minimum of d i.i.d. copies of τe, then
also
P

T(0, x) −ET(0, x) ≤−t
p
∥x∥1

≤e−C1t2 for all t ≥0.
Proof. Write T = T(0, x). We will show only the upper-tail inequality and we will assume
Eατe < ∞for all α > 0 for simplicity. The idea will be to set ψ(λ) = log Eeλ(T−ET) and
to note that it suﬃces to show for some C > 0 independent of x,
ψ(λ) ≤C∥x∥1λ2 for λ ∈(0, C).
(3.16)
49

<!-- page 50 -->

Indeed, one has for λ > 0
P(T −ET ≥t
p
∥x∥1) ≤P(eλ(T−ET) ≥etλ√
∥x∥1) ≤eψ(λ)−tλ√
∥x∥1
≤eC∥x∥1λ2−tλ√
∥x∥1.
Choosing λ = t/2C
p
∥x∥1 would then complete the proof.
To show the bound (3.16), we will use the Herbst argument. Setting
Ent X = EX log X −EX log EX
for a nonnegative random variable X, we will aim to show for some C > 0 independent
of x,
Ent eλT ≤C∥x∥1λ2EeλT for λ ∈(0, C).
(3.17)
This implies that
d
dλ
ψ(λ)
λ

≤C∥x∥1 for λ ∈(0, C),
and this gives (3.16).
So we focus on proving (3.17). Enumerating the edge variables (in any deterministic
order) as (τe1, τe2, . . .), then
Ent eλT ≤
∞
X
i=1
E Enti eλT ,
(3.18)
where Enti is entropy relative to only the edge-weight τei. This is known as “tensorization
of entropy”.
We now apply a modiﬁed Log-Sobolev inequality of Boucheron-Lugosi-
Massart [31, Theorem 6.15].
Lemma 3.14 (Symmetrized modiﬁed LSI). Let q(x) = x(ex −1). If X is a random
variable and X′ is an independent copy, then for all t ∈R,
Ent etX ≤E

etXq(λ(X′ −X)+)

.
Use the symmetrized LSI in (3.18):
Ent eλT ≤
∞
X
i=1
EeλTq(λ(T ′
i −T)+) .
(3.19)
Here T ′
i is the passage time from 0 to x in the edge-weight conﬁguration in which the
weight τei is replaced by an independent copy τ ′
ei, and all other edges remain the same
(that is, they are equal to τe). As in the FPP variance proof, we know that T ′
i −T is
only positive if ei is in a geodesic from 0 to x in the original edge-weights (τe). We claim
50

<!-- page 51 -->

more here: that if we deﬁne G(0, x) as the collection of edges in the intersection of all
geodesics from 0 to x (since there need not be a unique one), then
T ′
i −T > 0 ⇒ei ∈G(0, x) in the original weights (τe) .
To argue this, assume that T ′
i −T > 0. Then let γ be a geodesic from 0 to x in the edge-
weights (τe). If ei /∈γ then γ must have the same passage time in the new weights, since
we replace only τei by an independent copy. In other words, all edge weights for edges on
γ have the same value in both conﬁgurations. Therefore if T ′ denotes the passage time
in the new weights (with only τei replaced),
T ′
i = T ′(0, x) ≤T ′(γ) = T(γ) = T(0, x) = T ,
contradicting that T ′
i −T > 0.
Returning to (3.19), we can give the upper bound
Ent eλT ≤
∞
X
i=1
EeλTq(λ(T ′
i −T)+)1{ei∈G(0,x)} .
The function x 7→q(x) is monotone increasing for x ≥0 so using the bound T ′
i −T ≤τ ′
ei
(which we established during the FPP variance bound) and independence,
Ent eλT ≤
∞
X
i=1
EeλTq(λτ ′
ei)1{ei∈G(0,x)} = Eq(λτe) E

eλT#G(0, x)

.
(3.20)
To apply the Herbst argument we would love to decouple eλT from #G(0, x). Unfortu-
nately, the variable #G(0, x) is not bounded, so we cannot just pull it out. So we use a
variational characterization of entropy:
Ent X = sup{EXY : EeY ≤1},
which implies that for X ≥0 and any Y ,
EXY ≤Ent X + EX log EeY .
We therefore write for arbitrary a > 0,
EeλT#G(0, x) ≤a Ent eλT + a EeλT log E exp
#G(0, x)
a

.
Combining with (3.20), if aEq(λτe) < 1,
Ent eλT(1 −aEq(λτe)) ≤aEq(λτe)EeλT log E exp
#G(0, x)
a

,
or
Ent eλT ≤
aEq(λτe)
1 −aEq(λτe) log E exp
#G(0, x)
a

× EeλT .
To control these terms we will need a lemma from [56] which is a strong version of a
geodesic length bound.
51

<!-- page 52 -->

Lemma 3.15. Assuming Eeατe < ∞for all α > 0, there exist a, c1 > 0 such that
log E exp
#G(0, x)
a

≤c1∥x∥1 for all x ∈Zd .
Applying the lemma,
Ent eλT ≤c1∥x∥1
aEq(λτe)
1 −aEq(λτe) EeλT if aEq(λτe) < 1 ,
or
Ent eλT ≤2c1∥x∥1aEq(λτe)EeλT if aEq(λτe) < 1/2 .
Note that by dominated convergence,
lim
λ↓0
Eq(λτe)
λ2
= lim
λ↓0 E
τe(eλτe −1)
λ

= Eτ 2
e .
So we can ﬁnd c2 > 0 such that if λ ∈(0, c2) then Eq(λτe) ≤2λ2Eτ 2
e . For such λ, note
that aEq(λτe) < 1/2 if 2aλ2Eτ 2
e < 1/2, which occurs if λ is smaller than some positive
c3 < c2. So we obtain the desired equation (3.17): for some c4 > 0,
Ent eλT ≤4c1∥x∥1λ2Eτ 2
e EeλT ≤λ2c4∥x∥1EeλT for λ ∈(0, c3) .
3.5
Convergence of the mean for sub-additive ergodic processes
We return to the perspective adopted in Section 2.6. Let X be a subadditive ergodic
process satisfying the hypotheses of Theorem 2.2. As discussed previously, since the
leading order (shape theorem) limiting behavior of FPP is established using the general
framework of the subadditive ergodic theorem, one would hope that general abstract
arguments for subadditive sequences could be used to establish limit theorems and other
sharper results. Since the class of processes satisfying Theorem 2.2 is much too large
to characterize the exponents and limiting behavior of FPP, we will have to impose
additional axioms in order to derive useful results. Furthermore, many results will require
us to go beyond considering subadditive sequences X as above and instead consider the
d-dimensional structure of the model.
We will restrict our attention here to the convergence of the mean EX0,n/n of our
subadditive processes to their limiting µ (henceforth, we assume that the sequence is
ergodic).
Recall that in the FPP setting, the ﬂuctuations of T(0, x) −µ(x) can be
written as:
o(∥x∥1) = T(0, x) −ET(0, x)
|
{z
}
random fluctuations
+
ET(0, x) −µ(x)
|
{z
}
non−random fluctuations
.
Earlier in this section, we discussed bounds on the random ﬂuctuations in FPP. We
therefore shift to the problem of convergence of the mean in order to characterize the
other error term in the convergence to the FPP limit shape.
52

<!-- page 53 -->

3.5.1
Non-random ﬂuctuations for subadditive sequences
All of the existing methods for controlling the rate of convergence of the mean require also
some control of the corresponding non-random ﬂuctuations. As mentioned before, it is
expected that the random ﬂuctuations in FPP are governed by the ﬂuctuation exponent
χ. It is reasonable to postulate the existence of a similar exponent γ for the non-random
ﬂuctuations. Speciﬁcally, one expects that (for a suitable deﬁnition of “∼”)
T(0, x) −ET(0, x) ∼∥x∥χ
1,
ET(0, x) −µ(x) ∼∥x∥γ
1
for some exponents χ and γ. This motivates the following deﬁnitions from [18], made for
any subadditive ergodic sequence (Xm,n) satisfying the hypotheses of Theorem 2.2.
Deﬁnition 3.16. The exponents γ and γ are deﬁned as
γ = lim inf
n
log(EX0,n −nµ)
log n
and γ = lim sup
n
log(EX0,n −nµ)
log n
.
Here, recall µ is the almost-sure limit limn X0,n/n; we make the convention that
log 0 = −∞.
Deﬁnition 3.17. For p > 0, the ﬂuctuation exponents χp and χp are deﬁned as
χp = lim inf
n
log ∥X0,n −EX0,n∥p
log n
and χp = lim sup
n
log ∥X0,n −EX0,n∥p
log n
.
Note that χp ≤χp and by Jensen’s inequality,
χp ≤χq and χp ≤χq if p ≤q .
Remark 3.18. The deﬁnitions above only specify the ﬁrst order growth of the ﬂuc-
tuation terms. For instance, in FPP with Xm,n = T(me1, ne1), if there exist positive
constants c, C such that
cn2χ ≤Var T(0, ne1) ≤Cn2χ for all large n ,
then χ2 = χ2 = χ, and similarly for p̸ = 2 and γ, γ.
Our ﬁrst result will be a lower bound of the form “γ ≥χ” under some assumptions
on the sequence which are natural in the case of FPP. We explain these and the intuition.
Recall that the χ exponent in FPP is, except in the one-dimensional case, expected to
be strictly smaller than 1/2. Suppose that we had a sub-additive ergodic sequence X0,n
for which χ and γ exponents existed in a moderately strong sense (i.e. χp = χq for all
q, p ≥1 and γ = γ). Assume also that χ < 1/2 but γ < χ. Last, assume the sequence
(XKm,Kn) is weakly dependent for large K in a sense which we leave imprecise.
53

<!-- page 54 -->

For any positive integer K, we can write using subadditivity
[X0,n −EX0,n] ≤
n/K
X
i=1
[X(i−1)K,iK −EX0,K] + (n/K)EX0,K −EX0,n
≤
n/K
X
i=1
[X(i−1)K,iK −EX0,K] + (n/K)[EX0,K −Kµ] .
(3.21)
If the sequence Xm,n is suitably weakly dependent, then we could expect that the
sum on the right-hand side of (3.21) obeys a central limit theorem, and in particular
frequently takes values on the order of −n1/2Kχ−1/2. The second term is of order nKγ−1.
In particular, if we were to choose K on the order of n1−δ for 0 < δ ≪1, the ﬁrst term
would dominate, giving that X0,n −EX0,n typically has lower tail ﬂuctuations of order
−nχ+δ(1/2−χ). This is a contradiction to the deﬁnition of χ under the assumption that
χ < 1/2.
We are now ready to state some rigorous results following the line of reasoning above;
these results appeared in [18]. In order to avoid precisely describing the weak dependence
axioms needed, we restrict to the case of FPP on Zd (here, the “weak dependence” axiom
amounts to control of the diameter of geodesics – see Theorem 4.5). In the language of
FPP, this result says that given existence of the ﬂuctuation exponent χ, the non-random
ﬂuctuations are at least of the same order as the standard deviation of the passage time.
In the FPP results below, we consider exponents deﬁned as in Deﬁnitions 3.16 and
3.17 for some arbitrary ﬁxed direction x (that is, γ = lim sup[ET(0, nx) −nµ(x)]/ log n,
etc.).
Theorem 3.19. Consider the case of FPP on Zd for d ≥2. Assume Eτ 2+δ
e
< ∞for
some δ > 0.
1. If χ := χ2 = χ2+δ, then
(
γ ≥χ
if χ̸ = 1/2
γ ≥χ −1
2(2β−1 −1)
if Var T(0, nx) = O

n
(log n)β

for some β > 0 .
2. If χ := χ2 = χ2+δ, then if χ < 1/2, γ ≥χ.
Remark 3.20. The result above should be compared with the classical example of
a sum of i.i.d. random variables with ﬁnite mean and ﬁnite 2 + δ moment. In this case,
χ = 1/2 and γ = γ = −∞and these exponents do not satisfy the conclusions of the
theorem above. Note that this example is equivalent to FPP on Z with the appropriate
moment conditions.
Another theorem of [18] gives the bound γ ≥−1/2 under minimal assumptions. This
result should be compared to [115, Theorem 1], where it is shown that γ ≥−1.
54

<!-- page 55 -->

Theorem 3.21. Consider FPP on Zd for d ≥2. Assume F(0) < pc(d), that the distri-
bution of τe is not concentrated at a point, and that
Eeατe < ∞for some α > 0 .
(3.22)
One has the bound γ ≥−1/2: for any nonzero x ∈Zd and ϵ > 0,
ET(0, nx) −nµ(x) ≥n−1
2 −ϵ for inﬁnitely many n .
One may be interested in the sharpness of the result γ ≥χ. In the case of FPP on
Zd (d ≥2), it is reasonable to expect the result is sharp – that is, γ = χ – as explained
in the following remark.
Remark 3.22. Alexander [6] has remarked that if ˆχ is any number such that for
some a > 0,
P

|T(0, y) −ET(0, y)| ≥λ∥y∥ˆχ
1

≤e−aλ for all λ ≥0, y ∈Zd ,
then γ ≤ˆχ. This follows from an extension of his method for controlling nonrandom
ﬂuctuations, discussed in Section 3.5.3 below. Note that if this exponential inequality
holds for some ˆχ, then χp ≤ˆχ for all p > 0. Combining these observations with Theorem
3.19 above, if ˆχ can be taken to be χ := ˆχ = χ2, then
γ := γ = γ = χ
when χ < 1/2 and γ = χ under the assumption Var T(0, nx) = O(n/(log n)β) for every
β > 0.
Question 17. Show that χ2 = χ2+δ for some δ > 0.
A solution for the question above has the following consequence that improves The-
orem 3.21.
Remark 3.23. Combining Theorem 3.4 with Theorem 3.19, if χ2 = χ2+δ for some
δ > 0 and F(0) < pc(d) holds, then γ ≥0.
We now return to the general setting of subadditive ergodic sequences (Xm,n). Our
goal is to ﬁnd extensions to the axioms of Theorem 2.2 which could guarantee upper
bounds on the exponent γ, complementing the lower bounds given above.
Note that by subadditivity, for any n,
EX0,2n ≤EX0,n + EXn,2n = 2EX0,n.
The magnitude of non-random ﬂuctuations is related to the degree to which the above
inequality is strict. Let us abbreviate EX0,n =: H(n). As a ﬁrst step, assume the above
inequality were actually an equality; then
µ = lim
k H(2kn)/(2kn) = lim
k H(n)/n = H(n)/n ,
55

<!-- page 56 -->

in which case we would have H(n) = nµ.
In fact, we can give upper bounds on H(n) −nµ assuming just some weaker quan-
titative control on H(2n) −2H(n). To illustrate this, assume that for some function
j : N →[0, ∞), one has
H(2n) ≥2H(n) −j(n).
(3.23)
We then have the following lemma from [105].
Lemma 3.24. Suppose that j satisﬁes Ψ := lim supn j(2n)/j(n) < 2, j(n)/n →0, and
(3.23). Then for all large n, one has H(n) ≤nµ + cj(n).
Proof. For c > 1/(2 −Ψ), set a(n) = H(n) −cj(n). Then for large n,
a(2n) = H(2n) −cj(2n) ≥2H(n) −2cj(n)
 1
2c + j(2n)
2j(n)

≥2a(n).
So for m ≥1,
a(2mn)
2mn
≥a(n)
n .
Taking m →∞, we obtain
a(n)
n
≤µ,
or
H(n) ≤nµ + cj(n).
We would like to apply the above lemma to bound non-random ﬂuctuations in FPP,
taking H(n) = h(nx) (recall h(y) = ET(0, y)) and µ = µ(x). By the last lemma, if
we could show that h(nx) is close to being a linear function of n, we can upper bound
h(nx) −nµ(x). This approach is most successful in the axis direction e1 for symmetry
reasons, and a version of this method is the basis for several such bounds (see Zhang [173],
Rhee [140], and related reﬂection trick in Alexander [5]). One typically tries to use the
fact that the passage time along geodesics is additive. Deﬁne the events
A = {T(0, ne1) −h(ne1) ≥−j(n)/3} ∩{T(ne1, 2ne1) −h(ne1) ≥−j(n)/3}
B = {T(0, 2ne1) ≤2h(ne1) + j(n)/3}
C = {a geodesic from 0 to 2ne1 passes through ne1}.
Then if A ∩B ∩C has positive probability, on this event, one has
T(0, 2ne1) = T(0, ne1) + T(ne1, 2ne1)
and therefore
h(2ne1) ≥h(ne1) −j(n).
56

<!-- page 57 -->

It is reasonable that C has probability at least some power in n given the heuristics
for geodesic wandering, which would – along with strong bounds on the probabilities of A
and B, provided by concentration results on the passage time – complete the argument.
Unfortunately, not enough is known about geodesic wandering to make this argument
work in general directions. In practice, showing this positive probability statement has
been done using a reﬂection trick, and unfortunately this typically does not work in
non-axis directions.
3.5.2
Asymptotic geodesicity and rate of convergence of the expected ball
Since the axioms of Lemma 3.24 are not easy to verify in FPP, one is led to generalize from
subadditive sequences to some broader framework. Recalling that h(x) := ET(0, x), it is
obvious that h has more than the structure of a subadditive sequence: it is subadditive as
a function of Zd. Since the time constant and h itself depend on the entire d-dimensional
structure of the lattice, it is reasonable to hope that one could progress by considering
the geometry induced on the lattice by h.
We can think of h as a metric given by h(x, y) = h(y −x), and although the metric
space (Zd, h) does not have geodesics, we would like to show that it has approximate
midpoints. Recall that to implement the approach outlined following Lemma 3.24, we
would like to show that ne1 is nearly a midpoint between 0 and 2ne1. A realization of
Tessera is that we do not have to choose ne1: it actually suﬃces to ﬁnd skeletons of paths
0 = x0, x1, . . . , xk = 2ne1 which are approximately geodesics. That is, we would like to
give an upper bound for
h(2ne1) −
k−1
X
i=0
h(xi+1 −xi),
for some suitably chosen skeleton.
Deﬁnition 3.25. Let N
:
R+
→
R+ be an increasing function such that
limα→∞N(α) = ∞. A metric space (X, d) is called SAG(N) (strongly asymptotically
geodesic) if there exists α0 ≥0 such that for all integers m ≥1, and for all x, y ∈X
such that d(x, y)/m ≥α0, there exists a sequence x = x0, . . . , xm = y satisfying, for all
0 ≤i ≤m −1,
d(x, y)
m

1 −
1
N(α)

≤d(xi, xi+1) ≤d(x, y)
m

1 +
1
N(α)

,
where α = d(x, y)/m; and for all x ∈X and large enough r,
B

x,

1 +
1
N(r)r

⊂[B(x, r)]
6r
N(r).
Here, B(x, r) is the r-neighborhood of x and for any T > 0 and set S ⊂X, [S]T is the
T neighborhood of S.
57

<!-- page 58 -->

We brieﬂy try to motivate the above deﬁnition in the setting of (Zd, h). Although
h is not geodesic, the passage time metric T is. Given suitably strong concentration
results, one can ﬁx a realization of edge weights and hope to approximate the passage
times between vertices in a geodesic by the expected passage time between them, up to
a smaller order correction. Breaking this geodesic into a skeleton of m roughly equally
spaced vertices would then furnish the sequence (xi) in Deﬁnition 3.25.
The second
condition of Deﬁnition 3.25 guarantees that there are no large “jumps” in B(x, r) as one
increases r, which is also natural to expect in FPP.
The SAG condition guarantees that the metric d obeys a shape theorem with an
explicit rate of convergence estimate, which we now describe. We say that an increasing
function φ : R+ →R+ is sublinearly doubling if there exists a function η : R+ →R+ such
that limλ→∞η(λ)/λ = 0 and for all λ > 0, φ(λr) ≤η(λ)φ(r). The following theorem
(from [154]) gives an equivalence between asymptotic geodesicity in the form of SAG(N)
and a shape theorem with error term on the order N −1:
Theorem 3.26. Let δ be some translation-invariant metric on Zd and φ : R+ →[1, ∞)
be an increasing, sublinearly doubling function. The following are equivalent.
1. There exists c > 0 such that δ is SAG(N) with
N(α) ≥cφ(α),
for α large enough.
2. There exists a norm ∥· ∥on Rd and C > 0 such that for all large enough n,
B∥·∥

0, n −Cn
φ(n)

∩Zd ⊂Bδ(0, n) ⊂B∥·∥

0, n + Cn
φ(n)

.
A version of SAG with exponent 1/2 was shown for FPP in [154]. Namely, if the edge
weights have ﬁnite exponential moments and F(0) < pc(d) then h satisﬁes SAG(N) for
N(α) = c(α/(log α))1/2,
for some c > 0. Therefore there exist C > 0 and n0 ≥1 such that for any n ≥n0,
(n −C(n log n)1/2)Bν ∩Zd ⊆¯B(0, n) ⊆(n + C(n log n)1/2)Bν,
(3.24)
where ¯B(0, n) := {x ∈Zd : h(x) ≤n}.
Let us conclude this subsection by giving some account of why Theorem 3.26 should
hold; we will freely make simplifying assumptions to avoid getting bogged down. We
note that our simipliﬁcations obscure some of the technical triumph of [154] and direct
the reader there for a full account. We work in the metric δ. Let m and n be arbitrary
58

<!-- page 59 -->

positive integers, and recall that the convergence in Theorem 3.26 is convergence under
the Hausdorﬀmetric dH, where for two compact sets A and B
dH(A, B) = inf{r > 0 : A ⊆[B]r, B ⊆[A]r} .
We will try to show that dH(n−1B(0, n), (mn)−1B(0, mn)) is small; taking m to inﬁnity
gives a corresponding bound for dH(n−1B(0, n), B∥·∥(0, n).
A natural way to use the SAG assumption is to try to build the ball B(0, mn) out
of m copies of the ball B(0, n). This is because, given x at distance mn from 0, we can
ﬁnd a skeleton sequence of vertices (xi) as in Deﬁnition 3.25 each at distance ∼n from
each other. This reasoning may be made rigorous as in the following lemma, adapted
from [154, Lemma 5.1]:
Lemma 3.27. There exists a C such that, for all m and n positive integers,
dH
 (mn)−1B(0, n)m, (mn)−1B(0, mn)

≤C/φ(n) .
In the above, AB = {a + b : a ∈A, b ∈B} and A2 = AA, etc. Using the fact that
dH is a metric, we see
dH(n−1B(0, n), (mn)−1B(0, mn)) ≤dH(n−1B(0, n), (mn)−1B(0, n)m)
+dH((mn)−1B(0, n)m, (mn)−1B(0, mn)) ,
(3.25)
and the second term on the right-hand side is small by Lemma 3.27.
In other words, we see that B(0, mn) is approximately equal to m copies of B(0, n)
strung together, and it remains to show that mB(0, n) is also approximately equal to
this set. Note that if B(0, n) were a convex subset of Rn, then B(0, n)m = mB(0, n),
and we would be done. So the only way this line of attack could fail would be if B(0, n)
were irregular enough to be signiﬁcantly diﬀerent from its convex hull.
This complication is dealt with in [154] by inductively controlling the Hausdorﬀdis-
tance between balls B(0, n) and their convex hulls bB(0, n). A crucial fact used is that
in Rd, Carath´eodory’s Theorem for convex hulls guarantees that dH(B(0, n)m, bB(0, n)m)
can be bounded by a uniform constant multiple of dH(B(0, n), bB(0, n)).
3.5.3
Alexander’s method
The ﬁrst upper bounds on non-random ﬂuctuations go back to the work of Kesten [115,
Theorem 1], who established the bound
ET(0, ne1) ≤nµ(e1) + Cn5/6 log n.
The right side has since been improved to nµ(e1) + C(n log n)1/2; see (3.24) and (3.31).
An extraordinary contribution was the theory developed by Alexander in his famous
59

<!-- page 60 -->

paper [4]. It concerns the rate of convergence of a (deterministic) nonnegative function
h deﬁned on the lattice Zd to its limit g. That is, we assume
h(x + y) ≤h(x) + h(y),
x, y ∈Zd,
and set
g(x) := lim
n h(nx)/n.
Clearly, g(x) ≤h(x) for all x ∈Zd. Alexander’s theory has a very geometric ﬂavor and
explores the fact that the function h is deﬁned on the whole space Zd; h is not only a
(one-dimensional) subadditive sequence. In his point of view, one way to guarantee an
upper bound on the rate of convergence of h(nx)/n is to request h to have “good skeleton
paths” towards nx in Zd. It turns out that this condition is possible to verify in FPP
and other models and suitable to bound the rate of convergence.
To describe his approach, let 0 = v0, v1, v2, . . . , vk = x be points in Zd. We will soon
think of the v′
is as the skeleton of a path from 0 to x. From subadditivity, we have
h(x) = h

k
X
i=1
(vi −vi−1)

≤
k
X
i=1
h

vi −vi−1

.
Imagine for a moment that we can choose these vi’s in such a way that
h(vi −vi−1) = g(vi −vi−1) + o(|x|a),
(3.26)
for some a < 1. Well, this would lead to
h(x) −g(x) ≤
k
X
i=1
 g(vi −vi−1)

−g(x) + ko(|x|a).
(3.27)
Unfortunately, (3.27) is not enough for our purpose as we still need to control the diﬀer-
ence Pk
i=1
 g(vi −vi−1)

−g(x). Moreover, in a ﬁrst look, it seems hard to actually ﬁnd
vi’s satisfying (3.26).
To overcome these problems, we will do two modiﬁcations. First, instead of request-
ing (3.26), we will replace g by the projection gx of g on the line that goes through 0 and
x. This projection is linear and satisﬁes gx(nx) = g(nx) for any n ∈N. In Alexander’s
description, the value gx(y) is the amount of progress toward x made by a vector incre-
ment of y. The diﬀerence h(y) −gx(y) should be thought as the error or the ineﬃciency
associated with such an increment. An increment y will be called good if y satisﬁes
h(y) = gx(y) + C|x|a for a given C > 0.
This nice change from g to gx helps us to get rid of the sum in the right side of (3.27).
However, we still need to ﬁnd the skeleton vi’s of good increments. The main idea and
contribution of [4] is to come up with an assumption on h that allow us to roughly do so.
It will suﬃce to verify the weaker condition that x (divided by some constant larger than
60

<!-- page 61 -->

one) is in the convex hull of good increments. This weaker condition can be veriﬁed in
FPP and many other models. We will make this statement precise in the next deﬁnition.
Let Φ denote the set of all positive nondecreasing functions on (1, ∞) and ﬁx a < 1.
For φ ∈Φ, deﬁne Qx(a, φ) be the set of good increments with exponent a, correction
factor φ in direction x, that is,
Qx(a, φ) :=

y ∈Zd : h(y) ≤gx(y) + |x|aφ(|x|), gx(y) ≤g(x)

.
The condition gx(y) ≤g(x) is to avoid “overshooting” the increment. For a collection of
points in Rd, let Co(·) denote the convex hull.
Deﬁnition 3.28 (CHAP). The function h satisﬁes the convex hull approximation
property (or CHAP) with exponent a and correction φ if there exists L > 1 such that
x/α ∈Co(Qx(a, φ)) for some α ∈[1, L],
for all x ∈Qd with |x| suﬃciently large.
We will discuss in a second how one would verify that a subadditive function h satisﬁes
CHAP and why it is natural to work with Deﬁnition 3.28. The next deﬁnition is exactly
the control on the rate of convergence of h.
Deﬁnition 3.29 (GAP). For a > 0 and φ ∈Φ we say that h satisﬁes the general
approximation property (GAP) with exponent a and correction factor φ if there exist
M > 1 and C > 0 such that for all x ∈Zd and |x| ≥M,
g(x) ≤h(x) ≤g(x) + C|x|aφ(|x|).
Alexander’s main result is the following.
Theorem 3.30 (Theorem 1.8 [4]). Suppose h is a nonnegative subadditive function on
Zd with h(x) ≤r|x| for all x and some r > 0. If h satisﬁes CHAP with exponent a > 0
and correction φ, then it satisﬁes GAP with exponent a and correction φ.
We will now sketch the proof of Theorem 3.30. The following ﬁxed point argument
is the core of the proof.
Step 1. Suppose h satisﬁes CHAP(a, φ) and GAP(b, φ) with b ∈(a, 1]. Then h
satisﬁes GAP(b′, φ) with b′ = b/(1 + b −a) < b.
Sketch of the proof of Step 1. Suppose |x| is large. Take q ∈Q also large but smaller
than |x|. Applying CHAP(a, φ) to x/q we can express x/(αq) as a convex combination
of d + 1 good increments yq,i, 1 ≤i ≤d + 1. Now, use the increments yq,i’s to decompose
x =
d+1
X
i=1
αq,iyq,i + x(r)
61

<!-- page 62 -->

with αq,i ∈N ∪{0} and some remainder x(r). By subadditivity and the fact that g(x) =
gx(x), we have
h(x) −g(x) ≤
d+1
X
i=1
αq,ih(yq,i) + h(x(r)) −gx(x −x(r)) + gx(x(r))
≤
d+1
X
i=1
αq,i
 h(yq,i) −gx(yq,i)

+ h(x(r)) −gx(x(r)).
(3.28)
Our choice of the yq,i’s implies that the sum on the right side of (3.28) is bounded by
q1−a|x|aφ(|x|),
while the remainder term h(x(r)) −gx(x(r)) is controlled using GAP(b, φ). Optimizing
over q gives the desired result. (The optimal q turns out to be of the order of a small
power of |x|.) The details can be found in [4, Proposition 2.1].
Step 2. By hypothesis, h satisﬁes GAP(1, φ) with φ = 1. If a ≥1, there is nothing
to prove. If a < 1, we can iterate Step 1 to obtain GAP(a, φ) as the ﬁxed point of the
function x →x/(1 + x −a) is at a.
Let us now end this subsection by saying a few words on why it is natural to expect
CHAP and how one proves CHAP(a, φ) in FPP. The ﬁrst ingredient is the following
geometric lemma, that deﬁnes the skeleton {vi} discussed before. The proof is short, so
we include it here.
Lemma 3.31 (Lemma 1.6, [4]). Let h be a nonnegative subadditive function in Zd and
c > 1. Suppose that for each x ∈Qd, there exist n ≥1, a lattice path Γ from 0 to nx and
a sequence of sites 0 = v0, v1, . . . , vm = nx in Γ such that m ≤cn and vi−vi−1 ∈Qx(a, φ)
for all 1 ≤i ≤m. Then h satisﬁes CHAP(a, φ).
Proof. We need to ﬁnd α ≥1 so that
x/α ∈Co(Qx(a, φ)).
Let a(n, y) be the number of indices such that vi −vi−1 = y. Then, by hypothesis,
nx =
X
y∈Qx(a,φ)
a(n, y)y
and applying gx to the equation above
ng(x) =
X
y∈Qx(a,φ)
a(n, y)gx(y) ≤g(x)
X
y∈Qx(a,φ)
a(n, y)
so
n ≤
X
y∈Qx(a,φ)
a(n, y).
62

<!-- page 63 -->

On the other hand, by hypothesis,
n−1
X
y∈Qx(a,φ)
a(n, y) ≤c.
Taking α = n−1 P
y∈Qx(a,φ) a(n, y) and L = c gives us CHAP.
Now look at the m + 1 sites in the Lemma above. It is easy to ﬁnd inductively a
sequence of marked sites for any path Γ from 0 to nx. One starts with v0 = 0 and chooses
zi as the ﬁrst site in Γ after vi−1 such that zi −vi−1 /∈Q(a, φ). Now take vi as the last
vertex visited by Γ. We call the sequence of marked sites, obtained from a self-avoiding
path in this way, the Qx(a, φ) skeleton of Γ. The diﬃculty is to control the number of
vertices in a Qx skeleton. This is the role of the next proposition, which combined with
Lemma 3.31 ends the proof of CHAP(1/2, log |x|).
Proposition 3.32. Assume F(0) < pc(d) and Eeλτe < ∞for some λ > 0. Then there
exists C > 0 so that if |x| > C then for suﬃciently large n there exists a lattice path from
0 to nx with Qx(1/2, log |x|)-skeleton of 2n + 1 or fewer vertices.
We will take one geodesic from 0 to nx as our lattice path. To prove Proposition 3.32,
Alexander takes use of one very useful estimate of Kesten [113, (4.13)] and his skeleton
construction. The argument goes as follows where here the constant C will change from
line to line.
Consider a sequence of points vi ∈Zd so that the increments vi −vi−1 are the Qx
skeleton of some lattice path Γ. Let
Yi = ET(vi−1, vi) −T(vi−1, vi).
The hypothesis Eeλτe < ∞is suﬃcient (see Theorem 3.13) to ﬁnd a constant C > 0 so
that
P

|Yi| > λ|vi −vi−1|1/2

≤C exp

−Cλ

.
In particular, this combined with the fact that |vi −vi−1| ≤2d|x| implies that the
existence of a β > 0 so that E exp(βYi/|x|1/2)) < ∞. Now, imagine for a second that the
set of edges used by paths from vi to vi+1 were disjoint and deterministic. Then Yi would
independent. Moreover, taking exponentials and using Markov’s inequality we obtain:
P
 m−1
X
i=0
ET(vi, vi+1) −T(0, vm) ≥C|x|1/2 log |x|

≤P
 m−1
X
i=0
Yi ≥C|x|1/2 log |x|

≤C exp(−Cm log |x|).
(3.29)
The fact that we assumed the use of disjoint edges and independence (illegally) can
be handled by taking independent copies Y ′
i of Yi and using two extensions of the BK
63

<!-- page 64 -->

inequality (this is done in Kesten [113, (4.13)]), see also Theorem 3.41 below. So morally
(3.29) is correct. On the other hand, by geometric considerations, there are at most
C|x|md ways to choose vi, i = 1, . . . , m, with |vi −vi−1| ≤2d|x|. Thus, union bound
implies,
P
 m−1
X
i=0
ET(vi, vi+1) −T(0, vm) ≥C|x|1/2 log |x| for some m ≥1 and some Qx
skeleton with m vertices

≤C exp(−C log |x|).
(3.30)
Now take a Qx-skeleton with m vertices of a geodesic from 0 to nx. By deﬁnition, it
follows that n ≤m. The estimate (3.30) above implies that we can choose n large enough
so that with probability going to one
m−1
X
i=0
ET(vi, vi+1) −T(0, nx) ≤Cm|x|1/2 log |x|.
Hence, by the convergence to the time constant,
m−1
X
i=0
ET(vi, vi+1) ≤nµ(x) + Cm|x|1/2 log |x|.
Now the proof ends by arguing (see [4, (3.16)]) that the large majority of the terms in the
sum of the left-side are good increments, that is, they are of order Cµ(x), which leads
(when all constants are put in order) to m ≤2n for n large and |x| ≥C.
We ﬁnish this section with a summary of the current state of the art on (a) non-
random ﬂuctuation upper bounds and (b) convergence rate to the limit shape. For (a),
Alexander’s methods were used in the low moment case in [61] along with a concentration
inequality for the lower tail of T(0, x) to obtain the following.
Theorem 3.33. Assume F(0) < pc and EY 2 < ∞, where Y is the minimum of d i.i.d.
copies of τe. There exists C = C(F, d), such that for all x ∈Zd with |x| > 1,
µ(x) ≤ET(0, x) ≤µ(x) + C(|x| log |x|)1/2 .
(3.31)
For (b) above, the following version of convergence rate bounds under low moments
also appears in [61].
Theorem 3.34. Assume P(τe = 0) < pc. If EY 2 < ∞, where Y is the minimum of d
i.i.d. copies of τe, then there exists C > 0 such that almost surely,
B(t) ⊂(t + C(t log t)1/2)Bν for all large t.
If in addition Eτ α
e < ∞for some α > 1 + 1/d, then almost surely
(t −C(log t)4t1/2)Bν ⊂B(t) for all large t.
The inner bound given here can be improved under ﬁnite exponential moments, as
mentioned in (3.24).
64

<!-- page 65 -->

3.6
Large deviations
In this section we discuss some large deviation bounds for FPP. For given ϵ > 0, our ﬁrst
goal is to estimate the following probabilities
pu
n(ϵ) : = P

T(0, ne1) > n(µ(e1) + ϵ)

pℓ
n(ϵ) : = P

T(0, ne1) < n(µ(e1) −ϵ)

,
(3.32)
for n large. If we assume (2.1), Theorem 2.1 implies that both pu
n(ϵ) and pℓ
n(ϵ) go to
zero as n →∞. Stronger assumptions will allow us to derive exponentially small upper
bounds for these probabilities.
Typically, these bounds will exhibit a strong asymmetry: the behavior of pu
n is radi-
cally diﬀerent from pℓ
n. Indeed, If we want to increase the passage time T(0, ne1), we may
need to increase the passage time of every path from 0 to ne1. If the passage times are
bounded, this roughly amounts to increasing almost every edge in a box with diameter
of order n. Thus one should expect:
pu
n(ϵ) ≈exp(−ndIu(ϵ)),
(3.33)
where an ≈bn if log an/ log bn →1.
On the other hand, to decrease the passage time between two points, it suﬃces to
force the passage time of a single path to be very low. Thus,
pl
n(ϵ) ≈exp(−nIℓ(ϵ)),
(3.34)
as we should roughly decrease the passage time of order n edges.
The reader may have seen this asymmetry in large deviations in other areas of prob-
ability. In random matrix theory, for instance, the smallest eigenvalue of a n × n GOE
matrix satisﬁes a large deviation principle of speed n if forced out of the bulk; on the
other hand, the probability that the smallest eigenvalue is above the bottom edge of the
equilibrium measure is of order exp(−cn2).
Before stating the results, let’s make a remark.
The boundedness assumption is
important for (3.33). For example, suppose that τe is exponentially distributed with
mean 1. We know that if all 2d edges connected to the origin take values larger than
n(µ + ϵ), then T(0, ne1) ≥n(µ + ϵ). Therefore,
pu
n(ϵ) ≥P

τe ≥n(µ(e1) + ϵ)
2d
= exp

−2nd(µ(e1) + ϵ)

,
so (3.34) does not hold.
We now state the main theorems of this section. We start with the large deviation
lower bounds and we will assume that F(0) < pc(d) so that µ(e1) > 0. Let
β = sup

x : F(µ(e1) −x) > 0

.
65

<!-- page 66 -->

Theorem 3.35 (Theorem 5.2 [113]). For any ϵ > 0, there exist a positive constant A(ϵ)
and Iℓ(ϵ) ∈(0, ∞] such that for any n ≥1
pℓ
n(ϵ) ≤A(ϵ) exp(−nIℓ(ϵ)).
Furthermore,
lim
n
1
n log pℓ
n(ϵ) = −Iℓ(ϵ),
(3.35)
where 0 < Iℓ(ϵ) < ∞for any 0 < ϵ < β and Iℓ(ϵ) = ∞for β ≤ϵ. Setting Iℓ(ϵ) = 0 for
ϵ ≤0, the function Iℓ: (−∞, β) →R is convex on its domain and strictly increasing on
[0, β).
Remark 3.36. Theorem 3.35 requires no assumptions on the moments of the edge
passage times. In particular, it holds even when (2.6) is not satisﬁed and µ(e1) is not
deﬁned. In this case, µ(e1) in the deﬁnition of pℓ
n(ϵ) is replaced by the shell time constant
ˆµ(e1) := limn ˆT(0, ne1)/n from Section 2.4.
Remark 3.37. The function Iℓ(ϵ) should depend on d and on the underlying distri-
bution of passage times F.
We now turn to upper large deviation bounds. As remarked before, it does not suﬃce
to have ﬁnite exponential moments to derive (3.33). The most natural assumption is to
require the passage times to be bounded.
This was the assumption used by Kesten
in [113].
Theorem 3.38 (Theorem 5.9 [113]). Assume that
R
eλx dν < ∞for some λ > 0. Then
for each ϵ > 0 there exist constants A1 = A1(ϵ) and B1 = B1(ϵ) > 0 such that for all
n ≥1,
P
 T(0, ne1) > n(µ + ϵ)

≤A1e−B1n.
If τe is bounded with probability one, then for all ϵ > 0 there exist constants A2 = A2(ϵ)
and B2 = B2(ϵ) > 0 such that for all n ≥1
P
 T(0, ne1) > n(µ + ϵ)

≤A2e−B2nd .
(3.36)
Note that Kesten’s bound (3.36) for bounded passage times is not of the form (3.33).
Proving (3.33) is still an open question.
Question 18. Assume that the passage times are bounded and not concentrated at
a single point. Show that there exists a convex function Iu(ϵ) such that the following
limit holds:
lim
n
1
nd log pu
n(ϵ) = Iu(ϵ).
Theorem (3.38) was extended by Cranston-Gauthier-Mountford. In [55], the authors
found necessary and suﬃcient conditions on the distribution of edge weights so that
lim supn n−d log pu
n(ϵ) < 0. Their result remains the most recent contribution concerning
large deviation bounds of passage times.
66

<!-- page 67 -->

Theorem 3.39 (Theorem 1.3 in [55]). Assume that for positive M0 < ∞that there is a
positive increasing function f so that
log P
 τe > x

= −xdf(x),
x > M0.
Then for every suﬃciently small ϵ > 0,
lim sup
n
1
nd log P

T(0, ne1) ≥n(µ + ϵ)

< 0
if and only if
∞
X
n=1
1
f(2n)1/(d−1) < ∞.
If instead of considering the sequence of point to point passage times T(0, ne1), we
take the box-to-plane passage times
G(0, n) = inf

T(Γ) : Γ a path from {0} × [−n, n]d−1 to {n} × Zd−1

,
then [90]
lim
n
G(0, n)
n
= lim
n
T(0, ne1)
n
= µ(e1).
In this case, assuming F(0) < pc, Chow and Zhang [50] showed the existence of convex
functions c(ϵ) and c′(ϵ) so that
lim
n
1
nd log P

G(0, n) ≥n(µ + ϵ)

= c(ϵ) and lim
n
1
n log P

G(0, n) ≤n(µ −ϵ)

= c′(ϵ).
We end this section with a brief remark on the proofs of Theorems 3.35 and 3.38.
The most important step in the proof of these theorems is Proposition 3.40 below. It
estimates the distribution of T(0, ne1) in terms of sums of independent passage times
across large strips. These strips are still small with respect to n.
Let Hk be the hyperplane {x = (x1, . . . , xd) : x1 = k} ⊂Rd. For a path Γ we write
Hk < Γ < Hl
if all edges of Γ lie strictly between Hk and Hl. For N ≥1, we deﬁne the d−1 dimensional
box
HN
k :=

x = (x1, x2, . . . , xd) : x1 = k, 0 ≤xi ≤N, i ≥2

.
Note that HN
k ⊂Hk. For M ≥1 we also deﬁne the following cylinder passage time:
sN
0,n := inf

T(Γ) : Γ a path from some point in HN
0 to Hn for which H0 < Γ < Hn

.
67

<!-- page 68 -->

Proposition 3.40 (Proposition 5.23 in [113]). Let Xi(M, N) be independent random
variables, each with distribution of sN
0,M. Then for any n ≥M ≥N ≥1, t ≥0 one has
P

T(0, ne1) < t

≤
∞
X
ℓ=n/(M+N)−1
(2d(15M
N )d)ℓP

X1(M, N) + . . . + Xℓ(M, N) < t

.
(3.37)
The proof of Proposition 3.40 relies on the van den Berg-Kesten (BK) inequality.
In reality, it relies on the lemma below, one of the most useful illustrations of the BK
inequality in FPP. Given sets of vertices of Zd, A(i, j), i ≥1, 1 ≤j ≤n(i), with
A(i, 0) = {0} and A(i, n(i) + 1) = {x}, deﬁne
Πi(0, x, t) :=

Γ : Γ is a path from 0 to x that passes successively
through A(i, 1), . . . A(i, n(i)) and has T(Γ) < t

.
Lemma 3.41 ((4.13) in [113]). For any t > 0,
P
 [
i≥1
Πi(0, x, t)
̸
= ∅

≤
X
i≥1
P
 n(i)
X
j=0
T ′(i, j) < t

,
(3.38)
where T ′(i, j) are all independent with T ′(i, j)
d= T(A(i, j), A(i, j + 1)).
3.7
Cases where Gaussian ﬂuctuations appear
In this section, we discuss a few cases where Gaussian ﬂuctuations are known or expected
to appear. The ﬁrst example is Critical FPP, where the edge weights are taken with
F(0) = pc, exactly at the percolation threshold. The second case, is when we look at
FPP in a thin cylinder, where geodesics want to be thinner than usual.
3.7.1
Critical ﬁrst-passage percolation
From Kesten’s work, we know that
p := P(τe = 0) ≥pc ⇔g(x) = 0 for all x,
where pc is the critical threshold for bond percolation on Zd. Most theorems in FPP
have the assumption that the above probability is strictly less than pc for this reason.
Namely, if p > pc, then one can show that as x →∞, the family of passage times T(0, x)
is stochastically bounded above by a random variable Z that does not depend on x. The
reason is that to travel from 0 to x, one simply needs to go from 0 to the inﬁnite cluster
68

<!-- page 69 -->

of zero-weight edges, then from here to x. On the other hand, if p = pc, then very little
is known. This is due to the fact that one of the outstanding problems in probability
theory is to show that there is no inﬁnite cluster of open edges in bond percolation on
Zd at the critical point. Therefore, we do not know if vertices can take advantage of an
inﬁnite cluster of zero-weight edges when p = pc, and so we do not know even the growth
rate of T(0, x) as x →∞.
The ﬁrst result on the critical case came from Chayes [45], who showed in the Bernoulli
case:
Theorem 3.42 (Chayes [45]). Let d ≥3 and τe satisfy P(τe = 0) = pc = 1 −P(τe = 1).
If ϵ > 0, then
lim
n
T(0, ne1)
nϵ
= 0 almost surely.
See [113, Remark 3] also, where Kesten claims that Chayes’s argument can be extended
to
T(0, ne1) ≤exp(C
p
log n)
for large n almost surely.
Although this is the extent of knowledge on the critical case in general, much more
progress has been made in two dimensions. In this case, Zhang [173] showed that diﬀerent
critical edge-weight distributions can display completely diﬀerent asymptotic behavior:
deﬁning the passage time to inﬁnity as
ρ = lim
n T(0, ∂B(n)),
(3.39)
where B(n) = [−n, n]2, there are some critical distributions for which ρ = ∞almost
surely, and there are some for which ρ < ∞almost surely. He did this by comparing
critical FPP to a model called the CCD incipient inﬁnite cluster [47]. This is an inhomo-
geneous bond percolation model which is tuned to resemble a large percolation cluster at
criticality. A complete characterization was given recently in [62], using methods related
to 2d invasion percolation. For its statement, let F be the distribution function of τe and
set
F −1(y) = inf{x : F(x) ≥t} for t ≥0.
Theorem 3.43 (Damron-Lam-Wang [62]). For d = 2, one has ρ < ∞almost surely if
and only if P∞
n=2 F −1(pc + 2−n) < ∞.
The theorem implies that if F is concentrated enough near 0, then the passage time
to inﬁnity can be ﬁnite. Otherwise, if there are too many large edges, the passage time
to inﬁnity will be inﬁnite. For example, ρ < ∞when F decays to pc like a polynomial
near 0, but can be inﬁnity if F decays to pc like an exponential (with large enough rate)
near 0.
69

<!-- page 70 -->

In the case that ρ = ∞, one can ask about the rate of growth of T(0, ∂B(n)) as
n →∞.
The ﬁrst result of this type was due to Chayes-Chayes-Durrett [46], who
showed that in the critical Bernoulli case,
ET(0, ∂B(n)) ≍log n.
The general case was given in [62]. Let Y be the minimum of 4 i.i.d. copies of τe.
Theorem 3.44 (Damron-Lam-Wang [62]). For d = 2, suppose that F(0) = pc and
EY α < ∞for some α > 1/4. There exists C > 0 such that
1
C ET(0, ∂B(2n)) ≤
n
X
k=2
F −1(pc + 2−k) ≤CET(0, ∂B(2n)).
If EY α < ∞for some α > 1/2, then there exists C′ > 0 such that
1
C′Var T(0, ∂B(2n)) ≤
n
X
k=2
 F −1(pc + 2−k)
2 ≤C′Var T(0, ∂B(2n)).
For FPP on the sites of a triangular lattice, Yao [167] proved the existence of a
constant µ such that (log n)−1T(0, ∂B(n)) converges to µ almost surely.
Back in Zd, in the case F(0) < pc, establishing limiting laws for T(0, ne1) is far oﬀ,
except in some related exactly solvable models. Kesten and Zhang, however, showed that
in the Bernoulli critical case in two dimensions, one has a Gaussian central limit theorem
on scale √log n. Their theorem applied to all distributions with a gap near zero; that is,
ones which have F(0) = pc and F(a) = pc for some a > 0. This reason for this condition
is precisely the uncertain behavior of ρ in the critical case.
Theorem 3.45 (Kesten-Zhang [117]). Let d = 2. Suppose that for some C > 0, F(0) =
F(C) = 1/2. Also assume that Eτ 4+δ
e
< ∞for some δ > 0. Then there exist positive
constants C1, C2, γn with
C1(log n)1/2 ≤γn ≤C2(log n)1/2
such that for all t ∈R
lim
n→∞P
 γ−1
n

T(0, ne1) −ET(0, ne1)

≤t

= P(Z ≤t)
where Z is a standard Gaussian random variable.
This result was extended to general distributions in [62]. Note that if P
k F −1(pc +
2−k) = ∞but P
k(F −1(pc+2−k))2 < ∞, then the variance of the passage time converges,
whereas the mean diverges. In this case, we do not get a Gaussian limit, but rather
T(0, ∂B(n)) −ET(0, ∂B(n)) converges to some other variable.
70

<!-- page 71 -->

Theorem 3.46 (Damron-Lam-Wang [62]). Suppose that F(0) = pc, and P∞
k=2 F −1(pc +
2−k) = ∞. Suppose further that EY α < ∞for some α > 1/2.
1. If P∞
k=2(F −1(pc + 2−k))2 < ∞, then there is a random variable Z with EZ = 0 and
EZ2 < ∞such that as n →∞,
T(0, ∂B(n)) −ET(0, ∂B(n)) →Z almost surely and in L2.
2. If P∞
k=2(F −1(pc + 2−k))2 = ∞, then as n →∞,
T(0, ∂B(n)) −ET(0, ∂B(n))
p
Var T(0, ∂B(n))
⇒N(0, 1).
In the case of site FPP on the triangular lattice with Bernoulli weights (P(tv = 1) =
1/2 = P(tv = 0)), Yao [168] has recently used the conformal loop ensemble to ﬁnd
the exact values of the constants appearing the ﬁrst-order asymptotics of the mean and
variance of the passage time.
Theorem 1 (Yao [168]). Consider site FPP on the triangular lattice with weights
(tv) satisfying P(tv = 1) = 1/2 = P(tv = 0). The following statements hold:
lim
n
T(0, ∂B(n))
log n
=
1
2
√
3π almost surely,
lim
n
ET(0, ∂B(n))
log n
=
1
2
√
3π,
lim
n
Var T(0, ∂B(n))
log n
=
4
3
√
3π −1
π2.
There is a big diﬀerence between what is known in two dimensions and in three
dimensions. This is due mostly to the lack of knowledge about critical percolation in
three dimensions. A natural question is then:
Question 19. Let d = 3 and P(τe = 0) = pc = 1 −P(τe = 1). Is it true that
ET(0, ∂B(n)) ≍log n?
The two dimensional proof of this fact uses that with uniformly positive probability
in n, every path which crosses B(n) \ B(n/2) must contain at least one edge which does
not have zero weight. This is no longer expected to be true in high dimensions (above,
say, six [3, Theorem 4]), so it is conceivable that for large d, one has T(0, ∂B(n)) ≪log n
even in the Bernoulli case.
71

<!-- page 72 -->

3.7.2
Passage time in small cylinders
A second case where Gaussian ﬂuctuations are present is when we consider passage times
of paths constrained in a thin cylinder. Precisely, for h > 0 let an(h) be the ﬁrst passage
time from 0 to the point ne1 in the graph (Z × [−h, h]d−1) ∩Zd; that is, deﬁne
an(h) := inf{T(γ)|γ is a path from 0 to ne1 in Z × [−h, h]d−1}.
If h is small enough, the constrained geodesic from 0 to ne1 should not coincide with
the one in the full integer lattice. In this case, we have Gaussian ﬂuctuations after proper
centering and scaling.
Theorem 3.47 (Chatterjee-Dey [44]). Suppose that Eτ 2+δ
e
< ∞for some δ > 0. Assume
that h = h(n) satisﬁes h(n) = o(nα) with
α <
1
d + 1 + 2(d −1)/δ.
(3.40)
Then we have for all t ∈R
lim
n→∞P
an(hn) −Ean(hn)
p
Var(an(hn))
≤t

= P(Z ≤t)
(3.41)
where Z is a standard Gaussian random variable.
The theorem above gives rise to the following question, raised in [44].
Question 20. Is (3.40) optimal? More explicitly, let
γF(d) := sup

α : (3.41) holds for hn = nα

.
What is the value of γF(d)?
A partial answer was given by the same paper of Chatterjee and Dey. If one has
ﬁnite exponential moments for the passage time distribution and assumes (a) existence
of both ﬂuctuation exponent χ(d) and transversal exponent ξ(d) (See Section 4.2), (b)
strict convexity of the limiting shape and (c) d = 2 or d = 3 or ξ(d), χ(d) > 0 then
γF(d) = ξ(d).
72

<!-- page 73 -->

4
Geodesics
Recall the deﬁnition of a geodesic from x to y.
Deﬁnition 4.1. A path γ from x to y with T(γ) = T(x, y) is called a geodesic from
x to y.
Soon, we will deal with inﬁnite geodesics, so we will sometimes use the term ﬁnite
geodesic for a geodesic between two points. We will also require geodesics from points
to sets: if A is a subset of Zd, then a geodesic from x to A is a geodesic from x to y,
where y is the vertex of A minimizing T(x, y). If A is a subset of Rd, then we identify
each point of A with the nearest vertex of Zd, and deﬁne ˆA to be the union of all such
vertices; the set of all geodesics from x to A is deﬁned to be the set of all geodesics from
x to ˆA.
4.1
Existence of ﬁnite geodesics, their sizes, and the geodesic
tree
The very ﬁrst step in the study of geodesics is to determine if a geodesic from x to y
exists with probability one. Hammersley and Welsh [98] showed almost sure existence of
ﬁnite geodesics if the distribution of the passage times is bounded above and below away
from zero. They also conjectured that almost sure existence held for all distributions,
with no moment conditions. Smythe and Wierman [148,149] veriﬁed this conjecture in
two dimensions except for distributions with an atom at zero equal to the critical bond
percolation probability. In Z2, the conjecture was proven in full generality by Wierman
and Reh [161], and it is now the state-of-the-art. It was stated in [161, Corollary 1.3] for
geodesics between 0 and ne1, but the same proof applies to any pair of points.
Theorem 4.2 (Wierman-Reh [161]). For any passage time distribution F, With proba-
bility one, there exists a geodesic between any two points of Z2.
The theorem above was extended to connected, inﬁnite subsets of Z2 with connected
inﬁnite complement in [17]. Existence under no assumptions on the distribution of pas-
sage times is currently open in d dimensions if d > 2.
Question 21. Prove that under no assumptions on the distribution F of passage
times, geodesics exist in the d dimensional lattice, d ≥3.
Next we will now explain how, as noted in [113, (9.23)], the existence of geodesics
follows from standard estimates provided F(0) < pc. In the case F(0) > pc with Eτe < ∞,
existence of geodesics was proved in [174]. The case F(0) = pc is still unresolved in general
when d > 2. However, geodesics are known to exist for distributions F satisfying F(0) =
F(δ) = pc for some δ > 0 and in suﬃciently high dimensions. See [171, Theorem 8.1.7].
73

<!-- page 74 -->

Recall the deﬁnition of the passage time to inﬁnity (3.39), from [171]:
ρ = ρ(F) = lim
n [T(0, ∂B(n))] ,
where B(n) = [−n, n]d. By the Komogorov 0-1 law, either ρ = ∞almost surely, or
ρ < ∞almost surely.
Lemma 4.3. One has
ρ = inf{T(γ) : γ is an inﬁnite self-avoiding path from 0}.
Proof. Let γn be a (self-avoiding) geodesic from 0 to ∂B(n), with γ a subsequential limit
of (γn) along some subsequence (γnk). For a ﬁxed n, if ˆγn is the initial portion of γ from
0 to ∂B(n), then for all large k, the geodesic from 0 to ∂B(nk) contains ˆγn as its ﬁrst
segment, so
T(0, ∂B(nk)) ≥T(ˆγn) ≥T(0, ∂B(n)).
Letting k →∞, we obtain
ρ ≥T(ˆγn) ≥T(0, ∂B(n)).
As n →∞, T(ˆγn) →T(γ), so T(γ) = ρ, showing the inequality ≥in the lemma. The
other inequality holds because each inﬁnite self-avoiding path from 0 contains an initial
segment from 0 to ∂B(n).
From this lemma, we can show that F(0) < pc implies that geodesics exist.
Proposition 4.4. The following hold.
1. If ρ = ∞for some conﬁguration (τe), then for all x, y ∈Zd, a geodesic exists from
x to y.
2. If F(0) < pc, then almost surely, ρ = ∞.
Proof. If ρ = ∞in some passage time conﬁguration (τe), pick any x, y, and ﬁx a path γ
between them. Note that since ρ = ∞and
|T(x, ∂B(n)) −T(0, ∂B(n))| ≤T(0, x) < ∞,
one can choose n large enough so that (a) x, y and γ are in B(n) and (b) T(x, ∂B(n)) >
T(γ). If γ′ is any path from x to y that exits B(n), then
T(γ′) ≥T(x, ∂B(n)) > T(γ),
so the inﬁmum in the deﬁnition of T(x, y) is over a ﬁnite set of paths, and there must
be a minimizer.
If F(0) < pc, then choose δ > 0 such that F(δ) < pc. By the deﬁnition of pc, almost
surely there is no inﬁnite self-avoiding path of edges e with τe ≤δ, and in fact each such
path γ must have inﬁnitely many edges with weight at least δ. So by the characterization
in the previous lemma, ρ = ∞.
74

<!-- page 75 -->

The condition ρ = ∞is not necessary for existence of geodesics, as they exist when
F(0) > pc (and therefore ρ < ∞). It is unknown in general dimensions if ρ = ∞for a
given F, but in two dimensions, Theorem 3.43 gives necessary and suﬃcient conditions.
We conclude this section with estimates on the size of ﬁnite geodesics. Let Γ(0, x)
denote the set of geodesics from 0 to x. The diameter of the set Γ(0, x) is controlled by
the following result that requires no moment condition on the distribution of τe.
Theorem 4.5 (Theorem 6.2 [18]). Assume F(0) < pc(d).
Then there exist positive
constants M, C such that
P

diam Γ(0, x) ≥M∥x∥∞

≤exp(−C∥x∥∞),
for all x ∈Zd.
Furthermore, we have the following linear bound on the number of edges in a geodesic.
Theorem 4.6. Assume F(0) < pc(d). Then there exist positive constants M, C such
that
P (m(x) ≥M∥x∥∞) ≤exp
 −C∥x∥1/d
∞

for all x ∈Zd,
where m(x) = max{|σ| : σ a self-avoiding geodesic from 0 to x}.
Theorem 4.6 implies that if F(0) < pc, then for each α > 0, there is Cα such that for
all x̸ = 0,
Em(x)α ≤Cα∥x∥α
∞.
Under the additional condition Eeβτe < ∞for some β > 0, [57, Corollary 1.4(2)] shows
a similar bound for exponential moments: there exists α1 > 0 such that
sup
0̸=x∈Zd
log Eeα1m(x)
∥x∥∞
< ∞.
The proofs of the above two theorems use almost the same method, and since the
ﬁrst has appeared, we give the second.
Proof of Theorem 4.6. For p ∈[0, 1] let Pp be the product measure on Ω= {0, 1}Ed with
marginal Pp(ω(e) = 1) = p, where ω is a typical element of Ω. In a conﬁguration ω we
write x →y if there is a path from x to y with edges e satisfying ω(e) = 1. This gives
a connectivity equivalence relation and the equivalence classes are called open clusters.
It is known that for p > pc, there is almost surely a unique inﬁnite open cluster. Deﬁne
B(n) = {x ∈Zd : ∥x∥∞≤n} and ∂B(n) = {x ∈Zd : ∥x∥∞= n}.
Lemma 4.7. Let An be the event that every path from 0 to ∂B(n) intersects the inﬁnite
open cluster. There exists p0 ∈(pc, 1) such that if p ∈[p0, 1] then for some C1 > 0,
Pp(An) ≥1 −e−C1n for all n.
75

<!-- page 76 -->

Proof. The proof is a slight modiﬁcation of the result of Kesten [113, Theorem 2.24], and
is given as [18, Lemma 6.3].
We now recall the result of Antal-Pisztora and afterward prove the theorem.
Lemma 4.8 (Theorem 1.1 [12]). Let p > pc. Then there exists a constant ρ = ρ(p, d) ∈
[1, ∞) such that
lim sup
∥y∥∞→∞
1
∥y∥∞
log Pp(dI(0, y) > ρ∥y∥∞, 0, y ∈I) < 0,
where I is the inﬁnite open cluster and dI is the intrinsic distance in I.
We choose K > 0 such that p := P(τe ≤K) > p0 and deﬁne a percolation conﬁgura-
tion (ηe) from the weights (τe) by
ηe =
(
1
if τe ≤K
0
if τe > K .
For x ∈Zd and an integer M > 0, we ﬁrst write nx = ⌊∥x∥1/d
∞⌋and estimate
P(m(x) ≥M∥x∥∞) ≤2e−C1nx + P(m(x) ≥M∥x∥∞, Anx, Bnx),
(4.1)
where Anx is written for the event in Lemma 4.7 for the percolation conﬁguration (ηe)
and Bnx is the same event with 0 translated to x. From now on, we take ∥x∥∞≥C′
with C′ chosen so that nx ≥1 and 4nx < ∥x∥∞. On the event on the right we can select
a self-avoiding geodesic P from 0 to x which has length at least M∥x∥∞. Write P1 for
the portion of P from 0 to its ﬁrst intersection of ∂B(nx) and P2 for the portion from
its last intersection of ∂B(x, nx) to x (here B(x, n) is the translate of B(n) centered at
x). We can then choose u, v vertices of P1 and P2 respectively such that ∥u∥∞≤nx and
∥v −x∥∞≤nx and both u, v are in the inﬁnite open cluster of (ηe). By construction,
the portion P3 of P from u to v is a self-avoiding geodesic that has at least
M∥x∥∞−C2nd
x ≥(M −C2)∥x∥∞number of edges.
We now apply Lemma 4.8 to ﬁnd C3, C4 > 0 such that
Pp(dI(0, y) > ρ∥y∥∞, 0, y ∈I) ≤C4e−C3∥y∥∞for all y ∈Zd.
So for u ∈B(nx) and v ∈B(x, nx), as 2∥x∥∞≥∥u −v∥∞≥∥x∥∞/4,
Pp(dI(u, v) > 2ρ∥x∥∞, u, v ∈I) ≤C4e−(C3/4)∥x∥∞
and by a union bound,
Pp(dI(u, v) > 2ρ∥x∥∞for some u ∈B(nx) ∩I, v ∈B(x, nx) ∩I) ≤C5e−C6∥x∥∞
76

<!-- page 77 -->

for some C5, C6 > 0.
On the complement of this event, each u ∈B(nx) ∩I and v ∈B(x, nx) ∩I have
dI(u, v) ≤2ρ∥x∥∞and so T(u, v) ≤2Kρ∥x∥∞. Use this in the right side of (4.1) to
bound it above by
C7e−C8∥x∥1/d
∞+ P

∃u ∈B(nx), v ∈B(x, nx) with m(u, v) at
least (M −C2)∥x∥∞but T(u, v) ≤2Kρ∥x∥∞

.
(4.2)
Here, m(u, v) is the maximal number of edges in a self-avoiding geodesic from u to v.
Last, we appeal to Kesten’s result [113, Proposition 5.8], which states that under the
assumption F(0) < pc, there exist constants a, C9 > 0 such that for all n ≥1,
P (∃self-avoiding γ starting at 0 with #γ ≥n but with T(γ) < an) ≤e−C9n.
By a union bound, for all n ≥1 and x ∈Zd with ∥x∥∞≥C′,
P (∃self-avoiding γ from B(nx) to B(x, nx) with #γ ≥n but T(γ) < an) ≤C10nd
xe−C9n.
If there exist u ∈B(nx) and v ∈B(x, nx) with m(u, v) ≥(M −C2)∥x∥∞, then we may
select a self-avoiding geodesic between u and v that has at least this number of edges.
So ﬁxing any M with 2Kρ/a + C2 < M, the expression in (4.2) is bounded by
C7e−C8∥x∥1/d
∞+ C10nd
xe−C9(M−C2)∥x∥∞,
and this is bounded by C11e−C12∥x∥1/d
∞for some C11, C12 > 0. Increasing M if necessary,
we obtain e−C13∥x∥1/d
∞, completing the proof.
In several occasions we will be interested in the case where Γ(0, x) is a single path
from 0 to x.
Deﬁnition 4.9. A metric space is called uniquely geodesic if between any two points
there exists a unique geodesic.
A suﬃcient condition to be uniquely geodesic is to require the passage times to have
a continuous distribution. Given the existence of ﬁnite geodesics, one can deﬁne the
geodesic graph as follows. For each x ∈Zd we deﬁne the (directed) graph T (x) with
vertex set Zd and whose edge set is the collection of all directed edges e that belong
to some ﬁnite geodesic between x and another point y, oriented as they are traversed
moving from x to y. This graph has the following properties:
1. T (x) is almost surely connected.
2. Every ﬁnite directed path in T (x) is a ﬁnite geodesic.
3. If the space is uniquely geodesic, then T (x) is a spanning tree of Zd. In this case,
T (x) is called the geodesic tree or the tree of infection of x.
77

<!-- page 78 -->

Many of the geodesic questions discussed in this survey can be rephrased in the
language of geodesic trees. Many of the “large scale” questions – for instance, the number
of ends of the tree T (x) – will require the study of inﬁnite geodesics, which we will deﬁne
shortly. Before discussing these large scale questions, let us focus on a natural small scale
question: how straight is the geodesic between x and y?
4.2
The wandering exponent
Given a set A ⊆Rd, let Γ(x, A) denote the geodesic from x to A (in case A = {y} is a
singleton, write Γ(x, y)). If the ﬁrst-passage metric is not uniquely geodesic but geodesics
exist, let Γ(x, A) instead denote the union of all geodesics between x and A (similarly
for A = {y}).
A natural asymptotic question for geodesics is: what does Γ(0, x) “typically” look
like for ∥x∥large? Since the ﬁrst-passage metric is in some sense a random perturbation
of the ℓ1 metric, one could perhaps imagine that Γ(0, x) often stays order ∥x∥distance
from the midpoint x/2. On the other hand, as discussed in Section 2, it is reasonably
well-accepted that the limit shape should have uniformly positive curvature, suggesting
that Γ(0, x) should perhaps look more like the straight line between 0 and x, as in the
case of the ℓ2 metric.
Indeed, the ℓ2 intuition is believed to be the correct one, and Γ(0, x) is believed to lie
close to the line
Lx = {mx : m ∈R}.
It is easy to see very weak versions of this claim under the assumption that the limit
shape is strictly convex. For instance, if x = x(n) = (n, n, . . . , n), and we denote by αn
the path which ﬁrst takes n steps down the e1 axis to connect 0 to ne1, then moves n
steps up the e2 axis, etc., we see by the shape theorem that
lim
n→∞T(αn)/n = dµ(e1) .
On the other hand, limn→∞T(0, x(n))/n = µ(e1 + . . . + ed), which is strictly smaller
if the limit shape is strictly convex. In particular, αn is not only not the geodesic to
x(n), but actually pays order ∥x∥more in passage time than the geodesic does. This
simple argument only rules out a single path, however, and stronger arguments have been
developed to simultaneously control all candidate geodesics.
Set D(0, x) to be the maximal distance of Γ(0, x) to Lx. An early argument showing
that Γ(0, x) lies close to Lx (under assumptions on the limit shape) appeared in [113,
Lemma (9.10)], where it was shown that, assuming that the limit shape has an exposed
point in the e1 direction,
lim
ε→0 lim inf
n→∞P (dist(D(0, ne1) < εn) = 1 .
(4.3)
To discuss strengthening (4.3), we will introduce the “wandering exponent” ξ.
78

<!-- page 79 -->

Roughly speaking, ξ is the number such that the maximal distance between Γ(0, x)
and Lx is typically of order ∥x∥ξ.
A simple attempt at giving a concrete deﬁnition
could be to set ε = nα−1 in the left-hand side of (4.3) and deﬁne α to be the smallest
number such that (4.3) still holds; however, there is no general agreement on the “correct”
deﬁnition, and we will discuss several possible alternatives. The exponent ξ is believed
to be universal in the sense that, for all reasonable choices of continuous distributions
F, one should ﬁnd the same value of ξ, and ξ should govern geodesic ﬂuctuations in all
directions. In general, however, ξ should depend on the dimension d. It is widely believed
that, for d = 2, the value of ξ should be 2/3 [106,107,110,111]. (For distributions in Mp
with p >⃗pc, deﬁned in Section 2.5.1, inﬁnite oriented paths of 1-edges give ξ = 1 – and
χ = 0 – for directions in the percolation cone.)
n
nξ
Figure 8: Representation of the wandering exponent ξ. The deviation of a geodesic from
0 to ne1 from the straight line that connects these points is of order nξ.
4.2.1
Upper bounds on ξ
We now discuss the most robust current upper bounds, from [135,137].
Letting C(x, m) denote the cylinder with axis Lx = {mx : m ∈R} and radius m –
i.e.,
C(x, m) = {z ∈Rd : dist(z, Lx) ≤m} ,
(4.4)
then one expects C(x, ∥x∥ξ) is the “right scale” to contain Γ(0, x). To be more precise:
one expects
lim
∥x∥→∞P
 Γ(0, x) ⊆C(x, ∥x∥ξ+ε)

= 1 ,
(4.5)
lim
∥x∥→∞P
 Γ(0, x) ⊆C(x, ∥x∥ξ−ε)

= 0 ,
(4.6)
for any ε > 0.
Note that (4.5) and (4.6) carry an implicit isotropy; if these equations are true, then
the order of geodesic wandering must be the same in all directions. More should be
79

<!-- page 80 -->

true: if we replace Γ(0, x) by, for instance, a geodesic to a hyperplane at distance n,
then this geodesic should asymptotically be contained in cylinders of radius nξ+ε (and
not in cylinders of radius nξ−ε). Lastly, “near-geodesics”should have the same behavior:
if ΓM(x, y) is the set of all paths η connecting x to y satisfying T(η) ≤T(x, y)+M, then
ΓM should satisfy relations like (4.5) and (4.6).
The article [137], in conjunction with its studies of χ, provides a technique for bound-
ing ξ in terms of χ (a more detailed discussion of the relation between the exponents can
be found in the next subsection). For a unit vector x, let
ξx = inf

α > 0 : lim inf
n
P
 Γ(0, nx) ⊂C(x, nα)

> 0

.
We also deﬁne the ﬂuctuation exponent
χ′ = inf

κ : (t −tκ)Bν ⊆B(t) ⊆(t + tκ)Bν for large t, a.s.

The following inequality is a version of one of the inequalities of the scaling relation
χ = 2ξ −1, discussed in Section 4.3 below. In their paper, Newman-Piza explain that
their argument is essentially a rigorous version of the argument of Krug-Spohn [123] for
this scaling relation.
Theorem 4.10 (Theorem 6 [137]). Assume d ≥2 with Eτ 2
e < ∞and F(0)
<
pc.
Assume that x is a direction of curvature for the limit shape (see Section 2.8). Then,
ξx ≤1 + χ′
2
.
A sketch of the argument that leads to Theorem 4.10 is given below, after Re-
mark 4.12.
In the same paper, under the assumptions of Theorem 3.7, the authors
also obtained the following bound, which they refer to as the Wehr-Aizenman bound,
after the latter’s work on phase transitions in disordered systems [158]:
χx ≥1 −(d −1)ξx
2
,
(4.7)
where
χx = sup

α > 0 : for some C > 0, Var T(0, nx) ≥Cn2α for all n

.
In dimension 2, Theorem 4.10 can be strengthened when combined with the previous
bound χ′ ≤1/2, derived by Kesten [115] and Alexander [4]. One still assumes that x is
a direction of curvature of the limit shape. The result is the following.
80

<!-- page 81 -->

Theorem 4.11 (Theorem 7 [137]). Assume F(0) < pc, and that E exp(βτe) < ∞for
some β > 0. Then, if x is a direction of curvature for the limit shape Bν,
ξx ≤3/4,
χx ≥1/8,
the latter bound holding when d = 2.
Remark 4.12. Note that, for any limit shape Bν, there is at least one direction of
curvature. Indeed, for any nonempty bounded subset B of Rd, there is at least one point
z in ∂B such that, for some closed ball D with positive radius, z ∈∂D and ¯B ⊂D. If
B is not just the origin, just take a ball of radius r large enough that contains B and
decrease r to the ﬁrst time it intersects ∂B. This argument appeared in [137].
Without the curvature assumption on the limit shape the best current bounds are
still 0 ≤ξ ≤1, the upper bound following directly from the shape theorem.
We now give a sketch of the argument for Theorem 4.10. Rather than prove this
inequality, we show a more intuitive Newman-Piza type argument under stronger as-
sumptions, like those made in Chatterjee’s work [43]. The argument also explains the
inequality ξx ≤3/4 from Theorem 4.11. Suppose that χ is any number for which there
is c > 0 such that
P(|T(0, z) −µ(z)| ≥λ∥z∥χ) ≤e−cλc for all z ∈Zd and λ ≥0.
(4.8)
(In fact, it is only necessary here to make this assumption for deviations of T(0, z) about
its mean, not µ(z) because, as pointed out in [43], one can run Alexander’s argument
under exponential bounds for T(0, z) −ET(0, z) to obtain upper bounds for non-random
ﬂuctuations on scale ∥z∥χ, with logarithmic corrections.) We aim now to show that if x
is a direction of curvature, then
ξx ≤1 + χ
2
.
(4.9)
Since the above exponential decay is known to hold for χ = 1/2 + δ for any δ > 0 (see
Theorem 3.11), one obtains ξx ≤3/4.
Sketch of proof of (4.9). The version of the Newman-Piza argument we present here is
from Chatterjee [43]. For simplicity, take x = e1. Let ϵ > 0, ˆξ = 1+χ
2
+ ϵ, and note that
since χ can be chosen arbitrarily close to 1/2, for ϵ small one has ˆξ < 1. Let L be the e1
axis, and for n ≥1, let L′ the segment connecting 0 and ne1. Now deﬁne
V = {v ∈Zd : ∃w ∈L′ such that ∥w −v∥∈[n
ˆξ, 2n
ˆξ]},
Then we must show that the probability that a geodesic from 0 to ne1 enters V goes to
0 with n. This will show that ξx ≤ˆξ and, taking ϵ →0, we obtain the result.
By a union bound, the probability of this event is bounded by
X
v∈V
P(v in a geodesic from 0 to ne1) =
X
v∈V
P(T(0, ne1) = T(0, v) + T(v, ne1)).
81

<!-- page 82 -->

The goal now is to show that the sum on the right is itself summable in n. We claim
there is c > 0 such that for all n and all v ∈V ,
P(T(0, ne1) = T(0, v) + T(v, ne1)) ≤(1/c) exp (−cnc) .
(4.10)
After summing over v and then over n, we complete the proof.
The main step is to show that, due to the curvature assumption, there is C such that
for all v ∈V ,
µ(v) + µ(ne1 −v) −µ(ne1) ≥Cn2ˆξ−1 = Cnχ+2ϵ.
(4.11)
Note that this inequality implies that if we replace T by µ in the event {T(0, v) +
T(v, ne1) = T(ne1)}, then the event cannot occur and, furthermore, gives us large lower
bound on the diﬀerence of both sides. Assuming this inequality for the moment, on the
event {T(0, v) + T(v, ne1) = T(ne1)}, one has
0 ≥[T(0, v) −µ(v)] + [T(v, ne1) −µ(ne1 −v)] −[T(ne1) −µ(ne1)] + Cnχ+2ϵ,
and so at least one of the ﬁrst three terms on the right has absolute value at least
(C/3)nχ+2ϵ. Applying (4.8) gives (4.10) and would complete the proof.
To show (4.11) is actually quite simple in the case that v = (n/2)e1 + n
ˆξe2: for
illustration, in this case, we use curvature and symmetry for
µ(v) + µ(ne1 −v) −µ(ne1) = [µ(v) −µ(n/2e1)] + [µ(ne1 −v) −µ(n/2e1)]
= 2[µ(v) −µ(n/2e1)]
= n
"
µ
 
e1 + 2n
ˆξ
n e2
!
−µ(e1)
#
≥Cn
 
n
ˆξ
n
!2
= Cnχ+2ϵ.
To modify this argument for other v’s, we split into cases. Let H be the hyperplane
taken in the deﬁnition of curvature (which in the case x = e1 is just orthogonal to e1 by
symmetry), and set w to be the projection of v onto L along H. Note that convexity of
the limit shape implies that µ(v) ≥µ(w) and µ(ne1 −v) ≥µ(ne1 −w). Furthermore, at
least one of the two vectors v or ne1 −v has Euclidean norm at least n/2.
1. Case 1: w ∈L′. In this case, we may assume that v is the vector with norm at
least n/2, and use curvature as above:
µ(v) −µ(w) = ∥w∥

µ
 w
∥w∥+ v −w
∥w∥

−µ
 w
∥w∥

≥C∥w∥

v −w
∥w∥

2
≥Cn2ˆξ−1.
82

<!-- page 83 -->

Therefore, since w ∈L′, one has µ(ne1 −w) + µ(w) = µ(ne1), and so
µ(v) + µ(ne1 −v) −µ(ne1) = [µ(v) −µ(w)] + [µ(ne1 −v) −µ(ne1 −w)]
≥µ(v) −µ(w) ≥Cn2ˆξ−1.
2. Case 2: w ∈L \ L′, on the side closer to ne1. In this case we still have as above
µ(v) −µ(w) ≥C∥w∥

v −w
∥w∥

2
.
Now
n2ˆξ ≤∥v −ne1∥2 = (∥w∥−n)2 + ∥v −w∥2 = (∥w∥−n)2 + ∥w∥2

v −w
∥w∥

2
,
and n ≤∥w∥≤3n.
Therefore either ∥w∥2 v−w
∥w∥

2
> n2ˆξ/2, or ∥w −ne1∥=
∥w∥−n ≥n
ˆξ/
√
2. In the ﬁrst case,
µ(v) + µ(ne1 −v) −µ(ne1) ≥µ(v) −µ(w) ≥C∥w∥2

v −w
∥w∥

2
∥w∥−1
≥Cn2ˆξ−1,
and in the second case, since ˆξ < 1,
µ(v) + µ(ne1 −v) −µ(ne1) ≥µ(w) −µ(ne1) ≥C∥w −ne1∥≥Cn
ˆξ
is still at least Cn2ˆξ−1.
3. Case 3: w ∈L \ L′, on the side closer to 0. This case is similar to the previous one.
These three cases together prove (4.11) and complete the proof.
4.2.2
Licea-Newman-Piza lower bound on ξ
We turn to the question of lower bounds for the wandering exponent ξ, taken up in [128].
The lower bounds proved in [128] all depend on the particular deﬁnition of ξ given. We
note a basic guideline for the possible lower bounds expected. As discussed earlier, one
expects ξ = 2/3 when d = 2, and, in view of the scaling relation χ = 2ξ −1 discussed
below, one expects ξ ≥1/2 for all d.
We are now ready to present the results of [128]. They consider four diﬀerent deﬁ-
nitions of ξ (in their words, proving these deﬁnitions are equivalent is “one of the open
83

<!-- page 84 -->

foundational problems of the subject”). Their ﬁrst deﬁnition, ξ(0), is the only point-to-
point exponent considered:
ξ(0) := sup

α ≥0 : lim sup
∥x∥→∞
P (Γ(0, x) ⊆C(x, ∥x∥α)) < 1

,
where C is as in (4.4); compare (4.6). The ﬁrst result of [128] is
Theorem 4.13. Let d ≥2. Then ξ(0) ≥1/(d + 1) .
The bound in Theorem 4.13 is extended to ξ(1), the second deﬁnition of ξ made
in [128]; we omit it here for brevity.
All the remaining bounds in [128] are for point-to-hyperplane versions of ξ, so we will
need notation for hyperplanes. For ˆθ a unit vector of Rd, let
Λ(ˆθ, L) := {y ∈Zd : y · ˆθ < L}
be the corresponding family of half-spaces. Let ∂Λ(ˆθ, L) be the set of y ∈Zd \ Λ(ˆθ, L)
adjacent to some vertex of Λ(ˆθ, L); the ∂Λ are the families of hyperplanes to which we
will consider geodesics.
For the next deﬁnition of ξ, let ΓM(x, A) be the set of “near-geodesics” from x to a
set A as above. Let end[ΓM(x, A)] be the set of endpoints in A of these near-geodesics.
ξ(2)(M) = sup

α ≥0 : ∃(ˆθn)n, (Jn)n with Jn →∞such that
there is no deterministic An with diam(An) ≤Jα
n
such that P

end[ΓM(0, ∂Λ(ˆθn, Jn))] ⊆An

→1

.
Note that in some respects, this deﬁnition is considerably weaker than the deﬁnition
of ξ(0): it concerns only near-geodesics, and it concedes a great deal of uniformity over
direction and “scale” (Jn). On the other hand, it is stronger in the sense that we actually
know where the wandering occurs (at the endpoints; in the deﬁntion of ξ(0) it is not
speciﬁed). The second theorem of [128] is
Theorem 4.14. Let M > inf supp(τe) and d ≥2. Then ξ(2)(M) ≥1/2.
The last deﬁnition of ξ considered has in some ways all the weaknesses of both ξ(0)
and ξ(2):
ξ(3)(M) := sup

α ≥0 : ∃(ˆθn)n, (Jn)n with Jn →∞such that
there is no deterministic xn such that
P

ΓM(0, ∂Λ(ˆθn, Jn)) ⊆C(xn, Jα
n )

→1

.
84

<!-- page 85 -->

Theorem 4.15. Let d = 2 and M > inf supp(τe). Then ξ(3)(M) ≥3/5.
We will attempt to justify the above results without providing complete proofs. In
some cases, unessential details or complications will be elided for presentation; for com-
plete arguments, see [128].
Proof of Theorem 4.13. The idea is an extension of the variance estimate of Newman-
Piza [137] discussed in Section 3.3. This idea was later used in W¨uthrich [163], Johansson
[109], Auﬃnger-Damron [15]. Let α > ξ(0) and choose some sequence of vertices xn with
∥xn∥→∞such that
P(Γ(0, xn) ⊆C(xn, ∥xn∥α)) →1 .
(4.12)
We will need to deﬁne a shifted version of xn; let x⊥
n be a unit vector perpendicular
to xn/∥xn∥and deﬁne 0′ and x′
n to be the closest vertices to 0 + C∥xn∥αx⊥
n and xn +
C∥xn∥αx⊥
n . where C is some large constant independent of n.
Deﬁne ∆T := T(0, xn) −T(0′, x′
n); the proof is by analysis of ∆T. Note that ∆T
cannot be too large in magnitude. Indeed, there are paths Γ1 and Γ2 connecting 0 to 0′
and xn to x′
n respectively, such that #Γi ≤C∥xn∥α. Therefore,
|∆T| ≤T(Γ1) + T(Γ2) ⪅C∥xn∥α
(4.13)
where the precise strength of the bound is a function of the tail behavior of τe (if τe is
bounded, then the left-hand side is really ≤C∥xn∥α).
On the other hand, we can get a lower bound on Var ∆T using the machinery of
Newman-Piza. Deﬁning as before the sigma-algebra Σj generated by the ﬁrst j edge
weights in some enumeration (ej) of the edges, we have
Var ∆T =
∞
X
j=1
E [E [∆T | Σj] −E [∆T | Σj−1]]2 .
(4.14)
As before, one can provide a lower-bound for (4.14) in terms of inﬂuences of individual
edges. There is an extra complication: when decreasing an edge weight τei, we are only
guaranteed to decrease ∆T if ei simultaneously satisﬁes
ei ∈Γ(0, xn) and ei /∈Γ(0′, x′
n) .
(4.15)
Letting F i be the event described in (4.15) (compare to (3.8)), we have similarly to before
Var(∆T) ≥C
∞
X
j=1
P(F j)2 .
(4.16)
For a further lower bound, we can restrict the sum in (4.16) to j such that ej ∈
C(xn, ∥xn∥α). For this restricted set J of indices,
X
j∈J
P(F j) ≥E# (Γ(0, xn) ∩C(xn, ∥xn∥α)) −E# (Γ(0′, x′
n) ∩C(xn, ∥xn∥α))
≥c∥xn∥1
85

<!-- page 86 -->

Ai
Sj
θi
Figure 9: The construction of the polygon PL. The segments Ai (in blue) are close to
the midpoints of the sides Si and, with high probability, contain the end point of the
geodesic Γ(0, Si) (in red).
for some c > 0 and all n large (the last inequality holds using (4.12)). As before, using
the Cauchy-Schwarz inequality on (4.16), we see
Var ∆T ≥
X
j∈J
P(F j)2 ≥
1
#C(xn, ∥xn∥α)
 X
j∈J
P(Fj)
!2
(4.17)
≥C∥xn∥1−α(d−1) .
Comparing this with (4.13) (which gives a variance upper bound of ∼∥x∥2α) completes
the proof.
Proof of Theorem 4.14. The proof of this theorem relies on a geometrical construction.
Rather than trying to exactly give the parameters of this construction, we will outline the
main idea; we leave it to the reader to ﬁll in the details (and in particular to understand
why ξ(2) is deﬁned as it is). We ﬁrst outline the argument in case d = 2, then explain
the necessary changes for d > 2.
Assume ξ(2) < 1/2, and choose 1/2 > α > ξ(2) + ε. We choose a ﬁnite set of unit
vectors {ˆθi}N
i=0 and construct a polygon PL (having the symmetries of Z2) whose sides
are segments of lines perpendicular to ˆθi:
PL =
\
i
Λ(ˆθi, Li) , where Li ∼L .
Let Si be the side of PL corresponding to ˆθi–that is, the side of PL which is a segment of
the boundary of Λ(ˆθi, Li). We will choose ˆθi and Li such that the lengths |Si| ∼Lα. We
will also ensure our choice of ˆθi guarantees geodesics to the lines deﬁning PL are localized
on sides of PL.
Let Ai be a deterministic subset of ∂Λ(ˆθi, L) of diameter smaller than Lξ(2)+ε in
which end [Γ(0, ∂Λ(ˆθi, L))] is localized with high probability; such a subset exists by the
86

<!-- page 87 -->

deﬁnition of ξ(2). Our main restriction in choosing ˆθi is that Ai should be close to the
midpoint of Si, and in particular be at least distance Lξ(2)+ε from Sj, j̸ = i. This is done
inductively.
To start the induction: we can choose ˆθ0 = e1 and L0 = L, with |S0| = Lα. Then
by the symmetries of the lattice, we can choose A0 symmetric about the e1 axis (and
in particular localized near the midpoint of S0).
Our last vector ˆθN will be the “45
degree line” (
√
2)−1(e1 + e2); again, by the symmetries of the lattice, AN will be chosen
symmetric about ˆθN. The remaining vectors will only be chosen to ﬁll out the polygon
in the arc of angles [0, π/4]; once this is done, we ﬁll the rest by symmetry.
Assuming ˆθ0, . . . , ˆθj−1 and the corresponding {Li} have been chosen correctly, let
zj−1 be the counterclockwise endpoint of Sj−1. If zj−1 is within distance ≤CLα from
a point on the line {aˆθN : a ∈R}, then choose ˆθj = ˆθN and Lj is chosen so that
zj−1 ∈∂Λ(ˆθj, Lj). Otherwise, we will choose |Sj| = Lα, and we need to ﬁnd a direction
ˆθj with arg ˆθj ∈(arg ˆθj−1, π/4) such that Aj is appropriately localized on Sj (once we
choose ˆθj, then Lj is ﬁxed by the necessity of connecting to Sj−1).
The fact that such a choice is possible is somewhat subtle. If we were to try to choose
ˆθj extremely close to ˆθj−1, then we would fail. Indeed, Aj−1 is close to the midpoint
of Sj−1, and as ˆθj →ˆθj−1, the side Sj is essentially a subsegment of the boundary of
Λ(ˆθj−1, Lj−1) and so Aj should be identical to Aj−1, which lies order Lα distance below
the midpoint of Sj. Similarly, taking ˆθj →ˆθN would give an Aj lying too far above
the midpoint.
By a continuity argument, we can therefore ﬁnd a ˆθj that meets our
requirements.
After this somewhat lengthy construction, we can quickly conclude. Note that, be-
cause |Si| ≤CLα, we have
arg(ˆθi) −arg(ˆθi−1) ≤CLα−1.
In particular, ∂Λ(ˆθi−1, Li−1) comes within distance L2α−1 = o(1) of the midpoint of
Si.
But this means it takes only O(1) edges to extend a time-minimizing path for
T(0, ∂Λ(ˆθi, Li)) from Si to ∂Λ(ˆθi−1, Li−1)–since this extension costs very little in the way
of passage time, there is a “near-geodesic” to ∂Λ(ˆθi−1, Li−1) lying on Si, contradicting the
placement of Ai−1. We conclude that α cannot be chosen < 1/2; that is, ξ(2) ≥1/2. This
argument assumes that T(0, Si) < T(0, Si−1), which, without loss in generality, holds
with positive probability.
The case d ≥3 is very similar.
Now instead of using polygons PL, we consider
polygonal “barrels” of the form PL × Rd−2, where PL is constructed similarly to before.
The arguments above may be adapted to these barrels without substantially diﬀerent
ideas.
Proof of Theorem 4.15. The proof is in some ways heavily motivated by the proofs of
both of the preceding theorems. Consider α > ξ(3). We construct a polygon PL with
87

<!-- page 88 -->

sides {Sj} similarly to the arguments in the proof of Theorem 4.14. For ﬁxed j, we
consider
∆T = T(0, ∂Λ(ˆθj+1, Lj+1)) −T(0, ∂Λ(ˆθj, Lj)).
As before, ∂Λ(ˆθj, Lj) comes within distance ∼L2α−1 of every point on Sj+1, so
geodesics to ∂Λ(ˆθj, Lj) may be extended to ∂Λ(ˆθj+1, Lj+1) at an increase in passage time
⪅L2α−1. This (combined with a symmetric argument in the other direction) gives an a
priori bound
Var(∆T) ⪅L4α−2 .
(4.18)
As in the proof of Theorem 4.13, we ﬁnd an α-dependent lower bound which con-
tradicts this when α is too small. The geodesics to both hyperplanes are conﬁned with
high probability in cylinders with radius Lα. These cylinders are largely disjoint for L
large, and so there are typically order L edges which are on Γ(0, ∂Λ(ˆθj, Lj)) but not on
Γ(0, ∂Λ(ˆθj+1, Lj+1)).
With this in hand, we lower-bound Var(∆T) with a martingale expansion as in The-
orem 4.13; we ultimately arrive at an analogue of (4.17). As before, this yields
Var(∆T) ≥CL1−α .
Comparing this bound with (4.18) completes the contradiction.
4.3
The scaling relation χ = 2ξ −1
Consider a geodesic from the origin to a point x with passage time T(0, x).
One of
the central questions in ﬁrst-passage percolation (and in related models) is to establish
the following statement. There exists an intrinsic relation between the magnitude of
deviation of T(0, x) from its mean and the magnitude of deviation of the geodesic T(0, x)
from a straight line joining 0 and x. This relation is universal; that is, it is independent
of the dimension d and of the law of the weights (as long they satisfy certain moment
assumptions).
As discussed in Section 3.1 (see (3.1)), the ﬂuctuations of the passage time T(0, x)
about ET(0, x) should be of order |x|χ, where χ is the ﬂuctuation exponent. Analogously,
the transversal or wandering exponent ξ studied in Section 4.2 measures the maximal
Euclidean distance of a geodesic from 0 to x from the straight line that joins 0 to x. The
intrinsic relation described above should be given as
χ = 2ξ −1.
(4.19)
As previously discussed, the existence and the ‘correct’ deﬁnition of these exponents
is part of the problem. Before we discuss the mathematical history behind (4.19) and
current progress [15,43], let us stress that, as of today, there is not a single distribution
of passage times, where we know that (4.19) holds. Therefore proving equation (4.19) is
still an open question. (As we will see in the pargraphs below, one will probably need to
solve Question 4 ﬁrst.)
88

<!-- page 89 -->

Question 22. Find a distribution of passage times where (4.19) holds.
As shown in Theorem 4.10, a version of the inequality χ ≥2ξ −1 was proved and
understood in the 1995 work of Newman-Piza [137]. A proof of the other inequality
appeared ﬁrst in a model of a Brownian particle in a Poissonian potential [163]; see
also [109]. In FPP, Chatterjee [43] proposed a stronger deﬁnition of the exponents that
allows a proof of (4.19). Unlike the deﬁnitions of ξ and χ given in the previous sections,
Chatterjee’s exponents are not known to exist for any distribution of passage times. The
proof in [43] relies on a construction similar to that in [44]. One ﬁrst breaks a geodesic
into smaller segments and then uses an approximation scheme to compare the passage
time to a sum of nearly i.i.d. random variables. The proof is then a trade-oﬀbetween
minimizing the error and maximizing the variance of the passage time. Assuming that
the distribution is ‘nearly Gamma’ (see Section 3.2), the optimization can be achieved
by choosing diﬀerent parameters in the approximation.
The scaling relation should be valid only under the unproven assumption of uniform
positive curvature of the limit shape. This is one of the main reasons for the introduction
of a strong deﬁnition of exponents. It is reasonable to believe that given a diﬀerent
curvature exponent κ as deﬁned in (2.36), one would have the alternate scaling relation
χ = κξ −(κ−1); see [15, Section 3] for more details. The inequality ≤comes from (4.22)
below, and the case κ = 2 reduces to χ = 2ξ −1.
We will now sketch how one could derive (4.19) assuming positive curvature. Here
we follow the proof of (4.19) given in [15]. It starts with a standard fact: if X′ is an
independent copy of a random variable X then
Var X = 1
2E(X −X′)2.
Thus, if we want to estimate χ, it suﬃces to compare the diﬀerence of two independent
copies of T(0, ne1) as
n2χ ∼E(T(0, ne1) −T ′(0, ne1))2.
(4.20)
The exponent ξ tells us exactly how to build such an independent copy T ′. As in
the proof of Theorem 4.13, one just needs to consider passage times from starting points
0 and nξe2, as the geodesics from there to ne1 and nξe2 + ne1 will live in two disjoint
cylinders. Therefore, one should expect that
n2χ ≲E(T ′(0, ne1) −T(0, ne1))2 ∼E(T(nξe2, nξe2 + ne1) −T(0, ne1))2.
(4.21)
Now, recall the deﬁnition of the curvature exponent given in Deﬁnition 2.36. If we
are allowed to bound the passage time by its asymptotic value (as T(0, x) ∼µ(x)), we
would have
89

<!-- page 90 -->

T(nξe2, ne1) −T(0, ne1) ≤µ(ne1 −nξe2) −µ(ne1) + error
= n

µ(e1 −nξ−1e2) −µ(e1)

+ error
∼nnκ(e1)(ξ−1) + error
= n2ξ−1 + error,
(4.22)
if κ(e1) = 2. In the ﬁrst equality above, we used the fact that µ is a norm. The above
sequence of approximations tells us 2ξ −1 comes from the assumption that the curvature
exponent exists and it is equal to 2.
Now we compare the left side of (4.22) with the right side of (4.21). The diﬀerence
in passage times in both equations is almost the same except that the ending points are
not exactly equal. However, a crossing trick from [15] and an argument from [43] using
assumed exponential concentration on scale nχ from the strong deﬁnition of χ allows us
to take the same ending points, and thus we can replace (4.22) in (4.20) to get χ ≤2ξ−1.
The above argument is made rigorous in [15]. The assumptions on the exponents ξ
and χ are just to justify the approximations in (4.22).
4.4
Inﬁnite geodesics
In this section we investigate inﬁnite geodesics. We start with the following deﬁnitions.
Deﬁnition 4.16. We will say that an inﬁnite self-avoiding nearest-neighbor path Γ
is an inﬁnite geodesic (for a given edge-weight realization) if every ﬁnite subpath of Γ is
a ﬁnite geodesic.
Inﬁnite geodesics come in two varieties:
1. Indexed by N; that is, Γ’s edges in order are e1, e2, . . . . These are called unigeodesics,
singly inﬁnite geodesics, geodesic rays or (when it does not cause confusion) simply
geodesics.
2. Indexed by Z. These are called bigeodesics or geodesic lines.
Much work on geodesics in the model has focused on determining the number of
geodesics of either type which exist, as well as their properties and relationships.
4.4.1
Existence of geodesic rays
It is a simple argument to check that, given any x ∈Zd there exists (almost surely)
at least one geodesic ray originating from x. The argument goes as follows. Suppose
without loss of generality that x is the origin 0.
Consider the sequence of geodesics
90

<!-- page 91 -->

Γn = Γ(0, ne1). The ﬁrst vertex of each Γn is 0. There are 2d possible choices for the
second vertex of each Γn. Therefore, there must be some edge e incident to the origin
such that inﬁnitely many Γn have e as their ﬁrst edge. Repeating this argument on
subsequences yields a singly inﬁnite path which is seen to be a geodesic.
The natural question that follows is to determine if there exist multiple distinct
geodesic rays with positive probability (here, distinct means sharing at most ﬁnitely
many edges and vertices). A natural way to construct distinct geodesics would be to
consider the subsequential limit obtained in the last construction by taking a diﬀerent
sequence of endpoints, such as −ne1. Given the observations on the wandering exponent
given above, it is reasonable to believe that one should be able to construct many inﬁnite
geodesics this way. However, this intuition is hard to make rigorous.
In the two-dimensional case, with exponentially distributed edge weights, H¨aggstr¨om
and Pemantle [95] showed that there are at least two distinct geodesic rays with positive
probability. This result was extended to a wide range of ﬁrst-passage distributions by
Garet and Marchand [79] and Hoﬀman [101].
An important further advance was made in a paper of Hoﬀman [102]. He demon-
strated that it is in fact possible to derive the existence of more than two geodesic rays
directly from properties of the limit shape.
His theorem goes as follows.
We write
sides(Bν) for the number of sides of the limit shape Bν (that is, sides(Bν) is the number
of extreme points of Bν and is ﬁnite if and only if Bν is a polygon).
Theorem 4.17 (Hoﬀman [102]). Assume d = 2 and τe has continuous distribution with
Eτ 2+α
e
<
∞for some α > 0. Deﬁne G(x1, . . . , xk) to be the event that there exist
distinct geodesic rays beginning at vertices x1, . . . , xk. Then for any k ≤sides(Bν) and
any ϵ > 0, there exist x1, . . . , xk ∈Z2 such that
P(G(x1, . . . , xk)) > 1 −ϵ.
This result also holds for a wide range of non-i.i.d. passage time distributions and was
used by Hoﬀman to establish existence of 4 distinct geodesic rays. In light of Hoﬀman’s
results, the following question could be seen as a weak version of the uniform curvature
conjecture, since uniform curvature implies sides(Bν) = ∞. We will discuss this theorem
(and give some idea of its proof) in another guise below as Theorem 6.4.
Question 23. Show that for continuous distributions, there are inﬁnitely many
geodesic rays.
It is expected that sides(Bν) = ∞for any distribution that has more than one point
in its support. As we saw in Section 2, this fact is only known for distributions in Mp
with p ≥⃗pc.
Question 24. Find a distribution that is not in Mp that has an inﬁnite number of
geodesic rays.
91

<!-- page 92 -->

The last observation combined with Hoﬀman’s result allows us to conclude:
Theorem 4.18 (Auﬃnger-Damron [16]). Assume d = 2. For any measure ν ∈Mp,
p ≥⃗pc, almost surely, there exist inﬁnitely many geodesic rays.
It is important to note that when p >⃗pc, for ν ∈Mp, one can trivially ﬁnd inﬁnitely
many geodesic rays by choosing them to be oriented up-right and contain only edges e
with te = 1. In the above result, the geodesics constructed can be chosen to contain a
positive fraction of edges e with te > 1. Each of these rays will be directed in a sector
that is disjoint from the Durrett-Liggett ﬂat edge of the limit shape.
The argument that we provided at the beginning of this section for the existence of a
geodesic ray leads to another important question. If Γn is a sequence of ﬁnite geodesics,
we say that Γn →Γ for some path Γ if for all N > 0, the ﬁrst N steps of Γn equal those
of Γ for all n large.
Question 25. Is it true that the sequence of geodesics Γ(0, ne1) converges to some
Γ?
By the deﬁnition of convergence, such Γ is clearly a geodesic ray. The question above
was solved if one replaces Z2 by the upper-half plane or any other inﬁnite connected
subgraph of Z2 with inﬁnite connected complement in [17] (in the latter case, ne1 is
replaced by a sequence of points on the boundary of the domain). The following similar
open question (a version of the “BKS midpoint problem” posed in [23]) is related to the
existence of geodesic lines, which is the topic of Subsection 4.5.
Question 26. Assume that the space is uniquely geodesic. Show that in d = 2, there
is no ﬁxed vertex that is in inﬁnitely many of the geodesics between −ne1 and ne1.
Both of the last two questions are solved if one assumes that the limit shape boundary
is diﬀerentiable both in direction e1, and at θ1, θ2, where the θi’s are the endpoints of the
sector of angles of contact of the unique tangent line to the limit shape in direction e1.
See Theorem 5.20 and Remark 5.24.
4.4.2
Directions and coalescence
In this section, we investigate questions of whether there exist geodesic rays that have
an asymptotic direction, and if the geodesic in a given direction is unique. Much of the
work described here was done in the ’90s by Newman and co-authors; for a complete
pedagogical account of these results we refer to the book of Newman [136]. (Also see
Section 5 on Busemann functions, where geodesics are analyzed without a curvature
assumption.) For θ ∈Rd with ∥θ∥2 = 1 (identiﬁed with [0, 2π) when d = 2), we will
say that a self-avoiding singly inﬁnite path (v1, e1, v2, . . .) has direction θ if arg vn :=
vn/∥vn∥2 →θ as n goes to inﬁnity. Under the assumption of uniform positive curvature,
the following theorem shows that inﬁnite geodesics cannot wind around the origin; they
must have a direction.
92

<!-- page 93 -->

Theorem 4.19 (Newman [135]). For d ≥2, assume that F(0) < pc, there exists α > 0
such that E exp(ατe) < ∞and that Bν is uniformly curved (see Deﬁnition 2.37). Then
with probability one,
1. Every geodesic ray has a direction.
2. For every θ with ∥θ∥2 = 1, there exists a geodesic ray with direction θ.
Newman’s theorem from [135, Theorem 2.1] is stated for continuous distributions, in
which both items above are deduced not for geodesic rays, but for inﬁnite branches in
T (0), but the same proofs yield the above extension.
Proof. We will need to control geodesics from all points at once, and we will control them
by ensuring that they exit cones along certain parts of their boundaries. For this, we
need deﬁnitions of cone sectors: for x ∈Zd and ϵ > 0, deﬁne the set
Cx = {z ∈Zd : ∥z∥∈[∥x∥/2, 2∥x∥], ∥arg x −arg z∥≤∥x∥−1/4+ϵ},
its “forward boundary”
∂fCx = {z ∈Zd : ∃w ∈Cx such that ∥z −w∥= 1 and ∥z∥> 2∥x∥},
and the rest of the boundary is ∂′Cx = {z ∈Zd \ Cx : ∃w ∈Cx such that ∥z −w∥=
1} \ ∂fCx. Deﬁne the “out tree” of x as
out(x) = {z ∈Zd : T(0, z) = T(0, x) + T(x, z)}
(the set of vertices whose geodesics from 0 contain x) and the event
Gx = {out(x) ∩∂′Cx̸ = ∅}.
A union bound gives
P(Gx) ≤
X
z∈∂′Cx
P(T(0, z) = T(0, x) + T(x, z)).
Nearly the same argument as in the Newman-Piza bound (summing a similar bound to
that given in (4.10)) shows that there is c > 0 such that for all x ∈Zd,
P(Gx) ≤(1/c)e−c∥x∥c.
(Here it is important that the aperture of the cone Cx be at least ∥x∥−1/4+ϵ, so that
the distance between x and the side boundaries of Cx be larger order than ∥x∥3/4. The
exponent 3/4 is the bound for the wandering exponent from Theorem 4.11.)
By Borel-Cantelli, we can almost surely ﬁnd a random M such that for ∥x∥≥M, the
event Gc
x occurs. Now ﬁx any such conﬁguration; we claim that there is a C > 0 such
93

<!-- page 94 -->

that the following geodesic-direction estimate holds. Let γ be a geodesic from 0 with
vertices 0 = x0, x1, x2, . . . , xr. So long as n ∈[0, r] is such that ∥xn∥≥M, then
∥arg xr −arg xn∥≤C∥xn∥−1/4+ϵ,
(4.23)
where C = P∞
i=1 2−i(1/4−ϵ). To prove this inequality, deﬁne a sequence of points y1, . . .
inductively as follows. Let y1 = xn and for i = 2, . . ., deﬁne yi to be the ﬁrst point on
γ after yi−1 and after the ﬁrst exit from Cyi−1 (as long as there is one). Then there is I
such that xr ∈CyI, and we can bound
∥arg xr −arg xn∥≤
I−1
X
i=1
∥arg yi −arg yi+1∥+ ∥arg yI −arg xr∥
≤
I
X
i=1
∥yi∥−1/4+ϵ
≤
∞
X
i=1
(2i∥y1∥)−1/4+ϵ
≤C∥xn∥−1/4+ϵ.
In the second to last inequality, we have used that Gc
yi occurs, and thus the geodesic
exits Cyi along the forward boundary. This proves (4.23).
Now we can argue the two points of Newman’s theorem.
For any passage time
conﬁguration as above, let γ be an inﬁnite geodesic. We may assume it starts from 0 and
label its vertices x0, x1, . . .. Pick n such that ∥xn∥≥M and note now that by (4.23),
one has ∥arg xn −arg xm∥≤∥xn∥−1/4+ϵ for all m ≥n. This means (arg xn) is a Cauchy
sequence and thus converges, proving that γ has a direction.
For the second point, for our given passage time conﬁguration, pick any direction
θ and let (yn) be any sequence of vertices with arg yn →θ. Let γn be a self-avoiding
geodesic from 0 to yn and let (γnk) be a subsequence that converges to some γ. We claim
that γ, which is a self-avoiding inﬁnite geodesic, has direction θ. Indeed, if its vertices
are x0, x1, x2, . . ., then for any n, one must have for all large k that γnk has as its initial
segment the portion of γ from 0 to xn. Thus if ∥xn∥≥M, by (4.23),
∥arg ynk −arg xn∥≤∥xn∥−1/4+ϵ.
Taking k →∞, we obtain ∥θ −arg xn∥≤∥xn∥−1/4+ϵ, and this means arg xn →θ.
Under the unproved assumption of uniform curvature, one could now ask how many
geodesic rays with direction θ exist. Is it possible that there is more than one such ray,
and does the answer change if we insist that they share the same initial vertex?
In this direction, the following two results are due to Licea and Newman [127]. The
ﬁrst one roughly says that geodesic rays cannot bifurcate in the same deterministic direc-
tion. We note, however, that birfucations are expected in random directions, and indeed
94

<!-- page 95 -->

such a statement has been proved in a related last-passage model [49]. See Remark 4.21
below.
Theorem 4.20 (Licea-Newman [127], Newman [135]). Let d = 2 and assume that F is
continuous. Fix x ∈Z2. Then there exists some set D ⊆[0, 2π) of full Lebesgue measure
such that if θ ∈D, then there is zero probability that there exist distinct geodesic rays
starting at x with direction θ.
The proof of Theorem 4.20 is surprisingly simple, and we will provide a sketch of the
main idea below. Let us ﬁrst mention that not much is known about the set D, although
it is expected that D = [0, 2π). Establishing the question below would be a step towards
one of the main problems (absence of geodesic lines) discussed in Section 4.5.
Question 27. Show that D = [0, 2π).
In [136, Theorem 1.5], it was shown by Zerner that D has at most a countably inﬁnite
complement. But we do not even know that a particular angle - for instance, 0 - belongs
to D.
Question 28. Show that 0 ∈D.
These two questions have been recently solved by Damron-Hanson under additional
assumptions on the limit shape boundary [60]. Namely, the ﬁrst is solved if one assumes
that the limit shape boundary is diﬀerentiable everywhere. For the second, the condition
is diﬀerentiability of the boundary in the direction e1 and in directions θ1 and θ2, where
the θi’s are the endpoints of the sector of angles at which the unique tangent line to the
limit shape in direction e1 contacts the limit shape. See Theorem 5.20.
Proof of Theorem 4.20. Suppose there are two diﬀerent geodesic rays r1, r2 starting at
some x with the same direction θ. As they are diﬀerent, they will have to bifurcate at
a vertex u ∈Z2, one taking the directed edge (u, v1), the other the edge (u, v2), with
v2̸ = v1. As the space is uniquely geodesic, r1 and r2 do not intersect after u. The
construction that follows buids a geodesic between r1 and r2 by exploring the geodesic
tree T (x) between these two rays. By planarity, any other geodesic ray that sits between
these two geodesics must also have direction θ.
Assume without loss of generality that r1 is asymptotically counter-clockwise to r2.
Now construct the geodesic ray r+(u, v2) that starts with the edge (u, v2) and always
chooses the most counter-clockwise edge possible in T (x) (in particular this edge must
belong to at least one geodesic ray). Roughly speaking, the ray r+(u, v2) is the ‘closest’
geodesic ray to r1 that takes the edge (u, v2) in the region delimited by r1 and r2. As r2
is a geodesic ray, r+(u, v2) sits between r2 and r1 and must have direction θ. Through
this construction, it is not diﬃcult to see that to any direction θ with multiple geodesic
rays from x corresponds at least one edge e such that r+(e) has direction θ. Let G(e, θ)
be the event that r+(e) has direction θ.
95

<!-- page 96 -->

u
r1
r2
r+(u, v2)
Figure 10: Construction of the path r+(u, v2). Path r1 depicted in black, r2 in purple
and r+(u, v2) in blue.
We thus get the estimate
P
 there exist at least two geodesic rays from x with direction θ

≤
X
e
P(G(e, θ)).
Now, for each edge e, r+(e) cannot be a geodesic ray with direction θ for more than one
θ; thus Fubini’s theorem implies
Z
P(G(e, θ))dθ =
Z Z
1G(e,θ)(ω)dθdP = 0.
(One must check that 1G(e,θ)(ω) is jointly measurable in θ, ω.) Combining the last two
displays ends the proof of the theorem.
Remark 4.21. It is expected (and, as mentioned above, proved under the assumption
of diﬀerentiability of ∂Bν) that for any ﬁxed θ, one cannot ﬁnd disjoint geodesic rays
with direction θ. However, there should exist (and in fact, it is provable if one assumes
existence of geodesics in any deterministic direction) a random, dense, countably inﬁnite
set of θ’s for which there are disjoint geodesic rays with direction θ. Such a statement is
proved in related exactly solvable last-passage percolation model by Coupier [49], along
with the statement that almost surely, there exist no θ such that three geodesic rays
from 0 have direction θ. Although this last result does not yet have a counterpart in
non-exactly solvable models, progress toward the ﬁrst under limit shape assumptions was
made in [85].
The second theorem in this direction establishes that for all θ ∈D any two geodesic
rays with direction θ (starting from possibly diﬀerent points) must coalesce.
Theorem 4.22 (Licea-Newman [127], Newman [135]). For d = 2 and θ ∈D with F
continuous, there is zero probability that there are disjoint geodesic rays in the direction
θ.
The proof of the above theorem is two-dimensional. It is unclear if the result holds
for d > 2.
96

<!-- page 97 -->

Question 29. Decide whether Theorem 4.22 holds in arbitrary dimension.
Ending this subsection, we quote an apparently still open question raised in the end
of the original paper of Hammersley and Welsh.
Question 30. Let Γ(v) be a geodesic from 0 to v. Consider the collection of all edges
that belong to an inﬁnite number of Γ(v) and let f(r) be the number of such edges which
intersect the circle x2 + y2 = r2. Does f(r) →∞as r →∞; and, if so, how fast?
Sketch of the proof of Theorem 4.22. We sketch the proof in the case where the weights
are unbounded. The bounded weight case requires a little modiﬁcation of the arguments,
and the details can be found in the original paper [127, Theorem 1].
Consider all geodesic rays that have asymptotic direction θ ∈D.
Construct the
geodesic graph G composed of the union of all these geodesic rays. By Theorem 4.20, if
two such rays meet, they must coalesce. Thus, the geodesic graph must be a forest with
N ≥1 distinct trees. The proof in [127] has 3 steps.
Step 1.
If P(N ≥2) > 0 then P(N ≥3) > 0.
To see this, without loss of generality, assume that the direction θ has positive e1-
coordinate.
Also suppose that the two diﬀerent geodesics touch the e2 axis at two
points y1̸ = y2 and are contained in the half-plane x ≥0. By vertical translation in-
variance, it is possible to ﬁnd with positive probability two diﬀerent points y3, y4 that
also have two such disjoint geodesic rays. By planarity, at least three of these rays must
be distinct (see ﬁgure). In fact, if P(N ≥2) > 0 then there exist y1 < y2 < y3 in Z
such that P(E(y1, y2, y3)) > 0, where E(y1, y2, y3) is the event that there are three dis-
joint geodesic rays with direction θ starting from (0, yi) and contained in the half-plane
{(x, y) ∈Z2 : x > 0} except their initial points.
y1
y2
y4
y3
Figure 11: Construction of the three disjoint geodesics starting from y1, y2, y3. The
geodesic starting from y4 may intersect the geodesic emanating from y2 but all others
are disjoint.
97

<!-- page 98 -->

Let F be the subevent of E(y1, y2, y3) on which the tree in G containing the “middle”
point y2 does not intersect the half-plane {(x, y) ∈Z2 : x < 0}. The second step of the
proof is the following.
Step 2. If P(N ≥3) > 0 then P(F) > 0.
By the arguments in step 1, with positive probability there are three disjoint geodesic
rays r1, r2, r3 starting respectively from the points (0, y1), (0, y2), (0, y3) ∈L. We now
set y∗= max{|y1|, |y3|} and increase the values of the passage times on edges between
u = (−1, y) and v = (0, y), −y∗≤y ≤y∗. If we increase this ﬁnite collection of passage
times enough, we create a barrier through which no geodesic can pass. (Speciﬁcally, if
F is this set of edges, one needs that for all e ∈F, τe is greater than Te(x, y), where
x, y are the endpoints of e, and Te is the passage time among all paths which do not use
e.) Note that, by our construction, this change also does not alter the r′
is. Furthermore,
as we have unbounded passage times, one can construct such a barrier with positive
probability. Now, since there is zero probability that two geodesic rays with direction θ
meet and do not coalesce, this barrier attached to L together with r1 and r3 implies that
the tree that contains r2 is the desired tree.
Step 3. P(F) > 0 is impossible.
Consider all possible translates Lu of L by points of the form 3uy∗for u ∈Z2 and the
corresponding translated events F u of F. If F u and F v both occur, then, by deﬁnition,
the corresponding “middle” geodesics must be disjoint, as they are blocked from each
other by the “outer” geodesics constructed in the event F. Now, let CM be the number
of Lu’s contained in the box [0, M] × [0, M] and NM be the random number of the
corresponding F u’s which occur. By translation invariance, and the assumption that
P(F) > 0,
ENM = CMP(F) ≥cM 2,
(4.24)
for some c > 0. But NM is bounded above by the number of disjoint trees that touch
[0, M] × [0, M]. Since each tree is inﬁnite, this number cannot exceed the number of
boundary points in that box, which is at most c′M for some c′. As c′M < cM 2 for M
large, this contradicts (4.24), showing that P(F) = 0.
4.5
Absence of geodesic lines and connection to the two-
dimensional Ising ferromagnet
A ﬁnal (and one of the most) important question in the study of inﬁnite geodesics is the
following:
Question 31. Do geodesic lines exist?
4.5.1
Heuristic argument
We start with a heuristic argument that there should (a.s.)
exist no geodesic lines
(bigeodesics), at least in dimensions d where ξ > 1/2. In particular, the argument gives
98

<!-- page 99 -->

a plausible reason to disbelieve in bigeodesics in two dimensions. It is important to note
that this argument is far from rigorous even under usual unproven assumptions (e.g.
curvature) on the model.
Indeed, a main feature of the argument is that it uses various meanings of ξ in-
terchangably. For instance, it is assumed that the geodesics Γ(0, ne1) and Γ(mξ, ne1)
(m ≪n) merge after a distance of order m, in addition to Γ(0, ne1) having transversal
wandering typically of order nξ.
This argument was presented by C. Newman at the workshop “First-passage perco-
lation and related models” at AIM in Summer, 2015.
Setup. Consider the box [−N, N]d for N large; we denote the boundary of this box by
the symbol ∂N. Assume that there were positive probability that a bigeodesic existed
(and hence that this event had probability one). Then this bigeodesic would have to pass
through the origin, at least with positive probability.
Such a bigeodesic must pass through a vertex of ∂N before passing through 0 and
another vertex of ∂N afterwards. In particular, the event that a geodesic between a pair
of vertices of ∂N passes through 0 must have uniformly (in N) positive probability. To
rule out bigeodesics, it therefore suﬃces to show
lim
N→∞P
 ∃x, y ∈∂N such that 0 ∈Γ(x, y)

= 0 .
(4.25)
The key idea to prove (4.25) is to use the wandering property of individual geodesics to
show that they are likely to avoid 0.
Unfortunately, using a union bound on the events {0 ∈Γ(x, y)}x,y∈∂N is unlikely to
give (4.25), since there are many such pairs (x, y). To overcome this diﬃculty, we will
group geodesics together, tiling ∂N with (d −1)-dimensional cubes Bj of side length
∼N ξ. For a particular cube Bj, let B−j denote the corresponding cube on the opposite
side of ∂N.
Reducing to opposite blocks. We will need to make assumptions on the structure
of geodesics between vertices in blocks Bi, Bj.
The ﬁrst is that we can neglect the
contribution of the probability
P (∃x ∈Bi, y ∈Bj for some i̸ = −j such that 0 ∈Γ(x, y))
to (4.25).
This is a not unreasonable assumption, since the straight line segment between such
x and y will typically lie distance ≫∥x −y∥ξ from 0 (perhaps requiring the diameter
of Bj to be at least N ξ+ε). This of course is likely to necessitate a strong concentration
bound on the variable
D(0, ne1)/nξ
where D was deﬁned above (4.3). In any case, once such a result is established, it allows
us to reduce the problem of showing (4.25) to showing
lim
N→∞P
 ∃j, x ∈Bj, y ∈B−j such that 0 ∈Γ(x, y)

= 0 .
(4.26)
99

<!-- page 100 -->

Reducing within blocks. We have (at least assuming ξ < 1) succeeded in reducing
the number of pairs of geodesics under consideration by a power by reducing to opposite
pairs of blocks. Since geodesics emanating from vertices within the same Bj are likely
to be highly correlated, we can reduce further by grouping together such geodesics. In
particular, we will assume that all geodesics originating in one block Bj merge before
coming close to the origin.
For each block Bj, choose a particular vertex xj near the center of Bj. As discussed
earlier, it is likely that there is “geodesic merging on scale ξ”. In particular, for any y ∈
Bj, y′ ∈B−j, it is plausible that the event {0 ∈Γ(y, y′)} has extremely small symmetric
diﬀerence with the event {0 ∈Γ(xj, x−j)}. Thus, for ﬁxed j we can approximate
P
 ∃x ∈Bj, y ∈B−j such that 0 ∈Γ(x, y)

≈P
 0 ∈Γ(xj, x−j)

.
(4.27)
With this approximation established, we have reduced the problem suﬃciently to take
a union bound:
(4.26) ≤P
 ∃j such that 0 ∈Γ(xj, x−j)

≤#(boxes Bj) P
 0 ∈Γ(x1, x−1)

.
(4.28)
Individual geodesic wandering.
We now need to control the probability P(0 ∈
Γ(x1, x−1)), which will require a further assumption on geodesic wandering. In some
sense, what is required here is a lower bound on ﬂuctuations, that Γ(x1, x−1) is “spread
out” on scale N ξ far from its endpoints.
To be more precise, if L is the line between x1 and x−1 and H denotes the hyperplane
passing through 0 to which L is normal, we assume that Γ(x1, x−1) passes through ∼1
vertices of H and every vertex in a block of diameter N ξ is equally likely to be passed
through. In particular,
P(0 ∈Γ(x1, x−1)) ∼N −ξ(d−1) .
With this information, we can conclude by bounding the right-hand side of (4.28).
Note that the number of boxes Bi is of the order of
N d−1
N ξ(d−1) = N (1−ξ)(d−1) ,
giving
(4.28) ⪅N (1−ξ)(d−1)
N ξ(d−1)
=
 N 1−2ξd−1 .
In particular, if ξ > 1/2, then (4.25) holds, giving a.s. absence of bigeodesics.
100

<!-- page 101 -->

4.5.2
Rigorous results
In this section, we move to rigorous results, outlining the work of Licea-Newman and
Wehr from the ’90s. In Section 5.6, we give the most recent results, ruling out bigeodesics
in any deterministic direction under assumptions on the limit shape boundary.
The ﬁrst partial answer to the bigeodesic question comes from Theorem 4.23. We say
that a geodesic line has asymptotic directions (θ1, θ2) if its two rays have directions θ1,
θ2.
Theorem 4.23 (Licea-Newman [127]). For d = 2, θ1, θ2 ∈D, and F continuous, there
is zero probability that there exists a geodesic line that has directions (θ1, θ2).
It turns out that the proof of Theorem 4.23 is not diﬃcult once we have Theorem
4.22.
Proof of Theorem 4.23. If there are two distinct geodesic lines with directions (θ1, θ2),
then by Theorem 4.22 they must coalesce in both directions. This cannot happen, as we
assumed the model is uniquely geodesic. Thus, there must be at most one geodesic line
with directions (θ1, θ2). Let E be the event that there exists a unique geodesic line with
directions (θ1, θ2). By the ergodic theorem, the probability of the event E is either 0 or
1.
Suppose P(E) = 1. Then it is possible to ﬁnd L large enough such that the geodesic
line intersects the box BL := [0, L] × [0, L] with probability at least 2/3. Choosing a
sequence of points zn ∈Z2 with zn/|zn| →α̸ = θ1 or θ2, the geodesic line passes through
only ﬁnitely many shifted boxes zn + BL. Letting ES be the event that the geodesic line
intersects the set S ⊆Z2, this implies
P(EBL ∩Ezn+BL) →0 as n →∞.
On the other hand, by translation invariance,
P(EBL ∩Ezn+BL) ≥1 −P(Ec
BL) −P(Ec
zn+BL) ≥1/3,
a contradiction.
We note that it is possible to replace this proof with one that proceeds along the lines
of Wehr’s proof that the number of bigeodesics in all of Zd must be zero or inﬁnity:
Theorem 4.24 (Wehr [157]). Assume F is continuous, d ≥2, and τe has ﬁnite mean.
The number N of bigeodesics in a conﬁguration (τe) is almost surely constant, and N ∈
{0, ∞}.
The idea of the proof is that if there are only ﬁnitely many bigeodesics, then by
translation invariance, in a large box, at least one of the bigeodesics must take a positive
density of edges. So the passage time along this bigeodesic from its ﬁrst entry to its
101

<!-- page 102 -->

last exit of the box must be of the order of the volume of the box. This contradicts
the fact that the Euclidean distance between these two points is of the order of the
length of the box, and so the passage time between them is also of this same order, with
high probability. This theorem has been extended to some non-i.i.d. times by Boivin-
Derrien [28].
In that paper, they also construct non-i.i.d.
but translation invariant
distributions on (τe) for which there exist bigeodesics.
A positive answer to Question 27, combined with the above result, would exclude the
possibility of geodesic lines with deterministic asymptotic directions. It would not solve
Question 31, as one could still have geodesic lines in random directions.
Question 31 was answered negatively in the upper half-plane H of Z2 under passage
times with ﬁnite mean and continuous distribution by Wehr and Woo [159].
Theorem 4.25 (Theorem 1 [159]). Assume F is continuous and τe has ﬁnite mean.
With probability 1, there exist no geodesic lines in H.
We outline the proof of Theorem 4.25 at the end of this section. The assumption of
ﬁnite mean was relaxed in [17], where one just requires uniqueness of ﬁnite geodesics.
That work also included a range of non-i.i.d. times.
Before sketching Wehr-Woo’s argument, let us stress that the importance of determin-
ing non-existence of geodesic lines goes beyond the ﬁeld of FPP. Apparently, this problem
was ﬁrst posed by H. Furstenberg in a diﬀerent context. As far as we know, the ﬁrst math-
ematical reference is in [113]. Question 31 is equivalent to the problem of existence of
nontrivial (that is, non-constant) ground states in the two-dimensional Ising ferromagnet
with random exchange constants. Ground states of higher-dimensional random ferromag-
nets are similarly related to hypersurfaces with minimal random weights [114,157].
More precisely, consider the lattice dual to Z2, deﬁned by
(Z2
∗, E2
∗) = (Z2, E2) + 1
2(e1 + e2),
where E2 is the set of nearest-neighbor edges of Z2, and deﬁne a spin conﬁguration as an
element σ of {±1}Z2
∗. Let (Jx,y)(x,y)∈E2∗be a collection of i.i.d. positive random variables.
For any conﬁguration and any ﬁnite set S ⊂Z2
∗deﬁne the (random) energy
HS(σ) = −
X
{x,y}∈E2∗
x∈S
Jx,yσxσy.
We will call σ a ground-state for (Jx,y) if, for each conﬁguration ρ such that ρx = σx
for all x outside some ﬁnite set, we have
HS(σ) ≤HS(ρ) for all ﬁnite S ⊂Z2
∗.
It is an open problem to describe the set of ground states for this ferromagnetic model. In
particular, it is not known how many ground states there are for a given (Jx,y), although
102

<!-- page 103 -->

it is conjectured (see, for instance, the discussion in [136, Chapter 1]) that, if the law of
J is continuous, there should be only two almost surely.
These two are the constant conﬁgurations σ = ±1, which are clearly almost surely
ground states. If any non-constant ground states σ exist, they cannot have ﬁnite regions
of disagreement; that is, there can be no ﬁnite S such that σx = +1 for all x ∈S and
σy = −1 for all y ∈∂S, or vice-versa. Therefore, any non-constant ground state must
have a two-sided (and circuitless) inﬁnite (original lattice) path of edges dual to edges
{x, y} such that σx = −σy.
Now consider FPP with τe = J(x,y) if {x, y} is the dual of the edge e. If such a
ﬁrst-passage conﬁguration had a geodesic line, then the conﬁguration σ which takes the
value +1 on one side of the geodesic line and −1 on the other would be a non-constant
ground state for the associated spin model. Therefore in d = 2, based on the analogy
above with the Ising ferromagnet, it is believed that one should have a negative answer
for Question 31.
Partial results in this direction were obtained by Newman under the assumption of
uniform curvature. In [135], it was shown that there cannot exist geodesic lines other
that those that have directions on both ends. As we have seen also, Theorem 4.23 ruled
out almost all possible pairs of deterministic directions and Theorem 4.24 showed that
the number of geodesic lines is 0 or inﬁnity. There are also further arguments in the
physics literature [77] against the existence of non-constant ground states.
Let us end this subsection by presenting the argument of Wehr-Woo for absence of
geodesic lines in the upper half-plane H. For simplicity, we assume that the passage time
τe has ﬁnite mean and continuous distribution so that ﬁnite geodesics are unique almost
surely. Let
E = {there exists a geodesic line in H}.
By horizontal translation ergodicity, P(E) is zero or one; so let us assume for a contra-
diction that P(E) = 1.
Any geodesic line Γ divides R2 \ Γ into two components, say R+ = R+(Γ) and
R−= R−(Γ); that is,
R+(Γ) ∩R−(Γ) = ∅,
R+(Γ) ∪R−(Γ) = R2 \ Γ ,
∂R+ = ∂R−= Γ ,
where R−is a region that contains (0, −1) and where ∂A denotes the usual boundary of
a set A ⊂R2. Hence, by uniqueness of passage times, for any points x, y ∈R−(Γ), no
edge b belonging to the ﬁnite geodesic Γ(x, y) can intersect R+(Γ).
Step 1. The ﬁrst step of the proof is to establish the following fact [159, Proposi-
tion 4]. If a geodesic line exists with probability one, then there is a lowest one.
103

<!-- page 104 -->

Sn
Γ
Γn
Γm
R+(Γ)
R−(Γ)
Figure 12:
The geodesic line Γ0 (in red) and the regions R+(Γ0) and R−(Γ0). The inter-
section of Γ and the strip Sn is a union of bounded sets. The geodesics Γn monotonically
converge to Γ0.
Proposition 4.26. Assume that P(E) = 1. Consider the sequence of geodesics Γn from
(−n, 0) to (n, 0) for n ∈N. With probability 1, this sequence has a limit:
Γ0 = lim
n→∞Γn .
Moreover, Γ0 is a geodesic line, and for any geodesic line Γ,
Γ0 ⊂

R−(Γ) ∪Γ

.
Proof. Suppose that Γ is a geodesic line.
By uniqueness of passage times, for every
n, Γn must lie in R−(Γ) ∪Γ. Now, compactness (that is, the same argument used to
establish existence of geodesic rays) implies that (passing through a subsequence) Γnk
must converge.
For instance, there is only a ﬁnite number of possible edges for the
intersection of Γn with the e2, axis and this intersection is nonempty for every n ≥1.
However, the geodesics Γn are increasing in the sense that for n > m, Γn must lie above
Γm. Thus, they converge.
Step 2. If Γ0 intersects the strip Sn = {(x, y) ∈Z2 : 0 ≤y < n} then R+(Γ0) ∩Sn
is a non-empty union of bounded connected sets. See Figure 12.
This step is the main tool to prove:
Step 3. Uniquely geodesic implies that Γ0 must intersect any large box with proba-
bility bounded below uniformly in the position of the box.
Precisely, for l ∈N, write B = B(l) = [−l, l] × [0, 2l] and let K be the event that at
least one geodesic line intersects B. Deﬁne, for L ∈N, translations of B by
Bi,j = Bi,j(l, L) = B + (iL, jL) for (i, j) ∈Z2, j ≥0 .
For L > 2l, the Bi,j are mutually disjoint.
104

<!-- page 105 -->

Proposition 2. Let δ = 1 −P(K). Then
P

Bi,j ⊂R+(Γ0)

≤δ,
P

Bi,j ⊂R−(Γ0)

≤δ and
(4.29)
P

Γ0 ∩Bi,j̸ = ∅for all 0 ≤i, j < k

≥1 −2k2δ.
(4.30)
To prove (4.29), one uses Step 1 and Step 2. Consider the reﬂected half-planes
Hr
j =

(x, y) ∈R2 : y ≤jL + 2l
	
and set the event
Ki,j :=

at least one Hr
j −geodesic line intersects Bi,j
	
.
Since the model is invariant under translations and under rotations by 180 degrees,
P(K) = P(Ki,j).
(4.31)
If Bi,j ⊂R+(Γ0) then, by Step 2, Bi,j is contained in exactly one of the bounded
components of R+(Γ0) ∩Hr
j. Hence, there is no Hr
j-geodesic line intersecting Bi,j, since
otherwise two diﬀerent ﬁnite Hr
j ∩H-geodesics would meet more than once. Thus the
event {Bi,j ⊂R+(Γ0)} is contained in the event Kc
i,j, which combined with (4.31) leads
to the ﬁrst part of (4.29).
The second equation in (4.29) uses Step 1 and the fact that the box Bi,j would be
trapped by the union of the geodesics Γn for n large and the boundary of the upper-half
plane. Therefore no geodesic line in the shifted half-plane Hj = {(x, y) ∈R2 : y ≥jL}
can intersect Bi,j. To get (4.30), one uses (4.29) to obtain
P(Γ0 ∩Bi,j̸ = ∅) ≥1 −2δ
for any i, j. Thus,
P

\
0≤i,j<k
{Γ0 ∩Bi,j̸ = ∅}

≥1 −2k2δ.
Step 4. Equation (4.30) contradicts P(E) = 1. Thus, P(E) = 0.
To see this, assume P(E) = 1 and for any k ≥1, choose l = l(k) large enough so that
B = B(l) satisﬁes
P(Γ0 ∩B̸ = ∅) = P(K) > 1 −1/4k2.
By (4.30), this leads to
P

\
0≤i,j<k
{Γ0 ∩Bi,j̸ = ∅}

≥1/2.
(4.32)
105

<!-- page 106 -->

The proof now consists of showing that (4.32) cannot hold for k large if L, the distance
of the shift, is large enough.
To see this, choose L > 2l so all shifted boxes Bi,j are disjoint. If the event in (4.32)
happens, it is possible to ﬁnd k2 diﬀerent points in the geodesic line Γ0, each one inside
a diﬀerent box Bi,j. Let xm, 1 ≤m ≤k2 be these points, ordered as they appear in Γ0.
As the ℓ1-diameter of the union of the Bi,j is at most 2Lk, one has
L −2l < ∥x1 −xk2∥1 < 2Lk and |Γ(x1, xk2)| ≥k2L/4,
where |Γ(x1, xk2)| is the number of edges in the geodesic Γ(x1, xk2). There at most 4l2k2
choices for each x1 and xk2 so an union bound leads to
P

\
0≤i,j<k
{Γ0 ∩Bi,j̸ = ∅}

≤4l2k4
sup
∥x−y∥1≤2Lk
P

|Γ(x, y)| > k2L/4

.
If k (and therefore l) is ﬁxed large enough then, by Theorem 4.6, the supremum on the
right side is exponentially small in L. This contradicts (4.32) for L large, ending the
proof of Theorem 4.25.
106

<!-- page 107 -->

5
Busemann functions
The modern approach to inﬁnite geodesics involves Busemann functions. These functions
were ﬁrst exploited in FPP in the important papers of Hoﬀman [101,102], and they allow
one to obtain some information about geodesics without the unproven assumptions (for
example, positive curvature of the limiting shape) of Newman.
5.1
Basics of Busemann functions
Long before FPP existed, H. Busemann invented various tools to study the geometry
of geodesics in certain metric spaces [34]. Busemann functions grew out of the attempt
to understand parallelism between geodesics. We can give some basic properties in the
setting of FPP. Consider an inﬁnite geodesic R (which we will call a geodesic ray), and
list its vertices in order as
V (R) = {x0, x1, . . .}.
Here, x0 is taken as the initial vertex of the ray.
Deﬁnition 5.1. The Busemann function associated to the ray R is fR : Zd →R,
given by
fR(x) = lim
n [T(x, xn) −T(x0, xn)] .
This limit exists since the terms are bounded (|T(x, xn) −T(x0, xn)| ≤T(x0, x)) and
monotone:
T(x0, xn+1) −T(x, xn+1) = T(x0, xn) + T(xn, xn+1) −T(x, xn+1)
≥T(x0, xn) −T(x, xn).
The Busemann function for R measures asymptotically how far behind x0 is relative
to x when both points attempt to travel down the ray R. The next lemma collects a few
properties of this function.
Lemma 5.2. The following holds.
(a) For m < n, one has fR(xm) −fR(xn) = T(xm, xn).
(b) For all x, one has |fR(x)| ≤T(x, x0).
(c) If R1 and R2 are geodesic rays that coalesce (they have ﬁnite symmetric diﬀerence)
and have initial points x0 and y0, then
fR1(x) = fR2(x) −fR2(x0) for all x ∈Zd.
107

<!-- page 108 -->

Proof. Proofs of items (a) and (b) follow directly from the deﬁnition. As for (c), write
the vertices of R1 as x0, x1, . . . and the vertices of R2 as y0, y1, . . .. By coalescence, we
can ﬁnd N ∈Z such that for n large, one has xn = yn+N. For such n,
T(x0, xn) −T(x, xn) = T(x0, yn+N) −T(y0, yn+N) + T(y0, yn+N) −T(x, yn+N)
→fR2(x0) −fR2(x).
We will often use a Busemann function that does not depend on the initial point. If
R is a geodesic ray, then set
fR(x, y) = fR(x) −fR(y).
We note the following properties of this Busemann function:
1. fR is additive in the following sense:
fR(x, y) = fR(x, z) + fR(z, y) for x, y, z ∈Zd.
2. If R1 and R2 are geodesic rays that coalesce, then fR1(x, y) = fR2(x, y) for all
x, y ∈Zd.
Proof. Using item (c) in Lemma 5.2 and writing x0 for the initial point of R1, one
has
fR1(x) −fR1(y) = (fR2(x) −fR2(x0)) −(fR2(y) −fR2(x0))
= fR2(x) −fR2(y).
3. (Translation covariance) Let θ be a translation of the lattice by an integer vector.
Then
fR(x, y)(ω) = fθR(θ(x), θ(y))(θ(ω)).
Here, the translated weight-conﬁguration θ(ω) is deﬁned as te(θ(ω)) = tθ−1e(ω) for
any edge e. Furthermore θR is the translated geodesic ray.
5.2
Hoﬀman’s argument for multiple geodesics
In this section we present Hoﬀman’s [102] argument that there exist at least two inﬁnite
geodesics almost surely.
This result was also proved by Garet-Marchand [79] in less
generality, but Hoﬀman’s techniques, involving Busemann functions, led to many other
results. We will use the assumptions of Hoﬀman:
108

<!-- page 109 -->

1. P is ergodic under lattice translations.
2. P has unique geodesics.
3. Et2+δ
e
< ∞for some δ > 0.
4. The limit shape for P is bounded.
We say that two geodesic rays R1 and R2 are disjoint if the cardinality of both sets
R1 \ R2 and R2 \ R1 are inﬁnite. Let N = N(ω) be the maximal number of disjoint
geodesic rays in the edge-weight conﬁguration ω. It can be shown that N is a measurable
function (relative to the product sigma-algebra) and is invariant under lattice translation.
So by ergodicity, it is almost surely constant.
Theorem 5.3. The number N is at least two.
Proof. By taking a sub-sequential limit of the geodesics from 0 to ne1, we see that N > 0.
Assume for a contradiction that N = 1. Let R = R(ω) be any geodesic ray and deﬁne
f : Zd × Zd →R by
f(x, y) = fR(x, y) = fR(ω)(x, y)(ω).
On the probability 1 set on which N = 1, this function is well-deﬁned and is independent
of the choice of R. Indeed, if R′ = R′(ω) is another geodesic ray, then it cannot be disjoint
from R. Due to unique geodesics, R′ must coalesce with R and therefore fR = fR′.
One can also check that f is measurable. Due to translation covariance of fR, the
function f is also translation covariant:
f(x, y)(ω) = f
 θ(x), θ(y)

(θ(ω)) almost surely
(5.1)
for any lattice translation θ. This, combined with additivity and the ergodic theorem,
implies for each x ∈Zd,
1
nf(0, nx) = 1
n
n
X
k=1
f((k −1)x, kx) →Ef(0, x) almost surely and in L1.
To ﬁnd this limit, we use symmetry: for k = 1, . . . d,
0 = Ef(0, ek) + Ef(ek, 0) = Ef(0, ek) + Ef(−ek, 0)
and by translation invariance, this equals
Ef(0, ek) + Ef(0, ek).
We then obtain
Ef(0, mek) = mEf(0, ek) = 0
109

<!-- page 110 -->

for all integers m and then by additivity,
Ef(0, x) = 0 for all x ∈Zd.
Combining with (5.1), one has
1
nf(0, nx) →0 almost surely and in L1.
To upgrade this convergence to an omni-directional statement, one can use the con-
dition Eτ 2+δ
e
< ∞along with the bound |f(x, y)| ≤T(x, y) to prove a type of shape
theorem for the function f. That is, for each ϵ > 0,
P
 |f(0, x)| < ϵ∥x∥1 for all x ∈Zd with ∥x∥1 large enough

= 1.
In particular, if we deﬁne R to be the (almost surely unique) limit of geodesics from 0
to ne1, and we write its vertices in order as x0, x1, . . ., then one has
f(0, xn)/∥xn∥1 →0 almost surely.
On the other hand, by property 1 of fR, one has f(0, xn) = T(0, xn). Therefore
T(0, xn)/∥xn∥1 →0 almost surely.
This contradicts the fact that the limit shape is bounded and therefore shows N > 1.
5.3
Directions of geodesics via Busemann functions
In [135], Newman introduced the following limits, in an attempt to study the local
geometry of the boundary of the growing ball B(t):
lim
n [T(x, xn) −T(y, xn)] for x, y ∈Z2,
where (xn) is a deterministic sequence diverging to inﬁnity in a direction θ. He was able
to show that under a strong (unproven) curvature assumption on the limit shape, the
limit exists almost surely for Lebesgue-almost every direction.
In 2008, Hoﬀman [102] deﬁned Busemann-type functions similar to the above, and
introduced the idea of using Busemann functions to prove properties of geodesics. These
ideas were continued by Damron-Hanson in 2013 [59]. In the next two sections, we will
describe some of the ideas therein. For simplicity, we will assume that the τe’s are i.i.d.,
with continuous distribution having ﬁnite mean, although this is far from necessary.
Let Hn be the hyperplane {x ∈Rd : x · e1 = n} and deﬁne the function
Bn(x, y) = T(x, Hn) −T(y, Hn) for x, y ∈Zd.
For purposes of illustration, we will in this section assume that
B(x, y) := lim
n Bn(x, y) exists almost surely
(5.2)
for x, y ∈Zd. This assumption should be compared to assumption (2.14).
110

<!-- page 111 -->

Lemma 5.4. Assume (5.2). One has
EB(x, y) = (x −y) · ρ for x, y ∈Zd,
where ρ = e1g(e1).
Proof. By translation invariance,
EB(x, y) = EB(θ(x), θ(y))
for any lattice translation θ. As B is additive and satisﬁes EB(0, −x) = −EB(0, x), we
ﬁnd that EB(x, y) is a linear function of x −y. To ﬁnd this function, it suﬃces to ﬁnd
EB(0, ei) for i = 1, . . . , d. By symmetry, one has EB(0, ei) = 0 for i > 1. To ﬁnd the
other expectation, we use an averaging trick introduced by Garet-Marchand and Hoﬀman
(and also appears in some form in Kingman [118, Eq. (26)]). By translating by e1, we
obtain
ET(0, Ln) =
n
X
k=1
E[T((k −1)e1, Hn) −T(ke1, Hn)]
=
n
X
k=1
E[T(0, Hn−k+1) −T(e1, Hn−k+1)]
=
n
X
k=1
EBk(0, e1).
Using the bound |Bk(x, y)| ≤T(x, y) and dominated convergence, one has
EB(0, e1) = lim
n
1
nET(0, Hn) = µ(e1).
The last equality follows from the shape theorem.
Using the lemma, along with the ergodic theorem, one has almost surely and in L1,
1
nB(0, nx) = 1
n
n
X
k=1
B((k −1)x, kx) →x · ρ for x ∈Zd.
It is not diﬃcult to upgrade this to a sort of shape theorem for the Busemann function.
See, for example, Damron-Hanson [59, Section 4], where this is done for a “reconstructed
Busemann function.”
Lemma 5.5. Assume (5.2) and that ET(x, y) < ∞for all x, y. For each ϵ > 0,
P (|B(0, x) −x · ρ| > ϵ∥x∥1 for only ﬁnitely many x) = 1,
where ρ = e1g(e1).
111

<!-- page 112 -->

We now show how assumption (5.2), along with a mild curvature assumption, implies
existence of geodesic rays with asymptotic direction e1. Let S be the set
S = ∂B ∩{w : ρ · w = 1}.
Because the latter set is a supporting hyperplane for the limit shape at the point e1/g(e1),
the set S is a portion of the boundary ∂B containing e1/g(e1).
We will say that a
geodesic ray with vertices x0, x1, . . . is asymptotically directed in S if each limit point
of {xn/g(xn)} is contained in S. If S contains only one point, we say that this ray has
an asymptotic direction. If the point e1/g(e1) is exposed (there is a hyperplane for B
that touches B only at e1/g(e1)), then by symmetry, S = {e1/g(e1)}. In that case, the
following theorem implies that almost surely, there is an inﬁnite geodesic from 0 that has
asymptotic direction e1. It is an open problem to show without assumption (5.2) that
with positive probability, there is an geodesic ray with an asymptotic direction. The
following theorem is similar to [59, Theorem 5.3].
Theorem 5.6. Assume (5.2). With probability one, every subsequential limit of geodesics
from 0 to Hn is asymptotically directed in S.
Proof. Let γ be any subsequential limit of geodesics from 0 to Hn and label its vertices
0 = x0, x1, x2, . . .. First note that for m ≥0, one has
B(0, xn) = T(0, xn).
Indeed, choose (nk) a subsequence such that some geodesics γnk from 0 to Hnk converge
to γ. For k large, the ﬁrst n vertices of γnk coincide with those of γ, so
B(0, xn) = lim
N [T(0, HN) −T(xn, HN)] = lim
k [T(0, Hnk) −T(xn, Hnk)] = T(0, xn).
Choose z to be any limit point of {xn/g(xn)} so that for some subsequence (nj), one
has xnj/g(xnj) →z. Now by the Busemann shape theorem, almost surely for any choice
of γ and (nj),
lim
j
B(0, xnj)
g(xnj)
= lim
j ρ · (xnj/g(xnj)) = ρ · z.
On the other hand, by the shape theorem and the above remarks, this equals
lim
j
T(0, xnj)
g(xnj)
= 1.
Therefore ρ · z = 1. Since z ∈∂B, this means z ∈S.
Corollary 3. Assume (5.2) and that e1/g(e1) is an exposed point of B. With prob-
ability one, every subsequential limit of geodesics from 0 to Hn has asymptotic direction
e1.
112

<!-- page 113 -->

5.4
Busemann increment distributions and geodesic graphs
It is not known how to show that limits like (5.2) exist, and this is one of the ﬁrst and
main obstacles when dealing with Busemann functions.
Question 32. Prove that under some conditions on the edge weights the limit (5.2)
exists.
Remark 5.7. For the generalized Busemann functions introduced in Section 2.6,
the question above translates to Question 10. One may also want to consider Busemann
limits to points rather than lines.
Here we will explain one method of getting around Question 32 to obtain partial
results about directions of geodesic rays. This method shares some similarities with the
plan of attack used in the original paper of Kingman, described in Section 2.6. (We make
this relation precise below in Remark 5.15.) Instead of trying to establish the limit (5.2),
we will focus on a weak limit by considering translational averages. We follow Damron-
Hanson [59]. Their results were stated only for dimension d = 2, but the existence of
geodesics directed in sectors (the ﬁrst theorem below) also holds for general dimensions,
so we outline this more general argument.
We will focus on i.i.d. measures, although the theorems apply to a wide class of
translation-ergodic measures. Consider the following conditions on the distribution P of
the edge-weights:
I. P is i.i.d. with EY d < ∞, where Y is the minimum of 2d i.i.d. edge-weights.
Furthermore, P(τe = 0) < pc, the bond percolation threshold for d dimensions.
II. P has continuous marginals.
Recall that a supporting hyperplane H for B at a point z0 ∈∂B is a hyperplane that
contains H and such that B does not intersect both components of Hc. If there is only
one such hyperplane, we say that B is diﬀerentiable at z0. In this case, we write Hz0 for
this hyperplane. Write
Sz0 = Hz0 ∩∂B.
Theorem 5.8. Assume (I). If ∂B is diﬀerentiable at z0, then with probability one, there
is a geodesic ray γ containing the origin which is asymptotically directed in Sz0. This
means every limit point of {x/g(x) : x ∈γ} is contained in Sz0.
Note that if z0 is an exposed point of diﬀerentiability, then there is a geodesic ray
with asymptotic direction z0, since in that case, Sz0 = {z0}. This is an improvement
on Newman’s theorem because we have replaced the global curvature condition with a
local, directional condition.
The next result is related to coalescence of geodesic rays.
113

<!-- page 114 -->

Theorem 5.9. Take dimension d = 2, assume (I) and (II) and that ∂B is diﬀerentiable
at z0.
With probability one, there exists a collection of geodesic rays (Γx : x ∈Z2)
satisfying the following properties:
1. Each x ∈Z2 is a vertex of Γx.
2. Each Γx is asymptotically directed in Sz0.
3. For all x, y ∈Z2, the paths Γx and Γy coalesce.
4. Each x ∈Z2 is on Γy for only ﬁnitely many y ∈Z2.
As discussed in Section 4.2, Licea-Newman [127] (along with an improvement by
Zerner [136]) have shown (under a global curvature condition) that there exists a de-
terministic set D ⊂[0, 2π) with countable complement such that, for each θ ∈D, the
following holds almost surely.
1. (Existence) There exists a collection of inﬁnite geodesics (Γx : x ∈Z2) such that
each Γx starts from x, has asymptotic direction θ, and each Γx and Γy coalesce.
2. (Uniqueness) For each x, Γx above is the only inﬁnite geodesic starting from x with
asymptotic direction θ.
The above theorem shows that one can take D = [0, 2π) for existence, but does not
address uniqueness.
Idea of proof. The idea is to work with sub-sequential Busemann limits in distribution.
First we set Ω1 = [0, ∞)Ed to be a copy of our passage time space and let Ω2 = RZd×Zd
be our space for recording Busemann increments. Last, we have a space Ω3 = {0, 1}⃗Ed,
where⃗Ed is the set of oriented nearest-neighbor edges of Zd, which will record geodesic
graphs. Put eΩ= Q3
i=1 Ωi.
Let H be any supporting hyperplane for B and let ρ be the vector with H = {w :
w · ρ = 1}. Deﬁne Hα = {w : w · ρ = α} for α ∈R. We now deﬁne Busemann increments
and geodesic graphs toward Hα. For an edge-weight conﬁguration ω ∈Ω1, set
Bα(ω) =
 Bα(x, y) : x, y ∈Zd
∈Ω2,
where Bα(x, y) = T(x, Hα) −T(y, Hα). Furthermore deﬁne the geodesic graph conﬁgu-
ration ηα(ω) by
ηα(ω)(⟨x, y⟩) =
(
1
if {x, y} ∈Gα(z) for some z and Bα(x, y) ≥0
0
otherwise
.
Here, Gα(z) is the union of all edges in geodesics from z to Hα. Last, deﬁne Φα : Ω1 →eΩ
by
Φα(ω) = (ω, Bα(ω), ηα(ω))
114

<!-- page 115 -->

and µα the push-forward of µ through Φα. The measure µα is a Borel probability measure
on eΩ.
A generic element of eΩwe write as
eω = ((te), B, η).
B is called the reconstructed Busemann function. We would like to take α to inﬁnity
and sample a conﬁguration eω from the limit of (µα).
Such a measure would be the
distribution of Busemann increments and geodesic graphs “at inﬁnity” in a direction
related to H. However, we do not know how to show that this sequence converges, so
we settle for sub-sequential limits of an average of µα over α. The averaging is done to
ensure that any limit is translation invariant. That is, set
µ∗
k = 1
k
Z k
0
µα dα
and let µ∗be any sub-sequential limit of (µ∗
k). (Tightness is not diﬃcult to verify, due to
the condition |Bα(x, y)| ≤T(x, y).) The measure µ∗
k can be thought of as representing
ﬁrst sampling a hyperplane Hα uniformly at random for α ∈[0, k] and then sampling
the conﬁguration Φα(ω).
The idea is then to try to approximate the proof of last section to (a) ﬁnd asymptotics
of the function B(0, x) as x →∞and (b) use them to control geodesics. First, we give
the mean of B, writing Eµ∗for expectation relative to µ∗. Recall that ρ is the vector
with Hα = {w : w · ρ = α}.
Lemma 5.10. Assume (I). The mean of the reconstructed Busemann function is
Eµ∗B(x, y) = ρ · (y −x) for x, y ∈Zd.
The idea of this lemma is that the Busemann function on average measures distances
traveled perpendicular to H, but from a sort of stationary state, where the velocity is
given by the time constant in direction z0 (which is the same time constant corresponding
to the line H). The reader should compare to Lemma 5.4 in the previous section, where
the mean of the Busemann function is similarly given by a projection along the vertical
hyperplane H used to create the function.
Proof. The proof uses a version of the averaging trick (inspired by Gou´er´e [88]) from last
section. For k ≥1 we can use stationarity to write
Eµ∗
kB(−x, 0) = 1
k
Z k
0
ET(−x, Hα) dα −
Z k
0
ET(0, Hα) dα

= 1
k
Z k+ρ·x
k
ET(0, Hα) dα −
Z ρ·x
0
ET(0, Hα) dα

.
115

<!-- page 116 -->

Choosing kj so that µ∗
kj →µ∗, one can show using translation invariance and the bound
Eµ∗
kj B(−x, 0)2 ≤ET(−x, 0)2 that the left side converges to Eµ∗B(0, x). So we obtain the
formula
Eµ∗B(0, x) = lim
j→∞
Z ρ·x
0
ET(0, Hα+kj)
kj
dα.
By the shape theorem, the integrand converges to 1, so by dominated convergence, we
obtain ρ · x. By translation invariance, then, of µ∗, one has
Eµ∗B(x, y) = Eµ∗B(0, y −x) = ρ · (y −x).
As usual, we can prove asymptotics for the reconstructed Busemann function in the
form of a shape theorem. However, due to the extra randomness introduced through
α, we may have lost ergodicity for µ under translations, and therefore obtain a random
“shape.”
Lemma 5.11 (Shape theorem for B). Assume (I). There exists a random vector ϱ ∈R2
such that for any ϵ > 0,
µ∗ |B(0, x) −x · ϱ| > ϵ∥x∥1 for inﬁnitely many x ∈Zd
= 0.
The vector ϱ satisﬁes the following conditions:
1. µ∗-almost surely, the hyperplane
Hϱ := {w ∈Rd : w · ϱ = 1}
is a supporting line for B at z0.
2. The mean of ϱ under µ∗is ρ.
The proof of the above result establishes ﬁrst radial limits of the form limn B(0, nx)/n
and patches them together using shape theorem arguments. A key tool is that the vector
ϱ is invariant under translating the edge-weights, Busemann increments, and geodesic
graphs. The mean of ϱ follows directly from the previous lemma, whereas the fact that
Hϱ is a supporting line for B follows from the statements (a) x · ϱ ≤1 for all x ∈B and
(b) z0 · ϱ = 1. The claim (a) is shown by noting that for x ∈B, one has
x · ϱ = lim
n
B(0, nx)
n
≤lim
n
T(0, nx)
n
= µ(x) ≤1.
On the other hand, using the mean of ϱ, one has
µ∗(ϱ · z0 ≤1) = 1 and Eµ∗(ϱ · z0) = ρ · z0 = 1.
This implies that ϱ · z0 = 1 with µ∗-probability one.
At this point, we may identify the asymptotics of the reconstructed Busemann func-
tion if ∂B is diﬀerentiable at z0. In this case, there is a unique supporting hyperplane
for B at z0, so this must be Hϱ. In other words, one has
116

<!-- page 117 -->

Corollary 5.12. Assume (I). If ∂B is diﬀerentiable at z0, then
µ∗(ϱ = ρ) = 1.
In the presence of diﬀerentiability, the asymptotics of the Busemann function are
given exactly by projection onto the line H used to create the function. This result gives
us a major piece needed to establish directional properties of geodesic rays sampled from
µ∗. To complete the proofs of the directional results, we simply need to argue as in the
last section, combining asymptotics of B with the shape theorem for T.
To do this, we consider the conﬁguration η ∈Ω3 sampled from µ∗and build a directed
graph from it. Speciﬁcally, we set G = G(η) to be the directed graph with vertex set
equal to Zd and with edge set equal to {⟨x, y⟩: η(⟨x, y⟩) = 1}. Simple properties of
geodesics carry over to the weak limit. For example: with µ∗-probability one, for all
x, y ∈Zd,
• Each directed path in G is a geodesic and from each x there is an inﬁnite self-
avoiding directed path.
• If x →y in G (there is a directed path from x to y in G), then f(x, y) = T(x, y).
• Under assumption (II), viewed as an undirected graph, G has no circuits.
• Each vertex x has out-degree at least 1 in G. Under assumption (II), the out-degree
is exactly 1, and thus emanating from each x is a unique inﬁnite directed path Γx.
Given these properties, the next result follows as in the last section.
Theorem 5.13. Assume (I). With µ∗probability 1, each directed inﬁnite path from 0 in
G is asymptotically directed in Sϱ = ∂B ∩Hϱ.
The “asymptotic shape” for the Busemann function has boundary which is the line
Hϱ, whereas the asymptotic shape for T is simply the limit shape for the FPP model.
The picture here is that where these two make contact are angles in which geodesic rays
are created. Therefore if B has an exposed point in Hϱ, then we can produce geodesic
rays with asymptotic directions. To use this theorem, one should know the exact line
Hϱ, so we impose a condition of diﬀerentiability. By pulling this result back to the space
Ω1, we can obtain Theorem 5.8.
The proof of the coalescence statement of Theorem 5.9 follows the argument of Licea-
Newman outlined in Section 4.2. The main diﬀerence is that one must be careful to
construct barrier events on the original space Ω1 and port them over to eΩ, and this
creates a considerable headache. We refer the reader to [58, Section 6]. As for absence
of inﬁnite backward paths, we have
Theorem 5.14. Assume (I) and (II) in dimension d = 2. With µ∗probability 1, the set
{y : 0 ∈Γy} is ﬁnite.
117

<!-- page 118 -->

Proof. The idea is a modiﬁcation of that of Burton-Keane [38], from their proof of unique-
ness of the inﬁnite occupied cluster in Bernoulli percolation. For x ∈Z2, deﬁne the event
Ax that in G, x has an inﬁnite forward path and two disjoint inﬁnite backward paths. A
counting argument shows that in the box B(n) = [−n, n]2, one has surely
#{x ∈B(n) : Ax occurs} ≤Cn,
as the number of points for which Ax occurs is bounded above by the number of points on
the boundary of B(n). On the other hand, if we show that P(Ax) > 0, then translation
invariance gives
Eµ∗#{x ∈B(n) : Ax occurs} = #B(n)µ∗(Ax),
and this would be a contradiction for large n.
Assume then that with positive probability, the set {y : 0 ∈Γy} is inﬁnite; we will
show that this implies µ∗(Ax) > 0. By directedness of paths in G, we may ﬁnd a line
L intersecting Z2 such that each Γy eventually stays on one ﬁxed side of L. Choose y
on the other side of L and follow Γy on the event that {z : y ∈Γz} is inﬁnite until its
last intersection with L. Such a point has an inﬁnite forward path in G that does not
intersect L except at its initial point, and has an inﬁnite backward path. By the ergodic
theorem, we can then ﬁnd x1, x2 ∈L such that with positive µ∗probability, each of x1
and x2 has an inﬁnite forward path that does not touch L except at the initial point and
has an inﬁnite backward path. Because all forward paths in G coalesce, we may choose x
to be the coalescence point of Γx1 and Γx2 and so Ax occurs. By translation invariance,
µ∗(Ax) > 0.
Combining the previous theorems, this ends the proof of Theorem 5.9.
Remark 5.15. Here we explain the relation between the above method and King-
man’s approach. The analogue of Kingman’s variables deﬁned in (2.16) and (2.17) would
be, respectively
fk = 1
k
k
X
j=1
(T(0, Hj) −T(e1, Hj)) = 1
k
k
X
j=1
Bj(0, e1)
and the replacement for
fk + Tfk + . . . + T n−1fk
would be
1
k
k
X
j=1
Bj(0, ne1).
(These last two are not quite equal, but the average of Busemann increments is more
natural in our setup.) Note that fk is similar to the projection of our averaged measure
µ∗
k on the Busemann ﬁeld coordinate, but it is one level lower (at the level of random
118

<!-- page 119 -->

variables rather than measures). To make this precise, suppose that fk converges almost
surely to some random variable f and that µ∗
k converges to a measure µ. Then one would
have
f = Eµ [B(0, e1) | (te)] ,
where we have conditioned on the edge-weight conﬁguration (te).
In other words, f
would be an average of the reconstructed Busemann function B(0, e1) over the additional
randomness we introduced into µ (the uniformly random hyperplane Hn).
In fact, exactly this approach is taken in Liggett’s proof of his improved subaddi-
tive ergodic theorem [129]. In that result, Kingman’s assumptions are weakened, and
the question of building a decomposition into an additive and nonnegative subadditive
process is handled using weak limits, instead of employing weak-∗compactness to ﬁnd
a limit point for the sequence (fk) in (L1)∗∗, as Kingman did. Precisely, for a subaddi-
tive ergodic sequence (Xm,n), Liggett deﬁned an independent uniform {1, . . . , n} random
variable Un and set
Y (n)
i
= X0,i+Un −X0,i+Un−1.
Note the distribution of Y (n)
i
is the same as that given by averaging the distributions of
Xi,i+j −Xi+1,i+j for j = 1, . . . , n. Then he sets (Y1, Y2, . . .) to be any subsequential limit
in distribution of the sequences (Y (n)
1
, Y (n)
2
, . . .) and uses the distributional monotonicity
(Y1, Y1 + Y2, Y1 + Y2 + Y3 . . .) ≤st (X0,1, X0,2, X0,3 . . .)
in place of Kingman’s f ≤X0,1. If one wants to recover a pointwise monotonicity, one
can then simply couple these sequences together with a standard coupling.
One could ask whether the approach used in this section could work by only consid-
ering quantities at the level of the random variable fk rather than at the level of µ∗
k. This
may not be possible in general translation-ergodic environments. Suppose that geodesics
from 0 to Ln and e1 to Ln have exactly two pairs of subsequential limiting geodesics
γ(1)
0 , γ(2)
0
and γ(1)
e1 , γ(2)
e1 . Then if γ(i)
0
and γ(i)
e1 coalesce for each ﬁxed i, one would expect
the Busemann functions constructed from these limiting geodesics to be diﬀerent. How-
ever it is still possible that the limit f of fk exists, assuming the subsequences on which
we took limits are regular enough. Then f is simply an average of these two Busemann
functions corresponding to the diﬀerent geodesics and for this reason, it will no longer
necessarily have property (a) from Lemma 5.2. That is, if only one of the two geodesics
passes over the edge from 0 to e1, then it need not be true that f = τ{0,e1}. This property
is essential for deriving directional properties of geodesics from their Busemann functions.
5.5
Busemann functions along boundaries in Z2
In one case, Busemann limits such as (5.2) can be shown to exist surely, and this implies
that limits exist almost surely for sequences of ﬁnite geodesics. We explain in this section
the results of Auﬃnger-Damron-Hanson [17], where these statements are proved when
119

<!-- page 120 -->

Busemann limits are taken to points on boundaries of subsets in Z2.
The existence
of a boundary will allow for a “paths-crossing” trick due to Alm-Wierman [10] (ﬁrst
introduced in [8] in a special case), and existence of limits follow from this.
For simplicity, the results below will be stated for the half plane with vertices
VH = {(x, y) ∈Z2 : y ≥0}
and edges EH = {{x, y} : x, y ∈VH}, although we will remark about extending them
to general subsets. For x, y ∈VH, let TH(x, y) be the minimum passage time among all
paths with all vertices in VH from x to y.
Theorem 5.16. Let xn = (n, 0) and let (τe) be any edge-weight conﬁguration in [0, ∞)EH.
For all x, y ∈VH, the Busemann limit to xn exists:
BH(x, y) = lim
n [TH(x, xn) −TH(y, xn)] .
(5.3)
As usual, existence of Busemann limits gives us quite a bit of control on geodesics.
Using this theorem, we can prove existence of limiting geodesic graphs. Formally, we
represent geodesic graphs as elements of a directed graph space, as in the last section.
Let⃗EH be the set of directed edges of VH⃗
EH = {⟨x, y⟩: x, y ∈VH}
and write η for an arbitrary element of {0, 1}⃗EH. Build the graph G = G(η) as follows:
the vertices are all the vertices of VH and an edge e ∈⃗EH is present in the graph if and
only if η(e) = 1. For a sequence of vertices (xn), we let ηn(e) = 1 if e = ⟨x, y⟩is in a
geodesic from some point to xn and TH(x, xn) ≥TH(y, xn); we then set Gn = G(ηn).
The graphs Gn converge to a graph G = G(η) if for each e ∈⃗EH, one has ηn(e) →η(e).
Theorem 5.17. Let P be a probability measure on [0, ∞)EH such that
P(∃geodesic between x, y for all x, y ∈VH) = 1.
(5.4)
Then with probability one, (Gn) converges to a graph G. Each directed path in G is a
geodesic.
Remark 5.18. The above two theorems are valid in the following more general
setting. Let V be any subset of Z2 that is connected and inﬁnite and has inﬁnite con-
nected complement.
(The reader can think of slit planes, sectors, etc.)
There is a
unique doubly inﬁnite path in Z2 + (1/2, 1/2) which separates V and V c.
Enumer-
ate the vertices in V that are adjacent to this dual path as . . . , x−1, x0, x1, . . .. If we
deﬁne passage times between vertices x, y in V to be TV (x, y), the minimum passage
time among all paths with all vertices in V connecting x and y, then the Busemann
limits BV (x, y) = limn [TV (x, xn) −TV (y, xn)] exist for all edge-weight conﬁgurations in
[0, ∞)EV , and limiting geodesic graphs exist almost surely for any measure on [0, ∞)EV
that has the geodesic property (5.4). Here, EV is the set of nearest-neighbor edges with
both endpoints in V .
120

<!-- page 121 -->

In the case of the half-plane, one can say more about the structure of the limiting
graph G. We will state the theorem, but we refer to [18] for the complete proof.
Theorem 5.19. Let P be a product measure on [0, ∞)EH with continuous marginals.
The limiting geodesic graph G from Theorem 5.17 satisﬁes the following almost surely.
1. Each vertex in VH has out-degree 1. Therefore from each x ∈VH emanates a unique
inﬁnite directed path Γx.
2. Viewed as an undirected graph, G has no circuits.
3. For each x ∈VH, the backward cluster Bx = {y ∈VH : x ∈Γy} is ﬁnite.
4. For all x, y ∈VH, Γx and Γy coalesce. That is, their edge symmetric diﬀerence is
ﬁnite.
Proof. We give the ideas of the proofs of the ﬁrst two items only. Busemann limits are
proved as a consequence of the “paths crossing” trick of Alm and Wierman [10]. We
begin by showing the limit for x = m1e1, y = m2e1 with m1 < m2 ∈Z. For simplicity,
let us assume that for each w, z ∈VH, there is a geodesic between w and z. If this is not
the case, then the geodesics can be replaced by paths that have passage time within ϵ of
the inﬁmum.
Let n2 > n1 > m2 and let σ1 be a geodesic from x to n1e1, with σ2 a geodesic from y
to n2e1. Note that by planarity, the paths σ1 and σ2 must share a vertex z. Deﬁne the
path ˆσ1 by traversing σ1 from x to z and then σ2 from z to n2e1. Deﬁne ˆσ2 by traversing
σ2 from y to z and then σ1 from z to n1e1. Then (see Figure 13),
TH(x, n1e1) + TH(y, n2e1) = TH(σ1) + TH(σ2) = TH(ˆσ1) + TH(ˆσ2)
≥TH(x, n2e1) + TH(y, n1e1).
x
y
n1e1
n2e1
z
Figure 13: “Paths crossing” trick: the red path σ1 and the blue path σ2 must intersect
at a point z.
Rearranging this, we obtain
TH(x, n1e1) −TH(y, n1e1) ≥TH(x, n2e1) −TH(y, n2e1).
121

<!-- page 122 -->

Therefore the sequence in (5.3) is monotone and the limit exists.
To generalize this oﬀthe axis, we take x ∈VH be of the form (x1, x2) with x2 > 0 and
x′ ∈VH. Let B be a set of vertices of the type [−k, k] × [0, k] for some k large enough
so that B contains both x and x′. The set V ′ := VH \ B is connected and inﬁnite with
complement that is also connected and inﬁnite. So we deﬁne for y, z ∈V ′ the restricted
passage time T ′(y, z) to be the minimum of passage time of all paths from y to z which
have only vertices in V ′. One may then repeat the above paths crossing argument to see
that if y, z are on ∂V ′, the inner boundary of V ′ (they have a Z2-neighbor which is in
B), then the limit
B′(y, z) = lim
n [T ′(y, xn) −T ′(z, xn)]
exists. Furthermore, for all such y, z, the sequence deﬁning B′(y, z) is monotone.
Now the idea is to use existence of Busemann limits in the restricted region to prove
it for x and any other point. The crucial point is that for large n,
T(x, xn) = min{T(x, y) + T ′(y, xn) : y ∈∂V ′}.
For z ∈∂V ′ ﬁxed and y ∈∂V ′ variable, we set ψn(y) = T(x, y) + T ′(y, xn) −T ′(z, xn).
We ﬁnd then that ψn(y) has a ﬁnite limit ψ(y) for each y and furthermore
lim
n [T(x, xn) −T ′(z, xn)] = min{ψ(y) : y ∈∂V ′}.
The same argument works if we replace x by x′, so limn [T(x, xn) −T(x′, xn)] exists.
To prove existence of limiting geodesic graphs, we take e = ⟨x, y⟩∈⃗EH and again let
B be a box of the form [−k, k] × [0, k] which is large enough to contain x and y. If N is
large enough so that xn /∈B for all n ≥N, then if we again put V ′ = VH \ B one has for
all y ∈∂V ′,
y is on a geodesic from x to xn ⇔ψn(y) = min{ψn(u) : u ∈∂V ′},
where ψn(y) was deﬁned as above:
T(x, y) + T ′(y, xn) −T ′(z, xn)
for the ﬁxed z ∈∂V ′. However, due to monotonicity of the diﬀerences T ′(y, xn)−T ′(z, xn)
in n as y ranges over ∂V ′, the set of minimizers mn of ψn is eventually constant for n
large. Thus, the set of points in ∂V ′ which are in geodesics from x to xn is eventually
constant in n.
The above work readily implies that Gn converges. Indeed, one can check that ηn(e) =
1 if and only if e is in a geodesic from x to a vertex of mn. Since mn is eventually constant,
the value of ηn(e) must also be eventually constant.
122

<!-- page 123 -->

5.6
Nonexistence of Bigeodesics in ﬁxed directions
Recently, Busemann functions have been used in last-passage percolation models by
Sepp¨al¨ainen, Rassoul-Agha, and Georgiou to prove analogues (and improved versions)
of the directional results from [59] in FPP. This is a big advance because most work on
LPP has assumed exponential or geometric weights, where random matrix and queueing
theory analysis can be used.
A main tool in the general LPP models is the directedness of paths, which allows one
to use various forms of the “paths crossing” argument of Alm and Wierman. Applied to
Busemann functions, this gives a certain monotonicity of directional Busemann functions.
In FPP, paths are not directed, and this creates a fundamental obstacle to obtaining
similar improvements.
In this section, we describe recent work of Damron-Hanson [60] which shows that
the LPP results are also valid in FPP. These theorems address the “uniqueness” issue in
Newman’s conjectures about inﬁnite geodesics, and consequently rule out bigeodesics in
ﬁxed directions. In that paper, the following theorem is proved for general passage-time
distributions considered in the ﬁrst work of Damron-Hanson [59], but we will focus again
on i.i.d. weights, assuming I and II from Section 5.4.
We will again make some assumptions on the limiting shape: let θ ∈[0, 2π) and let
vθ be the point on the boundary of the limit shape in direction θ. When vθ is a point of
diﬀerentiability of the boundary of the limit shape, let Lθ be the unique tangent line and
let Sθ be the sector of angles of contact of Lθ with ∂B. Let θ1 and θ2 be the endpoints
of Sθ.
Theorem 5.20 (Damron-Hanson [60]). Assume I, II, and that the limit shape boundary
is diﬀerentiable at θ, θ1, and θ2. The following hold with probability one.
1. For each x ∈Z2, there is an inﬁnite geodesic Γx that is directed in Sθ such that for
any (random) sequence (xn) of points directed in Sθ,
Γx = lim
n Γ(x, xn).
2. For each x, y ∈Z2, the geodesics Γx and Γy coalesce.
3. There are no bigeodesics with one end directed in Sθ.
Notice that this result reduces Licea-Newman’s “exceptional set” Dc (see Theo-
rem 4.20) of directions to be empty.
Furthermore, their result rules out bigeodesics
with both ends directed oﬀthe exceptional set, whereas the above theorem only requires
directedness of one end (in a sector).
From this result we obtain the following progress on the bigeodesic conjecture.
Corollary 5.21. Assume I, II and that the limit shape boundary is diﬀerentiable. Then
for each θ,
P(there is a bigeodesic with one end in direction θ) = 0.
123

<!-- page 124 -->

Remark 5.22. Since geodesic lines with ﬁxed directions cannot exist, one should ask
if inﬁnite geodesics are even required to have directions. One can show using planarity
and the results of [59] that if ∂B is diﬀerentiable and I and II hold, then the following
statements are true with probability one:
1. for all θ, there is an inﬁnite geodesic starting from 0 directed in Sθ and
2. every inﬁnite geodesic is directed in Sθ for some θ.
As a consequence of the theorem and arguments like those presented elsewhere in this
section, one has
Theorem 5.23. Assume I, II, and that the limit shape boundary is diﬀerentiable at θ, θ1,
and θ2. With probability one, for each x, y ∈Z2 and sequence (xn) directed in Sθ, the
limit
B(x, y) = lim
n [T(x, xn) −T(y, xn)]
exists. Furthermore, letting ρ be the unique vector in R2 such that {r ∈R2 : r · ρ = 1} is
the tangent line to B in direction θ, one has:
1. EB(0, x) = ρ · x for all x ∈Z2.
2. For each ϵ > 0, the set of x ∈Z2 such that |B(0, x) −ρ · x| > ϵ∥x∥1 is ﬁnite.
In these theorems, one may ask the degree to which the diﬀerentiability assumption
is necessary. The main obstacle to removing diﬀerentiability is the result of H¨aggstrom-
Meester (Remark 2.23 and Theorem 7.6), which we recall states that given any compact,
convex subset C of Rd that has the symmetries of Zd that ﬁx the origin, there exists a
translation-ergodic distribution of edge-weights whose limit shape is C. In particular,
there are models of FPP whose limit shapes are polygons. In these situations, it is rea-
sonable to believe that one can construct models in which the only inﬁnite geodesics
are ones which do not have asymptotic directions – they wander across the sectors cor-
responding to the sides of the limit shape. This leads one to consider directedness in
sectors. However, it is reasonable to expect some FPP models which have polygonal
limit shapes, and for which there are inﬁnite geodesics directed along the corners. In
this case, uniqueness of inﬁnite geodesics in sectors corresponding to these sides will not
hold. This is prevented by making a diﬀerentiability assumption.
Remark 5.24. The above results give progress toward solving the “BKS midpoint
problem” introduced in [23]. The question is: is it true that the probability that (n/2)e1
is in a geodesic from 0 to ne1 goes to zero at n →∞? By a translation, if the answer is
no, then
P (0 is in a geodesic from (−n/2)e1 to (n/2)e2 for inﬁnitely many n) > 0.
On this event, 0 is in a bigeodesic, and under the diﬀerentiability assumptions of The-
orem 5.20, one of its ends is directed in the sector corresponding to θ = 0. This is
impossible by the theorem.
124

<!-- page 125 -->

Sketch of proof of Theorem 5.20. The proof we give is more similar to the proof in LPP
[85] because, for ease of exposition, we will omit the many applications of the Jordan
curve theorem needed to deal with undirected paths.
For simplicity, we will take θ = π/2 and π > θ1 ≥θ ≥θ2 > 0. We will ﬁrst reduce
to the half-plane, where it is easier to order inﬁnite geodesics. We will skip this step,
although it is the most work, because it does not involve Busemann functions, but just
topological arguments. However we mention that this reduction is done using the fact
that every path that is directed in S must have at most ﬁnitely many intersections with
L0 = {(x1, 0) : x1 ∈Z}.
So consider the upper half-plane H with vertices VH and
nearest-neighbor edges EH, and deﬁne passage times TH(x, y) for vertices x, y ∈VH by
considering only paths in the upper half-plane. We will content ourselves with proving
the theorem for geodesics constructed in the upper half-plane using TH.
Let GH(x) be the set of inﬁnite geodesics starting from x in the half-plane which
are directed in S = Sθ. By ordering in the half-plane, for each x ∈L0, there is a left-
most inﬁnite geodesic in GH(x), written as ΓL
x,H. Similarly there is a right-most inﬁnite
geodesic written as ΓR
x,H. The main step is to argue that that
ΓL
0,H = ΓR
0,H with probability one.
(5.5)
Once we show this, then part 1 of the theorem follows from the following lemma, whose
proof uses the assumption that ∂B is diﬀerentiable at the endpoints θ1, θ2.
Lemma 5.25. With probability one, the following holds. If (xn) is any (random) sequence
of points directed in S, and Γ is a subsequential limit of half-plane geodesics from 0 to
xn, then Γ is directed in S.
The idea of the proof of this lemma is to use trapping. Since ∂B is diﬀerentiable at
θ1, one can ﬁnd inﬁnitely many extreme points of B which converge to vθ1 outside of
S. One then derives a half-plane version of the existence of inﬁnite geodesics directed in
sectors from [59], and uses this to construct inﬁnitely many geodesics which are directed
in distinct sectors whose endpoints converge to vθ1 from outside S. The subsequential
limit Γ cannot intersect any of these geodesics more than once, due to unique passage
times, and so it must be directed “to the right” of θ1. A similar argument works for θ2
and traps Γ in S.
If we assume (5.5) and use the previous lemma, we ﬁnd that any subsequential limit
of geodesics from 0 to a (random) sequence (xn) directed in S must itself be directed in
S, and must then sit between ΓL
0,H and ΓR
0,H. This means all three geodesics are equal
and this proves item 1.
Now we argue for (5.5). Similarly to the deﬁnition of left and right-most geodesics in
the half-plane, one can deﬁne left-most and right-most inﬁnite geodesics from x in the
full-plane; we write these as ΓL
x and ΓR
x . One can show that with positive probability,
ΓL
x,H = ΓL
x and similarly for ΓR
x,H = ΓR
x . This is argued using the fact that any path
that is directed in S and starts on L0 has a last intersection with L0. One can then use
125

<!-- page 126 -->

a Licea-Newman style argument to show that each ΓL
x and ΓL
y coalesce (and similarly
for right-most) and then deduce this same statement for the half-plane geodesics. This
allows one to deﬁne Busemann functions
B∗
H(x, y) = lim
n [TH(x, xn) −TH(y, xn)]
for ∗= L or R, where (xn) is the sequence of vertices on Γ∗
0,H. A similar statement holds
for full-plane geodesics and in this case, one may deﬁne the Busemann functions on all
of Z2.
By the ergodic theorem, if we set ∆H(x, y) to be BL
H(x, y) −BR
H(x, y) (and similarly
for ∆(x, y) as the diﬀerence of full-plane Busemann functions), then there is a c such
that
1
n∆H(0, ne1) →c and 1
n∆H(0, −ne1) →−c
almost surely. Furthermore, one can use equality of full-plane and half-plane geodesics
(with positive probability) to show that these limits exist for ∆(the full-plane diﬀerence)
and are equal to c and −c, respectively.
Now (5.5) follows from two results:
Lemma 5.26. Under our diﬀerentiability assumption on ∂B, c = 0.
Proof. For ∗= L, R the function x 7→EB∗(0, x) is linear, and so there is a vector ρ∗
such that
EB∗(x, y) = ρ∗· (y −x) for x, y ∈Z2.
Once again, one may upgrade this to a shape theorem: for each ϵ > 0,
P
 |B∗(0, x) −ρ∗· x| > ϵ∥x∥1 for ﬁnitely many x ∈Z2
= 1.
Label the last intersections of Γ∗
0 with L0, L1, L2, . . . as x∗
0, x∗
1, x∗
2, . . .. (Here, Lk = L0 +
ke2.) Then choose a subsequence (x∗
nk) such that x∗
nk/∥x∗
nk∥1 converges to some vector
z. By the shape theorem,
B∗(0, x∗
n)/∥x∗
n∥1 = T(0, x∗
n)/∥x∗
n∥1 →µ(z).
On the other hand, this converges to ρ∗· z, so
ρ∗· z = µ(z).
Furthermore, as usual, ρ∗·r ≤µ(r) for all r ∈R2, using the fact that Busemann functions
are bounded by the passage time. As z must have argument in S, one concludes that
L∗= {r ∈R2 : r · ρ∗= 1}
is a supporting line for the limit shape at z. By diﬀerentiability, there is only one such
supporting line, and so we conclude that EB∗(0, e1) = 0 for both ∗= L, R. The ergodic
theorem completes the proof.
126

<!-- page 127 -->

Lemma 5.27. With probability one, ∆H(0, e1) ≥0. If ΓL
0,H̸ = ΓR
0,H with positive proba-
bility, then c > 0.
Proof. If ΓL
0,H̸ = ΓR
0,H, then by coalescence, the same must be true for geodesics starting
from any point x on L0. Since the edge-weights are continuous, one then deduces that
∆(0, e1)̸ = 0. So if we can show that ∆(0, e1) ≥0 almost surely, one then has E∆(0, e1) >
0, and the ergodic theorem will prove the lemma.
To show the almost sure inequality, we use the “paths crossing” trick. We claim that
ΓR
0,H and ΓL
e1,H must share a vertex. If ΓL
0,H = ΓR
0,H, then by coalescence, this is obvious.
Otherwise if follows from the fact that ΓR
0,H coalesces with ΓR
e1,H, which is to the right
of ΓL
e1,H. Let y be a vertex on ΓL
e1,H beyond z and let y′ be a vertex on ΓR
0,H beyond z.
Then
TH(0, y′) + TH(e1, y) = TH(0, z) + TH(z, y′) + TH(e1, z) + TH(y, z)
≥TH(0, y) + TH(e1, y′).
Rearranging and taking y, y′ →∞along their respective geodesics, we obtain BL
H(0, e1) ≥
BR
H(0, e1).
Now that ΓL
0,H = ΓR
0,H almost surely, coalescence of left and right-most geodesics
implies item 2 of the theorem. For the last item, we use the argument of Theorem 5.14.
To recall, this result states that in general cases in which from each point x, there is an
inﬁnite geodesic Γx such that (a) for x, y ∈Z2, Γx and Γy coalesce, and (b) there is a line
L such that almost surely, each Γx intersects this line ﬁnitely often, the tree of geodesics
(Γx : x ∈Z2) cannot have inﬁnite backward paths. In other words, with probability one,
the set {y : 0 ∈Γy} is ﬁnite. If there were a bigeodesic with one end directed in S,
then by uniqueness, it would have to be the forward inﬁnite path of the tree of geodesics
directed in S. However this tree cannot have an inﬁnite backward path, which would be
a contradiction.
127

<!-- page 128 -->

6
Growth and infection models, competition inter-
face, and more
6.1
Eden Model and the limit shape in high dimensions
The Eden model was introduced by Eden [71] in 1961 with the goal of analyzing a simple
question in cell reproduction:
‘Starting from a single cell which may divide, its daughters cells divide again and
again, what are the structural properties of the resulting colony of cells and how do
various possible constraints eﬀect the architecture? ’
We assume that each cell is identical to every other cell, that each cell is connected
to at least one other cell and that each of their locations is speciﬁed by a node in a
lattice. Eden considered the lattice to be two-dimensional, and reasoned that the close
biological counterpart to this process was the growth of bacterial cells or tissue cultures
that are constrained from moving. One example given is Ulva latuca [37], the common sea
lettuce, that grows as a sheet of two cells thick and apparently grows only at its periphery.
Although sometimes (wrongly) attributed as a possible model for cancer, it is important
to point out that there is strong in vitro and in vivo evidence that Eden model does not
explain any possible tumor growth (see for instance [33] and the references therein).
The model can be described as follows. We start with a cell at the origin. Next, we
adjoin one of the four neighbor vertices with equal probability. This two-celled conﬁgu-
ration has six adjacent nodes. We attach one of these six neighbor vertices with equal
probability 1/6. We continue this procedure indeﬁnitely, each time adjoining a single
cell uniformly from the collection of neighbor vertices. The construction is similar in any
dimension.
The model just described can be seen as a time change of a site ﬁrst passage perco-
lation model in Zd with exponential mean one passage times (on the vertices). Indeed,
let A(0) = {0} and, for n ≥1, deﬁne A(n) as the subset of Zd after we attach n cells
(vertices) to A(0) according to the above construction. Let σ0 = 0 and
σk = inf

t ≥0 : B(t) contains k + 1 points of Zd

,
where B(t) := {x ∈Zd : T(0, x) ≤t} and the vertex weights have distribution F(t) =
P(τe ≤t) = 1 −e−t. Then, Richardson [141] noticed that

A(k) : k ≥1
	
and

B(σk) : k ≥1
	
have the same distribution.
(6.1)
One of the most important consequences of (6.1) is that, contrary to other vertex
weight distributions, the evolution of the ball B(t) when the weights are exponentially
distributed is Markovian. Precisely, in the Eden model, the process
B(σ1), B(σ2), B(σ3), . . .
128

<!-- page 129 -->

B(σ14)
B(σ15)
Figure 14: The evolution of the Eden model. On the left, B(σ14) in yellow and the pos-
sible 16 neighbor squares that can be added. We choose one (shaded in pink) uniformly
with probability 1/16. On the right, B(σ15) is shown after the attachment. Note that
the order in which we attached previous boxes does not matter for the decision, only the
knowledge of the boundary is needed to establish the next transition. The vertices in the
picture are vertices of the dual lattice Z2
∗and a box represents a vertex of the lattice Z2.
is Markovian on the space of possible shape conﬁgurations. Indeed, one only needs to
know the (boundary) of the ball at time σi to know its possible transitions at time σi+1.
The Markovian property of the growing ball allows diﬀerent tools that are not avail-
able for general passage times. This has been explored since [141]. The reader may think
that this powerful property would open doors to a better knowledge of the limit shape
and its properties. The truth is, however, distressing. As of today, there are almost
no results that successfully explored the Markovian nature of B to describe limit shape
and geodesic results. The state of the art of the Eden model is essentially the same as
for any FPP model. In Section 7.2.3, we discuss the last-passage analogue of the Eden
model, where through combinatorics and explicitly formulae (achieved through the RSK
correspondence), one indeed can say much more. We stress however that the methods
of that section are not available in FPP. Connections between ﬁrst-passage percolation
and growth processes are further explored in Section 6.2.
We end this section with the few theorems that explicitly use the Markovian nature
of the Eden model. The ﬁrst one is the determination of the growth of the time constant
as the dimension diverges.
129

<!-- page 130 -->

Theorem 6.1 (Dhar [67]). For exponential passage times,
lim
d→∞µ(d)
d
log d = 1
2
The next result is the only one which rules out certain shapes as possible limit shapes.
However it was only proved in astronomically large dimensions.
Theorem 6.2 (Kesten [113]). If d > 650, 000 then the limit shape of the Eden model is
not an Euclidean ball.
The theorem above was later improved by Couronn´e, Enriquez, Gerin [51] to d ≥
35. Eden and Richardson observed by simulation that the limit shape looks circular in
dimension two. However, one of the most fascinating long-standing questions remains:
Question 33. Show that for d ≥2, the limit shape of the Eden model is not an
Euclidean ball.
6.2
Growth and competition models
First-passage percolation is closely related to certain growth and competition models.
In fact, the original version of the shape theorem was proved by Richardson [141] for
a simpliﬁed growth model now known as the original 1-type Richardson model. In this
model on Zd, we suppose that the origin houses an infection at time t = 0 and all other
sites are healthy. At each subsequent time t = 1, 2, . . ., any healthy site with at least
one infected neighbor becomes infected independently of all other sites with probability
p ∈(0, 1). Richardson proved a shape theorem for the infected region at time t as t →∞
and believed, on the basis of computer simulations, that as p varies from 1 to 0, the limit
shape varies from a “diamond to a disk.” This model was shown to be equivalent to
an FPP model with i.i.d. weights on sites with a geometric distribution [70]. It is also
common to describe the infection as species that are trying to colonize all sites of the
space Zd.
We can build an infection/species model based on edge FPP in an analogous manner.
Let (τe) be a realization of passage times. We infect the origin at time 0 and the infection
spreads at unit speed across edges, taking time τe to cross the edge e. In the case that
the edge weights are exponential, the memoryless property of the distribution implies
that the growth of the infected/colonized region is equal in distribution to a time change
of a type of Eden process (discussed in the previous section). For general distributions,
this growth process is called a 1-type Richardson model.
Building on the above deﬁnition, we may describe ﬁrst-passage competition models,
ﬁrst introduced in the plane by H¨aggstr¨om and Pemantle [95] with two competing in-
fections (species). We will consider the same model in Zd with k ∈N ∪{∞} diﬀerent
species. Fix k vertices x1, . . . , xk ∈Zd and at time t = 0, infect site xi by an infection
of type i. Each species spreads at unit speed as before, taking time τe to cross an edge
130

<!-- page 131 -->

e. An uninhabited site is exclusively and permanently colonized by the ﬁrst species that
reaches it; that is, y ∈Zd is occupied at time t by the i-th species if τ(y, xi) ≤t and
τ(y, xi) < τ(y, xj) for all j̸ = i. If multiple xi’s try to colonize a site at the same time (an
event that happens with probability 0 if the passage time distribution is continuous), we
use some rule to break ties.
Now, consider the set colonized by the i-th species starting at vertex xi:
C(xi) = {y ∈Zd : y is eventually occupied by i}.
We are interested in describing the geometry of these random sets. In particular, much
of the research so far has been driven to answer the following questions:
1. Can multiple species simultaneously succeed in invading an inﬁnite subset of the
lattice? In other words, does coexistence occur?
2. Assuming coexistence of at least two species, what can we say about the boundary
interface? Does it have an asymptotic direction? A possible scaling limit?
3. What happens if one species moves with faster speed than the others? Can we
show that coexistence is impossible?
These questions will be discussed in the following sections.
6.3
Competition with same speed
Before describing previous work, let’s be precise and put the above questions on solid
ground.
We start with the ﬁrst question on coexistence.
Recall that ν denotes the
probability measure of the passage times. One says that ν admits coexistence of k species
if for some choice of initial sites x1, . . . , xk,
P

|C(xi)| = ∞for all i = 1, . . . , k

> 0.
(6.2)
The ﬁrst point to make is that the choice of the initial sites is not important. One just
needs to make sure that with positive probability, a new set of initial sites y1, y2, . . . , yk are
respectively infected by the initial sites x1, x2, . . . , xk. When the passage times are expo-
nentially distributed, this can be easily done using the Markov property (see [95, Proposi-
ton 1.1], when d = k = 2). When the passage times are not exponentially distributed,
this modiﬁcation argument follows the same lines as those described in Section 3 and
can be found in [79, Section 5 and 6]. One therefore has the following proposition:
Proposition 6.3. Assume that the measure ν is non-atomic and the supp ν = R+.
Suppose that x1, . . . , xk and y1, . . . , yk are points in Zd such that there exists disjoint
paths from xi to yi. Then
P

|C(xi)| = ∞for all i = 1, . . . , k

> 0 ⇔P

|C(yi)| = ∞for all i = 1, . . . , k

> 0.
131

<!-- page 132 -->

The proposition above deals with choosing the initial sites. One has a similar result
if starting with k diﬀerent infected sets instead of single points [63]. However, the key
and diﬃcult step to prove coexistence is to establish that there indeed exist sources such
that (6.2) occurs. Heuristically, the shape theorem gives the intuition that the larger the
distance between the sources, the harder it is for one infection to surround the other one.
If one hopes to prove that (6.2) holds, then the best strategy is to choose the sources far
from each other and hope that their infections will grow at least in diﬀerent directions.
This gives an indication that coexistence is strongly related to the existence of distinct
geodesic rays with diﬀerent asymptotic directions. We will explore this avenue in a few
moments. First, let us brieﬂy say what is the state of the art on Question 1 above.
When ν is the exponential distribution, H¨aggstr¨om and Pemantle [95] proved coexis-
tence of 2 species (one could also see [25] for a review of results on Richardson models,
focused on exponential passage times). Shortly thereafter, Garet and Marchand [79] and
Hoﬀman [101] independently extended these results to prove coexistence of 2 species for
a broad class of translation-invariant measures, including some non-i.i.d. ones. Later,
Hoﬀman [102] demonstrated coexistence of 4 species in d = 2 for a similarly broad class
of measures by establishing a relation with the number of sides of the limit shape in the
associated FPP. Tying the question of coexistence with the number of sides of the limit
shape is a very neat argument that we try to sketch below. Recall that the number of
sides of the limit shape Bν is equal to k if ∂Bν is a polygon with k sides. If ∂Bν is not a
polygon, the number of sides is deﬁned to be inﬁnity.
Theorem 6.4 (Hoﬀman [102]). Suppose that Bν is bounded, Eτ 2+η
e
< ∞for some η > 0,
and that ν is uniquely geodesic. If the number of sides of Bν is at least k then for any
ϵ > 0 there exist x1, . . . , xk such that
P

|C(xi)| = ∞for all i = 1, . . . , k

> 1 −ϵ.
(6.3)
Remark 6.5. J.-B. Gou´er´e has extended these methods to give a geometric descrip-
tion of invaded regions; see [88].
We now sketch the proof of the above theorem.
Step 1 - Moving to Busemann functions.
We start the proof by some geometric considerations. Let the set V consist of all
v ∈∂Bν such that there is a unique line Lv which is tangent to Bν through v. For such
a v, let w(v) be a unit vector parallel to Lv. Let Ln,v be the line through nv in the
direction of w(v). It is not diﬃcult to see that if the number of sides of ∂Bν is at least
k then one can ﬁnd points v1, . . . , vk ∈V such that the lines Lvi are distinct for all i.
Fix such vi’s from now on. The ﬁrst step is to establish the following lemma. Deﬁne the
Busemann-type function BS(x, y) for x, y ∈Z2 and S ⊂R2 as
BS(x, y) = T(x, S) −T(y, S).
132

<!-- page 133 -->

From the subadditivity of T and the deﬁnition of BS, we have BS(x, y) ≤T(x, y) and
BS(x, y) + BS(y, z) = BS(x, z).
(6.4)
These Busemann functions are useful because of the following fact.
Lemma 6.6. Assume that the number of sides of Bν is at least k and let vi, 1 ≤i ≤k
be deﬁned as above. If there exist ϵ > 0, and x1, . . . , xk such that
P

BLn,vi(xj, xi) > 0 ∀i̸ = j

≥1 −ϵ
(6.5)
for inﬁnitely many n, then (6.3) holds.
Proof. It is easy to be convinced that such lemma must be true. Recall that geometrically
BLn,vi(x, y) is roughly given by the orthogonal projection of the vector y −x on the line
orthogonal to Ln,vi. If for all j̸ = 1, BLn,v1(xj, x1) > 0 then roughly speaking x1 is closer
than xj is to the line Ln,v1 and the species 1 will infect in the direction orthogonal to
Ln,v1.
Precisely, if for all j̸ = 1
BLn,v1(xj, x1) > 0,
then one can ﬁnd zn ∈Ln,v1 such that T(zn, x1) < T(zn, xj) for all j̸ = 1. Thus the
species 1 will infect inﬁnitely many sites. Repeating the same argument for the other
species, we obtain exactly what we want.
Step 2 - The density argument. Now, to end the proof of Theorem 6.4, one needs to
construct points xi, 1 ≤i ≤k such that (6.5) holds. The idea is to choose xi = Mvi with
M large enough. The proof then boils down to proving the following stronger statement.
Lemma 6.7. There exists c > 0 such that for all ϵ > 0 there exists M such that the
density of n with
P

BLn,vi(Mvj, Mvi) > cM
∀i̸ = j

> 1 −ϵ
is at least 1 −ϵ.
Proof. We sketch the argument here. The idea is that BLn,vi(Mvj, Mvi) for large M
behaves as f(M(vj −vi)) where f is a linear functional with null space in direction w(vi)
such that f(vi) = µ(vi). Precisely, for each pair (vi, vj) with i̸ = j we will ﬁnd w ∈Z2
and c > 0 (depending on the angle between the points vi and vj) such that for all large
enough M the events:
BLn,vi(Mvj, Mw) > Mc(1 −ϵ)
(6.6)
and
BLn,vi(Mw, Mvi) > −Mϵ
(6.7)
133

<!-- page 134 -->

happen with probability 1 −ϵ in a set of values of n with density at least 1 −ϵ. Then
(6.4) implies
BLn,vi(Mvj, Mvi) > Mc(1 −2ϵ).
Establishing (6.6) and (6.7) was the central part Hoﬀman’s paper (see Lemma 4.4, 4.5
in [102]) and these have much in common with ﬁnding the mean of the Busemann function
(see Lemma 5.10). For (6.6), Hoﬀman developed his averaging trick, which was explained
in Lemma 5.4 in Section 5. The proof of (6.7) relies on the assumption that the limit
shape is diﬀerentiable at vi.
Coexistence of inﬁnitely many species is deﬁned similarly and clearly implies coexis-
tence of k species for any k. In [58], Damron and Hochman extended the above Theorem
to include the case k = ∞and used it to construct a ﬁrst example of a passage time
distribution that admits inﬁnite coexistence. The current state of the art is given by a
combination of Theorem 6.4 and Theorem 2.25. It establishes that inﬁnite coexistence
is not only possible but is necessary for all measures in Mp:
Theorem 6.8. If ν ∈Mp for p ∈[⃗pc, 1) then ν admits coexistence of inﬁnitely many
species.
6.4
Competition with diﬀerent growth speeds
If two species have two diﬀerent growth rules, with one faster than the other, one may
expect that non-coexistence occurs with probability 1. We try to make this precise in
this section. We ﬁrst consider the original case treated in [95], where the passage times
are exponentially distributed. We restrict ourselves to the case of 2 infections.
The model now has two parameters, denoted by λ1 and λ2, indicating the intensities
of the infections. We assume that λ1 < λ2 so species 1 spreads faster than species 2. The
infection rule is the same as in the last section. We start with one infected site of each
species. An uninfected site becomes type i infected at a rate proportional to the number
of type i infected nearest neighbors and then stays type i for the rest of its life.
The best result to date in this case is still from the original paper [96]. With one of
the intensities ﬁxed, coexistence has probability 0 for all but at most countably many
values of the other intensity. It is however expected that this should be true for all values
of λ1̸ = λ2.
Question 34. Show that for all values of λ1̸ = λ2 and all x1̸ = x2,
P

|C(x1)| = |C(x2)| = ∞

= 0.
134

<!-- page 135 -->

It is worth mentioning that at this point, we are still left in the embarrassing situation
where we know that for almost every value of λ2, coexistence does not occur, but we are
unable to exhibit any single value of λ2 such that coexistence actually does not happen.
A possible way to answer this conjecture is to show:
Question 35. Fix λ1. Show that the event {|C(x1)| = |C(x2)| = ∞} is monotone
in λ2.
Monotonicity is believed to hold on Zd, although there are examples of graphs where
such a statement does not hold [64]. As in Proposition 6.3, the initial placement of the
infections does not matter and one could start with any bounded initial sets. Certain
special cases where one of the initial sets is unbounded (a half-space against a point, for
instance) are treated in [65] and extended in [13].
As in Section 6.3, one can deﬁne competition with diﬀerent speeds using more general
passage times. This was studied by Garet and Marchand in [81], where they extended
the results of [95] when the passage time distributions for the two species are stochas-
tically comparable. Under certain conditions on the passage times, in this setting, they
established in any dimension d that if the slow species survives, the fast species cannot
occupy a very high density of space. For d = 2, they show that almost surely, one species
must ﬁnally occupy a set of full density in the plane while the other species occupies only
a set of null density.
6.5
The competition interface
Before ending the chapter, it is worth mentioning that in the past twenty years, there
has also been a growing interest in questions related to the asymptotic shape of infected
regions.
In the physics literature, Derrida and Dickman [66] describe simulations for the Eden
model in which two clusters grow into a vacant sector of the plane with angle θ. They
obtain values for the roughening exponents of the competition interface depending on the
angle. Two species competition with mutation and selection was also studied numerically
in the Eden model in a paper by Kuhr, Leisner and Frey [125] (see the references therein
for other variations of the model).
However, precise results are scarce, except for solvable models of Last Passage Perco-
lation, where there is a coupling between the competition interface and the second class
particle in TASEP [74,75,139]. One of the main diﬃculties (again!) in the case of FPP
is the lack of knowledge of the limit shape Bν.
One way to deﬁne the competition interface is as follows (see Figure 15). For a given
θ ∈(0, 2π), consider the sets
Blue = {(x, y) ∈Z2 : θ
2 < arctan(x/y) < π}
135

<!-- page 136 -->

1
1
2
2
Figure 15: Competition interface between two species. On the left, the starting conﬁgu-
ration, with initial angle θ > 180◦. In blue, the region infected by species 1 while, in red,
the region infected by species 2. In black, the random curve deﬁning the competition
interface.
and
Red = {(x, y) ∈Z2 : π < arctan(x/y) < 2π −θ
2}.
Blue and Red are the initial infected sites. We now color a vertex v in Z2 \ (Blue ∪Red)
blue if
T(v, Blue) < T(v, Red),
and red otherwise. If the passage time distribution is continuous, almost surely there will
be no ties and a color will be assigned to each vertex of Z2. The competition interface is
the unique path on the dual lattice that sits between the blue and red regions. Now for
v ∈Z2 \ (Blue ∪Red) deﬁne
p(v) = P(v is blue)
and let the width of the interface at absissa x be deﬁned as
ω(x) =
X
v:v=(x,y)
p(v)(1 −p(v)).
For x large, [66] predicted that the function ω(x) has power law behavior with exponent
ω(x) ∼xζ where
ζ =





1
3, if θ < 180◦
2
3, if θ = 180◦
1, if θ > 180◦.
136

<!-- page 137 -->

As far as we know, none of these questions were investigated except in the case of LPP
with exponential weights.
Some recent result in this exactly solvable case were also
obtained in [85].
137

<!-- page 138 -->

7
Related models and questions
In this section, we discuss a sequence of related questions that were not addressed in
the previous sections. They all share something in common: they are all variants of a
ﬁrst-passage percolation model. Here we try to summarize the most well-studied con-
nections and variants for the reader; however, we refer to the papers for their proofs and
methodologies.
7.1
The maximum ﬂow
The maximum ﬂow problem can be considered a higher-dimensional version of ﬁrst-
passage percolation, where the passage times τ(e) play the role of the maximal amount
of ﬂuid that can go through an edge e per unit of time. Consider an open bounded
connected subset Υ of Rd such that the boundary ∂Υ is piecewise of class C1. Let Λ1, Λ2
be two open disjoint subsets of ∂Υ. For simplicity, let’s consider the case where Υ is a
d-dimensional cube of diameter n and Λ1 and Λ2 are two opposite faces. Λ1 will play the
role of the source of the network, while Λ2 denotes the sink. A ﬂow g from Λ1 to Λ2 is
an assignment of non-negative values g(e) and a direction to each edge e in Υ such that
0 ≤g(e) ≤τe,
for all e
and such that for each vertex v ∈Υ \ (Λ1 ∪Λ2), the total inﬂow equals the total outﬂow;
that is,
X
e∈v+
g(e) =
X
e∈v−
g(e)
where v+ (v−) is the set of edges incident to v and directed towards v (away from v).
Furthermore, “inﬂow” should be equal to the “outﬂow”; that is,
X
v∈Λ1
X
x/∈Λ1,x∼v
g(⟨v, x⟩) =
X
z∈Λ2
X
x/∈Λ2,x∼z
g(⟨x, z⟩),
where ⟨v, x⟩denotes the oriented edge that starts at v and ends at x and x ∼v means
x and v are neighboring vertices. The strength of the ﬂow g is now deﬁned as
S(g) =
X
v∈Λ1
X
x/∈Λ1,x∼v
g(⟨v, x⟩).
Clearly, if we choose g(e) = 0 at all edges, we have a ﬂow with S(g) = 0. Intuitively, as
we increase the strength, it becomes harder to ﬁnd ﬂows. The main question is therefore:
Given the passage times, the source Λ1, the sink Λ2 and the network Υ, what is the
maximum strength ﬂow from Λ1 to Λ2? How does its strength scale as Λ1 and Λ2 are
moved far apart?
138

<!-- page 139 -->

Precisely, we deﬁne the max-ﬂow from Λ1 to Λ2 as
Φ(Λ1, Λ2) =
max
g:g ﬂow from Λ1 to Λ2 S(g).
Our goal is to study the behavior of Φ(Λ1, Λ2) as n (the diameter of the d dimensional
cube) goes to inﬁnity.
The connection with ﬁrst-passage percolation is not only in the construction of the
network with the passage times τe. The maximum ﬂow can be represented in a diﬀerent
way, as the next theorem shows. We need a few deﬁnitions. We say that a set of edges
E separates Λ1 from Λ2 if there is no path in Υ \ E from Λ1 to Λ2. We call E a cut of
(Λ1, Λ2) if E separates Λ1 from Λ2 and if no proper subset of E separates Λ1 from Λ2.
Now the max-ﬂow and min-cut theorem (see [29]) states that
Theorem 7.1. The following equality holds:
Φ(Λ1, Λ2) = min
 X
e∈E
τe : E is a cut of (Λ1, Λ2)
	
.
Thus, the max-ﬂow can be viewed as the ﬁrst-passage time from Λ1 to Λ2 over all
possible cuts. In two dimensions, a cut E reduces to a path and we can clearly see the
connection. In higher dimensions, E is a hypersurface and thus the claim is that the
maximum ﬂow is the analogue of FPP in higher dimensions.
Now, let d = 3 and write Λ = [0, k] × [0, l] × [0, m]. The following result is due to
Kesten, [114].
Theorem 7.2 (Kesten [114]). Assume that m(k, l) →∞as k ≥l →∞and that for
some δ > 0,
k−1+δ log m(k, l) →∞.
(7.1)
Then there exists p∗, with 1/27 ≤p∗≤pc such that for any distribution F with F(0) < p∗
and ﬁnite exponential moments
lim
k,l,m→∞
1
klΦ(Λ1, Λ2) = ♥
a.s. and in L1
(7.2)
for some ♥> 0.
In the same paper, Kesten showed that if F(0) > 1 −pc the limit above exists and it
is equal to zero. The critical case when F(0) = 1−pc was also considered by Zhang [170].
Theorem 7.3 (Zhang [170]). If F(0) = 1−pc and F(0−) = 0 and if in addition Eτ < ∞
then for any l > 0,
lim
k,m→∞
1
klΦ(Λ1, Λ2) = 0
and
lim
k,l,m→∞
1
klΦ(Λ1, Λ2) = 0,
where k, l, m go to inﬁnity without any restriction (in particular wiithout assuming (7.1)).
139

<!-- page 140 -->

Remark 7.4. The arguments in [170] are also valid for dimensions d ≥2.
Assuming ﬁnite exponential moments, Zhang [170] improved the above results and
proved the following theorem.
Theorem 7.5 (Zhang [170]). Let d ≥2. Assume ﬁnite exponential moments and that
for all m = m(k) that goes to inﬁnity when all the length of the sides ki, i = 1, . . . , d −1
go to inﬁnity in such a way that
∃δ ∈(0, 1),
log m(k) ≤
max
i=1,...,d−1 k1−δ
i
.
Then we have
lim
k1,...,kd−1→∞
1
Qd−1
i=1 ki
Φ(Λ1, Λ2) = ♥(d)
a.s. and in L1.
Moreover, this limit is positive if and only if F(0) < 1 −pc(d).
After the work of Kesten and Zhang, remarkable progress has been made over a period
of years. Flows between general subsets of Zd were studied in the papers of Garet [78] and
Rossignol-Theret [144,145] in dimension 2. The study culminates in a series of papers of
Cerf and Th´eret [40–42]. Under ﬁnite exponential moments for the passage times, they
establish the analogue of Theorem 7.2 for very general networks and derive lower and
upper large deviations (of surface order nd−1, and volume order nd, respectively).
Question 36. Determine if the limit in (7.2) exists without the exponential moment
assumption.
7.2
Variants
7.2.1
Diﬀerent edge weights
Over the years, some attention has been devoted to cases of FPP with non-i.i.d edge
weights. Boivin [27] considers the general stationary and ergodic case and proved the
analogue of Theorem 2.1 and Theorem 2.16. (See Bjorklund [26] for the most general
known version.)
In the special case where the τe’s are also bounded, this was ﬁrst
noted by Derrienic (see [116, page 259]). A remarkable contribution came in the work
of H¨aggstr¨om and Meester [94], who classify all possible limit shapes in the stationary,
ergodic case:
Theorem 7.6 (Haggstrom-Meester [94]). Let C be the set of all subsets of Rd that are
compact, convex, and symmetric about the axes with nonempty interior. Let C∗be the
set of all compact subsets of Rd with nonempty interior that can arise as limiting shapes
for stationary ﬁrst-passage percolation. Then C = C∗.
140

<!-- page 141 -->

Figure 16: Representation of the limit shape in the left-right inhomogeneous model. On
the left, the edge distribution F+ is more variable than F−. On the right, the case where
no such relation is present.
The proof of the theorem above given by H¨aggstr¨om and Meester is constructive.
That is, given a set B ∈C∗, they construct a stationary measure P whose limit shape is
B. The measure P has bounded support and trivial tail σ-ﬁelds. For any set B, their
construction does not lead to an i.i.d. environment. The following open question is a
rephrasing of Question 1.
Question 37. Find a convex, compact set B with non-empty interior and symmetry
about the axes, not equal to the ℓ1 unit ball, and an i.i.d. ﬁrst passage model that has
B as its limit shape.
Other types of dependent passage times were considered in the literature. Fontes and
Newman [76] and Chayes and Winﬁeld [48] considered dependent ﬁrst-passage weights
coming from random colorings of vertices of Zd. Van den Berg and Kiss [165] extended
the BKS sublinear variance theorem to a class of dependent weights that includes passage
times that correspond with the minimal number of sign changes in a subcritical Ising
landscape, a model investigated by Higuchi and Zhang [100].
The following variation was also considered recently in [2]. Consider an inhomoge-
neous version of FPP where edges in the left and right half-planes of the Z2 lattice are
assigned weights independently according to distributions F−and F+. Complementing
the classical approach of the proof of Theorem 2.16 with large deviation estimates for
half-plane passage times as in [1], the authors proved a shape theorem. They also show
how the limit shape behaves as a function of the limit shapes of regular (homogeneous)
FPP with weights from F−and F+. If either of F−and F+ is more variable than the other,
then the asymptotic shape is the convex hull of the restriction to respective half-planes of
the asymptotic shapes for F−and F+. When no such relation is present, the asymptotic
shape equals the convex hull of the two half-shapes and a potentially wider additional
141

<!-- page 142 -->

Figure 17: Examples of directed paths from the origin to (n, n). These paths are used in
oriented ﬁrst passage and last passage percolation models.
line segment along the vertical axis (see Figure 7.2.1). The question introduced in that
work is:
Question 38. Find necessary and suﬃcient conditions on F−and F+ so that the
asymptotic speed in the vertical direction is strictly greater than that given by either F−
or F+ in the standard homogeneous model.
7.2.2
FPP on diﬀerent graphs
FPP on graphs diﬀerent from Zd has also appeared in the recent literature. In the past
decade, several authors were interested in FPP on complete graphs, on the hypercube
and also on random graphs such as Erdos-Renyi and the conﬁguration model. As most of
these examples are ﬁnite graphs, the questions are slightly diﬀerent from those considered
here. We describe them below.
Given a ﬁnite graph G = (V, E) and weights (τe)e∈E, the passage time between two
vertices u and v is deﬁned analogously to (1.1) as the minimum passage time among
paths between u and v. A few new and old quantities emerge as the main object of
study. First, one would like to understand the typical passage time and the number of
edges of a geodesic between two vertices. That is, choosing U and V uniformly at random
from the vertices of G, we want to understand the behavior of the passage time T(U, V )
and the length of a geodesic |Γ(U, V )|. Further objects of interest include the ﬂooding
time maxv∈V T(u, v) and the the weight diameter maxu,v∈V T(u, v).
For results in this direction we refer the reader to the beautiful notes of Remco van
der Hofstad [166, Chapter 8] and the references therein.
142

<!-- page 143 -->

7.2.3
Last Passage Percolation
In last-passage percolation (LPP), paths are constrained to have non-decreasing coordi-
nates. Usually the weights are placed on vertices instead of the edges and the passage
time between two points is given by
L(u, v) = max
Γ:u→v T(Γ) .
Existence of time constants and the limit shape theorem are proved as in Section 2. The
proof follows the same strategy with the diﬀerence that the process is now superadditive
and Bν is concave instead of convex. The directness of paths, however, requires that
vertices u, v, w be ordered in a non-decreasing fashion so that L(u, w) ≥L(u, v)+L(v, w).
The details for the proof of the limit shape can be found in [133], including a necessary
extra step to handle the continuity of Bν up to the boundary of the upper right quadrant.
Here, we will discuss the model in two dimensions and in this subsection, we will use the
notation µ(x, y) = µ(xe1 + ye2) for the time constant.
The importance of this small variation of FPP is that in dimension 2 there are very
natural correspondences between LPP models and certain queueing networks, systems
of queues in tandem, and interacting particle systems. This connection reaches a deeper
level as very precise scaling laws have been obtained for very particular distributions. If
the passage times are exponentially distributed with mean 1, then the time constant can
be explicitly computed:
µ(x, y) = (√x + √y)2,
as it was ﬁrst shown by Rost [143]. If ν is the probability measure of a geometric random
variable with parameter p then [108]
µ(x, y) = 1
p(x + y + 2
p
xy(1 −p)) .
In both cases, ﬁner asymptotics are available [108], as the distribution of
L(0, n(x, y)) −nµ(x, y)
n1/3
(7.3)
converges to a non-degenerate limit as n goes to inﬁnity. The proof of (7.3) goes through
the following special identity that identiﬁes the law of the passage time with the law of
the largest eigenvalue of the Laguerre Unitary ensemble. Let A be an n × n matrix with
entries that are complex Gaussian random variables with mean zero and variance 1/2.
Theorem 7.7 (Johansson [108]). If ν is exponentially distributed with mean 1 then for
all t > 0,
P(L(0, (n, n)) ≤t) = P(λn ≤t),
(7.4)
where λn is the largest eigenvalue of the n × n matrix AA∗.
143

<!-- page 144 -->

The identity above is achieved by considering independent geometric passage times
with parameter p ﬁrst. In this case, one can condition on the sum of the weights in the
box [0, N −1] × [0, N −1] to discover:
1. The conditional law is uniform on the space MN(s) of all possible non-negative
integer valued N × N matrices with ﬁxed sum s.
Indeed, letting A(i, j) be the passage time attached to the vertex (i, j), we have for the
matrix A = (A(i, j)),
P

A = (aij)

=
Y
1≤i,j≤N
(1 −p)aijp = pN2(1 −p)
P
ij aij.
2. The space MN(s) is in one-to-one correspondence with the space of generalized
permutations of length s on the set {1, . . . , N}. The last-passage time is exactly
the length of the longest non-decreasing subsequence.
This allowed Johansson to use the RSK correspondence. At the end of the day, it is
possible to connect the probability of the event L(0, (n, n)) ≤t to the number of pairs of
semi-standard Young tableaux with the size of the ﬁrst row. Determining the distribution
of the passage times is thus a combinatorial problem and it has an explicit expression
for any ﬁnite N. By approximating an exponential random variable by geometric, one
derives Theorem 7.7. We refer the reader to the very readable Section 2 of [108] for the
details.
The law of λn is well-known to be explicit and amenable to asymptotic analysis
through the asymptotics of Laguerre orthogonal polynomials (see [11] and the references
therein). In particular, (7.3) is a simple combination of the theorem above and the fact
that [108]
Zn := λn −4n
24/3n1/3
D→W2 ∼F2,
where
F2(s) = exp

−
Z ∞
s
(x −s)q(x)2dx

,
with q the solution of the Painlev´e II diﬀerential equation
q′′(x) = xq(x) + 2q3(x)
q(x) ∼Ai(x) as x →+∞,
and Ai(x) denotes the Airy function. The F2 distribution is called the Tracy-Widom
distribution due to the famous paper [156].
In the case of geometric passage times,
the asymptotic analysis was also carried out in [108] using the asymptotics of Meixner
polynomials. Once proper rescaled the passage times ﬂuctuate according to F2.
144

<!-- page 145 -->

As mentioned before, the explicit formula for the limit shape and the connection with
TASEP was also explored to understand the asymptotic behavior of inﬁnite geodesics
and the competition interface. The reader interested in this direction should look at the
work of Ferrari, Martin and Pimentel [74] and the references therein.
For passage times that are not geometric or exponentials, the results in LPP are
somewhat limited to same stage as those in FPP. Some recent progress was obtained in
the sequence of papers [85–87] where, using results from queueing theory, the authors
manage to (1) derive variational formulas for the time constant, similar to those in Section
2.2, (2) obtain asymptotic results for inﬁnite geodesics and Busemann functions under
mild assumptions on the limit shape.
145

<!-- page 146 -->

8
Summary of open questions
Here, just for convenience, we catalogue all questions that appear in the manuscript.
They are certainly not exhaustive, but we believe that this set provides a good scope of
the open and important areas of the ﬁeld.
Question 1. Find a non-trivial explicit distribution for which we can actually de-
termine µ(e1).
Question 2. Suppose that the support of the distribution of τe equals R+. Let
Xn = max{τe : e is in a geodesic from 0 to ne1}.
How does Xn scale with n?
Question 3. Extend the comparision theorem to the case d > 2, r > 0 and⃗pc(d) ≤
F(r).
Question 4. Show that if F is a continuous distribution then the limit shape is
strictly convex.
Question 5. Show that the d-dimensional cube is not a possible limit shape for a
FPP model with independent, identically distributed passage times.
Question 6. Show that for any measure ν ∈Mp the boundary of the limit shape is
not ﬂat outside the percolation cone.
Question 7. Show that if the limit shape of a measure ν has a ﬂat piece then ν ∈Mp,
for p > pc and the ﬂat piece is delimited by the percolation cone.
Question 8. Show that for any measure ν, in direction e1 the boundary of the
boundary of the limit shape does not contain any segment parallel to the e2 axis.
Question 9. Find a non-trivial example of a measure ν such that in direction e1 the
limit shape boundary is not parallel to the e2 axis.
Question 10. Find conditions that guarantee the existence of the limit (2.14).
Question 11. Show that for continuous distributions of passage times, the limit
shape is uniformly curved.
Question 12. Assume F(0) = 0. Show that µ ≥E4τe.
Question 13. Show that for any d ≥2, under suitable conditions on F, χ < 1/2. If
d = 2, show that χ = 1/3.
Question 14. Show that
lim
d→∞χ(d) = 0.
146

<!-- page 147 -->

Question 15. For suitable d, show that χ > 0.
Question 16. Show that in d = 2 there exists a > 0 such that the expected overlap
of two geodesics from 0 to nx sampled independently is at least na.
Question 17. Show that χ2 = χ2+δ for some δ > 0.
Question 18. Assume that the passage times are bounded and not concentrated at
a single point. Show that there exists a convex function Iu(ϵ) such that the following
limit holds:
lim
n
1
nd log pu
n(ϵ) = Iu(ϵ).
Question 19. Let d = 3 and P(τe = 0) = pc = 1 −P(τe = 1). Is it true that
ET(0, ∂B(n)) ≍log n?
Question 20. Is (3.40) optimal? More explicitly, let
γF(d) := sup

α : (3.41) holds for hn = nα

.
What is the value of γF(d)?
Question 21. Prove that under no assumptions on the distribution F of passage
times, geodesics exist in the d dimensional lattice, d ≥3.
Question 22. Find a distribution of passage times where (4.19) holds.
Question 23. Show that for continuous distributions, there are inﬁnitely many
geodesic rays.
Question 24. Find a distribution that is not in Mp that has an inﬁnite number of
geodesic rays.
Question 25. Is it true that the sequence of geodesics Γ(0, ne1) converges to some
Γ?
Question 26. Assume that the space is uniquely geodesic. Show that in d = 2, there
is no ﬁxed vertex that is in inﬁnitely many of the geodesics between −ne1 and ne1.
Question 27. Show that D = [0, 2π).
Question 28. Show that 0 ∈D.
Question 29. Decide whether Theorem 4.22 holds in arbitrary dimension.
147

<!-- page 148 -->

Question 30. Let Γ(v) be a geodesic from 0 to v. Consider the collection of all edges
that belong to an inﬁnite number of Γ(v) and let f(r) be the number of such edges which
intersect the circle x2 + y2 = r2. Does f(r) →∞as r →∞; and, if so, how fast?
Question 31. Do geodesic lines exist?
Question 32. Prove that under some conditions on the edge weights the limit (5.2)
exists.
Question 33. Show that for d ≥2, the limit shape of the Eden model is not an
Euclidean ball.
Question 34. Show that for all values of λ1̸ = λ2 and all x1̸ = x2,
P

|C(x1)| = |C(x2)| = ∞

= 0.
Question 35. Fix λ1. Show that the event {|C(x1)| = |C(x2)| = ∞} is monotone
in λ2.
Question 36. Determine if the limit in (7.2) exists without the exponential moment
assumption.
Question 37. Find a convex, compact set B with non-empty interior and symmetry
about the axes, not equal to the ℓ1 unit ball, and an i.i.d. ﬁrst passage model that has
B as its limit shape.
Question 38. Find necessary and suﬃcient conditions on F−and F+ so that the
asymptotic speed in the vertical direction is strictly greater than that given by either F−
or F+ in the standard homogeneous model.
148

<!-- page 149 -->

References
[1] Ahlberg, D. (2013). Convergence towards an asymptotic shape in ﬁrst-passage perco-
lation on cone-like subgraphs of the integer lattice. To appear in J. Theoret. Probab..
[2] Ahlberg, D., Damron, M. and Sidoravicius, V. (2013). Inhomogeneous ﬁrst-passage
percolation. arXiv:1311.4058.
[3] Aizenman, M. (1997). On the Number of Incipient Spanning Clusters. Nuclear Physics
B, 485, 551.
[4] Alexander, K. (1997). Approximation of subadditive functions and convergence rates
in limiting-shape results. Annals of Probab., 1, 30–55.
[5] Alexander, K. (1993). A Note on Some Rates of Convergence in First-Passage Per-
colation. Ann. Appl. Probab. ,3 (1), 81–90.
[6] Alexander, K. (2011). Subgaussian rates of convergence of means in directed ﬁrst
passage percolation. arXiv: 1101.1549.
[7] Alexander, K. and Zygouras, N. (2013). Subgaussian concentration and rates of con-
vergence in directed polymers. Electronic Journal of Probability 18, 1–28.
[8] Alm, S. E. (1998). A Note on a Problem by Welsh in First-Passage Percolation.
Combin. Probab. Comput.7 (1), 11–15.
[9] Alm, S. E. and Deifjen, M. (2014). First passage percolation on Z2 - a simulation
study. arXiv:1412.5924.
[10] Alm, S. E. and Wierman, J. C. (1999). Inequalities for means of restricted ﬁrst-
passage times in percolation theory. Combin. Probab. Comput., 8, 307–315.
[11] Anderson, G., Guionnet, A. and Zeitouni, O. (2010) An introduction to random
matrices. Cambridge studies in advanced mathematics, 118.
[12] Antal, P. and Pisztora, A. (1996). On the chemical distance for supercritical
Bernoulli percolation. Annals of Probab., 24, 1036–1048.
[13] Antunovic, T. and Procaccia, E. (2015). Stationary Eden model on groups. arXiv
1410.4944
[14] Auﬃnger, A. and Damron, M. (2013). The universal scaling relation for polymers in
a random environment and related models ALEA, Lat. Am. J. Probab. Math. Stat.,
10 (2), 857–880.
149

<!-- page 150 -->

[15] Auﬃnger, A. and Damron, M. (2014). A simpliﬁed proof of the relation between
scaling exponents in ﬁrst-passage percolation. Annals of Probab., 42, No. 3, 1197–
1211.
[16] Auﬃnger, A. and Damron, M. (2013). Diﬀerentiability at the edge of the limit
shape and related results in ﬁrst passage percoaltion. Probab. Theory and Rel. Fields,
156,193–227.
[17] Auﬃnger, A., Damron, M. and Hanson, J. (2013). Limiting geodesics for ﬁrst-
passage percolation on subsets of Z2, To appear in Annals of Applied Probability.
[18] Auﬃnger, A., Damron, M. and Hanson, J. (2014). Rate of convergence of the mean
for sub-additive ergodic sequences. To appear at Adv. in Math.
[19] Auﬃnger, A. and Louidor, O. (2011). Directed polymers in random environment
with heavy tails. Comm. on Pure and Applied Mathematics, 64, 183–204.
[20] Basdevant, A. Enriquez, N., Gerin, L., Gou´er´e, J.-B. (2014). The shape of large balls
in highly supercritical percolation. Elect. Journal Probab. 19, 1–14.
[21] Basdevant, A. Enriquez, N., Gerin, L. (2013). Distances in the highly supercritical
percolation cluster. Annals of Probab., 41, no. 6, 4342–4358.
[22] Benaim, M. and Rossignol, R. (2008). Exponential concentration for ﬁrst passage
percolation through modiﬁed Poincar´e inequalities. : Ann. Inst. H. Poincar´e Probab.
Statist. , 44, 544–573.
[23] Benjamini, I. Kalai, G. and Schramm, O. (2003). First passage percolation has
sublinear distance variance. Annals of Probab., 31, 1970–1978.
[24] Benjamini, I. and Tessera, R. (2014). First passage percolation on nilpotent Cayley
graphs and beyond. arXiv:1410.3292.
[25] Blair-Stahn, N. D. (2010). First passage percolation and competition models.
arXiv:1005.0649.
[26] Bjorklund, M. (2010). The asymptotic shape theorem for generalized ﬁrst passage
percolation. Annals of Probab., 38, 632–660.
[27] Boivin, D. (1990). First passage percolation: the stationary case. Probab. Theory
Relat. Fields. 86 491–499.
[28] Boivin, D. and Derrien, J.-M. (2002). Geodesics and recurrent of random walks in
disordered systems. Ann. Probab. 7, 101-115.
[29] Bollob´as, B. (1979). Graph theory. Graduate Texts in Mathematics. Springer- Ver-
lag, New York. An introductory course.
150

<!-- page 151 -->

[30] Bollob´as, B. and Riordan, O. Percolation. Cambridge University Press, New York,
2006. x+323 pp.
[31] Boucheron, S, Lugosi, G., and Massart, P. Concentration inequalities. A nonasymp-
totic theory of independence. With a foreword by Michel Ledoux. Oxford University
Press, Oxford, 2013. x+481 pp.
[32] Bridson, M. and Haeﬂiger, A. Metric Spaces of Non-Positive Curvature, Grundlder
Math. Wiss. 319, Springer Verlag, 1999.
[33] Bru, A., Albertos, S. Subiza, J., Garcia, J., and Bru, I. The Universal Dynamics of
Tumor Growth. (2003). Biophys J., 85, 2948–2961.
[34] Busemann, H. The geometry of geodesics. Vol. 6. Dover Publications, 1985.
[35] Burago, D., Burago, Y. and Ivanov, S. (2001). A Course in Metric Geometry. Amer-
ican Mathematical Society.
[36] Burkholder, D. L. (1973). Discussion of Prof. Kingman’s paper. Annals of Probab.,
1, 900–902.
[37] Burrows, E.M. (1991). Seaweeds of the British Isles. London: Natural History Mu-
seum. ISBN 0-565-00981-8.
[38] Burton, R.M. and Keane, M. (1989). Density and uniqueness in percolation. Comm.
Math. Phys., 121, Number 3, 501–505.
[39] Cerf, R. and Th´eret, M. (2014). Weak shape theorem in ﬁrst passage percolation
with inﬁnite passage times. arXiv:1404.4539.
[40] Cerf, R. and Th´eret, M. (2011). Law of the large numbers for the maximal ﬂow
through a domain of Rd in ﬁrst passage percolation. Trans. Amer. Math. Soc., 363,
3665–3702.
[41] Cerf, R. and Th´eret, M. (2011). Lower large deviations for the maximal ﬂow through
a domain of Rd in ﬁrst passage percolation. Probab. Theory and Relat. Fields, 150
no. 3, 635–661.
[42] Cerf, R. and Th´eret, M. (2011). Upper large deviations for the maximal ﬂow through
a domain of Rd in ﬁrst passage percolation. Ann. of Appl. Probab., 21, No. 6, 2075–
2108.
[43] Chatterjee, S. (2013). The universal relation between scaling exponents in ﬁrst-
passage percolation. Ann. of Math. (2), 177, no. 2, 663–697.
[44] Chatterjee, S. and Dey, P. (2009). Central limit theorem for ﬁrst-passage percolation
time across thin cylinders. arXiv:0911.5702.
151

<!-- page 152 -->

[45] Chayes, L. (1991). On the critical behavior of the ﬁrst passage time in d ≥3.
Helvetica Phys. Acta, 67, 1055–1071.
[46] Chayes, J.T., Chayes, L. and Durrett, R. (1986). Critical behavior of the two-
dimensional ﬁrst passage time. Journal of Stat. Physics, 45, 933–951.
[47] Chayes, J. T., Chayes, L. and Durrett, R. (1987). Inhomogeneous percolation prob-
lems and incipient inﬁnite clusters. J. Phys. A, 20, 1521–1530.
[48] Chayes, L. and Winﬁeld, C. (1993). The density of interfaces: A new ﬁrst passage
problem. J. Appl. Prob., 30, 851–862.
[49] Coupier, D. (2011). Multiple geodesics with the same direction. Elect. Com. Probab.,
16, 517–527.
[50] Chow, Y. and Zhang, Y. (2003). Large deviations in ﬁrst-passage percolation. Annals
of Appl. Probab., 13, 1601–1614.
[51] Couronn´e, O., Enriquez, N. and Gerin, L. (2011) Construction of a short path in
high-dimensional ﬁrst passage percolation. Elect. Comm. in Probab., 16, 22–28.
[52] Cox, J. T. (1980). The time constant of ﬁrst-passage percolation on the square
lattice. Adv. Appl. Prob., 12, 864–879.
[53] Cox, J. T. and Durrett, R. (1981). Some limit theorems for percolation with neces-
sary and suﬃcient conditions. Annals of Probab., 9, 583–603.
[54] Cox, J. T. and Kesten, H. (1981). On the continuity of the time constant of ﬁrst-
passage percolation. J. Applied Prob., 18, 809–819.
[55] Cranston, M. Gauthier, D. and Mountford, T. S. (2009). On large deviation regimes
for random media models. Annals of Appl. Probab., 19, 826–862.
[56] Damron, M., Hanson, J. and Sosoe, P. Subdiﬀusive concentration in ﬁrst-passage
percolation. arXiv:1401.0917.
[57] Damron, M., Hanson, J. and Sosoe, P. Sublinear variance in ﬁrst-passage percolation
for general distributions. To appear in Probab. Theory Rel. Fields.
[58] Damron, M. and Hochman, M. (2013) Examples of non-polygonal limit shapes in
i.i.d. ﬁrst-passage percolation and inﬁnite coexistence in spatial growth models. An-
nals of Applied Probability, 23, No. 3, 1074–1085.
[59] Damron, M. and Hanson, J. (2014). Busemann functions and inﬁnite geodesics in
two-dimensional ﬁrst-passage percolation. Comm. Math. Phys., 325, 917–963.
152

<!-- page 153 -->

[60] Damron, M. and Hanson, J. (2015). Bigeodesics in ﬁrst-passage percolation. In
preparation.
[61] Damron, M. and Kubota, N. (2014). Gaussian concentration for the lower tail in
ﬁrst-passage percolation under low moments. arXiv:1406.3105.
[62] Damron, M., Lam, W.-K. and Wang, X. (2015). Asymptotics for 2D Critical First
Passage Percolation. arXiv:1505.07544.
[63] Deijfen, M. and H¨aggstr¨om, O. (2006). The initial conﬁguration is irrelevant for the
possibility of mutual unbounded growth in the two-type Richardson model. Combin.
Probab. Comput. 15 345–353.
[64] Deijfen, M. and H¨aggstr¨om, O. (2006). Nonmonotonic coexistence regions for the
two- type Richardson model on graphs. Electron. J. Probab., 11, 331–344.
[65] Deijfen, M. and H¨aggstr¨om, O. (2007). The two-type Richardson model with un-
bounded initial conﬁgurations. Ann. Appl. Probab., 17, 1639–1656.
[66] Derrida, B. and Dickman, R. (1991). On the interface between two growing Eden
clusters. Journal of Physics A: Mathematical and General, 24, 1385–1422.
[67] Dhar, D. (1986). Asymptotic shape of Eden clusters, On growth and form. Ed. H.
E. Stanley and N. Ostrowsky, Martinus Nijhoﬀ., 288–292.
[68] Durrett, R. (2005). Probability: theory and examples. Duxbury Advanced Series.
Third Edition.
[69] Durrett, R. (1984). Oriented percolation in two dimensions. Annals of Probab. 12,
999–1040.
[70] Durrett, R. and Liggett, T. (1981). The shape of the limit set in Richardson’s growth
model. Annals of Probab., 9, 186–193.
[71] Eden, M. A two-dimensional growth process. 1961 Proc. 4th Berkeley Sympos. Math.
Statist. and Prob., Vol. IV, 223–239, Univ. California Press, Berkeley, Calif.
[72] Efron, B. and Stein, C. (1981). The jackknife estimate of variance. Annals of Statis-
tics, 9, 586–596.
[73] Falik, D. and Samorodnitsky, A. (2007). Edge-isoperimetric inequalities and inﬂu-
ences. Combin. Probab. Comput., 16, 693–712.
[74] Ferrari, P. A., Martin, J. B, and Pimentel, L. (2009). A phase transition for compe-
tition interfaces. Ann. Appl. Probab., 19, 281–317.
153

<!-- page 154 -->

[75] Ferrari, P. A. and Pimentel, L. (2005). Competition interfaces and second class
particles. Annals of Probab., 33, 1235–1254.
[76] Fontes, L. and Newman, C.M. (1993). First passage percolation for random colorings
of Zd. Ann. Appl. Probab., 3, 746–762.
[77] Forgacs, G., Lipowsky R. and Nieuwenhuizen, T. (1991). The behavior of interfaces
in ordered and disordered systems. In Phase transitions and critical phenomena (C.
Domb and J. Lebowitz, eds.) 14 135–363. Academic, London.
[78] Garet, O. (2009). Capacitive Flows on a 2D random net. Ann. of Appl. Probab., 19,
641–660.
[79] Garet, O. and Marchand, R. (2005). Coexistence in two-type ﬁrst-passage percola-
tion models. Ann. Appl. Probab., 15(1A), 298–330.
[80] Garet, O. and Marchand, R. (2004). Asymptotic shape for the chemical distance
and ﬁrst-passage percolation on the inﬁnite Bernoulli cluster. ESAIM Probab. Stat.
8, 169–199.
[81] Garet, O. and Marchand, R. (2008). First-passage competition with diﬀerent speeds:
positive density for both species is impossible. Elect. Journal. of Probab., 13, 2118–
2159.
[82] Garet, O. and Marchand, R. (2014). Large deviations for the contact process in
random environment. Annals of Probab., 42, 1438–1479.
[83] Garet, O., Marchand, R., Procaccia, E. and Th´eret, M. Continuity of the time and
isoperimetric constants in supercritical percolation. arXiv:1512.00742.
[84] Garet, O., Gouere, J-B. and Marchand, R. (2013). The number of open paths in
oriented percolation. arXiv:1312:2571.
[85] Georgiou, N., Rassoul-Agha, F. and Sepp¨al¨ainen, T. (2015). Geodesics and the
competition interface for the corner growth model. arXiv:1510.00860.
[86] Georgiou, N., Rassoul-Agha, F. and Sepp¨al¨ainen, T. (2015). Stationary cocycles and
Busemann functions for the corner growth model. arXiv:1510.00859.
[87] Georgiou, N., Rassoul-Agha, F. and Sepp¨al¨ainen, T. (2013). Variational formulas
and cocycle solutions for directed polymer and percolation models. arXiv:1311.3016.
[88] Gou´er´e, J.-B. (2007). Shape of territories in some competing growth models. Ann.
Appl. Probab., 17(4), 1273–1305.
[89] Grimmett, G. (1999). Percolation. 2nd edition. Springer, Berlin.
154

<!-- page 155 -->

[90] Grimmett, G. and Kesten, H. (1984). First-passage percolation, network ﬂow and
electrical resistances Z..Wahrsch. Verw. Gebiete, 66, 335–366.
[91] Grimmett, G. and Kesten, H. Percolation since Saint-Flour. arXiv:1207.0373.
[92] Grimmett, G. and Marstrand, J. M. (1990). The supercritical phase of percolation
is well behaved. Proc. Roy. Soc. London Ser. A, 430, 439–457.
[93] Gross, L. (1975). Logarithmic Sobolev inequalities. Amer. J. Math., 97, 1061–1083.
[94] H¨aggstr¨om, O. and Meester, R. (1995). Asymptotic shapes for stationary ﬁrst pas-
sage percolation. Annals of Probab., 23, 1511–1522.
[95] H¨aggstr¨om, O. and Pemantle, R. (1998). First passage percolation and a model for
competing spatial growth. J. Appl. Probab., 35, 683–692.
[96] H¨aggstr¨om, O. and Pemantle, R. (2000). Absence of mutual unbounded growth for
almost all parameter values in the two-type Richardson model. Stochastic Process.
Appl., 90, 207–222.
[97] Hambly, B. and Martin, J. B. (2007). Heavy tails in last-passage percolation.Probab.
Th. Rel. Fields, 137, 227–275.
[98] Hammersley, J. and Welsh, D. First-passage percolation, subadditive processes,
stochastic networks, and generalized renewal theory. 1965 Proc. Internat. Res.
Semin., Statist. Lab., Univ. California, Berkeley, Calif., 61–110, Springer-Verlag,
New York.
[99] Hanson, J. (2013). Ph.D. thesis.
[100] Higuchi, Y. and Zhang, Y. (2000). On the speed of convergence for two-dimensional
ﬁrst passage Ising percolation. Annals of Probab., 28, 353–378.
[101] Hoﬀman, C. (2005). Coexistence for Richardson type competing spatial growth
models. Ann. Appl. Probab., 15(1B), 739–747.
[102] Hoﬀman, C. (2008). Geodesics in ﬁrst passage percolation. Ann. Appl. Probab., 18,
1944-1969.
[103] Howard, C. D. Models of ﬁrst-passage percolation. Probability on discrete struc-
tures, 125–173, Encyclopaedia Math. Sci., 110, Springer, Berlin, 2004.
[104] Howard, C. D. and Newman, C. M., From greedy lattice animals to Euclidean ﬁrst-
passage percolation. Perplexing Problems in Probability, (R. T. Durrett, ed.), Progr.
Probab., 44, Birkh¨auser, Boston, 107–119.
155

<!-- page 156 -->

[105] Howard, C. D. and Newman, C. M. (2001). Geodesics and spanning trees for Eu-
clidean ﬁrst-passage percolation. Annals of Probab., 29, 577–623.
[106] Huse, D. and Henley, C. (1985). Pinning and roughening of domain walls in Ising
systems due to random impurities. Phys. Rev. Lett., 54, 2708–2711.
[107] Huse, D. A., Henley, C. L., Fisher, D. S.: (1985). Huse, Henley, and Fisher respond.
Phys. Rev. Lett., 55, 2924–2924.
[108] Johansson, K. (2000). Shape ﬂuctuations and random matrices. Comm. Math.
Phys., 209, 437–476.
[109] Johansson, K. (2000). Transversal ﬂuctuations for increasing subsequences on the
plane. Probab. Theory and Rel. Fields, 116, 445–456.
[110] Kardar, M. (1985). Roughening by impurities at ﬁnite temperatures. Phys. Rev.
Lett., 55, 2923–2923.
[111] Kardar, M. and Parisi, G. and Zhang, Y. (1986). Dynamic scaling of growing
interfaces. Phys. Rev. Lett. 56 889–892.
[112] Kardar, M. and Zhang, Y.-C. (1987). Scaling of directed polymers in random media.
Phys. Rev. Lett. , 56, 2087–2090.
[113] Kesten, H. (1986). Aspects of ﬁrst passage percolation. ´Ecole d’´Et´e de Probabilit´es
de Saint Flour XIV, Lecture Notes in Mathematics, 1180, 125–264.
[114] Kesten, H. (1987). Surfaces with minimal random weight and maximal ﬂows: A
higher dimensional version of ﬁrst passage percolation, Illinois J. Math., 31, 99–166.
[115] Kesten, H. (1993). On the speed of convergence in ﬁrst-passage percolation. Ann.
Appl. Probab., 3, 296–338.
[116] Kesten, H. First-passage percolation. From classical to modern probability, 93–143,
Progr. Probab., 54, Birkh¨auser, Basel, 2003.
[117] Kesten, H. and Zhang, Y. (1996). A central limit theorem for critical ﬁrst-passage
percolation in two dimensions. Probab. Theory Relat. Fields, 107, 137–160.
[118] Kingman, J. F. C. (1968). The ergodic theory of sub additive stochastic processes.
J. Roy. Statist. Soc. Ser., 30, 499–510.
[119] Kingman, J. F. C. (1973). Subadditive Ergodic Theory. Annals of Probab., 1, 883–
899.
[120] Krishnan, A. (2013) Variational formula for the time-constant of ﬁrst-passage per-
colation I: Homogenization. arXiv:1311.0316
156

<!-- page 157 -->

[121] Krishnan, A. (2013). Ph.D. Thesis.
[122] Krug, J. (1987). Scaling relation for a growing surface. Phys. Rev. A, 36, 5465–
5466.
[123] Krug, J. and Spohn, H. Kinetic roughening of growing interfaces. in Solids far from
equilibrium. Edited by C. Godr`eche. Collection Al´ea-Saclay: Monographs and Texts
in Statistical Physics, 1. Cambridge University Press, Cambridge, 1992. xvi+588 pp.
ISBN: 0-521-41170-X
[124] Kubota, N. (2013). Rates of convergence in ﬁrst passage percolation with low mo-
ment conditions. arXiv:1306.5917.
[125] Kuhr, J.-T., Leisner, M. and Frey, E. (2011). Range expansion with mutation and
selection: dynamical phase transition in a two-species Eden model. New Journal of
Physics, 13.
[126] Lalley, S. (2003). Strict convexity of the limit shape in ﬁrst-passage percolation.
Electron. Comm. Probab., 8, 135–141.
[127] Licea, C. and Newman, C. M. (1996). Geodesics in two-dimensional ﬁrst-passage
percolation. Annals of Probab., 24, 399–410.
[128] Licea, C. and Newman, C. and Piza, M. (1996). Superdiﬀusivity in ﬁrst-passage
percolation. Probab. Theory Related Fields 106, 559–591.
[129] Liggett, T. (1985). An Improved Subadditive Ergodic Theorem. Annals of Probab.,
13, 4, 1279–1285.
[130] Lions, P.-L. and P. E. Souganidis. (2005). Hamilton-Jacobi equations in stationary
ergodic media. Communications in Partial Diﬀerential Equations, 30(1-3), 335–375.
[131] Lions, P.-L. and P. E. Souganidis. (2010). Stochastic homogenization of Hamilton-
Jacobi and viscous Hamilton-Jacobi equations with convex nonlinearities - revisited.
Communications in Mathematical Sciences, 8 (2), 627– 637.
[132] Marchand, R. (2002). Strict inequalities for the time constant in ﬁrst passage per-
colation. Ann. Appl. Probab., 12, 1001–1038.
[133] Martin, J. (2004). Limiting shape for directed percolation models. Annals of
Probab., 32, 2908–2937.
[134] Mourrat, J.-C. (2012). Lyapunov exponents, shape theorems and large deviations
for the random walk in random potential ALEA, Lat. Am. J. Probab. Math. Stat., 9,
165–209.
157

<!-- page 158 -->

[135] Newman, C. M. A surface view of ﬁrst-passage percolation. Proceedings of the
International Congress of Mathematicians, Vol. 1, 2 (Z¨urich, 1994), 1017–1023,
Birkh¨auser, Basel, 1995.
[136] Newman, C. M. Topics in disordered systems. Springer, 1997.
[137] Newman, C. M. and Piza, M. (1995). Divergence of shape ﬂuctuations in two
dimensions. Annals of Probab., 23, 977–1005.
[138] Pemantle, R. and Peres, Y. (1994). Planar ﬁrst-passage percolation times are not
tight. Probability and phase transition (Cambridge, 1993), 261-264, NATO Adv. Sci.
Inst. Ser. C Math. Phys. Sci., 420, Kluwer Acad. Publ., Dordrecht.
[139] Pimentel, L. P. R. (2007). Multitype shape theorems for ﬁrst passage percolation
models. Adv. in Appl. Probab., 39(1), 53–76.
[140] Rhee, W. (1995). On Rates of Convergence for Common Subsequences and First
Passage Time. Ann. Appl. Probab., 5, 44–48.
[141] Richardson, D. (1973). Random growth in a tessellation. Proc. Cambridge Philos.
Soc., 74, 515–528.
[142] Rossignol, R. (2006). Threshold for monotone symmetric properties through a log-
arithmic Sobolev inequality. Annals of Probab., 35, 1707–1725.
[143] Rost, H. (1981). Nonequilibrium behaviour of a many particle process: density
proﬁle and local equilibria. Z. Wahrsch. Verw. Gebiete, 58, 41–53.
[144] Rossignol, R. and Th´eret. (2010). Law of large numbers for the maximal ﬂow
through tilted cylinders in two-dimensional ﬁrst passage percolation. Stoch. Proc.
Appl., 120, 873–900.
[145] Rossignol, R. and Th´eret. (2013). Lower large deviations for the maximal ﬂow
through tilted cylinders in two-dimensional ﬁrst passage percolation. ESAIM: Prob-
ability and Statistics, 17, 70–104.
[146] Sepp¨al¨ainen,
T.
(2009).
Lecture
Notes
on
the
Corner
Growth
Model.
http://www.math.wisc.edu/ seppalai/cornergrowth-book/ajo.pdf.
[147] Sepp¨al¨ainen, T. (1998). Exact limiting shape for a simpliﬁed model of ﬁrst- passage
percolation on the plane. Annals of Probab., 26(3), 1232–1250.
[148] Smythe, R. and Wierman, J. C. (1977). First-Passage Percolation on the Square
Lattice. I. Advances in Applied Probability, 9, No. 1, 38–54.
[149] Smythe, R. and Wierman, J. C. (1978). First-Passage Percolation on the Square
Lattice. III. Advances in Applied Probability, 10, No. 1, 155–171.
158

<!-- page 159 -->

[150] Sodin, S. (2014). Positive temperature versions of two theorems on ﬁrst-passage
percolation. Geometric Aspects of Functional Analysis, 2116, Lecture Notes in Math-
ematics, 441–453.
[151] Steele, J. (1986). An Efron-Stein inequality for nonsymmetric statistics. Annals of
Statistics, 14, 753–758.
[152] Talagrand, M. (1994). On Russo’s approximate zero-one law. Annals of Probab.,
22, 1576–1587.
[153] Talagrand, M. (1995). Concentration of measure and isoperimetric inequalities in
product spaces. Publ. Math. I.H.E.S., 81, 73–205.
[154] Tessera, R. (2014). Speed of convergence in ﬁrst passage percolation and geodesicity
of the average distance. arXiv:1410.1701.
[155] Tim´ar, A. (2013). Boundary-connectivity via graph theory. Proceedings of the AMS,
141, 475–480.
[156] Tracy, C. A.; Widom, H. (1994). Level-spacing distributions and the Airy kernel,
Comm. in Math. Phys., 159, 151–174
[157] Wehr, J. (1997). On the number of inﬁnite geodesics and ground states in disordered
systems. J. Statist. Phys., 87, 439–447.
[158] Wehr, J. and Aizenman, M. (1987). Fluctuations of extensive functions of quenched
random couplings. J. Stat. Phys., 60, 287–306.
[159] Wehr, J. and Woo, J. (1998). Absence of geodesics in ﬁrst-passage percolation on
a half-plane. Annals of Probab., 26, 358–367.
[160] Wierman, J. C., (1980). Weak moment conditions for time coordinates in ﬁrst-
passage percolation models. Journal of Applied Probability, 17, 968–978.
[161] Wierman, J. C. and Reh, W. (1978). On conjectures in ﬁrst passage percolation
theory. Annals of Probab., 6, 388–397.
[162] Wolf, D. E. and Kertesz, J. (1987). Noise reduction in Eden models: I. J. Phys. A,
20, L257-L261.
[163] Wuthrich, M. V. (1998). Scaling identity for crossing Brownian motion in a Posso-
nian potential. Probab. Theory Related Fields, 112, 299–319.
[164] Van den Berg, J. and Kesten, H. (1993). Inequalities for the time constant in
ﬁrst-passage percolation. Ann. Appl. Probab., 3, 56–80.
159

<!-- page 160 -->

[165] Van den Berg, J. and Kiss, D. (2012). Sublinearity of the travel-time variance for
dependent ﬁrst-passage percolation. Annals of Probab., 40, 743–764.
[166] Van der Hofstad, R. (2014). Random Graphs and Complex Networks. Vol. II. On-
line at http://www.win.tue.nl/˜rhofstad/NotesRGCNII.pdf.
[167] Yao, C.-L. (2014). Law of large numbers for critical ﬁrst-passage percolation on the
triangular lattice. Electronic Communications in Probability,19 1–14.
[168] Yao, C.-L. (2015). Limit theorems for critical ﬁrst-passage percolation on the tri-
angular lattice. arXiv: 1602.00065.
[169] Zabolitsky, J. G. and Stauﬀer, D. (1986). Simulation of large Eden clusters. Phys.
Rev. A 34 1523–1530.
[170] Zhang, Y. (2000). Critical behavior for maximal ﬂows on the cubic lattice. J. Stat.
Phys. 98, 799–811.
[171] Zhang, Y. (1999). Double behavior of critical ﬁrst-passage percolation. Perplexing
Problems in Probability, 143–158.
[172] Zhang, Y. (2010). On the concentration and the convergence rate with a moment
condition in ﬁrst-passage percolation. Stoch. Proc. Appl. 120 1317-1341.
[173] Zhang, Y. (2008). Shape ﬂuctuations are diﬀerent in diﬀerent directions. Annals of
Probab., 36, 331–362.
[174] Zhang, Y. (1995). Supercritical behaviors in ﬁrst-passage percolation. Stochastic
Proc. Appl. 59, 251–266.
160
