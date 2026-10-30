# Nonequilibrium Critical Phenomena and Phase Transitions into Absorbing States

> Adv. Phys. 49, 815 (2000) · [arXiv](https://arxiv.org/abs/cond-mat/0001070)
> Local copy: `docs/ref/Hinrichsen_2000_AdvPhys_absorbing_states.pdf` (153 pages)

<!-- page 1 -->

arXiv:cond-mat/0001070v2  [cond-mat.stat-mech]  22 May 2000
Nonequilibrium Critical Phenomena and
Phase Transitions into Absorbing States
Haye Hinrichsen 1
Theoretische Physik, Fachbereich 10
Gerhard-Mercator-Universit¨at Duisburg
D-47048 Duisburg, Germany
May 2000
Abstract:
This review addresses recent developments in nonequilibrium statistical physics. Focusing
on phase transitions from ﬂuctuating phases into absorbing states, the universality class
of directed percolation is investigated in detail. The survey gives a general introduction
to various lattice models of directed percolation and studies their scaling properties, ﬁeld-
theoretic aspects, numerical techniques, as well as possible experimental realizations. In
addition, several examples of absorbing-state transitions which do not belong to the di-
rected percolation universality class will be discussed. As a closely related technique, we
investigate the concept of damage spreading. It is shown that this technique is ambiguous
to some extent, making it impossible to deﬁne chaotic and regular phases in stochastic
nonequilibrium systems. Finally, we discuss various classes of depinning transitions in
models for interface growth which are related to phase transitions into absorbing states.
Keywords:
nonequilibrium phase transitions, stochastic lattice models, directed percolation, contact process, absor-
bing state transitions, parity-conserving class, damage spreading, interface growth, depinning transitions,
roughening transitions, nonequilibrium wetting
PACS: 05.70.Ln, 64.60.Ht, 64.60.Ak
submitted to Advances in Physics
1email: hinrichs@comphys.uni-duisburg.de

<!-- page 2 -->

Contents
1
Introduction
4
2
Stochastic many-particle systems
8
2.1
The one-dimensional random walk
. . . . . . . . . . . . . . . . . . . . . . .
8
2.2
The master equation . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
9
2.3
Diﬀusion of many particles: The asymmetric exclusion process
. . . . . . .
11
2.4
Reaction-diﬀusion processes . . . . . . . . . . . . . . . . . . . . . . . . . . .
13
2.5
Mean ﬁeld approximation . . . . . . . . . . . . . . . . . . . . . . . . . . . .
15
2.6
The inﬂuence of ﬂuctuations . . . . . . . . . . . . . . . . . . . . . . . . . . .
16
2.7
Numerical simulations . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
18
2.8
Exact results
. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
20
2.9
Experimental veriﬁcation of ﬂuctuation eﬀects . . . . . . . . . . . . . . . . .
24
2.10 Dynamic processes approaching thermal equilibrium . . . . . . . . . . . . .
25
2.11 Matrix product states
. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
27
3
Directed percolation
33
3.1
Directed percolation as a spreading process
. . . . . . . . . . . . . . . . . .
34
3.2
Lattice models for directed percolation . . . . . . . . . . . . . . . . . . . . .
37
3.3
Phenomenological scaling theory
. . . . . . . . . . . . . . . . . . . . . . . .
43
3.4
Estimation of the critical exponents
. . . . . . . . . . . . . . . . . . . . . .
50
3.5
Field-theoretic formulation of directed percolation
. . . . . . . . . . . . . .
57
3.6
Surface critical behavior . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
64
3.7
The inﬂuence of quenched disorder . . . . . . . . . . . . . . . . . . . . . . .
69
3.8
Related models . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
72
3.9
Experimental realizations of directed percolation . . . . . . . . . . . . . . .
79
4
Other classes of spreading transitions
87
4.1
Long-range spreading processes . . . . . . . . . . . . . . . . . . . . . . . . .
87
4.2
Absorbing-state transitions in systems with additional symmetries
. . . . .
92
4.3
Activated random walks . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
97
4.4
Absorbing phase transitions and self-organized criticality . . . . . . . . . . .
98
4.5
The annilation/ﬁssion process . . . . . . . . . . . . . . . . . . . . . . . . . .
99
4.6
The Lipowski model . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
99

<!-- page 3 -->

5
Damage spreading
101
5.1
Damage spreading phases . . . . . . . . . . . . . . . . . . . . . . . . . . . . 102
5.2
Universality of damage spreading transitions . . . . . . . . . . . . . . . . . . 107
5.3
Applications of damage spreading . . . . . . . . . . . . . . . . . . . . . . . . 109
5.4
Damage spreading and experiments . . . . . . . . . . . . . . . . . . . . . . . 111
6
Interface growth
112
6.1
Roughening transitions – a brief introduction . . . . . . . . . . . . . . . . . 112
6.2
Depinning transitions of driven interfaces
. . . . . . . . . . . . . . . . . . . 114
6.3
Polynuclear growth . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 116
6.4
Growth with evaporation at the edges of plateaus . . . . . . . . . . . . . . . 118
6.5
Unidirectionally coupled directed percolation processes . . . . . . . . . . . . 121
6.6
Parity-conserving roughening transitions . . . . . . . . . . . . . . . . . . . . 125
6.7
Nonequilibrium wetting transitions . . . . . . . . . . . . . . . . . . . . . . . 127
Appendices
132
A
Vector space notation and tensor products
. . . . . . . . . . . . . . . . . . 132
B
Derivation of the eﬀective action . . . . . . . . . . . . . . . . . . . . . . . . 132
C
Shell integration . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 133
D
One-loop integrals for directed percolation
. . . . . . . . . . . . . . . . . . 135
E
Notations . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 136
References
138

<!-- page 4 -->

1
Introduction
Random behavior is a common feature of complex physical systems. Although systems in
nature generally evolve according to well-known physical laws, it is in most cases impossible
to describe them by means of ab initio methods since details of the microscopic dynamics
are not fully known. Instead, it is often a good approximation to assume that the individual
degrees of freedom behave randomly according to certain probabilistic rules.
For this
reason methods of statistical mechanics become essential in order to study the physical
properties of complex systems.
In this approach a physical system is described by a
reduced set of dynamical variables while the remaining degrees of freedom are considered
as an eﬀective noise with a certain postulated distribution. The actual origin of the noise,
which may be related to chaotic motion, thermal interactions or even quantum-mechanical
ﬂuctuations, is usually ignored. Thus, statistical mechanics deals with stochastic models
of systems that are much more complicated in reality.
A complete description of a stochastic model is provided by the probability distri-
bution Pt(s) to ﬁnd the system at time t in a certain conﬁguration s. For systems at
thermal equilibrium this probability distribution is given by the stationary Gibbs ensem-
ble Peq(s) ∼e−H(s)/kBT , where H(s) denotes the microscopic Hamiltonian [1]. In principle,
the Gibbs ensemble allows us to compute the expectation value of any time-independent
observable by summing over all accessible conﬁgurations of the system. However, in most
cases it is very diﬃcult to perform the conﬁgurational sum. In fact, although numerous
exact solutions have been found [2], the vast majority of stochastic models cannot yet be
solved exactly. In order to investigate such nonintegrable systems, powerful approxima-
tion techniques such as series expansions [3] and renormalization group methods [4] have
been developed. Thus, in equilibrium statistical mechanics, we have a well-established
theoretical framework at our disposal.
From the physical point of view it is particularly interesting to investigate stochas-
tic systems in which the microscopic degrees of freedom behave collectively over large
scales [5,6]. Collective behavior of this kind is usually observed when the system undergoes
a continuous phase transition. The best known example is the order-disorder transition in
the two-dimensional Ising model, where the typical size of ordered domains diverges when
the critical temperature Tc is approached [7]. In most cases the emerging long-range cor-
relations are fully speciﬁed by the symmetry properties of the model under consideration
and do not depend on details of the microscopic interactions. This allows phase transi-
tions to be categorized into diﬀerent universality classes. The notion of universality was
originally introduced by experimentalists in order to describe the observation that several
apparently unrelated physical systems may be characterized by the same type of singular
behavior near the transition. Since then universality became a paradigm of the theory of
equilibrium critical phenomena. As the number of possible universality classes seems to
be limited, it would be an important theoretical task to provide a complete classiﬁcation
scheme, similar to the periodic table of elements.
The most remarkable breakthrough
in this respect was the application of conformal ﬁeld theory to equilibrium critical phe-
nomena [8–10], leading to a classiﬁcation scheme of continuous phase transitions in two
dimensions.
4

<!-- page 5 -->

In nature, however, thermal equilibrium is rather an exception than a rule. In most
cases the temporal evolution starts out from an initial state which is far away from equi-
librium.
The relaxation of such a system towards its stationary state depends on the
speciﬁc dynamical properties and cannot be described within the framework of equilib-
rium statistical mechanics. Instead it is necessary to deal with a probabilistic model for
the microscopic dynamics of the system. Assuming certain transition probabilities, the
time-dependent probability distribution Pt(s) has to be derived from a diﬀerential equa-
tion, the so-called Fokker-Planck or master equation. Nonequilibrium phenomena are also
encountered if an external current runs through the system, keeping it away from thermal
equilibrium [11]. A simple example of such a driven system is a resistor in an electric
circuit. Although the resistor eventually reaches a stationary state, its probability dis-
tribution will no longer be given by the Gibbs ensemble. As a physical consequence the
thermal noise produced by the resistor is no longer characterized by a Gaussian distri-
bution. Similar nonequilibrium phenomena are observed in catalytic reactions, surface
growth, and many other phenomena with a ﬂow of energy or particles through the system.
Since nonequilibrium systems do not require detailed balance, they exhibit a potentially
richer behavior than equilibrium systems. However, as their probability distribution can-
not be expressed solely in terms of an energy functional H(s), the master equation has to
be solved, being usually a much more diﬃcult task. Therefore, compared to equilibrium
statistical mechanics, the theoretical understanding of nonequilibrium processes is still at
its beginning.
The simplest nonequilibrium situation is encountered if a single or several particles in
a potential are subjected to a random force. Important examples are the Kramers and
the Smoluchowski equations describing the evolution of the probability distribution for
Brownian motion of classical particles in an external ﬁeld (for a review see [12]). But even
more complicated systems, for example one-dimensional tight-binding fermion systems as
well as electrical lines of random conductances or capacitances, can be described in terms
of discrete single-particle equations [13]. A more complex situation, on which we will focus
in the present work, emerges in stochastic lattice models with many interacting degrees
of freedom. A well-known example is the Glauber model [14] which describes the spin
relaxation of an Ising system towards the stationary state.
The corresponding master
equation in one dimension was solved exactly by Felderhof [15], who mapped the time
evolution operator onto a quantum spin chain Hamiltonian that can be treated by similar
methods as in Ref. [16].
One motivation for today’s interest in particle hopping models originates in the study
of superionic conductors in the 70’s [17,18]. In the superionic conductor AgI, for example,
the Ag+ ions may be viewed as particles moving stochastically through a lattice of I−
ions. Each lattice site can be occupied by at most one Ag+ ion, i.e., the particles obey an
exclusion principle. It was observed experimentally that the conductivity of AgI changes
abruptly when the temperature is increased, indicating an underlying order-disorder phase
transition of Ag+ ions. Assuming short-range interactions, this phase transition was ex-
plained in terms of a model for diﬀusing particles on a lattice [19]. Subsequently, particle
hopping models have been generalized to so-called reaction-diﬀusion models by including
various types of particle reactions or external driving forces [20]. It should be noted that
particles in a reaction-diﬀusion model do not always represent physical particles. More-
over, the reactions are not always of chemical nature. For example, in models for traﬃc
ﬂow individual cars are considered as interacting particles [21]. Similarly, electronic ex-
5

<!-- page 6 -->

citations of certain polymer chains may be viewed as particles subjected to a stochastic
temporal evolution [22].
The dynamic properties of a reaction-diﬀusion model on a lattice are fully speciﬁed
by its master equation [23–25]. In a few cases it is possible to solve the master equa-
tion exactly. During the last decade there has been an enormous progress in the ﬁeld
of exactly solvable nonequilibrium processes. This development was mainly triggered by
the observation that the Liouville operator of certain (1+1)-dimensional reaction-diﬀusion
models is related to Hamiltonians of previously known quantum spin systems. For exam-
ple, as ﬁrst realized by Alexander and Holstein [26], the symmetric exclusion process can
be mapped exactly onto the Schr¨odinger equation of a Heisenberg ferromagnet. This type
of mapping was extended to various other one-dimensional reaction-diﬀusion processes by
Alcaraz et al. [27], allowing exact methods of many-body quantum mechanics such as the
Bethe ansatz and free-fermion techniques to be applied in nonequilibrium physics [28].
Moreover, novel algebraic techniques have been developed in which the stationary state
of certain reaction-diﬀusion models is expressed in terms of products of non-commuting
algebraic objects [29].
In spite of this remarkable progress, the majority of reaction-diﬀusion models cannot
be solved exactly. It is therefore necessary to use approximation techniques in order to
describe their essential properties. The oldest approximation method is the law of mass ac-
tion, where the reaction rate of two reactants is assumed to be proportional to the product
of their concentrations. This mean ﬁeld approach is justiﬁed if diﬀusive mixing of particles
is much stronger than the inﬂuence of correlations produced by the reactions. Mean-ﬁeld
techniques have been applied successfully to a large variety of reaction-diﬀusion systems.
The study of pattern formation in nonlinear reaction-diﬀusion models, for example, is
essentially based on a mean-ﬁeld approach [30]. However, as has already been realized
by Smoluchowski [31], ﬂuctuations may be extremely important in low-dimensional sys-
tems where the diﬀusive mixing is not strong enough [32]. For example, if particles of
one species diﬀuse and annihilate by the reaction A + A →Ø, the standard mean-ﬁeld
approximation predicts an asymptotic decay of the particle concentration as ρ(t) ∼t−1.
In one dimension, however, the density is found to decay as ρ(t) ∼t−1/2. This slow decay
is due to ﬂuctuations produced by the dynamics, leading to spatial anticorrelations of the
particles. The existence of such ﬂuctuation eﬀects has been conﬁrmed experimentally by
measuring the luminescence of annihilating excitons on polymer chains [22].
As in equilibrium statistical mechanics, nonequilibrium phenomena are particularly
interesting if the system undergoes a phase transition, leading to a collective behavior of the
particles over long distances. There is a large variety of phenomenological nonequilibrium
phase transitions in nature, ranging from morphological transitions of growing surfaces [33]
to traﬃc jams [34]. It turns out that the concept of universality, which has been very
successful in the ﬁeld of equilibrium critical phenomena, can be applied to nonequilibrium
phase transitions as well.
However, the universality classes of nonequilibrium critical
phenomena are expected to be even more diverse as they are governed by various symmetry
properties of the evolution dynamics. On the other hand, the experimental evidence for
universality of nonequilibrium phase transitions is still very poor, calling for intensiﬁed
experimental eﬀorts.
In the present work we will focus on nonequilibrium phase transitions in models with
so-called absorbing states, i.e., conﬁgurations that can be reached by the dynamics but
6

<!-- page 7 -->

cannot be left. The most important universality class of absorbing-state transitions is
directed percolation (DP) [35]. This type of transition occurs, for example, in models for
the spreading of an infectious disease. In these models the lattice sites are considered as
individuals which can be healthy or infected. Infected individuals may either recover by
themselves or infect their nearest neighbors. Depending on the infection rate, the spreading
process may either survive or evolve into a passive state where the infection is completely
eliminated. In the limit of large system sizes the two regimes of survival and extinction are
separated by a continuous phase transition. As in equilibrium statistical mechanics, the
critical behavior close to the transition is characterized by diverging correlation lengths
associated with certain critical exponents.
Similar spreading processes with the same
exponents can be observed in models for catalytic reactions, percolation in porous media,
and even in certain hadronic interactions. It turns out that all these phase transitions
belong generically to a single universality class, irrespective of microscopic details of their
dynamic rules. In view of its robustness, the DP class may therefore be as important
as the Ising universality class in equilibrium statistical mechanics. Amazingly, directed
percolation is one of very few critical phenomena which cannot be solved exactly in one
spatial dimension. Although DP is easy to deﬁne, its critical behavior is highly nontrivial.
This is probably one of the reasons why DP continues to fascinate theoretical physicists.
The present review addresses several aspects of nonequilibrium phase transition2. In
the following Section we introduce elementary concepts of nonequilibrium statistical me-
chanics such as the master equation, reaction diﬀusion processes, Monte Carlo simulations,
as well as the most important analytical methods and approximation techniques. The
third Section discusses the problem of directed percolation, including a comprehensive in-
troduction to DP lattice models, basic scaling concepts, approximation techniques, as well
as ﬁeld-theoretic methods. In view of the robustness of DP, it is particularly interesting
to search for non-DP phase transitions which usually emerge in presence of additional
symmetries. These exceptional universality classes, which have attracted considerable at-
tention during the last few years, will be reviewed in Sec. 4. Sec. 5 discusses a simulation
technique, called damage spreading, which has been used in the past to search for chaotic
behavior in random processes. It is shown that this technique suﬀers of severe conceptual
problems, making it impossible to deﬁne chaotic phases. We also discuss the critical be-
havior of damage spreading transitions which are closely related to phase transitions into
absorbing states. Finally, we turn to depinning transitions in models of growing interfaces
which are related to nonequilibrium phase transitions into absorbing states as well. As
it is not intended to cover the whole ﬁeld of nonequilibrium critical phenomena, we will
not address various related topics such as self-organized critical phenomena [36], modiﬁed
reaction-diﬀusion processes [32], the dynamics of reacting fronts [37], driven diﬀusive sys-
tems [11], and recent results on spontaneous symmetry breaking and phase separation in
one-dimensional systems [38]. For further reading we will give references to related ﬁelds.
Supplementary information concerning the deﬁnition of tensor products, the derivation of
the eﬀective action of Reggeon ﬁeld theory, Wilson’s shell integration, and the one-loop
integrals for DP are given in the appendices A-D. For easy reference we also append a list
of frequently used symbols and abbreviations.
2The review is based on a Habilitation thesis submitted by the author to the Free University of Berlin
in June 1999.
7

<!-- page 8 -->

2
Stochastic many-particle systems
In this Section we discuss elementary concepts of nonequilibrium statistical mechanics. In
order to introduce basic notions, we ﬁrst consider the example of a simple random walk.
Turning to many-particle systems we introduce the asymmetric exclusion process which
is a model for biased diﬀusion of many particles on a one-dimensional line. Moreover,
we explain the standard mean ﬁeld approach to reaction-diﬀusion processes. In order to
demonstrate the importance of ﬂuctuations, two simple lattice models with particle reac-
tions will be discussed, namely coagulation 2A →A and pair annihilation 2A →Ø. It
turns out that in one dimension the temporal evolution of these systems diﬀers signiﬁ-
cantly from the mean ﬁeld prediction, proving that ﬂuctuations may play an important
role. Furthermore, we review basic concepts of numerical simulation techniques comparing
diﬀerent update schemes. Finally we turn to certain analytical methods by which reaction
diﬀusion models can be solved exactly.
In particular we discuss a recently introduced
algebraic technique which allows the stationary state of certain nonequilibrium models to
be expressed in terms of products of noncommutative operators.
2.1
The one-dimensional random walk
In order to introduce basic concepts of nonequilibrium statistical physics, let us ﬁrst con-
sider a simple symmetric random walk on a one-dimensional line. The ‘conﬁguration’ of
this dynamical system at time t is characterized by the position of the walker s(t). A
random walk may be deﬁned either on a continuous manifold or on a lattice.
If both
position s and time t are discrete variables, an unbiased random walk may be realized by
the random process
s(t + 1) = s(t) + X(t) .
(1)
In this expression X(t) = ±1 is a ﬂuctuating random variable with correlations
⟨X(t)⟩= 0 ,
⟨X(t)X(t′)⟩= δt,t′ ,
(2)
where ⟨. . . ⟩denotes the average over many realizations of randomness. While the individ-
ual space-time trajectory of a random walker is not predictable, the probability distribution
Pt(s) to ﬁnd the walker after t time steps at position s evolves deterministically according
to the so-called master equation
Pt+1(s) = 1
2

Pt(s −1) + Pt(s + 1)

.
(3)
Assuming the particle to be initially located at the origin P0(s) = δs,0, this diﬀerence
equation is solved by
Pt(s) =
( 1
2t
 t
(t+s)/2

if t + s even,
0
if t + s odd.
(4)
If space and time are continuous, the motion of a random walker may be described by a
stochastic Langevin equation
∂ts(t) = ζ(t) ,
(5)
8

<!-- page 9 -->

where, according to the central limit theorem, ζ(t) is a Gaussian white noise with zero
mean and correlations ⟨ζ(t)ζ(t′)⟩= Γδ(t −t′). The Langevin equation may be regarded
as a continuum version of Eq. (2). Starting from the origin s(0) = 0 the mean square
displacement of the random walker grows as ⟨s2(t)⟩=
R t
0 dt1
R t
0 dt2⟨ζ(t1)ζ(t2)⟩= Γt. The
resulting probability distribution
Pt(s) =
1
2πΓte−s2/(2Γt)
(6)
is a solution of the Fokker-Planck equation [12]
∂tPt(s) = Γ
2 ∂2
s Pt(s)
(7)
which can be seen as a variant of the master equation in a continuum (3).
The example of a random walk is particularly simple as it involves only one degree
of freedom.
In order to describe systems with many particles, it would seem natural
to introduce several degrees of freedom s1, s2, . . . , where sn denotes the position of the
n-th particle. However, this approach is restricted to systems with a conserved number
of particles. For systems with non-conserved particle number it is more convenient to
introduce local degrees of freedom for the number of particles located at certain positions
in space.
2.2
The master equation
Stochastic systems with many particles are usually deﬁned on a d-dimensional Euclidean
manifold representing the physical ‘space’. Attached to this manifold are local degrees
of freedom characterizing the conﬁguration of the system.
Depending on whether the
spatial manifold is continuous or discrete, the local degrees of freedom are introduced
as continuous ﬁelds or local variables residing at the lattice sites. Furthermore, a time
coordinate t is introduced which may be interpreted as an additional dimension of the
system. Therefore, stochastic models are said to be deﬁned in d+1 dimensions. Since t may
be continuous or discontinuous, we have to distinguish between models with asynchronous
and synchronous dynamics.
Asynchronous dynamics
Stochastic models with continuous time evolve by asynchronous dynamics, i.e., transitions
from a state s into another state s′ occur spontaneously at a given rate ws→s′ ≥0 per
unit time. It can be shown that in the limit of very large systems sizes the temporal
evolution of the probability distribution Pt(s) evolves deterministically according to a
master equation with appropriate initial conditions [23–25].
The master equation is a
linear partial diﬀerential equation describing the ﬂow of probability into and away from a
conﬁguration s:
∂
∂tPt(s) =
X
s′
ws′→sPt(s′)
|
{z
}
gain
−
X
s′
ws→s′Pt(s)
|
{z
}
loss
.
(8)
The gain and loss terms balance one another so that the normalization P
s Pt(s) = 1 is
conserved. Since the temporal change of Pt(s) is fully determined by the actual probability
9

<!-- page 10 -->

distribution at time t, the master equation describes a Markov process, i.e., it has no
intrinsic memory. Moreover, it is important to note that the coeﬃcients ws→s′ are rates
rather than probabilities. Thus, they may be larger than 1 and can be rescaled by changing
the time scale.
Using a vector notation (see Appendix A) the master equation (8) may be written as
∂t|Pt⟩= −L|Pt⟩,
(9)
where |Pt⟩denotes a vector whose components are the probabilities Pt(s). The Liouville
operator L generates the temporal evolution and is deﬁned through the matrix elements
⟨s′|L|s⟩= −ws→s′ + δs,s′
X
s′′
ws→s′′ .
(10)
A formal solution of the master equation is given by |Pt⟩= exp(−Lt)|P0⟩, where |P0⟩
denotes the initial probability distribution.
Therefore, in order to determine |Pt⟩, the
Liouville operator has to be diagonalized which is usually a nontrivial task.
Apart from very few exceptions, stochastic processes are irreversible and therefore not
invariant under time reversal. Hence the Liouville operator L is generally non-hermitean.
Moreover, it may have complex conjugate eigenvalues, indicating oscillatory behavior.
Oscillating modes are not only a mathematical artifact, but can be observed experimentally
in certain chemical reactions such as the Belousov-Zhabotinski reaction [39]. Due to the
positivity of rates, the real part of all eigenvalues is nonnegative, i.e., the amplitude of
excited eigenmodes decays exponentially in time. The spectrum of the Liouville operator
includes at least one zero mode L|Ps⟩= 0, representing the stationary state of the system.
Moreover, probability conservation can be expressed as ⟨1|Pt⟩= 1, where ⟨1| denotes the
sum vector over all conﬁgurations (cf. Appendix A). Consequently the Liouville operator
obeys the equation ⟨1|L = 0, i.e., the sum over each column of L vanishes.
Synchronous dynamics
If the time variable t is a discrete quantity, the model evolves by synchronous dynamics,
i.e., all lattice sites are simultaneously updated according to certain transition probabilities
ps→s′ ∈[0, 1]. The corresponding master equation is a linear recurrence relation
Pt+1(s) = Pt(s) +
X
s′
ps′→sPt(s′)
|
{z
}
gain
−
X
s′
ps→s′Pt(s)
|
{z
}
loss
,
(11)
which can be written in a compact form as a linear map
|Pt+1⟩= T |Pt⟩,
(12)
where T is the so-called transfer matrix. A formal solution is given by |Pt⟩= T t|P0⟩. As
can be veriﬁed easily, the conservation of probability ⟨1|Pt⟩= 1 implies that ⟨1|T = ⟨1|,
i.e., the sum over each column of the transfer matrix is equal to 1.
There has been a long debate which of the two update schemes is the more ‘realistic’
one. For many researchers models with uncorrelated spontaneous updates appear to be
more ‘natural’ than models with synchronous dynamics where all particles move simul-
taneously according to an artiﬁcial clock cycle. On the other hand, many computational
10

<!-- page 11 -->

t
i+1
s
q-1
i
q
t
t+dt
t
t+1/2
t+1
asynchronous:
synchronous:
s
Fig. 1:
The partially asymmetric exclusion process in one dimension with synchronous (sublattice-parallel)
and asynchronous (random-sequential) dynamics.
physicists prefer stochastic cellular automata with synchronous dynamics since they can be
implemented eﬃciently on parallel computers. However, as a matter of fact, in both cases
the dynamic rules are simpliﬁed descriptions of a much more complex physical process.
Therefore, it would be misleading to consider one of the two variants as being more ‘nat-
ural’ than the other. Instead, the choice of the dynamic procedure should depend on the
speciﬁc physical system under consideration. Very often both variants display essentially
the same physical properties. In some cases, however, they lead to diﬀerent results. For
example, models for traﬃc ﬂow with synchronous updates turn out to be more realistic
than random sequential ones. Another example is polynuclear growth (see Sec. 6.3), where
a roughening transition occurs only when synchronous updates are used.
2.3
Diﬀusion of many particles: The asymmetric exclusion process
One of the simplest stochastic many-particle models is the partially asymmetric exclusion
process on a one-dimensional chain with N sites [5]. In this model hard-core particles
move randomly to the right (left) at rate q (q−1). An exclusion principle is imposed, i.e.,
each lattice site may be occupied by at most one particle. Therefore, attempted moves are
rejected if the target site is already occupied. In the following we assume closed boundary
conditions, i.e., particles cannot leave or enter the system at the ends of the chain. The
conﬁguration s = {s1, s2, . . . , sN} of the system is given in terms of local variables si,
indicating presence (si = 1) or absence (si = 0) of a particle at site i.
Exclusion process with asynchronous dynamics
Let us ﬁrst consider the exclusion process with asynchronous dynamics. In this case a
pair of sites i and i + 1 is randomly selected. If only one of the two sites is occupied, the
particle moves with probability q/(q+q−1) to the right and with probability q−1/(q+q−1)
to the left, as shown in Fig. 1. Each update attempt corresponds to a time increment of
1/(N(q + q−1)). Thus, the transition rates are deﬁned by
ws→s′ =
N−1
X
i=1
i−1
Y
j=1
δsj,s′
j

N
Y
j=i+2
δsj,s′
j

q δsi,1 δsi+1,0 δs′
i,0 δs′
i+1,1 +
q−1 δsi,0 δsi+1,1 δs′
i,1 δs′
i+1,0

.
(13)
11

<!-- page 12 -->

The corresponding Liouville operator can be written as
L =
N−1
X
i=1
1 ⊗1 ⊗. . .
⊗Li⊗
| {z }
i-th position
. . . ⊗1 =:
N−1
X
i=1
Li ,
(14)
where 1 denotes a 2 × 2 unit matrix and Li is a 4 × 4 matrix generating particle hopping
between sites i and i + 1. In the standard basis (see Appendix A) this matrix is given by
Li =




0
0
0
0
0
q−1
−q
0
0
−q−1
q
0
0
0
0
0



.
(15)
As can be veriﬁed easily, a stationary state of the system is given (up to normalization)
by the tensor product
|Ps⟩=
1
1

⊗
q−1
q

⊗. . . ⊗
q−N
qN

=
N
O
j=1
q−j
qj

.
(16)
Since the vector |Ps⟩can be written as a tensor product, the local variables si are com-
pletely uncorrelated. Such a state is said to have a product measure.
As the total number of particles is conserved in the asymmetric exclusion process, the
dynamics decomposes into independent sectors. In fact, the Liouville operator commutes
with the particle number operator
M =
N
X
i=1
mi ,
mi =
1
0
0
0

.
(17)
Obviously the vector |Ps⟩is a superposition of solutions belonging to diﬀerent sectors,
i.e., it represents a whole ensemble of stationary states. To obtain a physically mean-
ingful solution for a given number of particles, the vector |Ps⟩has to be projected onto
the corresponding sector. In this sector the stationary system evolves through certain
conﬁgurations with speciﬁc weights given by the normalized components of the projected
vector.
Exclusion process with synchronous dynamics
The asymmetric exclusion process with synchronous updates may be realized by introduc-
ing two half time steps. In the ﬁrst half time step the odd sublattice is updated whereas
the even sublattice is updated in the second half time step (see Fig. 1). Note that the use
of sublattice-parallel updates admits local dynamic rules3. Assuming the number of sites
N to be odd, the corresponding transfer matrix reads
T = (T2 ⊗T4 ⊗T6 ⊗. . . ⊗TN−1) (T1 ⊗T3 ⊗T5 ⊗. . . ⊗TN−2) ,
(18)
where
Ti =
1
q + q−1




q + q−1
0
0
0
0
q
q
0
0
q−1
q−1
0
0
0
0
q + q−1




(19)
3In the exclusion process with fully parallel updates local moves may overlap, leading to subtle long-
range correlations (see Ref. [40]).
12

<!-- page 13 -->

diffusion
t
t+dt
coagulation
decoagulation
annihilation
Fig. 2:
Nearest-neighbor processes in reaction-diﬀusion models with only one particle species.
is the local hopping matrix. Again the product state (16) is a stationary eigenvector of
the transfer matrix. Thus, both the asynchronous and the synchronous exclusion process
have exactly the same stationary properties.
Asymmetric diﬀusion in a continuum
Let us ﬁnally turn to asymmetric diﬀusion on a continuous manifold.
In principle it
would be possible to trace trajectories of individual particles. However, it is much more
convenient to characterize the state of the system by a density ﬁeld ρ(x, t), rendering the
coarse-grained density of particles at position x at time t. The Langevin equation of such
a system may be written as
∂tρ(x, t) = D[ρ(x, t)] + U[ρ(x, t)] + ζ(x, t) ,
(20)
where D is a sum of linear diﬀerential operators describing spatial couplings, U a potential
for on-site particle interactions, and ζ(x, t) a noise term taking the stochastic nature of
Brownian motion into account. Since particles do not interact in the present case, the
potential U vanishes. Moreover, as will be shown below, the noise ζ(x, t) is irrelevant on
large scales and can be neglected. Thus, the resulting Langevin equation reads
∂tρ(x, t) = q + q−1
2
∂2
xρ(x, t) + q −q−1
2
∂xρ(x, t) .
(21)
The second term describes the bias of the diﬀusive motion which may be eliminated in a
co-moving frame. Notice that this equation is linear and does not incorporate the exclusion
principle. In a co-moving frame it reduces to the ordinary diﬀusion equation.
The diﬀusion equation provides a simple example of dynamic scaling invariance. As
can be veriﬁed easily, the equation ∂tρ(x, t) = D∇2
xρ(x, t) is invariant under rescaling of
space and time
x →Λx ,
t →Λzt ,
(22)
where z is the so-called dynamic exponent. Since space and time are diﬀerent in nature,
the exponent z is usually larger than 1. The value z = 2 indicates diﬀusive behavior.
2.4
Reaction-diﬀusion processes
Reaction-diﬀusion processes are stochastic models for chemical reactions in which parti-
cles are predominantly transported by thermal diﬀusion. Usually a chemical reaction in
a solvent or on a catalytic surface consists of a complex sequence of intermediate steps.
In reaction-diﬀusion models these intermediate steps are ignored and the reaction chain
is replaced by simpliﬁed probabilistic transition rules. The involved atoms and molecules
are interpreted as particles of several species, represented by capital letters A, B, C, . . . .
These particles neither carry a mass nor an internal momentum, instead the conﬁguration
13

<!-- page 14 -->

of a reaction-diﬀusion model is completely speciﬁed by the position of the particles. On a
lattice with exclusion principle such a conﬁguration can be expressed in terms of local vari-
ables si = 0, 1, 2, 3, . . . , representing a vacancy Ø and particles A, B, C, . . . , respectively.
Sometimes it is even not necessary to keep track of all substances involved in a chemical
reaction. For example, if a molecule of a gas phase is adsorbed at a catalytic surface,
this process may be eﬀectively described by spontaneous particle creation Ø →A without
modeling the explicit dynamics in the gas phase. Therefore, the number of particles in
reaction-diﬀusion models is generally not conserved.
Apart from spontaneous particle creation Ø →A many other reactions are possible.
Unary reactions are spontaneous transitions of individual particles, the most important
examples being
A →Ø
self-destruction,
A →2A
oﬀspring production or decoagulation,
A →B
transmutation,
A →A + B
induced creation of particles.
On the other hand, binary reactions require two particles to meet at the same place (or
at neighboring sites). Here the most important examples are
2A →Ø
pair annihilation,
2A →A
coagulation (coalescence),
A + B →Ø
two-species annihilation,
A + B →B
induced adsorption of particles.
In addition, the particles may diﬀuse with certain rates in the same way as in the previ-
ously discussed exclusion process. A process is called diﬀusion-limited if diﬀusion becomes
dominant in the long-time limit, i.e., the diﬀusive moves become much more frequent than
reactions. This happens, for example, in reaction-diﬀusion models with binary reactions
when the particle density is very low. On the other hand, if particle reactions become
dominant after very long time, the process is called reaction-limited.
In the following we will focus on simple reaction-diﬀusion models with only one type
of particles (see Fig. 2). They can be considered as two-state models since each site can
either be occupied by a particle (A) or be empty (Ø). Examples include the so-called
coagulation model in which particles diﬀuse at rate D, coagulate at rate λ and decoagulate
at rate κ:
AØ ↔
D ØA ,
AA →
λ AØ, ØA ,
AØ, ØA →
κ AA .
(23)
Using the same notation as in Eq. (15), the corresponding nearest-neighbor transition
matrix is given by
Lcoag
i
=




0
0
0
0
0
D + κ
−D
−λ
0
−D
D + κ
−λ
0
−κ
−κ
2λ



.
(24)
The special property of this model lies in the fact that the empty state cannot be reached
by the dynamic processes. An exact solution of the coagulation model will be discussed
in Sec. 2.8.
14

<!-- page 15 -->

Another important example is the annihilation model where particles diﬀuse and an-
nihilate:
AØ ↔
D ØA ,
AA →
α ØØ .
(25)
The corresponding interaction matrix is given by
Lannh
i
=




0
0
0
−α
0
D
−D
0
0
−D
D
0
0
0
0
α



.
(26)
In the annihilation model the number of particles is conserved modulo 2. As will be shown
in Sec. 2.8, both models are equivalent and can be related by a similarity transformation.
2.5
Mean ﬁeld approximation
In many cases the macroscopic properties of a reaction-diﬀusion process can be predicted
by solving the corresponding mean ﬁeld theory. In chemistry the simplest mean ﬁeld ap-
proximation is known as the ‘law of mass action’: For a given temperature the rate of
a reaction is assumed to be proportional to the product of concentrations of the react-
ing substances. This approach assumes that the particles are homogeneously distributed.
It therefore ignores any spatial correlations as well as instabilities with respect to inho-
mogeneous perturbations. Thus, the homogeneous mean ﬁeld approximation is expected
to hold on scales where diﬀusive mixing is strong enough to wipe out spatial structures.
Especially in higher dimensions, where diﬀusive mixing is more eﬃcient, the mean ﬁeld ap-
proximation provides a good description. It becomes exact in inﬁnitely many dimensions,
where all particles can be considered as being neighbored.
The mean ﬁeld equations can be constructed directly by translating the reaction scheme
into a diﬀerential equation for gain and loss of the particle density ρ(t). For example,
in the mean ﬁeld approximation of the coagulation model (23) the process A
κ→2A
takes place with a frequency proportional to κρ(t), leading to an increase of the particle
density. Similarly, the coagulation process 2A λ→A decreases the number of particles with
a frequency proportional to λρ2(t). Ignoring diﬀusion, the resulting mean ﬁeld equation
reads
∂tρ(t) = κρ(t) −λρ2(t) ,
(27)
where κ and λ are the rates for decoagulation and coagulation, respectively. In contrast
to the master equation this diﬀerential equation is nonlinear. For κ > 0 it has two ﬁxed
points, namely an unstable ﬁxed point at ρ = 0 and a stable ﬁxed point at ρ = κ/λ. The
physical meaning of the two ﬁxed points is easy to understand. The empty system remains
empty, but as soon as we perturb the system by adding a few particles, it quickly evolves
towards a stationary active state with a certain average concentration ρ > 0. This active
state is then stable against perturbations (such as adding or removing particles).
The mean ﬁeld equation (27) can also be used to predict dynamic properties of the
system. Starting from a fully occupied lattice ρ(0) = 1 the time-dependent solution is
15

<!-- page 16 -->

given by
ρ(t) =
κ
λ −(λ −κ)e−κt .
(28)
In the limit of a vanishing decoagulation rate κ →0, the two ﬁxed points merge into a
marginal one. As in many physical systems this leads to a much slower dynamics. In fact,
for κ = 0 Eq. (28) turns into
lim
κ→0 ρ(t) =
1
1 + λt ,
(29)
i.e., the particle density decays asymptotically according to a power law as
ρ(t) ∼t−1 .
(30)
The stability of the mean ﬁeld solution with respect to inhomogeneous perturbations may
be studied by adding a term for diﬀusion
∂tρ(x, t) = κρ(x, t) −λρ2(x, t) + D∇2ρ(x, t) .
(31)
In the present case the diﬀusive term suppresses perturbations with short wavelength and
therefore stabilizes the homogeneous solutions. However, in certain chemical reactions with
several particle species such a diﬀusive term may have a destabilizing inﬂuence. The study
of mean-ﬁeld instabilities is the starting point for the theory of pattern formation which
has become an important ﬁeld of statistical physics [30]. A very interesting application is
the Belousov-Zhabotinski reaction [39] that produces rotating spirals in a Petri dish.
It may be surprising that even simple reaction-diﬀusion processes are described by
nonlinear mean ﬁeld rate equations, whereas the corresponding master equation is always
linear. However, mean ﬁeld and Langevin equations are always deﬁned in terms of coarse-
grained particle densities involving many local degrees of freedom. These coarse-grained
densities, which can be thought of as observables in conﬁguration space, may evolve ac-
cording to a nonlinear laws. A similar paradox occurs in quantum physics: Although the
Schr¨odinger equation is strictly linear, most observables evolve in a highly nonlinear way.
2.6
The inﬂuence of ﬂuctuations
Although the mean ﬁeld equation (31) includes a term for particle diﬀusion, it still ignores
ﬂuctuation eﬀects and spatial correlations. However, especially in low-dimensional sys-
tems, ﬂuctuations may play an important role and are able to entirely change the physical
properties of a reaction-diﬀusion process.
In order to demonstrate the inﬂuence of ﬂuctuations, let us consider the coagulation
process 2A →A. The full Langevin equation for this process reads
∂tρ(x, t) = −λρ2(x, t) + D∇2ρ(x, t) + ζ(x, t) ,
(32)
where ζ(x, t) is a noise term which accounts for the ﬂuctuations of the particle density at
position x at time t. Clearly, the noise amplitude depends on the magnitude of the density
ﬁeld ρ(x, t). In particular, without any particles present, there will be no ﬂuctuations.
According to the central limit theorem, the noise is expected to be Gaussian with a
16

<!-- page 17 -->

squared amplitude proportional to the frequency of events leading to a change of the
particle number. Since the particle number only ﬂuctuates when two particles coagulate,
this frequency should be proportional to ρ2(x, t). Following these naive arguments, the
noise correlations should be given by
⟨ζ(x, t)⟩
=
0 ,
(33)
⟨ζ(x, t)ζ(x′, t′)⟩
=
Γ ρ2(x, t) δd(x −x′) δ(t −t′) ,
where Γ denotes the noise amplitude and d the spatial dimension. The next question
is to what extent the macroscopic behavior of the system will be aﬀected by the noise.
Typically there are three possible answers:
1. The noise is irrelevant on large scales so that the macroscopic behavior is correctly
described by the mean ﬁeld solution.
2. The noise is relevant on large scales, leading to a macroscopic behavior that is dif-
ferent from the mean-ﬁeld prediction.
3. The noise is marginal, producing (typically logarithmic) deviations from the mean-
ﬁeld solution.
In order to ﬁnd out whether the noise is relevant on large scale we need to introduce the
concept of renormalization [41]. The term ‘renormalization’ refers to various theoretical
methods investigating the scaling behavior of physical systems under coarse-graining of
space and time. Roughly speaking, it describes how the parameters of a system have to be
adjusted under coarse-graining of lengths scales without changing its physical properties.
A ﬁxed point of the renormalization ﬂow is then associated with certain universal scaling
laws of the system. The simplest renormalization group (RG) scheme ignores the inﬂuence
of ﬂuctuations. This approach is referred to as ‘mean ﬁeld renormalization’. Approaching
the ﬁxed point, the noise amplitude may diverge, vanish or stay ﬁnite, corresponding to
the classiﬁcation given above.
Hence, by studying mean ﬁeld renormalization, we can
predict whether ﬂuctuations are relevant or not.
In the mean ﬁeld approximation the Langevin equation (32) may be renormalized by
a scaling transformation
x →Λ x ,
t →Λz t ,
ρ(x, t) →Λχ ρ(Λx, Λzt) ,
(34)
where z denotes the dynamic exponent. The exponent χ describes the scaling properties of
the density ﬁeld itself. If the particles were distributed homogeneously, the ﬁeld would scale
as an ordinary density, that is, with the exponent χ = −d. However, in the coagulation
process nontrivial correlations between particles lead to a diﬀerent scaling dimension of
the particle distribution. In fact, invariance of Eqs. (32)-(33) under rescaling implies that
z = 2 and χ = −2. Therefore, the noise amplitude scales as
Γ →Λ1−d/2Γ ,
(35)
where d is the spatial dimension. Hence in one spatial dimension ﬂuctuations are relevant
whereas they are marginal in two and irrelevant in d > 2 dimensions. The value of d where
the noise becomes marginal is denoted as the upper critical dimension dc. For the coagula-
tion model the upper critical dimension is dc = 2. Above the critical dimension the mean
17

<!-- page 18 -->

t
Fig. 3:
Monte Carlo simulation of the (1+1)-dimensional coagulation-diﬀusion model A + A →A with
random-sequential updates. The ﬁgure shows an individual run starting with a fully occupied lattice of
200 sites. As can be seen, the coagulation process is diﬀusion-limited.
ﬁeld approximation provides a correct description, whereas for d < dc ﬂuctuation eﬀects
have to be taken into account. This can be done by using improved mean ﬁeld approaches,
exact solutions, as well as ﬁeld-theoretic renormalization group techniques [42].
A systematic ﬁeld-theoretic analysis of the coagulation process 2A →A leads to an
unexpected result: The noise amplitude Γ in Eq. (33) turns out to be negative [43].
Consequently, the noise ζ(x, t) is imaginary. This result is rather counterintuitive as we
expect the noise to describe density ﬂuctuations which, by deﬁnition, are real. However,
since the noise amplitude is a measure of annihilation events, it is subjected to correlations
that are produced by the annihilation process itself. In one dimension these correlations
are negative, i.e., particles avoid each other. This simple example demonstrates that it can
be dangerous to set up a Langevin equation by considering the mean ﬁeld equation and
adding a physically reasonable noise ﬁeld. Instead it is necessary to derive the Langevin
equation directly from the microscopic dynamics, as explained in Ref. [44].
2.7
Numerical simulations
To verify analytical results, it is often helpful to perform Monte Carlo simulations. In order
to demonstrate this numerical technique, let us again consider the coagulation process
2A →A on a one-dimensional chain. For simplicity we assume the rates for diﬀusion and
coagulation to be equal. This ensures that particles can move at constant rate irrespective
of the state of the target site. If the target site is empty, it will be occupied by the moving
particle. On the other hand, if the target site is already occupied, the two particles will
coagulate into a single one. Such a move from site i to site j may be realized by the pseudo
code instruction
Move(i,j) { if (s[i]==1) { s[i]=0; s[j]=1; }};
where s[i] denotes the occupation variable si = 0, 1 at site i. In one dimension particles
move randomly to the left and to the right. Thus a local update at sites (i, i + 1) may be
realized by the instruction
18

<!-- page 19 -->

10
0
10
1
10
2
10
3
10
4
t
10
−4
10
−3
10
−2
10
−1
10
0
ρ(t)
d=1
d=2
d=3
10
0
10
1
10
2
10
3
10
4
t
0.1
1.0
10.0
 t ρ(t)
d=3
d=2
Fig. 4:
Monte Carlo simulation of the coagulation-diﬀusion model. Left: Decay of the particle density
ρ(t) measured in a system with 105 sites in one, two, and three spatial dimensions, averaged over 104
independent runs. The dashed lines indicate the slopes −1/2 and −1, respectively. Right: tρ(t) versus t,
illustrating logarithmic corrections in d = 2 dimensions.
Update(i) {
if (rnd(0,1)<0.5) Move(i,i+1);
else Move(i+1,i); };
where rnd(0,1) returns a real random number from a ﬂat distribution between 0 and 1.
Since the coagulation model evolves by asynchronous dynamics it uses so-called random-
sequential updates, i.e., the update attempts take place at randomly selected pairs of sites.
A Monte Carlo sweep consists of N such update attempts:
for (i=1; i<=N; i++) Udpate(rndint(1,N-1));
where N denotes the lattice size and rndint(1,N-1) returns an integer random number
between 1 and N −1. Since on average each site is updated once, such a sweep corresponds
to a unit time step. It can be proven that the statistical ensemble of space-time trajec-
tories generated by random-sequential updates converges to the solution of the master
equation (8) in the limit N →∞. The above update algorithm can easily be generalized
to more complicated reaction schemes and higher dimensional lattices.
The coagulation process with synchronous updates may be simulated by using parallel
updates on alternating sublattices:
for (i=1; i<=N-1; i+=2) Update(i);
for (i=2; i<=N-2; i+=2) Update(i);
In Monte Carlo simulations most of the CPU time is consumed for generating random
numbers. Therefore, models with parallel updates are usually more eﬃcient since it is
not necessary to determine random positions for the updates. In addition, models with
parallel updates can be implemented easily on computers with parallel architecture [45].
Fig. 4 shows the particle density as a function of time for the coagulation model
with random-sequential updates and closed boundary conditions in various dimensions.
The particle concentration is averaged over 104 independent runs and plotted in a double-
logarithmic representation, where straight lines indicate power-law behavior. As expected,
19

<!-- page 20 -->

the mean-ﬁeld prediction ρ(t) ∼1/t is reproduced in d > 2 dimensions. In one dimension,
however, the graph suggests the density to decay as
ρ(t) ∼t−1/2
(d = 1) .
(36)
Thus, the simulation result demonstrates that ﬂuctuation eﬀects can change the asymp-
totic behavior (an exact solution will be discussed in Sec. 2.8). At the critical dimension
d = dc = 2 the density ρ(t) deviates slightly from the mean-ﬁeld prediction, indicating
logarithmic corrections. In fact, as can be shown by a ﬁeld-theoretic analysis [44], the
density decays asymptotically as
ρ(t) ∼



t−d/2
for d < 2 ,
t−1 ln t
for d = dc = 2 ,
t−1
for d > 2 .
(37)
2.8
Exact results
Equivalence of annihilation and coagulation processes
Sometimes it is possible to relate diﬀerent stochastic processes by an exact similarity
transformation [46–48]. For example, the coagulation process 2A →A and the annihi-
lation process 2A →Ø deﬁned in Sec. 2.4 are fully equivalent if their rates are tuned
appropriately. More precisely, for a particular choice of the rates it is possible to ﬁnd a
similarity transformation U such that
Lcoag = ULannhU−1 ,
(38)
where we assume the chains to have closed ends, i.e.,
Lcoag =
N−1
X
i=1
Lcoag
i
,
Lannh =
N−1
X
i=1
Lannh
i
.
(39)
Since Lcoag and Lannh are non-hermitean operators, the similarity transformation U is
not orthogonal.
However, if U exists, the two operators will have the same spectrum
of eigenvalues. As can be veriﬁed easily, the local operators Lcoag
i
and Lannh
i
have the
eigenvalues {0, 0, 2D + κ, 2λ + κ} and {0, 0, 2D, a}, respectively. Therefore, choosing the
rates
D = λ = 1,
α = 2,
κ = 0 ,
(40)
both operators obtain the same spectrum {0, 0, 1, 1}. Moreover, it can be shown that they
both obey the same commutation relations, namely the so-called Hecke algebra [49,50]
L2
i = Li ,
LiLi+1Li −Li+1LiLi+1 = 2(Li −Li+1) ,
(41)
[Li, Lj] = 0 for (|i −j| ≥2) .
Since this algebra generates the spectrum of L, we can conclude that the spectra of Lcoag
and Lannh coincide for an arbitrary number of sites. Obviously, the equivalence of the
spectra is a necessary condition for the existence of a similarity transformation between
20

<!-- page 21 -->

q
0
1
0
1
α
β
MAX
HD
LD
α
β
2
1
N
-1
q
-1
q
q
Fig. 5:
The partially asymmetric exclusion process with particle adsorption and desorption at the bound-
aries. Left: Particles adsorb at the left boundary at rate α, perform a biased random walk in the bulk
of the chain, and ﬁnally desorb from site N at rate β. Right: Phase diagram for the totally asymmetric
exclusion process, denoting the low (LD), high (HD) and maximal (MAX) density phases. Analogous
phase diagrams for partial asymmetry were obtained in Ref. [57].
the two systems. In the present case it is even possible to compute the similarity transfor-
mation explicitly. It turns out that U can be expressed in terms of local tensor products
(see Appendix A)
U = u ⊗u ⊗. . . ⊗u =
N
O
i=1
u = u⊗N ,
(42)
where
u =
1
−1
0
2

,
u−1 =
1
1/2
0
1/2

.
(43)
Consequently the n-point density correlation functions of both models are related by
⟨sj1sj2 · · · sjn⟩coag = 2n⟨sj1sj2 · · · sjn⟩annh .
(44)
In particular, the particle densities in both models diﬀer by a factor of 2:
ρcoag(t) = 2ρannh(t).
It should be noted that only a subset of initial conditions in the coagulation model can
be mapped onto physically meaningful initial conditions in the annihilation model (see
Ref. [51]).
Exact mapping between equilibrium and nonequilibrium systems
A remarkable progress has been achieved by realizing that certain nonequilibrium models
can be mapped onto well-studied integrable equilibrium models. More speciﬁcally, it has
been shown that the Liouville operator L of a nonequilibrium models may be related to the
Hamiltonian H of an integrable quantum spin systems by similarity transformation [15].
This allows the nonequilibrium model to be solved by exact techniques of equilibrium
statistical mechanics such as free-fermion diagonalization, the Bethe ansatz, or other al-
gebraic methods [2]. For example, the exclusion process can be mapped onto the XXZ
quantum chain [18,26,52], whereas the coagulation-decoagulation process is related to the
XY chain in a magnetic ﬁeld [53–55]. Exact mappings were also found for higher spin
analogs [56]. A complete summary of the known results can be found in Ref. [28].
21

<!-- page 22 -->

In order to demonstrate this technique, let us consider the partially asymmetric exclu-
sion process in one spatial dimension (see Eq. 15). As will be shown in the following, this
model is related to the quantum spin Hamiltonian of the ferromagnetic XXZ Heisenberg
quantum chain with open boundary conditions H = PN−1
i=1 Hi, where
Hi = −1
2

σx
i σx
i+1 + σy
i σy
i+1 + q + q−1
2
(σz
i σz
i+1 −1) + q −q−1
2
(σz
i −σz
i+1)

.
(45)
This quantum chain Hamiltonian generates translations in the corresponding two-dimen-
sional XXZ model in a strongly anisotropic scaling limit. As can be veriﬁed easily, the
Hamiltonian is non-hermitean for q̸ = 1.
Using the standard basis of Pauli matrices
σx = ( 0 1
1 0 ), σy =
  0 −i
i 0

, and σz =
  1 0
0 −1

the interaction matrix is given by
Hi =




0
0
0
0
0
q−1
−1
0
0
−1
q
0
0
0
0
0



.
(46)
The XXZ Heisenberg chain is integrable by means of Bethe ansatz methods [28]. The
integrability is closely related to two diﬀerent algebraic structures. On the one hand, the
Hamiltonian (45) commutes with the generators K, S± of the quantum algebra Uq[SU(2)]
(see Ref. [58])
KS±K−1 = qS± ,
[S+, S−] = K2 −K−2
q −q−1
,
(47)
where
K = (qσz/2)⊗N ,
S± =
N
X
k=1
(qσz/2)⊗(k−1) ⊗σ± ⊗(q−σz/2)⊗(N−k) .
(48)
Roughly speaking, the quantum group symmetry determines the degeneracies of the eigen-
values of H. On the other hand, the generators Hi are a representation of the Temperley-
Lieb algebra [59]
H2
i = (q + q−1)Hi ,
Hi = HiHi±1Hi ,
(49)
[Hi, Hj] = 0 for |i −j| ≥2 .
This algebra determines the actual numerical value of the energy levels. As realized by
Alcaraz and Rittenberg [56], the same commutation relations are satisﬁed by the transi-
tion matrix Li of the asymmetric exclusion process in Eq. (15). In fact, it can be shown
that the XXZ chain and the exclusion process with α = β = 0 are related by a simi-
larity transformation L = UHU−1, where U can be written as a tensor product of local
transformations
U =
1
0
0
q

⊗
1
0
0
q2

⊗. . . ⊗
1
0
0
qN

=
N
O
i=1
1
0
0
qi

.
(50)
In order to illustrate how symmetries of the equilibrium model translate into physical
properties of the stochastic process, let us consider the quantum group symmetry of the
22

<!-- page 23 -->

gain:
loss:
Fig. 6:
Interparticle distribution functions. The ﬁgure illustrates gain and loss processes for the empty-
interval probability Iℓ(t) by diﬀusion (left) and coagulation (right). The same processes take place at the
left boundary of the interval.
XXZ model. Regarding the exclusion process this symmetry emerges as a conservation
of the total number of particles n. The generators S± act as ladder operators between
diﬀerent sectors with a ﬁxed number of particles. The diagonal operator K is proportional
to q−n, weighting the sectors as in a grand-canonical ensemble, where q plays the role of
a fugacity. The partially asymmetric exclusion process is therefore a physical realization
of a quantum group symmetry with a real-valued deformation parameter q.
A similar mapping relates the coagulation model 2A →A and the ferromagnetic XY
quantum chain in a magnetic ﬁeld, which is deﬁned by the interaction matrix
Hi = −1
2

σx
i σx
i+1 + σy
i σy
i+1 + σz
i + σz
i+1 −2

.
(51)
In fact, it can be easily veriﬁed that the operators Hi satisfy the Hecke algebra (41). The
XY chain is exactly solvable in terms of free fermions [16]. The integrability of the model
is closely related to a SUq(1|1) quantum group symmetry. In the XY chain this symmetry
shows up as a fermionic zero mode, leading to two-fold degenerate energy levels. In the
coagulation model this symmetry emerges as a state without particles that can neither be
reached nor left by the dynamics. In the (equivalent) annihilation model the symmetry
appears as a parity conservation law.
It should be noted that reaction-diﬀusion models are usually related to ferromagnetic
quantum chains, the reason being that the diﬀusion process always corresponds to a fer-
romagnetic interaction in the quantum spin model.
Interparticle distribution functions
Even if a stochastic model can be mapped onto a known equilibrium system by a similarity
transformation, it is often technically diﬃcult to derive physical quantities such as density
proﬁles and correlation functions [60]. For models with an underlying fermionic symmetry
an alternative approach has been developed which does not explicitly use a similarity
transformation. Instead it expresses the state of a model in terms of so-called interparticle
distribution functions (IPDF) [61–63]. Consider, for example, the coagulation-diﬀusion
process with asynchronous dynamics on an inﬁnite chain where particles coagulate (A +
A →A) and diﬀuse at unit rates. Let Iℓbe the probability that an arbitrarily chosen
interval of ℓsites contains no particles.
In terms of these empty-interval probabilities
the master equation can be written in a particularly simple form. Since Iℓ−Iℓ+1 is the
probability to ﬁnd a particle at a neighboring site next to the interval of length ℓ, it is
possible to rewrite diﬀusion and coagulation in terms of gain and loss processes (see Fig. 6)
∂tIℓ(t) = Iℓ−1(t) −Iℓ(t)
|
{z
}
gain
−Iℓ(t) + Iℓ+1(t)
|
{z
}
loss
(52)
23

<!-- page 24 -->

with I0(t) ≡1. It is important to note that this particularly simple form requires the
rates for diﬀusion and coagulation to be identical. This ensures that the gain processes do
not depend on whether the target site is already occupied by a particle. If the two rates
are diﬀerent, higher-order probabilities for several adjacent intervals have to be included,
resulting in a coupled hierarchy of equations. The IPDF method exploits the fact that
this complicated hierarchy of equations decouples for a particular choice of the rates.
By solving the above equation we can compute the particle density ρ(t) which is given
by the probability for an empty interval of length 1 to be absent, i.e.,
ρ(t) = 1 −I1(t) .
(53)
In order to determine the asymptotic behavior of ρ(t) let us consider the continuum limit
of Eq. (52)
(∂t −∂2
ℓ) I(ℓ, t) = 0 ,
I(0, t) = 1 ,
(54)
where ρ(t) = ∂ℓI(ℓ, t)|ℓ=0. This equation has the solution I(ℓ, t) = 1 −erf(x/
√
t). In the
long time limit, the particle density therefore decays algebraically as
ρ(t) ∼t−1/2 ,
(55)
conﬁrming the numerical result of Sec. 2.7. It is interesting to compare this result with the
mean ﬁeld approximation (29) which lead to the incorrect result ρ(t) ∼t−1. Therefore,
the above exact solution demonstrates that ﬂuctuations may inﬂuence the entire temporal
evolution of a stochastic process.
The IPDF technique [61,63] was extended to the coagulation-decoagulation model by
including the inverse reaction A →A + A. Other exact solutions revealed phenomena
such as anomalous kinetics, critical ordering, nonequilibrium dynamic phase transitions,
as well as the existence of Fisher waves [64–68]. The IPDF technique was also used to
study the ﬁnite-size scaling behavior of coagulation processes [54, 69]. Even anisotropic
systems [55] and models with homogeneous or localized particle input [65, 70] have been
solved. Nevertheless the IPDF method is a rather special technique which seems to be
restricted to models with an underlying fermionic symmetry.
2.9
Experimental veriﬁcation of ﬂuctuation eﬀects
The preceding exact calculation proves that the particle concentration in a one-dimensional
coagulation process decays as ρ(t) ∼t−1/2. This result diﬀers signiﬁcantly from the mean-
ﬁeld prediction ρ(t) ∼1/t. Therefore, the coagulation model provides one of the simplest
examples where ﬂuctuation eﬀects change the entire temporal behavior of a reaction-
diﬀusion process.
It is quite remarkable that this result could be veriﬁed experimentally by analyzing
the kinetics of laser-induced excitons on tetramethylammonium manganese trichloride
(TMMC) [22, 71]. TMMC is a crystal consisting of parallel manganese chloride chains.
Laser-induced electronic excitations of the Mn2+ ions, so-called excitons, migrate along
the chain and may be interpreted as quasi-particles. The chains are separated by large
tetramethylammonium ions so that the exchange of excitons between diﬀerent chains is
24

<!-- page 25 -->

Fig. 7:
Decay of the luminescence of TMMC after excitation by a laser pulse for various energies (ﬁgure
reprinted from Ref. [71]). The straight lines are best ﬁts to Eq. (56).
suppressed by a factor of 104. Therefore, the polymer chains can be considered as one-
dimensional systems. Because of exciton-phonon induced lattice distortions the motion
of excitons is diﬀusive. Moreover, when two excitons meet at the same lattice site, the
Mn2+ ion is excited to twice the excitation energy. Subsequently, the ion relaxes back
to a simply exited state by the emission of phonons. Thus, the fusion of excitons can be
viewed as a coagulation process 2A →A+heat on a one-dimensional lattice.
The concentration of quasiparticles can be measured indirectly by detecting the lumi-
nescence intensity I(t) which is proportional to the number of excitons. Eq. (52) predicts
that I(t) ≃I(0)/(1 + αtδ), where α ﬁxes the time scale and δ = 1/2. This equation can
be rewritten as
log
I(0)
I(t) −1

= δ log t + log α .
(56)
The experimental results are shown in Fig. 7. The best ﬁts according to Eq. (56) yield
estimates of about δ = 0.48(3), being in perfect agreement with the theoretical prediction
δ = 1/2. In other experiments the polymers chains are conﬁned to small pores. Here
the excitons perform both annihilation and coagulation processes.
The estimates δ =
0.55(4) [72] and δ = 0.47(3) [73] are again in agreement with the theoretical result.
To summarize, the experimental investigation of excitons on polymer chains conﬁrms
that the concept of stochastic reaction-diﬀusion processes is well justiﬁed in order to quan-
titatively predict the behavior of certain complex systems. In addition, these experiments
prove that ﬂuctuations eﬀects do exist in nature and may change the physical properties
of the system in agreement with the theoretical prediction.
2.10
Dynamic processes approaching thermal equilibrium
Stochastic dynamic processes also play an important role in the context of equilibrium
models.
As outlined in the Introduction, equilibrium statistical mechanics deals with
many-particle systems in contact with a thermal reservoir (heat bath) of temperature T.
In the long-time limit such a system approaches a statistically stationary state where it
25

<!-- page 26 -->

evolves through certain conﬁgurations according to a well-deﬁned probability distribution
Peq(s). The key property of equilibrium models is the existence of an energy functional H
associating each conﬁguration s with a certain energy H(s). The equilibrium distribution
Peq(s) is then given by the canonical ensemble [1]
Peq(s) = 1
Z exp(−H(s)/kBT) ,
(57)
where T is the temperature, kB the Boltzmann constant, and Z the partition sum. This
probability distribution can be used to determine averages of certain macroscopic observ-
ables by summing over all accessible states.
It is important to note that equilibrium
statistical mechanics does not involve any dynamical aspect. In other words, it is irrele-
vant how the system evolves through diﬀerent conﬁgurations, one is only interested in the
relative frequency of certain conﬁgurations to be visited in the stationary state.
Although there is no ‘time’ in equilibrium statistical physics, one may use dynamic
random processes as a tool to generate the equilibrium ensemble Peq(s) of a particular
equilibrium model. More precisely, such a dynamic process evolves into a stationary state
P∞(s) := limt→∞Pt(s) that coincides with the equilibrium ensemble Peq(s). Generally
there is a large variety of dynamic random processes that can be used to generate the
stationary ensemble of a particular equilibrium model. Let us, for example, consider the
ferromagnetic Ising model on a d-dimensional square lattice [7]. Its energy functional is
given by
H(σ) = −J
X
<i,j>
σiσj ,
(58)
where the sum runs over pairs of adjacent sites, J is a coupling constant, and σi = ±1
denotes the local spin at site i. The equilibrium ensemble of the Ising model may be
generated by a dynamic process with synchronous dynamics mimicking the contact of the
system with a thermal reservoir. These dynamic rules – usually referred to as heat bath
dynamics – are deﬁned through the transition probabilities
pσ→σ′ =
Y
i
pi(σ) ,
pi(σ) =
ehi(σ)
ehi(σ) + e−hi(σ) ,
hi(σ) =
1
kBT
X
j
σj ,
(59)
where σ denotes the actual state of the model and j runs over the nearest neighbors of i.
In order to verify the coincidence of the stationary distribution P∞(σ) and the equilibrium
ensemble of the Ising model, it is suﬃcient to prove that the dynamic processes are ergodic
and obey detailed balance
Peq(σ)pσ→σ′ = Peq(σ′)pσ′→σ .
(60)
Detailed balance means that the probability currents between two states are exactly equal
in both directions, i.e., the currents cancel each other in the stationary state. Heat bath
dynamics is only one out of inﬁnitely many dynamic processes generating the Ising equilib-
rium ensemble. Examples include Glauber, Metropolis, and Kawasaki dynamics, as well
as the Swendsen-Wang and Wolf cluster algorithms. Although these stochastic models
have very diﬀerent dynamic properties, they all evolve towards the same stationary state
which is just the equilibrium state of the Ising model.
26

<!-- page 27 -->

2.11
Matrix product states
For the majority of reaction-diﬀusion models it is quite diﬃcult or even impossible to solve
the master equation analytically. In some cases, however, it is still possible to compute
the stationary state of the system.
In recent years a powerful algebraic approach has
been developed by which n-point correlation functions of certain nonequilibrium systems
can be computed exactly (see Ref. [29] for a general review). This approach generalizes
states with product measure to so-called matrix product states (MPS) by replacing real-
valued probabilities with non-commutative operators.
Representing these operators in
terms of matrices it is possible to compute correlation functions by evaluating certain
matrix products.
MPS for the exclusion process
In order to introduce the matrix product technique, let us consider the partially asymmet-
ric exclusion process in one spatial dimension with asynchronous dynamics and particle
adsorption (desorption) at the left (right) boundary. The Liouville operator of this system
is given by
L = S1 + SN +
N−1
X
i=1
Li ,
(61)
where Li describes the hopping in the bulk while S1 and SN are the surface contributions
for adsorption and desorption, respectively. In the standard basis (see Appendix A), these
operators read
Li =




0
0
0
0
0
q−1
−q
0
0
−q−1
q
0
0
0
0
0



,
S1 =
 α
0
−α
0

,
SN =
0
−β
0
β

.
(62)
The partially asymmetic exclusion process with particle input and output is particularly
interesting in presence of a current, that is, for q̸ = 1. Depending on α and β, the system
is in a phase of low, high, or maximal density (see Fig. 5). It seems to be quite surprising
that this simple model could be solved only a few years ago by recursion techniques [74,75]
and, at the same time, by using the matrix product method [57,76].
The solution is particularly simple along the coexistence line between the high and low
density phases α + β = q −q−1. Here the stationary state |Ps⟩may be written as
|Ps⟩= 1
Z
e
d

⊗. . . ⊗
e
d

=
e
d
⊗N
,
(63)
where e = 1/α and d = 1/β. This state has a product measure, i.e., all spatial corre-
lations vanish. The constant Z = (e + d)N normalizes the probability distribution. The
stationarity L|Ps⟩= 0 can be veriﬁed by proving the relations
S1
e
d

=
 1
−1

,
SN
e
d

= −
 1
−1

,
Li
e
d

⊗
e
d

= −
 1
−1

⊗
e
d

+
e
d

⊗
 1
−1

,
(64)
27

<!-- page 28 -->

which provide the following cancellation mechanism: First the action of S1 generates a
factor
  1
−1

at the leftmost position in the product (63). This factor is then commuted to
the right by successive action of L1, . . . , LN−1 until it is adsorbed at the right boundary
by the action of SN. This cancellation mechanism proves that the homogeneous product
state (63) is stationary along the coexistence line α + β = q −q−1.
For α+β̸ = q−q−1 the stationary state has no longer the form of a simple product state.
However, as will be shown below, it can be expressed as a matrix product state. To this end
we replace the probabilities e and d in Eq. (63) by noncommutative operators E and D.
These operators act in an auxiliary space which is diﬀerent from the conﬁgurational vector
space of the lattice model. Introducing boundary vectors ⟨W| and |V ⟩with ⟨W|V ⟩̸ = 0 a
matrix product state may be written as
|Ps⟩= 1
Z ⟨W|
E
D
⊗N
|V ⟩.
(65)
Notice that ⟨W| and |V ⟩are vectors in the auxiliary space while |Ps⟩denotes the sta-
tionary probability distribution in conﬁguration space. By selecting the matrix element
⟨W| . . . |V ⟩, products of the operators can be mapped to real-valued probabilities4. The
normalization constant is given by Z = ⟨W|CN|V ⟩, where C = D + E. The cancellation
mechanism
⟨W|S1
E
D

= ⟨W|
 1
−1

,
SN
E
D

|V ⟩= −
 1
−1

|V ⟩,
Li
E
D

⊗
E
D

= −
 1
−1

⊗
E
D

+
E
D

⊗
 1
−1

(66)
is the same as in Eq. (64), leading to the algebra
DE = D + E
(67)
and the boundary conditions
⟨W|E = α−1⟨W| ,
D|V ⟩= β−1|V ⟩.
(68)
This ansatz is not restricted to the exclusion process but can be applied to any reaction-
diﬀusion process with random-sequential updates. It converts the dynamic rules of the
model into a set of algebraic relations and boundary conditions. The problem of calculating
the stationary state is then shifted to the problem of ﬁnding a matrix representation of
the algebra. In the present case the quadratic algebra (67) can be mapped onto a bosonic
algebra for which an inﬁnite-dimensional matrix representation exists [76]. For particular
values of α and β, however, there are also ﬁnite-dimensional representations [77,78]. For
example, along the coexistence line α + β = q −q−1 the product state (63) is nothing
but a one-dimensional representation of the algebra. Moreover, it can be shown that for
α + β + αβq = q + q−1 a two-dimensional matrix representation is given by
E =
 
1
α
0
1
1+αq
αq2
!
,
D =
 1
β
−1
0
1+βq
βq2
!
,
⟨W| = (1 0) ,
|V ⟩=
1
0

.
(69)
4The matrix product technique can also be applied to models with periodic boundary conditions. In
this case the matrix elements ⟨W | . . . |V ⟩have to be replaced by a trace operation Tr[. . . ].
28

<!-- page 29 -->

0
20
40
60
80
100
i
0
0.5
1
ρi
stat
α=0.55
α=0.5
α=0.45
Fig. 8:
Exact stationary density proﬁle for the partially asymmetric exclusion process with particle
adsorption and desorption at the boundaries for q = 2 and N = 100 in the low density phase α = 0.45,
along the coexistence line α = 0.5, and in the high density phase α = 0.55.
For a given matrix representation the stationary particle concentration proﬁle ρstat
i
can be
computed by evaluating the matrix product
ρstat
i
= ⟨W|Ci−1DCN−i|V ⟩
⟨W|CN|V ⟩
.
(70)
In the case of the above 2 × 2 representation we obtain the exact result
ρstat
i
=
1
αqλN−1
1
−
αq
β2(1+αq)λN−1
2
+ λi−1
1
λN−i
2
β(1+αq)
αq
λN
1 −
αq
β(1+αq)λN
2
,
(71)
where
λ1 =
q −q−1
αβ(1 + αq) ,
λ2 =
q −q−1
αβ(1 + βq)
(72)
are the eigenvalues of C. The density proﬁle is shown in Fig. 8 for various values of α,
visualizing the transition between the low and the high density phase. Similarly, one can
compute two-point correlation functions
⟨sisj⟩= ⟨W|Ci−1DCj−i−1DCN−j|V ⟩
⟨W|CN|V ⟩
.
(73)
In this expression the matrix C plays the role of a transfer matrix between the sites i and j.
Therefore, the long-distance behavior of correlation functions in the bulk will be governed
by the largest eigenvalue of C. In particular, for any ﬁnite-dimensional representation
of the algebra, the correlation functions in the stationary state will decay exponentially.
Only inﬁnite-dimensional representations can lead to long-range correlations with power-
law decay in the stationary state.
It is also possible to apply the matrix product technique to the totally asymmetric
exclusion process with sublattice-parallel updates [79–81], where a diﬀerent cancellation
29

<!-- page 30 -->

mechanism is needed [82–85]. Recently, the matrix product method could even be extended
to the case of fully parallel updates [40].
MPS for models with particle reactions
Most models which have been solved so far by using the matrix product method are
diﬀusive systems, i.e., they describe stochastic transport of particles.
Up to now only
one exception is known where particles react with each other, namely the anisotropic
decoagulation model with closed boundary conditions and random sequential updates. As
will be explained in the following, this model exhibits a boundary-induced phase transition
and can be solved by using a generalized matrix product ansatz [86].
The anisotropic decoagulation model is deﬁned by the following dynamic rules:
diﬀusion:
ØA →
q AØ ,
AØ →
q−1 ØA ,
coagulation:
AA →
q AØ ,
AA →
q−1 ØA ,
decoagulation:
ØA →
κq AA ,
AØ →
κq−1 AA .
The corresponding Liouville operator reads
Li =




0
0
0
0
0
(κ + 1)q
−q−1
−q−1
0
−q
(κ + 1)q−1
−q
0
−κq
−κq−1
q + q−1



.
(74)
Thus, the model is controlled by two parameters, namely the anisotropy q and the decoag-
ulation rate κ. The phase diagram displays two phases, a low-density phase for κ < q2 −1
and a high-density phase for κ > q2 −1. From the physical point of view these phases are
diﬀerent from those observed in the asymmetric exclusion process since here the number
of particles is not conserved. Notice that the rates for diﬀusion and coagulation coincide.
Moreover, all reactions have the same bias. This special choice ensures that the model
is integrable. Various exact results have been obtained by using IPDF and free-fermion
techniques [27, 54, 64, 87]. At the critical point κc = q2 −1, the relaxational spectrum
becomes massless and algebraic long-range correlations can be observed [55].
In order to express the stationary state of the model as a matrix product state of
the form (65), the cancellation mechanism of Eq. (66) has to be generalized. This can
be achieved by replacing the vector
  1
−1

by

¯E¯D

, where ¯E and ¯D are two additional
operators acting in the auxiliary space. The generalized cancellation mechanism reads
⟨W|S1
E
D

= ⟨W|
 ¯E
¯D

,
SN
E
D

|V ⟩= −
 ¯E
¯D

|V ⟩,
Li
E
D

⊗
E
D

= −
 ¯E
¯D

⊗
E
D

+
E
D

⊗
 ¯E
¯D

.
(75)
This ansatz leads to the bulk algebra
0 = E ¯E −¯EE ,
(κ + 1)q ED −q−1DE −q−1DD = E ¯D −¯ED ,
−q ED + (κ + 1)q−1DE −q DD = D ¯E −¯DE ,
−κq ED −κq−1DE + (q + q−1) DD = D ¯D −¯DD ,
(76)
30

<!-- page 31 -->

2
2.5
3
3.5
4
κ
0
0.2
0.4
0.6
0.8
ρ=N
−1Σiρi
N=100
N=250
N=20
N=50
0
50
100
150
200
250
i
0
0.2
0.4
0.6
0.8
1
ρi
κ=3.0
κ=3.2
κ=10
κ=2.8
κ=0
Fig. 9:
The anisotropic coagulation-decoagulation model with bias q = 2. Left: Density proﬁles in a
system with N = 250 sites for various values of the decoagulation rate κ. Right: Density of particles
averaged over the entire system as a function of κ for various system sizes, illustrating the ﬁrst-order phase
transition.
and the boundary conditions
⟨W| ¯E = ⟨W| ¯D = ¯E|V ⟩= ¯D|V ⟩= 0 .
(77)
A trivial one-dimensional representation of this algebra is given by E = C = 1, ¯E = ¯C = 0,
describing a system without particles. In the symmetric case q = 1 there is another one-
dimensional representation E = 1, C = γ2, ¯E = ¯C = 0, corresponding to a homogeneous
product state with particle density ρ = κ/(1 + κ). For q̸ = 1 and κ̸ = q2 −1 we ﬁnd the
four-dimensional representation
E =




q−2
q2 −γ−2
q2 −1
q2(1 −γ2)
0
γ−2
0
γ2 −q2
0
0
1
γ2(q2 −1)
0
0
0
q2



,
C =




q−2
0
0
0
0
1
0
0
0
0
γ2
0
0
0
0
q2



,
⟨W| =

1
1 −q2γ2 , 0 ,
q2
q2γ2 −1 ,
a (q2 −q−2)(γ2 −q2)γ2 −q2γ2
(γ2 −1)(q2 + 1)

,
(78)
|V ⟩=
b (q4 −1)(q2γ2 −1) + q4
q2 + 1
, 0 ,
q2(γ2 −1)
γ2 −q2
,
(γ2 −1)q2
γ4 −γ2q2

.
Using this representation, it is easy to compute the stationary particle density
ρstat
j
=
γ2N

(γ2−1)+(q2−1)γ2(qγ)−2j

−q2N

(γ2−1)q2−4j+(q2−1)(q/γ)−2j

γ2 (γ2N+γ−2N −q2N−q−2N)
(79)
which is in agreement with the result obtained by using the IPDF method [55].
For
κ̸ = q2 −1 a similar four-dimensional matrix representation can be constructed.
In the thermodynamic limit the anisotropic coagulation model exhibits a ﬁrst-order
phase transition (see Fig. 9). If the decoagulation rate κ is small enough, the particles
are swept towards one of the boundaries where they coagulate. The stationary particle
density is therefore zero in the thermodynamic limit. Increasing κ this region grows until
its size diverges at a critical value κc = q2 −1. Above κc the decoagulation process is
31

<!-- page 32 -->

strong enough to maintain a non-vanishing density of particles in the bulk. It should be
emphasized that this type of phase transition is induced by the boundaries. In particular,
there is no such transition if periodic boundary conditions are used.
The matrix product technique has also been applied to various other systems such
as valence-bond-state models [88], spin-one quantum antiferromagnets [89,90], hard-core
diﬀusion of oppositely charged particles [91], systems with ﬁxed [92] or moving impuri-
ties [93,94], as well as n-state diﬀusion processes [95,96]. Furthermore, a dynamic matrix
product ansatz has been introduced by which time-dependent properties of the exclusion
process can be described [97,98]. Although the full range of possible applications is not yet
known, the matrix product technique seems to be limited to a few classes of models. By
deﬁnition the method is restricted to one-dimensional systems. Moreover, there seems to
be a subtle connection between matrix product states and integrability. This connection
is not yet fully understood. In Ref. [99] it was shown that the stationary state of any
reaction-diﬀusion model can be expressed in terms of a MPS. However, in this generic
case the corresponding matrix representation depends on the system size and is there-
fore useless from the practical point of view. A systematic classiﬁcation scheme for matrix
product states is not yet known. In this context it is interesting to note that the method of
dynamic density matrix renormalization allows ﬁnite-dimensional matrix representation to
be detected. As shown in Ref. [100], the existence of a ﬁnite-dimensional MPS is indicated
by the fact that the density matrix has only a ﬁnite number of nonvanishing eigenvalues.
Thus, by scanning the spectrum over a certain range of the systems parameter space, it
is possible to search systematically for ﬁnite-dimensional matrix product representations.
32

<!-- page 33 -->

3
Directed percolation
Spreading processes are encountered in many diﬀerent situations in nature as diverse as
epidemics [101], forest ﬁres [102], and transport in random media [103, 104]. Spreading
phenomena are usually characterized by two competing processes. For example, in an
infectious disease the spreading agent (bacteria) may multiply and infect neighboring
individuals. On the other hand, infected individuals may recover, decreasing the total
amount of the spreading agent. Depending on the relative rates for infection and recovery,
two diﬀerent situations may emerge.
If the infection process dominates, the epidemic
disease will spread over the entire population, approaching a stationary state in which
infection and recovery balance one another. However, if recovery dominates, the total
amount of the spreading agent continues to decrease and eventually vanishes.
Theoretical interest in models for spreading stems mainly from the emerging phase
transition between survival and extinction. The simplest model exhibiting such a transi-
tion is directed percolation (DP)5. In DP sites of a lattice can either be active (infected)
or inactive (healthy). Depending on a parameter controlling the balance between infection
and recovery, activity may either spread over the entire system or die out after some time.
In the latter case the system becomes trapped in a completely inactive state, the so-called
absorbing state of the model. Since the absorbing state can only be reached but not be
left, it is impossible to obey detailed balance, i.e., DP is a nonequilibrium process. The
transition between the active and the inactive phase is continuous and characterized by
universal critical behavior.
In many respects, the nonequilibrium critical behavior of DP is similar to that of
equilibrium models. In particular, it is possible to use the concept of scale invariance and to
identify various critical exponents. As in equilibrium statistical mechanics, these exponents
allow phase transitions of diﬀerent lattice models to be categorized into universality classes.
As we will see below, DP is a fundamental class of nonequilibrium phase transitions,
playing a similar role as the Ising universality class in equilibrium statistical mechanics.
Although DP is very robust and easy to simulate, its critical behavior turns out to be
highly nontrivial. This is what makes DP so fascinating.
In this Section we give a general introduction to DP, discussing the most important lat-
tice models, basic scaling concepts, ﬁnite-size properties, as well as mean-ﬁeld approaches.
We also summarize various approximation techniques such Monte Carlo simulations, series
expansions, and numerical diagonalization. Introducing basic ﬁeld-theoretic methods we
discuss the critical behavior of DP at surfaces, the early-time behavior, the inﬂuence of
fractal initial conditions, persistence probabilities, and the inﬂuence of quenched disorder.
Finally we review possible experimental realizations of DP and discuss the question why
it is so diﬃcult to verify the critical exponents in experiments.
5 For further reading in this ﬁeld we recommend the review on directed percolation by Kinzel [35,105],
a summary of open problems by Grassberger [106], and the relevant chapters in a recent book by Marro
and Dickman [107].
33

<!-- page 34 -->

directed bond percolation
isotropic bond percolation
Fig. 10:
Isotropic and directed bond percolation on a diagonal square lattice with free boundary conditions.
Open (closed) bonds are represented by solid (dotted) lines. In both cases a cluster, indicated by bold
lines, is generated from the lattice site in the center. In the directed case the spreading agent is restricted
to follow sense of the arrows, leading to a directed cluster of connected sites.
3.1
Directed percolation as a spreading process
From isotropic to directed percolation
Although DP is often regarded as a dynamic process, it was originally deﬁned as a geo-
metric model for connectivity in random media which generalizes isotropic (undirected)
percolation [108,109]. Such a random medium could be a porous rock in which neighbor-
ing pores are connected by channels of varying permeability. An important question in
geology would be how deep the water can penetrate into the rock.
In ordinary percolation the water propagates isotropically in all directions of space.
One of the simplest models for isotropic percolation is bond percolation on a d-dimensional
square lattice, as shown in the left part of Fig. 10. In this model the channels of the porous
medium are represented by bonds between adjacent sites of a square lattice which are open
with probability p and otherwise closed6. For simplicity it is assumed that the states of
diﬀerent bonds are uncorrelated. Clearly, if p is suﬃciently large, the water will percolate
through the medium over arbitrarily long distances. However, if p is small enough, the
penetration depth is expected to be ﬁnite so that large volumes of the material will be
impermeable. Both regimes are separated by a continuous phase transition.
The left part in Fig. 10 shows a typical conﬁguration of open and closed bonds in
two dimensions. As can be seen, each site generates a certain cluster of connected sites
corresponding to the maximal spreading range if the water was injected into a single pore.
The site from where the cluster is generated is called the origin of the cluster. The size
(or mass) of a cluster is deﬁned as the number of connected sites. Notice that diﬀerent
origins may generate the same cluster. Consequently the whole lattice decomposes into a
set of disjoint clusters.
6Alternatively, we could have blocked sites instead of bonds with a certain probability. The resulting
model, called site percolation, exhibits the same type of universal critical behavior at the transition.
34

<!-- page 35 -->

Directed percolation, introduced in 1957 by Broadbent and Hammersley [110], is an
anisotropic variant of isotropic percolation. As shown in the right part of Fig. 10, this
variant introduces a speciﬁc direction in space. The channels (bonds) function as ‘valves’
in a way that the spreading agent can only percolate along the given direction, as indicated
by the arrows. For example, we may think of a porous medium in a gravitational ﬁeld
that forces the water to propagate downwards7. Thus, ﬁlling in the spreading agent at a
particular site, the resulting cluster of wet sites is a subset of the corresponding cluster in
the isotropic case (see Fig. 10). Since in DP each site generates an individual cluster, a
decomposition of the lattice into disjoint clusters is no longer possible. As in the case of
isotropic percolation, DP exhibits a continuous phase transition.
The phase transitions of isotropic and directed percolation are similar in many re-
spects. They both can be characterized by an order parameter P∞which is deﬁned as the
probability that a randomly selected site generates an inﬁnite cluster. If P∞is ﬁnite, the
spreading agent is able to percolate over arbitrarily long distances wherefore the system
is said to be in the wet phase. If P∞vanishes, the system is in the so-called dry phase
where the spreading range is ﬁnite. Both isotropic and directed percolation are trivial in
one dimension: Since an inﬁnite cluster on a line requires all bonds to be open, the wet
phase consists of a single point p = 1. Another trivial case is the limit of inﬁnitely many
dimensions. Since each site is connected with inﬁnitely many neighbors, an inﬁnite cluster
will be generated for any p > 0. Consequently the inactive phase consists of single point
in phase space, namely p = 0. In ﬁnite dimensions 2 ≤d < ∞there is a continuous phase
transition separating the wet phase from the dry phase at some critical value 0 < pc < 1.
In the supercritical phase p > pc the medium is permeable (P∞> 0) while in the subcrit-
ical phase p < pc the medium becomes impermeable (P∞= 0). As expected, the critical
threshold pc for directed percolation is larger than in the isotropic case.
Although isotropic and directed percolation have several common features, their critical
behavior near the phase transition turns out to be diﬀerent. In the isotropic case, the
critical properties are in all directions the same (apart from lattice eﬀects which are usually
irrelevant on large scales) and hence the emerging long-range correlations are rotationally
invariant. Because of a duality symmetry, the critical point of isotropic bond percolation
is pc = 1/2 [109]. Moreover, the critical exponents are given by simple rational numbers.
Contrarily, the critical properties of DP reﬂect the anisotropy in space, leading to diﬀerent
critical exponents. In contrast to the isotropic case, the numerical values of the critical
point and the exponents of DP are not yet known analytically and seem to be given by
irrational numbers (see Sec. 3.4).
Interpretation of directed percolation as a dynamic process
Regarding the given direction as ‘time’, directed percolation may be interpreted as a d+1-
dimensional dynamic process that describes the spreading of some non-conserved agent8.
For example, as already mentioned before, DP may be viewed as a simple model for
epidemic spreading of some infectious disease without immunization [111]. In recent years
the dynamic interpretation has become increasingly popular, partly because the time-
dependent formulation is the natural realization of DP on a computer. In what follows
7This assumption is highly idealized since water is a conserved quantity. Moreover, the water can even
ﬂow against the gravity ﬁeld (see Sec. 3.9).
8 Note that DP diﬀers from ‘dynamic percolation’ which is used as a epidemic processes with immu-
nization (see Sec. 3.8.
35

<!-- page 36 -->

...
1
2
N
initial configuration
final configuration
t=3
t=2
t=1
t=0
i=
...
Fig. 11:
Directed bond percolation in 1+1 dimensions interpreted as a time-dependent stochastic process.
Open (closed) bonds are indicated by solid (dashed) lines. Filled (hollow) circles denote active (inactive)
sites. The conﬁguration of the horizontal row at t = 0 is the initial state. Starting from a fully occupied
initial state the model ‘evolves’ through intermediate conﬁgurations according to the dynamic rules of
Eq. (80) and reaches a ﬁnal state at t = 3.
we will adopt the dynamic interpretation. However, one should keep in mind that the
geometric deﬁnition in terms of directed paths is fully equivalent.
The dynamic interpretation of DP is illustrated in Fig. 11, where the lattice sites of
a (1+1)-dimensional directed bond percolation process are enumerated horizontally by a
spatial coordinate i and vertically by a discrete time variable t. Since we use a diagonal
square lattice, odd and even time steps have to be distinguished. A local binary variable
si(t) is attached to each site. si = 1 means that the site is active (occupied, wet) while
si = 0 denotes an inactive (empty, dry) site. The set s = {si} at a given time t speciﬁes
the conﬁguration of the system.
For a given conﬁguration at time t, the next conﬁguration at time t + 1 can be de-
termined as follows.
For each pair of bonds between the sites (i ± 1, t) and (i, t + 1)
two random number z±
i
∈(0, 1) are generated. A bond is considered to be open (with
probability p) if z±
i < p, leading to the update rule
si(t + 1) =



1
if si−1(t) = 1 and z−
i < p ,
1
if si+1(t) = 1 and z+
i < p ,
0
otherwise .
(80)
Thus, directed bond percolation can be considered as a Markov process with parallel
dynamics. As in any dynamic system, we have to specify the initial state. Common initial
states are the fully occupied lattice, random initial conditions, and conﬁgurations with a
single particle at the origin (also called ‘active seed’).
Even very simple numerical simulations demonstrate that the temporal evolution of a
DP process changes signiﬁcantly at the phase transition. Typical space-time histories for
random initial conditions are shown in the upper part of Fig. 12. For p < pc the number
of particles decreases exponentially until the system reaches the absorbing state, whereas
for p > pc the average particle number saturates at some constant value. At the critical
point the mean particle number decays very slowly and the emerging clusters of active
sites remind of fractal structures. A similar behavior can be observed if the DP process
starts from a single seed (see lower part of Fig. 12). For p < pc the average number of
particles ﬁrst grows for a short time and then decays exponentially. For p > pc there is a
ﬁnite probability that the resulting cluster is inﬁnite. In this case activity spreads within
a certain triangular region, the so-called spreading cone. At p = pc a critical cluster is
generated from a single seed, whose scaling properties will be discussed in Sec. 3.3.
It is often helpful to regard DP as a reaction-diﬀusion process of interacting particles.
36

<!-- page 37 -->

time t
position i
pc
pc
c
p
p>
p=
p<
Fig. 12:
Directed bond percolation in 1+1 dimensions starting from random initial conditions (top) and
from a single active seed (bottom). Each horizontal row of pixels represents four updates. As can be seen,
critical DP is a reaction-limited process.
Associating active sites with particles A and inactive sites with vacancies Ø, a DP process
corresponds to the reaction-diﬀusion scheme
self-destruction:
A →Ø ,
diﬀusion:
Ø + A →A + Ø ,
oﬀspring production:
A →2A ,
coagulation:
2A →A .
(81)
To understand this reaction-diﬀusion scheme, let us again consider the example of directed
bond percolation. Depending on the conﬁguration of the bonds, each active site (particle)
may activate two neighboring sites of the subsequent row (next time step). If both bonds
are closed, the particle self-destructs. If only one bond is open, the particle will diﬀuse
to the left or to the right with equal probability, whereas an oﬀspring is produced when
both bonds are open. On the other hand, if two particles reach the same target site, they
coalesce into a single particle, giving rise to the reaction 2A →A. This process limits the
maximal density of active sites. In fact, as will be shown below, the coagulation process is
the essential nonlinear ingredient of DP. In ‘fermionic’ models with an exclusion principle
it is automatically included. However, in ‘bosonic’ models allowing for an inﬁnite number
of particles per site one would have to add this process explicitly.
3.2
Lattice models for directed percolation
In the literature there is a vast variety of DP models following the spirit of the above
reaction-diﬀusion scheme. As we will see below, they all exhibit the same type of critical
behavior at the transition. The common feature of all these models is the existence of an
absorbing state, i.e., a conﬁguration that the model can reach but from where it cannot
37

<!-- page 38 -->

p
2
p
1
s
p
1
i
t+1
t
i-1
s
i+1
s
Fig. 13:
Transition probabilities in the (1+1)-dimensional Domany-Kinzel cellular automaton.
escape. In most cases, the absorbing state is just the empty lattice. The existence of
an absorbing state implies that certain microscopic processes are forbidden (for example,
spontaneous creation of particles Ø →A). In the sequel we will discuss three examples,
namely the Domany-Kinzel cellular automaton, the contact process, and the Ziﬀ-Gulari-
Barshad model for heterogeneous catalysis.
The Domany-Kinzel cellular automaton
Cellular automata are algorithms that map a conﬁguration of a lattice onto a new conﬁgu-
ration. The automaton evolves in time by iteration of the map. Thus, the time variable t is
discrete. Usually the map can be decomposed into independent local updates. Since these
updates can be processed simultaneously, cellular automata can eﬃciently be implemented
on parallel computers. Depending on the type of updates, we distinguish between deter-
ministic and stochastic cellular automata. A general classiﬁcation of stochastic cellular
automata was presented by Wolfram [112].
Various stochastic cellular automata are known to exhibit a DP transition from a
ﬂuctuating phase into an absorbing state.
One of the simplest models in this class is
the (1+1)-dimensional Domany-Kinzel (DK) model [105,113]. It is deﬁned on a diagonal
square lattice and evolves by parallel updates according to certain conditional transition
probabilities P[si(t + 1)|si−1(t), si+1(t)]. These probabilities depend on two parameters
and are deﬁned by
P[1|0, 0] = 0 ,
P[1|0, 1] = P[1|1, 0] = p1 ,
(82)
P[1|1, 1] = p2 ,
where P[0|·, ·] = 1 −P[1|·, ·]. The corresponding update scheme may be realized by the
following algorithm (see Fig. 13): For each site i we generate a uniformly distributed
random number zi(t) ∈(0, 1) and set
si(t + 1) =





1
if si−1(t)̸ = si+1(t)
and zi(t) < p1 ,
1
if si−1(t) = si+1(t) = 1 and zi(t) < p2 ,
0
otherwise .
(83)
In contrast to directed bond percolation, the DK model depends on two percolation prob-
abilities p1 and p2. The corresponding phase diagram is shown in Fig. 14. It comprises
an active and an inactive phase, separated by a phase transition line (the solid line in the
ﬁgure). In the active phase a ﬂuctuating steady state exists on the inﬁnite lattice, whereas
in the inactive phase the model always reaches the absorbing state. The DK model in-
cludes three special cases.
The previously discussed case of directed bond percolation
corresponds to the choice p1 = p and p2 = p(2 −p). Another special case is directed site
38

<!-- page 39 -->

0
1
p1
0
1
p2
INACTIVE PHASE
ACTIVE
PHASE
compact DP
bond DP
site DP
W18
Fig. 14:
Phase diagram of the Domany-Kinzel model.
transition point
p1,c
p2,c
Ref.
Wolfram rule 18
0.801(2)
0
[114]
site DP
0.705489(4)
0.705489(4)
[115]
bond DP
0.6447001(1)
0.8737620(2)
[116]
compact DP
1/2
1
[35]
Table 1:
Special transition points in the (1+1)-dimensional Domany-Kinzel model.
percolation [35], corresponding to the choice p1 = p2 = p. The third special case p2 = 0 is
equivalent to the rule ‘W18’ of Wolfram’s classiﬁcation scheme [112]. Numerical estimates
for the corresponding critical points are summarized in Table 1.
There is strong numerical evidence that the critical behavior along the whole phase
transition line (except for its upper terminal point) is that of DP. This means that all
these transition points exhibit the same type of long-range correlations. The short-range
correlations, however, are non-universal and may change when moving along the phase
transition line. In order to understand the signiﬁcance of short-range correlations from
the physical point of view, let us consider spatial conﬁgurations (snapshots) of a critical
DP cluster. Such conﬁgurations are typically characterized by localized spots of activity
separated by large voids in between. Approaching the phase transition line the average
size of inactive voids diverges, whereas the mean size of active spots remains ﬁnite and
converges to a certain value ⟨Sact⟩. By moving along the transition line of the DK model
the asymptotic average size ⟨Sact⟩of active spots varies. As shown in Fig. 15, it is minimal
for p2 = 0, grows monotonically with p2, and ﬁnally diverges at the terminal point p2 = 1,
where the system crosses over to a diﬀerent type of critical behavior. Thus, by moving
along the phase transition line the nonuniversal short-range properties change while the
long-range properties remain unaﬀected.
The exceptional behavior at the upper terminal point of the phase transition line is
due to an additional symmetry between active and inactive sites along the line p2 = 1 [35].
Here the DK model has two symmetric absorbing states, namely the empty and the fully
39

<!-- page 40 -->

0.5
0.6
0.7
0.8
p1
 0
 5
 10
 15
 20
<sact>
bond DP
compact DP
Wolfram rule 18
Fig. 15:
Numerical estimates for the average size of active spots ⟨Sact⟩in the Domany-Kinzel model
measured along the phase transition line. The insets show typical clusters for three special cases discussed
in the text.
occupied lattice. As the corresponding symmetry transformation maps p1 to 1 −p1, the
phase transition line must end in the terminal point (p1, p2) = (1/2, 1). As shown in the
corresponding inset of Fig. 15, the resulting clusters of active sites are compact. There-
fore, this special case is referred to as compact directed percolation (CDP) [117]. Unfor-
tunately this expression is misleading since CDP stands for a universality class which is
completely diﬀerent from DP. As can be veriﬁed easily, the DK model at the terminal
point is equivalent to the (1+1)-dimensional voter model [5] or the Glauber-Ising model
at zero temperature. Alternatively, one may describe CDP as a pair-annihilation process
of diﬀusing kinks separating inactive and active domains (cf. Sec. 2.4). We may there-
fore identify the critical behavior of CDP with the exactly solvable universality class of
diﬀusing-annihilating random walks [44]. We will come back to CDP in Sec. 3.8.
In more than one spatial dimension the DK model may be deﬁned by local updates with
conditional probabilities P(si(t + 1)|ni(t)) depending on the number ni(t) = P
j∈<i> sj(t)
of active neighboring sites. Thus, the model is controlled by 2d parameters p1, . . . , p2d:
P[1|0] = 0 ,
P[1|n] = pn .
(1 ≤n ≤2d)
(84)
Notice that for d = 1 this deﬁnition is compatible with the usual deﬁnition of the DK
model in Eq. (82). The special case of directed bond percolation corresponds to the choice
pn = 1 −(1 −p)n while for equal parameters pn = p one obtains directed site percolation
in d+1 dimensions.
The contact process
Another important lattice model for DP is the contact process. In contrast to cellular
automata this model uses asynchronous updates. The contact process was ﬁrst introduced
by Harris [111] as a model for epidemic spreading without immunization (for a review see
Ref. [107]).
Here the lattice sites represent infected and healthy individuals. Infected
40

<!-- page 41 -->

λ/2
λ/2
1/2
1/2
1/2
1/2
is
is
is
i+1
s
λ/2
λ/2
λ
1
i-1
s
i+1
s
t
t+dt
t
t+dt
Fig. 16:
Stochastic processes in the (1+1)-dimensional contact process. Infected sites (bold dots) infect
their neighbors at rate λ/2 and recover at rate 1. The upper part shows the deﬁnition of the rules according
to Eq. (85). The lower part shows an equivalent deﬁnition as a two-site process (see text).
individuals can either heal themselves or infect their nearest neighbors. Depending on
the relative rates of infection and recovery, the epidemic disease may either spread over
the whole population or vanish after some time. In contrast to the DK model, infection
and healing processes are assumed to occur spontaneously without correlation in space
and time, i.e. spatially separated processes are not synchronized. To mimic this kind
of asynchronous dynamics, the contact process uses random sequential instead of parallel
updates (cf. Sec. 2.2).
The contact process is deﬁned on a d-dimensional square lattice whose sites can be
either active (si(t) = 1) or inactive (si(t) = 0).
For each attempted update a site i
is selected at random. Depending on its state si(t) and the number of active neighbors
ni(t) = P
j∈<i> sj(t) a new value si(t+dt) = 0, 1 is assigned according to certain transition
rates w[si(t) →si(t + dt), ni(t)]. In the standard contact process these rates are deﬁned
by
w[0 →1, n] = λn/2d ,
w[1 →0, n] = 1 .
(85)
Here the parameter λ controls the infection rate and plays the role of the percolation
probability. For the (1+1)-dimensional case the dynamic processes are sketched in the
upper row of Fig. 16. Monte Carlo simulations and series expansions suggest that the phase
transition in 1+1 dimensions takes place at the critical point λc ≃3.29785(8) [118–120].
Whereas computational physicist often prefer the DK model for simulations, the con-
tact process is more popular in the mathematical community because it is easy to write
down the corresponding master equation. Following the notation of Sec. 2.2, the master
equation of the 1+1 dimensional contact process with periodic boundary conditions is
given by
∂tPt(s1, . . . , sN) =
N
X
i=1
(2si −1)

λsi−1 Pt(s1, . . . , si−2, 1, 0, si+1, . . . , sN)
+ λsi+1 Pt(s1, . . . , si−1, 0, 1, si+2, . . . , sN)
−Pt(s1, . . . , si−1, 1, si+1, . . . , sN)
	
,
(86)
where Pt(s1, . . . , sN) denotes the probability to ﬁnd the system at time t in the conﬁgura-
tion {s1, . . . , sN}. As can be seen, the master equation involves only two-site interactions,
41

<!-- page 42 -->

as illustrated in the lower row of Fig. 16. Using the vector notation of Eq. (9) the corre-
sponding Liouville operator LCP = P
i Li is given by
Li(λ) = 1
2




0
−1
−1
0
0
1 + λ
0
−1
0
0
1 + λ
−1
0
−λ
−λ
2



.
(87)
Notice that this operator does not satisfy simple algebraic relations as in Eq. (41), indicat-
ing that DP is a highly nonintegrable process. Finite-size spectra of LCP will be analyzed
in Sec. 3.4.
The Ziﬀ-Gulari-Barshad model for heterogeneous catalysis
Many catalytic reactions such as the oxidation of carbon monoxide on a platinum surface
mimic the reaction scheme of directed percolation. The key property of these reactions
is the existence of catalytically poisoned states where the system becomes trapped in a
frozen state. Thus, poisoned states play the role of absorbing conﬁgurations. A simple
model for surface catalysis of the chemical reaction CO + O →CO2 was introduced in
1986 by Ziﬀ, Gulari, and Barshad (ZGB) [121]. The model describes a gas composed of
CO and O2 molecules with ﬁxed concentrations y and 1−y, respectively, which is brought
into contact with a catalytic material. The catalytic surface is represented by a square
lattice whose sites can be either vacant (Ø), occupied by a CO molecule, or occupied by an
O atom. The ZGB model evolves by random sequential updates according to the following
probabilistic rules:
1. CO molecules ﬁll any vacant site at rate y.
2. O2 molecules dissociate on the surface into two O atoms and ﬁll pairs of adjacent
vacant sites at rate 1 −y.
3. Neighboring CO molecules and O atoms recombine instantaneously to CO2 and
desorb from the surface, leaving two vacancies behind.
On the lattice the three processes correspond to the reaction scheme
Ø →CO
at rate y ,
Ø + Ø →O + O
at rate 1 −y ,
(88)
O + CO →Ø + Ø
at rate ∞.
Clearly, this reaction is irreversible and thus the dynamic processes do not obey detailed
balance. Moreover, if the whole lattice is covered either with pure CO or O, the system
is trapped in a poisoned absorbing state. As shown in Fig. 17, the ZGB model can be
in three diﬀerent phases. For y < y1 ≃0.389 the system evolves into the O-poisoned
state whereas for y > y2 ≃0.525 it always reaches the CO-poisoned state. In between
the model is catalytically active. The model exhibits two diﬀerent phase transitions, a
continuous one at y = y1 and a discontinuous one at y = y2.
Grinstein et al. [122]
expected the continuous transition to belong to the DP universality class. In order to
verify this hypothesis, extensive numerical simulations were performed. Initially it was
believed that the critical exponents were diﬀerent from those of DP [123], while later
42

<!-- page 43 -->

y1
y2
0
1
y
0.0
0.5
1.0
0.0
0.0
0.0
0.0
[O]
[CO]
DP transition
first order
transition
Fig. 17:
Schematic phase diagram of the Ziﬀ-Gulari-Barshad model. The concentrations of oxygen (solid
line) and carbonmonoxide (dashed line) are plotted versus the CO adsorption rate y.
the transition at y = y1 was found to belong to DP [124]. Very precise estimates of the
critical exponents were recently obtained in Ref. [125], conﬁrming the existence of a DP
transition in the ZGB model. DP exponents were also obtained in a simpliﬁed version of
the ZGB model [126]. However, so far it has been impossible to observe DP exponents in
experiments. We will come back to this problem in Sec. 3.9.
3.3
Phenomenological scaling theory
In equilibrium statistical physics continuous phase transitions are usually characterized
by universal scaling laws. For example, the magnetization order parameter in the ordered
phase of the two-dimensional Ising models vanishes close to the critical point as |T −
Tc|β, where β is a universal exponent. Similarly the correlation length ξ, which is the
characteristic macroscopic length scale of the model, diverges as ξ ∼|T −Tc|−ν.
At
the critical point the correlation length is inﬁnite, i.e., there is no macroscopic length
scale in the system.
As a consequence, the system is invariant under suitable scaling
transformations. It turns out that a very similar picture emerges in the nonequilibrium
case. In the following we introduce a phenomenological scaling theory that can be applied
to DP and other types of phase transitions into absorbing states.
The critical exponents β, ν⊥, and ν∥
The order parameter of a spreading process is the density of active sites
ρ(t) = ⟨1
N
X
i
si(t)⟩,
(89)
where ⟨. . . ⟩denotes the ensemble average. Let us ﬁrst consider the case of an inﬁnite
system. In the active phase ρ(t) decays and eventually saturates at some stationary value
ρstat. The stationary density varies continuously with p −pc and vanishes at the critical
point (see Fig. 18). Close to the transition the order parameter varies according to a power
43

<!-- page 44 -->

0
0.02 0.04 0.06 0.08
0.1
p−pc
0
0.2
0.4
0.6
0.8
ρ
stat
10
−4
10
−3
10
−2
10
−1
p−pc
10
−1
10
0
ρ
stat
d=1
d=2
d=1
d=2
Fig. 18:
Stationary density ρstat in the active phase of directed bond percolation in 1+1 and 2+1
dimensions, plotted in a linear and a double-logarithmic representation.
law
ρstat ∼(p −pc)β ,
(90)
where β is the critical exponent associated with the particle density. In a double-logarith-
mic representation the power law behavior manifests itself as a straight line with slope β.
As can be seen in Fig. 18, the value of β depends on the dimensionality of the system.
The numerical value β ≃0.277 in 1+1 dimensions is comparatively small, indicating a
signiﬁcant change of ρstat near the transition. In 2+1 dimensions a larger value β ≃0.58
is observed.
In addition, spreading processes are characterized by certain correlation lengths. In
contrast to equilibrium models without any dynamical aspect, nonequilibrium critical
phenomena involve ‘time’ as an additional dimension. Since ‘time’ and ‘space’ are diﬀerent
in character, we have to distinguish spatial and temporal properties, denoting them by
the indices ⊥and ∥, respectively. In fact, nonequilibrium phase transitions are usually
characterized by two independent correlation lengths, namely a spatial length scale ξ⊥and
a temporal length scale ξ∥. Close to the transition, these length scales are expected to
diverge as
ξ⊥∼|p −pc|−ν⊥,
ξ∥∼|p −pc|−ν∥
(91)
with generally diﬀerent critical exponents ν⊥and ν∥. In the scaling regime the two corre-
lation lengths are related by ξ∥∼ξz
⊥, where z = ν∥/ν⊥is the so-called dynamic exponent.
In many models the triplet (β, ν⊥, ν∥) is the fundamental set of bulk exponents that
labels the universality class. Other critical exponents are usually related to these three
exponents by simple scaling relations (see below). Nonequilibrium phase transitions in
diﬀerent physical systems are believed to belong to the same universality class if their
critical exponents coincide9. In fact, the DK model, the contact process, the ZGB model,
and a vast variety of other DP models are characterized by the same triplet of exponents.
Fig. 19 illustrates the physical meaning of the correlation lengths ξ⊥and ξ∥. As in
equilibrium statistical mechanics, they are present below and above the critical point. In
9In addition, it should be proven that universal scaling functions coincide as well.
44

<!-- page 45 -->

ξ
||ξ
||ξ
ξ
ξ||
ξ
D
||
C
x
p>pc
p<pc
p>pc
p<pc
ξ
A
B
t
Fig. 19:
Interpretation of the correlation lengths ξ⊥and ξ∥in an almost critical (1+1)-dimensional DP
process below (left) and above criticality (right). In panels A and B a cluster is grown from a single active
seed while in panel C a fully occupied lattice is used as initial state. Panel D shows a stationary DP
process in the active phase. The indicated length scales ξ⊥and ξ∥must be interpreted as averages over
many independent realizations.
the inactive phase clusters originating from a single seed have the typical form of a droplet
(panel A). Averaging over many independent realizations the lateral size and the lifetime
of such droplets are proportional to ξ⊥and ξ∥, respectively. Above criticality the surviving
clusters grow within a spreading cone (panel B) whose opening angle is determined by the
ratio ξ⊥/ξ∥. The correlation lengths can also be seen if homogeneous initial conditions
are used. In the inactive phase the scaling length ξ∥plays the role of a typical decay
time (panel C), while in the stationary state of the active phase the correlation lengths
appear as the average sizes of inactive islands (panel D). This interpretation can be easily
generalized to higher dimensions.
As suggested by the scaling properties of the density (90) and the correlation lengths
(91), a spreading process should be invariant under dilatation x →Λx accompanied by
an appropriate rescaling of time and the deviation from criticality ∆= p −pc:
x →Λx ,
t →Λzt ,
∆→Λ−1/ν⊥∆,
ρ →Λ−β/ν⊥ρ .
(92)
This allows scale-invariant combinations to be constructed such as t/xz, ∆t1/ν∥and
∆x1/ν⊥. As we will see below, universal scaling functions can only depend on such scale-
invariant ratios.
Scaling theory for phase transitions into absorbing states
So far we have seen that the stationary density in the active phase scales as ρstat ∼∆β,
where ∆= p −pc denotes the distance from the critical point. A very similar quantity
is the ultimate survival probability P∞that a randomly chosen site belongs to an inﬁnite
cluster (cf. Sec. 3.1). In the active phase this probability is ﬁnite and scales as
P∞∼∆β′
(93)
with some critical exponent β′.
Although β and β′ coincide in the case of DP, they
may be diﬀerent in more general contexts, for example, in models with many absorbing
states. Therefore, phase transitions into absorbing states are generally described by four
exponents β, β′, ν⊥, and ν∥. The diﬀerent roles of β and β′ become apparent in a ﬁeld-
theoretic formulation (see Sec. 3.5). It turns out that β is associated with the particle
annihilation operator. Therefore, this exponent emerges whenever a particle density is
45

<!-- page 46 -->

‘measured’ in some ﬁnal state. The exponent β′, on the other hand, is associated with
the creation of particles and thus plays a role whenever particles are ‘introduced’. This
happens, for example, if an initial conﬁguration is speciﬁed.
In correlation functions,
which involve creation as well as annihilation operators, both exponents are expected to
appear.
Turning to time-dependent scaling properties in an inﬁnitely large system, there are
two important complementary quantities, namely the particle density ρ(t) starting from a
fully occupied lattice, and the survival probability P(t) that a cluster grown from a single
seed is still active after t time steps. Following the usual scaling concept of equilibrium
statistical mechanics, both quantities are expected to scale as
ρ(t) ≃t−αf(∆t1/ν∥) ,
P(t) ≃t−δg(∆t1/ν∥) ,
(94)
where α and δ are certain critical exponents for decay and survival, respectively [127,128].
f and g are universal scaling functions, i.e., they have the same functional form in all DP
models. For small arguments they both tend to a constant, whereas for large arguments
they scale in a way that the time dependence drops out:
f(ζ) ∼ζαν∥,
g(ζ) ∼ζδν∥.
(ζ →∞)
(95)
In the active phase the two quantities therefore saturate at ρstat = ρ(∞) ∼∆αν∥and
P∞= P(∞) ∼∆δν∥. Comparison with Eqs. (90) and (93) yields
α = β/ν∥,
δ = β′/ν∥.
(96)
An important quantity that combines both creation and annihilation of particles is the pair
connectedness function c(x′, t′, x, t), which is deﬁned as the probability that the sites (x′, t′)
and (x, t) are connected by a directed path of open bonds. Since the pair connectedness
function is translationally invariant in space and time, we may also write c(x′, t′, x, t) ≡
c(x−x′, t−t′). Starting from an initial condition with a single active site at the origin x′ =
0, the pair connectedness function c(x, t) is just the density of active sites in the resulting
clusters averaged over many realizations of randomness. Because of scaling invariance, the
pair connectedness function c(x, t) obeys the scaling form [127]
c(x, t) ∼tθ−d/z F(x/t1/z, ∆t1/ν∥) ,
(97)
where d denotes the spatial dimension and z = ν∥/ν⊥. The so-called critical initial slip
exponent θ describes the growth of the average number of particles as a function of time
(see Sec. 3.6). In order to determine θ we note that in the active phase ∆> 0 surviving
clusters will create an average density ρstat ∼∆β in the interior of the spreading cone (cf.
panel B of Fig. 19). Thus the autocorrelation function c(0, t) should saturate at the value
c(0, ∞) = lim
t→∞c(0, t) ∼∆β+β′ .
(98)
On the other hand, the scaling form (97) implies that c(0, t) saturates in the active phase
at a constant with the scaling behavior10
c(0, ∞) ∼∆ν∥(d/z−θ) .
(99)
10To prove this relation, notice that F(0, ζ) ∼ζ−(θ−d/z)ν∥for large values of ζ.
46

<!-- page 47 -->

Comparing the two expressions we obtain the generalized hyperscaling relation [129] for
phase transitions into absorbing states
θ −d
z = −β + β′
ν∥
.
(100)
It should be noted that the scaling argument (98) relies on the assumption that the cluster
spreads around the origin, i.e., the spreading cone surrounds the origin. In suﬃciently high
spatial dimensions, however, the cone becomes sparse and diﬀuses away from the origin
so that the autocorrelation function c(0, ∞) vanishes. For example, in a DP process this
happens above the upper critical dimension dc = 4. In fact, the generalized hyperscaling
relation (100) turns out to be valid only below the upper critical dimension of the spreading
process under consideration.
The scaling theory outlined above assumes the system size to be inﬁnite. For ﬁnite
systems sizes the scaling functions also depend on the invariant ratio ξd
⊥/N = td/z/N,
where N = Ld is the total number of sites. The generalized scaling forms read
ρ(t)
∼
t−β/ν∥f(∆t1/ν∥, td/z/N) ,
(101)
P(t)
∼
t−β′/ν∥g(∆t1/ν∥, td/z/N) ,
(102)
c(x, t)
∼
t−(β+β′)/ν∥F(x/t1/z, ∆t1/ν∥, td/z/N) .
(103)
Derived scaling properties
The scaling behavior of various other quantities can be derived directly from the scaling
relations (101)-(103). For example, the mean cluster mass M is given by the total integral
of the pair connectedness function
M =
Z
ddx
Z ∞
0
dt c(x, t) .
(104)
Inserting the scaling relation (97) and substituting the scaling variables we obtain a scaling
law for the average cluster mass measured in an inﬁnite system below criticality:
M ∼
Z
ddx
Z ∞
0
dt tθ−d/z F(x/t1/z, ∆t1/ν∥) ∼|∆|−ν∥(1+θ) .
(105)
Similarly, the mean survival time T, the mean spatial volume V , and the mean size S of a
cluster in the inactive phase are given by
T =
Z
dtP(t) =
Z
dt t−δ G(∆t1/ν∥) ∼|∆|−ν∥(1−δ) ,
(106)
V =
Z
dt P(t)td/z−1 =
Z
dt td/z−δ−1 G(∆t1/ν∥) ∼|∆|−ν∥(d/z−δ) ,
(107)
S =
Z
dt P(t)td/z =
Z
dt td/z−δ G(∆t1/ν∥) ∼|∆|−ν∥(d/z+1−δ) .
(108)
For these quantities, we obtain the following scaling relations:
M ∼|∆|−γ ,
γ = ν∥(1 + θ) = ν∥+ dν⊥−β −β′ ,
(109)
T ∼|∆|−τ ,
τ = ν∥(1 −δ) = ν∥−β′ ,
(110)
V ∼|∆|−v ,
v = ν∥(d/z −δ) = dν⊥−β′ ,
(111)
S ∼|∆|−σ ,
σ = ν∥(d/z + 1 −δ) = ν∥+ dν⊥−β′ .
(112)
47

<!-- page 48 -->

Spreading processes in an external ﬁeld
Let us ﬁnally consider a spreading process in an external ﬁeld h.
Using the particle
interpretation, such a ﬁeld may be realized by spontaneous creation of particles Ø →A
at rate h during the temporal evolution. Clearly, spontaneous particle creation destroys
the absorbing state and therefore the transition itself. That is, the external ﬁeld drives
the system away from criticality. For small h the resulting distance from criticality obeys
certain scaling laws.
In principle the presence of an external ﬁeld requires to introduce another independent
critical exponent for the coupling constant. In the case of DP, however, this exponent is not
independent, it is rather identical with the mean cluster size exponent γ. To understand
this relation, let us consider the stationary state of a subcritical DP process in presence of
a weak ﬁeld. Obviously, a site can only become active if it is connected with at least one
other site backwards in time where a particle was spontaneously created by the external
ﬁeld. Since the number of such sites is equal to the cluster size, the probability to become
active is given by ρstat ∼1 −(1 −h)M(∆).
For weak ﬁelds, the stationary density is
therefore linear in h:
ρstat ∼h M(∆) .
(∆< 0)
(113)
Consequently, the susceptibility of a supercritical DP process scales as
χ = ∂
∂hρstat ∼|∆|−γ .
(114)
Invariance under rescaling (92) requires the external ﬁeld to change as
h →Λβ′/ν⊥−z−dh = Λ−σ/ν⊥h .
(115)
Thus, at criticality the stationary response of a DP process is given by
ρstat ∼hβ/σ .
(116)
More generally, we may extend the scaling forms (101)-(103) by including the scale-
invariant argument htσ/ν∥. For example, the density ρ(t) evolves as
ρ(t) ∼t−β/ν∥f(∆t1/ν∥, td/z/N, htσ/ν∥) .
(117)
Time reversal symmetry of directed percolation
As shown in Ref. [127], a special symmetry of DP under time reversal implies the additional
scaling relation
β = β′ .
(118)
This is the reason why DP is characterized by only three instead of four critical exponents.
In order to understand this duality symmetry from the physical point of view, let us
consider the special case of directed bond percolation. By reversing the arrows shown
in the right part of Fig. 10 one obtains a directed bond percolation process that evolves
‘backwards’ in time. Obviously, the reversed process follows exactly the same probabilistic
rules as the original one. Moreover, if two sites were connected by a directed path in the
original process, they will also be connected in the reversed process. Hence, if the reversed
48

<!-- page 49 -->

process was started from a fully occupied lattice at t > 0, the resulting active sites at t = 0
would be precisely those sites which – in the original process – would generate clusters
that are still alive at time t. Therefore, in the case of directed bond percolation we obtain
P(t) = ρ(t) ,
(119)
i.e., the survival probability of a single seed P(t) is exactly equal to the density of ac-
tive sites ρ(t) in a DP process starting with a fully occupied lattice. Thus, in the active
phase the two quantities saturate at the same value P∞= ρstat wherefore the correspond-
ing critical exponents β and β′ have to be identical. It should be emphasized that this
time reversal symmetry of DP is nontrivial and does not hold for other systems such as
models with several absorbing states [129] or spreading processes with a ﬂuctuating back-
ground [130]. Together with Eq. (118) the generalized hyperscaling relation (100) reduces
to the DP hyperscaling relation
θ = d/z −2δ = dν⊥−2β
ν∥
.
(120)
Consequently, the autocorrelation function c(0, t) for DP saturates in the active phase at
the value c(0, ∞) ∼∆2β.
The DP conjecture
One of the most fascinating properties of DP models is their robustness with respect to the
microscopic dynamic rules. In fact, the DP class covers a wide range of models. It includes,
for example, the vast majority of spreading models such as the contact process [5, 131],
epidemic spreading without immunization [132], and forest ﬁre models [102, 133, 134].
Moreover, the DP class includes models for catalytic reactions [121, 135, 136], interact-
ing particles [137], as well as branching-annihilating random walks with odd number of
oﬀspring [138–140]. Furthermore, certain growth processes [141, 142] and coupled map
lattices with asynchronous updates [143] display DP behavior. In fact, this list is far from
being complete.
The variety and robustness of DP models led Janssen and Grassberger to the conjec-
ture that a model should belong to the DP universality class if the following conditions
hold [144,145]:
1. The model displays a continuous phase transition from a ﬂuctuating active phase
into a unique absorbing state.
2. The transition is characterized by a positive one-component order parameter.
3. The dynamic rules involve only short-range processes.
4. The system has no special attributes such as additional symmetries or quenched
randomness.
Although this conjecture has not yet been proven rigorously, it is highly supported by
numerical evidence. In fact, DP seems to be even more general and may be identiﬁed in
systems that violate some of the four conditions, for example in certain models with non-
unique [146–149] or ﬂuctuating passive states [129]. Even complicated spreading processes
with several spreading agents and multicomponent order parameters were shown to exhibit
49

<!-- page 50 -->

DP behavior [122,124,146,150–153]. Some models with inﬁnitely many absorbing states,
that were initially thought to belong to diﬀerent universality classes, were later found to be
in the DP class as well (see Sec. 3.8). Not only the bulk exponents β, ν⊥, ν∥are universal
but also other quantities as, for example, scaling functions and higher moments of the
order parameter [120].
It is remarkable that DP is one of very few critical phenomena in 1+1 dimensions
which has not yet been solved exactly. Despite of its simplicity and robustness it seems to
be impossible to compute the critical exponents exactly. In fact, the numerical estimates
suggest that the critical exponents are given by irrational numbers rather than simple
rational values. The lack of analytical results may be related to the fact that DP – in
contrast to ordinary (isotropic) percolation – is not conformally invariant since there is
no symmetry between ‘space’ and ‘time’. Attempts to replace conformal invariance by an
anisotropic scaling theory have not yet been successfully applied to DP [154].
Only few exceptions of DP are known so far.
They all violate at least one of the
four conditions listed above.
For example, a diﬀerent universality class emerges when
the system has two or more symmetric absorbing states (see Sec. 4.2). Another example
is the activated random walk model with a conserved order parameter (see Sec. 4.3)
Diﬀerent universal properties are also encountered in models where activity spreads over
long distances by L´evy ﬂights (see Sec. 4.1).
3.4
Estimation of the critical exponents
Directed percolation is one of very few critical phenomena whose critical exponents in 1+1
dimensions are not known exactly. However, thanks to extensive numerical simulations,
transfer matrix techniques, series expansions, and ﬁeld-theoretical calculations the critical
exponents have been estimated in various dimensions to an extremely high accuracy. This
subsection brieﬂy summarizes the available methods and the most precise estimates.
Mean-ﬁeld approximation
In order to estimate the critical exponents β, ν⊥, and ν∥of directed percolation, let us
ﬁrst consider a simple mean-ﬁeld (MF) approximation. Denoting by ρ(t) the density of
active sites at time t averaged over the entire system, the MF rate equation for the contact
process (85) reads
∂tρ(t) = (λ −1)ρ(t) −λρ2(t) .
(121)
This equation has two stationary solutions, namely ρstat = 0 and ρstat = (λ−1)/λ. Hence
the mean-ﬁeld critical point is λc = 1. The solution ρstat = 0 represents the absorbing
state from where the system cannot escape. In the inactive phase λ < λc the absorbing
state is stable while the other solution with negative density is unphysical. In the active
phase λ > λc the absorbing state becomes unstable against small perturbations while the
second solution represents a stable active state. Near criticality the stationary density
vanishes linearly as ρstat ∼λ −λc. Therefore, the mean ﬁeld density exponent is given
by βMF = 1. On the other hand, Eq. (121) implies the density to decay in the inactive
phase asymptotically as ρ(t) ∼exp(−|λ −λc|t) ∼exp(−t/ξ∥), hence ξ∥∼|λ −λc|−1 and
νMF
∥
= 1.
50

<!-- page 51 -->

In order to determine the spatial scaling exponent νMF
⊥
, the mean-ﬁeld rate equa-
tion (121) has to be extended by a term for particle diﬀusion
∂tρ(x, t) = ∆ρ(x, t) −λρ2(x, t) + D∇2ρ(x, t) ,
(122)
where D is the diﬀusion constant and ∆= λ −λc is the deviation from criticality. In a
lattice model this term corresponds to nearest-neighbor interactions. According to Eq. (94)
the density ρ(x, t) changes under rescaling (92) as
ρ(x, t) →Λ−β/ν⊥ρ(Λx, Λzt) .
(123)
Simple dimensional analysis shows that Eq. (122) is invariant under rescaling if
βMF = 1 ,
νMF
⊥
= 1
2 ,
νMF
∥
= 1 .
(124)
The above mean ﬁeld approximation becomes exact in the limit of inﬁnitely many di-
mensions. In a ﬁnite-dimensional contact process it is not clear whether the mean ﬁeld
approximation still applies, it is even not obvious that a continuous phase transition still
exists. However, Liggett [5] was able to rigorously prove the existence of a phase transi-
tion for a contact process in d ≥1 dimensions. As will be shown below, the mean ﬁeld
exponents turn out to be exact for d ≥4, where dc = 4 is the upper critical dimension of
DP. Note that the mean ﬁeld exponents satisfy the hyperscaling relation (120) precisely
in d = 4 dimensions.
In order to go beyond the standard mean ﬁeld approximation in low-dimensional sys-
tems spatial correlations have to be taken into account.
An improved mean ﬁeld ap-
proximation for the contact process in 1+1 dimensions was developed by Ben-Naim and
Krapivsky [155], who expressed the temporal evolution of empty intervals on the lattice by
an inﬁnite hierarchy of diﬀerential equations. This approach is very similar to the IPDF
technique introduced in Sec. 2.8 Approximating the probability to ﬁnd pairs of neighbor-
ing empty intervals by the product of single-interval probabilities, they derived a set of
equations which can be solved exactly. In this approximation the critical exponents are
given by β = 1
2, ν⊥= 1, and ν∥= z = 3
2. ´Odor could improve these estimates by using a
generalized mean ﬁeld approximation combined with coherent anomaly techniques [156],
reaching an accuracy of almost 1%.
Monte Carlo simulations with homogeneous initial conditions
In order to study nonequilibrium phase transitions quantitatively, numerical techniques
such as Monte Carlo simulations have become an important tool. Because of the steadily
growing computer capacity the critical exponents can nowadays be estimated within a few
per cent, in some cases even up to four digits. Further progress is expected as reaction-
diﬀusion models can easily be simulated on parallel computers with a large number of
simple processors [45].
The simplest numerical method that allows the critical exponents to be estimated is
a Monte Carlo (MC) simulation starting with a fully occupied lattice. This technique is
based on the scaling properties of Eq. (101). In a ﬁrst series of simulations the critical
percolation threshold has to be determined by measuring deviations from the asymptotic
power-law decay ρ(t) ∼t−δ in a suﬃciently large system.
To this end ρ(t) is plotted
versus t in a double-logarithmic graph (see Fig. 20a). Positive (negative) curvature for
51

<!-- page 52 -->

10
0
10
1
10
2
10
3
10
4
t
10
−1
10
0
ρ(t)
10
−6 10
−4 10
−2
10
0
10
2
t (p−pc)
ν||
10
−6
10
−4
10
−2
10
0
ρ(t)t
δ
10
−6 10
−4 10
−2
10
0
10
2
t/N
z
10
−4
10
−3
10
−2
10
−1
10
0
ρ(t)t
δ
p<pc
p=pc
p>pc
(a)
(b)
(c)
Fig. 20:
Ordinary Monte Carlo simulations of a (1+1)-dimensional directed bond percolation process
starting from a fully occupied lattice.
Part (a) shows the particle density as a function of time and
demonstrates the determination of the critical point. Parts (b) and (c) show data collapses for oﬀ-critical
and ﬁnite-size simulations, respectively (see text).
large t indicates that the system is still in the active (inactive) phase. It should be carefully
analyzed to what extent the estimate depends on the system size used in the simulation. If
ﬁnite-size eﬀects play a role, extrapolation techniques should be used in order to improve
the estimate [157].
Having determined pc and δ the exponent β may be estimated by measuring the
stationary density of active sites ρstat ∼∆β in the active phase. However, this type of
estimate is known to be quite inaccurate since the equilibration time to reach the stationary
state grows rapidly as the critical point is approached. This critical slowing down can be
controlled by plotting ρ(t) tδ versus t ∆ν∥for diﬀerent values of ∆and tuning ν∥in a way
that all curves collapse (see Fig. 20b). The exponent β is then given by β = δν∥. Finally,
the exponent ν⊥can be determined by ﬁnite-size simulations. According to Eq. (101),
ρ(t) tδ has to be plotted against t/N z/d for various system sizes (see Fig. 20c). By tuning z
the data points collapse onto a single curve which gives an estimate for ν⊥= ν∥/z.
Monte Carlo simulations with localized initial conditions
More accurate estimates for the critical exponents can be obtained by dynamic simulations
starting from a single particle (active seed) [127].
This technique exploits the scaling
properties of the pair-connectedness function c(x, t). Starting from a single particle, one
measures the survival probability P(t), the number of active sites N(t), and the mean
square of spreading from the origin R2(t) averaged over surviving runs. According to
Eqs. (102) and (103) these quantities obey the scaling forms
P(t) ≃t−δ g(∆t1/ν∥, td/z/N) ,
(125)
N(t) =
Z
ddx c(x, t) ≃tθ ¯F(∆t1/ν∥, td/z/N) ,
(126)
R2(t) = ⟨x2(t)⟩=
1
N(t)
Z
ddx x2 c(x, t) ≃t2/z ˜F(∆t1/ν∥, td/z/N) .
(127)
At criticality, they are expected to display asymptotic power laws
P(t) ∼t−δ ,
N(t) ∼tθ ,
R2(t) ∼t2/z ,
(128)
52

<!-- page 53 -->

0
0.01
0.02
1/t
−0.19
−0.17
−0.15
−0.13
0
0.01
0.02
1/t
0.3
0.31
0.32
0.33
0.34
0.35
0
0.01
0.02
1/t
1.25
1.3
1.35
1.4
δ(t)
θ(t)
2/z(t)
Fig. 21:
Dynamic Monte Carlo simulations of a directed bond percolation process. The eﬀective exponents
δ(t), θ(t), and 2/z(t) are shown as a function of time for p = 0.6446, 0.6447, 0.6448. The arrows mark the
actual values of the critical exponents.
i.e., they show straight lines in double logarithmic plots.
Oﬀcriticality, the lines are
curved, allowing a precise determination of the percolation threshold pc. Technically it is
often useful to consider local slopes of these curves by introducing eﬀective exponents
−δ(t) = log10
 P(t)/P(t/b)

log10 b
(129)
and similarly θ(t) and 2/z(t), where log10 b is the distance used for estimating the slope.
Plotting the local slopes as functions of 1/t, the curves may be extrapolated to t →∞, as
illustrated in Fig. 21. The same method works also in higher dimensional systems [158].
In order to improve the estimates, it is useful to eliminate the curvature of the data points
at criticality by plotting the quantities (128) against 1/tδ instead of 1/t.
Since the spatial size of the growing cluster at a given time is ﬁnite, the simulation
can be accelerated considerably by storing the coordinates of active particles in a dynam-
ically generated list. Especially at criticality, where the density of particles is low, such
algorithms are much more eﬃcient. Moreover, ﬁnite-size eﬀects are eliminated completely.
Numerical diagonalization
The critical exponents may also be approximated by numerical diagonalization of the
evolution operator. Although the resulting estimates are usually less accurate than those
obtained by other methods, this technique is of conceptual interest. Let us, for example,
consider the (1+1)-dimensional contact process on a ﬁnite lattice with N sites and periodic
boundary conditions which is deﬁned by the Liouville operator (87). Solving the eigenvalue
problem
LCP (λ) |ψi⟩= µi |ψi⟩
(130)
we obtain a spectrum of eigenvalues {µi}, as shown in the left panel of Fig. 22. As in all
reaction-diﬀusion models, the lowest eigenvalue µ0 vanishes. The corresponding stationary
state |ψ0⟩is the absorbing state of the contact process. The other eigenvectors represent
the relaxational modes of the system.
As can be seen in Fig. 22, all of them have a
short lifetime except for the ﬁrst excited state |ψ1⟩whose eigenvalue µ1 tends to zero as
λ increases. This eigenvector represents the active state of the system. In ﬁnite systems
53

<!-- page 54 -->

0
1
2
3
4
5
λ
0
1
2
3
Re[µi]
−10
−5
0
5
10
(λ−λc)N
1/ν
0
1
2
3
4
5
µ1 N
−z
 N=4
 N=6
 N=8
λc
µ0
µ1
T
Fig. 22:
Numerical diagonalization of the Liouville operator of the contact process. Left panel: Low-lying
part of the relaxational spectrum of the (1+1)-dimensional contact process on a lattice with 8 sites and
periodic boundary conditions. Right panel: Data collapse according to the scaling form (131).
there is always a ﬁnite probability to reach the absorbing state, hence µ1 > 0. In inﬁnite
systems, however, this eigenvalue decreases with λ and vanishes at the critical point. Since
the amplitude of |ψ1⟩decays in the inactive phase as e−µ1t, we may identify µ−1
1
with the
temporal scaling length ξ∥. In an inﬁnite system we therefore expect µ1 to decrease as
µ1 ∼|λ −λc|ν∥for λ < λc and to vanish for λ > λc. The corresponding scaling form reads
µ1 ∼N −z/d h(∆N 1/dν⊥) ,
(131)
where ∆= λ −λc. Thus, by plotting µ1N z/d against ∆N 1/dν⊥, the exponents z and ν⊥
can be determined by data collapse, as demonstrated in the right panel of Fig. 22. In
order to determine the exponent β, it would be necessary to analyze the components of
the eigenvector |ψ1⟩with respect to the particle density in the active phase. A similar
analysis of DP models with parallel updates, which are deﬁned by transfer matrices instead
of Liouville operators, can be found in Ref. [159].
Density matrix renormalization group methods
The method of numerical diagonalization can be improved considerably by using density
matrix renormalization group (DMRG) techniques. The concept of DMRG was introduced
in 1992 by White [160] in the context of equilibrium statistical physics as a tool for the
diagonalization of quantum spin chains. The main idea is to prolongate a given spin chain
by inserting additional spins and to reduce the resulting conﬁguration space by a suitable
projection mechanism, keeping only the most relevant eigenstates. This renormalization
procedure is then repeated many times and the spectrum of the iterated Hamiltonian is
analyzed. Recently DMRG techniques have also been applied to various (1+1)-dimensional
nonequilibrium systems [161–163] (see [100] for a general overview). The method yields
surprisingly accurate results. For example, Carlon et al. [162] were able to estimate the
critical exponents of DP by β/ν⊥= 0.249(3), ν⊥= 1.08(2), and z = 1.580(1), deviating
from the currently accepted values by less than 1.5%.
Series expansions
54

<!-- page 55 -->

The most precise estimates of the DP exponents in 1+1 dimensions have been obtained
by series expansions [119].
This technique is very similar to low- or high-temperature
expansions in equilibrium statistical physics. As an example let us consider the (1+1)-
dimensional contact process. Its Liouville operator (87) may be separated into two parts
L(λ) = L0 + λL1, where L0 and L1 describe spontaneous self-destruction A →Ø and
oﬀspring production A →2A, respectively.
The basic idea is to regard λ as a small
perturbation and to express physical quantities as power series in λ. To this end it is useful
to introduce the Laplace transform of the probability distribution |Pt⟩and to expand it
in powers of λ:
| ˜P(s)⟩=
Z ∞
0
dt e−st |Pt⟩=
∞
X
n=0
λn | ˜Pn(s)⟩.
(132)
By applying L(λ) from the left one can easily derive the recursion relation
(s −L0) | ˜Pn(s)⟩=

|P0⟩
if n = 0
L1| ˜Pn−1(s)⟩
if n ≥1
,
(133)
where |P0⟩denotes the initial particle conﬁguration. Hence, if the process started from a
conﬁguration with a single particle, the vector | ˜Pn(s)⟩describes an ensemble of conﬁgura-
tions with at most n particles. It is therefore possible to explicitly construct the vectors
| ˜Pn(s)⟩, as described in detail in Ref. [119].
The above expansion allows the temporal integral of any observable X(t) = ⟨1|X|Pt⟩
to be expressed as a power series in λ (for notations see Appendix A):
Z ∞
0
dt X(t) = lim
s→0 ⟨1|X| ˜P(s)⟩=
∞
X
n=0
λn lim
s→0 ⟨1|X| ˜Pn(s)⟩.
(134)
Let us, for example, consider the survival probability P(t) that the system has not yet
reached the absorbing state at time t [cf.
Eq. (94)].
Using the vector formalism this
quantity may be written as
P(t) = 1 −⟨0|Pt⟩= ⟨1|Pt⟩−⟨0|Pt⟩,
(135)
where ⟨0| = (1, 0, 0, . . . , 0) denotes the absorbing state. The critical exponents can be
estimated as follows.
On the one hand, the mean survival time T of clusters in the
inactive phase can be expanded in powers of λ by
T =
Z ∞
0
dt P(t) =
∞
X
n=0
λn lim
s→0

⟨1| ˜Pn(s)⟩−⟨0| ˜Pn(s)⟩

.
(136)
On the other hand, according to Eq. (110) we have T ∼(−∆)β′−ν∥so that
d
dλ ln T ≃ν∥−β′
λc −λ + const .
(137)
Therefore, in order to estimate λc and β −ν∥, three steps have to be taken. At ﬁrst, the
vectors | ˜Pn(s)⟩have to be determined by iterating Eq. (133) up to order nmax. Although
this recursion relation is quite complicated, it is still simple enough to be implemented on a
computer (for example, in Ref. [119] the iteration was carried out up to order nmax = 24).
55

<!-- page 56 -->

critical
MF
IMF
d = 1
d = 2
d = 3
d = 4 −ǫ
exponent
[155]
[168]
[125]
[170]
[171]
β
1
1/2
0.276486(8)
0.584(4)
0.81(1)
1 −ǫ/6 −0.01128 ǫ2
ν⊥
1/2
1
1.096854(4)
0.734(4)
0.581(5)
1/2 + ǫ/16 + 0.02110 ǫ2
ν∥
1
3/2
1.733847(6)
1.295(6)
1.105(5)
1 + ǫ/12 + 0.02238 ǫ2
z
2
3/2
1.580745(10)
1.76(3)
1.90(1)
2 −ǫ/12 −0.02921 ǫ2
δ
1
1/2
0.159464(6)
0.451
0.73
1 −ǫ/4 −0.01283 ǫ2
θ
0
1/2
0.313686(8)
0.230
0.12
ǫ/12 + 0.03751 ǫ2
γ
1
3/2
2.277730(5)
1.60
1.25
1 + ǫ/6 + 0.06683 ǫ2
v
1
1/2
0.82037(1)
0.88
0.94
1 −ǫ/12 + 0.03317 ǫ2
σ
2
2
2.554216(13)
2.18
2.04
2 + ǫ2/18
Table 2:
Estimates for the critical exponents of directed percolation obtained by mean ﬁeld (MF),
improved mean ﬁeld (IMF), numerical, as well as ﬁeld-theoretical methods.
Next, one has to express T as a power series in λ. Finally, λc and β −ν∥can be estimated
by determining the location and the amplitude of the singularity in Eq. (137) and by using
a Pad´e approximation [3]. Since the singularity is approached from the inactive phase, we
are dealing with a subcritical expansion. Similarly one may also consider the supercritical
case by expanding L(µ) = µL0 + L1 in powers of µ.
A general review on series expansion can be found in Ref. [3]. Series expansions were
applied to (1+1)-dimensional DP ﬁrstly in Ref. [164], where the critical exponents could
be determined with a relative accuracy of about 10−3. Reﬁned simulations [165] led the
authors to the conjecture that the DP exponents should be given by the rational values
β = 199/720, ν⊥= 26/15, and ν∥= 79/72. In a sequence of papers [116,119,166,167] the
error margins could be further reduced down to 10−4 . . . 10−5. These improved estimates
showed that the conjectured rational values were incorrect, indicating that the critical
exponents of DP could be given by irrational numbers. This should be taken as a warning
that critical exponents of non-integrable systems are usually not given by simple rational
values. Currently, the most precise estimates are given in Ref. [168]. Series expansions for
DP were also performed in two spatial dimensions [169]. In addition, the exponents were
found to be independent of the type of lattice under consideration. For easy reference we
listed the most precise estimates in Table 2.
Field-theoretical approximations
By a ﬁeld-theoretical renormalization group calculation (see Sec. 3.5) it is possible to com-
pute ﬂuctuation corrections of the critical exponents close to the upper critical dimension
dc = 4 in powers of ǫ = dc −d. In a two-loop approximation [144, 171] these corrections
are given by
β
=
1 −ǫ/6 +
11
12 −53
6 ln 4
3

(ǫ/12)2 + 0(ǫ3) ,
ν⊥
=
1
2 + ǫ/16 +
107
32 −17
16 ln 4
3

(ǫ/12)2 + 0(ǫ3) ,
(138)
ν∥
=
1 + ǫ/12 +
109
24 −55
12 ln 4
3

(ǫ/12)2 + 0(ǫ3) .
For d ≤2 these approximations are quite inaccurate. However, in three spatial dimensions,
where numerical simulations and series expansions are diﬃcult to perform, the two-loop
approximations are regarded as the most precise estimates available.
56

<!-- page 57 -->

3.5
Field-theoretic formulation of directed percolation
The robustness of the DP universality class can be partly understood by studying the
corresponding ﬁeld theory. It is interesting to note that the ﬁeld theory of DP was ﬁrst
discovered in a quite diﬀerent ﬁeld of physics, namely in the context of hadronic interac-
tions at ultra-relativistic energies. In order to predict the cross sections of such particles
at high energies quantitatively, a ﬁeld-theoretic approach, called Reggeon ﬁeld theory, was
developed in the 70’s (see Refs. [172–177], a general review is given in [178]). Surprisingly
it took almost another ten years to realize that Reggeon ﬁeld theory was nothing but a
ﬁeld-theoretic realization of the contact process [179–181], sometimes also called Gribov
process [182,183]. In the following we sketch the main ideas of a ﬁeld-theoretic approach
to DP.
The DP Langevin equation
The Langevin equation of motion for directed percolation can be derived directly from the
master equation for the contact process [144] and reads
∂tρ(x, t) = κρ(x, t) −λρ2(x, t) + D∇2ρ(x, t) + ζ(x, t) .
(139)
It diﬀers from the mean ﬁeld equation (122) by a density-dependent Gaussian noise ﬁeld
ζ(x, t), which is deﬁned by its correlations
⟨ζ(x, t)⟩= 0 ,
⟨ζ(x, t)ζ(x′, t′)⟩= Γ ρ(x, t) δd(x −x′) δ(t −t′) .
(140)
Since the amplitude of ζ(x, t) is proportional to
p
ρ(x, t), the noise is said to be multiplica-
tive. This ensures that the absorbing state ρ(x, t) = 0 does not ﬂuctuate. The square-root
behavior stems from the deﬁnition of ρ(x, t) as a coarse-grained density of active sites
averaged over some mesoscopic box size. Only active sites in this box give rise to ﬂuctua-
tions of the density, generating a bounded uncorrelated noise. The noise ﬁeld ζ(x, t) can
be viewed as the sum of all these noise contributions in the box. According to the central
limit theorem, if the number of particles in the box is suﬃciently high, ζ(x, t) tends to a
Gaussian distribution with an amplitude proportional to the square root of the number
of active sites in the box. This type of noise has to be distinguished from other nonequi-
librium systems with multiplicative noise where the noise amplitude is proportional to
the ﬁeld ρ(x, t) itself without square root [184,185]. These systems do not belong to the
DP class, rather they are related to the KPZ universality class [186]. The DP Langevin
equation was also tested numerically in Ref. [187], conﬁrming that the critical exponents
are in agreement with those of ordinary DP lattice models. Note that in contrast to the
annihilation process discussed in Sec. 2.6, the noise (140) is real due to positive density
correlations in the bulk.
The Langevin equation (139) can be seen as a minimal equation needed to describe
DP. It may also include higher order terms such as ρ3(x, t) or ∇4ρ(x, t), but these contri-
butions turn out to be irrelevant under renormalization group transformations. The same
applies to higher-order contributions to the noise.
These additional terms account for
(nonuniversal) short-range correlations while they are irrelevant on large scales. In fact,
the robustness of DP originates in the irrelevance of higher-order terms in the Langevin
equation.
57

<!-- page 58 -->

Relation to Reggeon ﬁeld theory
In ﬁeld-theoretic calculations it is often more convenient to characterize the dynamic
system by a partition sum Z. The sum is carried out over all realizations of the ﬁeld
φ(x, t) = ρ(x, t) and the noise ζ(x, t), weighted by an appropriate eﬀective action. More
precisely, the partition sum is deﬁned as the integral over all realizations of the ﬁeld φ(x, t)
and the noise ζ(x, t) which satisfy the Langevin equation. Therefore, we may write the
integrand as a δ-function with Eq. (139) as its argument:
Z ∼
Z
Dζ P[ζ]
Z
Dφ I[φ] δ

∂tφ −D∇2φ −κφ + λφ2 −ζ

.
(141)
Here Dζ and Dφ denote functional integration, P[ζ] is the probability distribution of the
noise ﬁeld, and I[φ] stands for an appropriate Jacobian which turns out to be irrelevant
in the present problem. As shown in Appendix B, it is possible to integrate the noise by
introducing a Martin-Siggia-Rosen response ﬁeld ˜φ(x, t). The resulting action S = S0+Sint
with
S0[φ, ˜φ]
=
Z
ddx dt ˜φ(x, t)

τ∂t −D∇2 −κ

φ(x, t) ,
(142)
Sint[φ, ˜φ]
=
Γ
2
Z
ddx dt ˜φ(x, t)

φ(x, t) −˜φ(x, t)

φ(x, t) ,
(143)
is the eﬀective action of Reggeon ﬁeld theory [178]. In momentum space it may also be
written as
S0[φ, ˜φ]
=
Z
dkω ˜φ(−k, −ω) (−iτω + Dk2 −κ)φ(k, ω)
(144)
Sint[φ, ˜φ]
=
Γ
2
Z
dkω
Z
dk′ω′ φ(−k, −ω) ˜φ(−k′, −ω′)
(145)
×
h
φ(k + k′, ω + ω′) −˜φ(k + k′, ω + ω′)
i
,
where dkω = (2π)−d−1ddk dω. Formally the free part of the action can be expressed as
S0[φ, ˜φ] = 1
2
Z
dkω⃗Φ(−k, −ω) S0(k, ω)⃗Φ(k, ω) ,
(146)
S0(k, ω) =

0
Dk2 −κ + iτω
Dk2 −κ −iτω
0

,
where⃗Φ = (φ, ˜φ). Introducing external currents J(k, ω) and ˜J(k, ω) we can rewrite the
partition sum as
Z[J, ˜J] ∼
Z
Dφ D ˜φ I′[φ, ˜φ] exp

−S0[φ, ˜φ] −Sint[φ, ˜φ] +
Z
ddx dt (Jφ + ˜J ˜φ)

∼exp

−Lint
h δ
δJ , δ
δ ˜J
i
exp
1
2
Z
dkω⃗J(−k, −ω) G0(k, ω)⃗J(k, ω)

,
(147)
where the free propagator G0 = S−1
0
is given by
G0(k, ω) =

0
G0(k, ω)
G0(−k, −ω)
0

(148)
58

<!-- page 59 -->

(x,t)
time
(x’,t’)
Fig. 23:
Critical DP cluster (grey) and the backbone of the pair connectedness function c(x′, t′, x, t),
illustrating the physical meaning of Feynman diagrams in directed percolation.
with G0(k, ω) = (Dk2 −κ −iτω)−1. Because of G0(q, ω)̸ = G0(q, −ω) DP is an irreversible
process.
Cluster backbone and Feynman diagrams
Before turning to ﬁeld-theoretic renormalization group techniques let discuss the physical
meaning of Feynman diagrams in directed percolation. The full propagator of the ﬁeld
theory is the pair connectedness function c(x′, t′, x, t) which is deﬁned as the probability
that two sites (x′, t′) and (x, t) are connected by a directed path of open bonds (see
Sec. 3.3). In a given realization of open and closed bonds there may be several possible
directed paths connecting the two sites, as illustrated in Fig. 23. The union of all possible
paths constitutes the so-called backbone of the pair connectedness function [188]. More
precisely, the backbone consist of all sites that are connected with the sites (x′, t′) and (x, t)
by a directed walk, i.e., we cut oﬀall dangling ends of the cluster. From the topological
point of view the backbone is a directed graph consisting of branching and merging lines.
Because of the duality symmetry (119), it is statistically invariant under time reversal.
In principle the full propagator is given by a weighted sum over all possible backbone
conﬁgurations. Fluctuation eﬀects are mainly due to the inﬂuence of closed loops. Above
the upper critical dimension dc = 4 the degree of spatial freedom for propagating lines
is so high that the probability to merge tends to zero. Thus, the contribution of closed
loops can be neglected and hence the full propagator is eﬀectively described by the free
propagator. Below the critical dimension, however, loops occur more frequently and begin
to play a signiﬁcant role (see Fig. 23).
The loops of the backbone may be associated with the Feynman diagrams of Reggeon
ﬁeld theory.
As shown in Fig. 24(a)-(c), the backbone can be decomposed into three
elementary components. The arrow stands for the free propagator (148) while the diagrams
for branching and merging represent the cubic vertices in Eq. (143), associated with the
weights ±Γ/2. Because of self-destruction A →Ø, free paths have a ﬁnite lifetime, as
expressed by the bare mass κ of the free propagator. Consequently, the paths in a given
conﬁguration of the backbone have to be weighted by their length. Moreover, each closed
loop carries a weight −Γ2/4.
The negative sign for the weight of closed loops can be explained as follows. The pair
connectedness function is the sum over all backbone conﬁgurations b connecting the sites
59

<!-- page 60 -->

(f)
(d)
(b)
(c)
(a)
(e)
Fig. 24:
Feynman diagrams of directed percolation.
Left: Elementary components of the backbone:
(a) free propagator, (b) branching vertex, and (c) merging vertex.
Right: One-loop diagrams for (d)
propagator renormalization and (e,f) vertex renormalization.
(x′, t′) and (x, t) weighted by their probability Pb:
c(x′, t′, x, t) =
X
b
Pb .
(149)
In order to ﬁnd an recurrence relation for the probability Pb, let us consider directed bond
percolation on a lattice. Obviously, Pb is the weighted sum over all lattice conﬁgurations
compatible with the backbone conﬁguration b. The backbone itself contributes with a
factor pnb, where p is the percolation probability and nb denotes the number of bonds
occupied by the backbone b.
Another factor comes from the bonds outside the back-
bone. This factor can be expressed as the probability that b is not contained in a larger
backbone b′. Thus, Pb satisﬁes the recurrence relation
Pb = pnb(1 −p−nb X
b′,b⊂b′
Pb′) .
(150)
Since b′ contains always more loops than b, this relation can be used to expand the pair
connectedness function (149) in the number of loops. As can be easily veriﬁed, the one-loop
correction carries a negative sign. More generally, it is possible to show by an inclusion-
exclusion argument that each closed loop contributes with a negative weight. An analogous
proof for isotropic percolation is explained in detail in the review article by Essam [108].
Thus, apart from the negative weight of closed loops, the backbone may be interpreted
as a graph consisting of Feynman diagrams. This interpretation is possible because in a
DP process the ﬁeld φ(x, t) represents the local density of particles. It should be noted
that this is not always true. Moreover, various authors prefer to shift the ﬁeld φ by its
mean ﬁeld expectation value so that the interpretation as a density is no longer valid. For
this reason we will continue to use the unshifted ﬁelds φ and ˜φ.
One-loop approximation
Slightly below the upper critical dimension one-loop diagrams start to contribute to the full
propagator while higher-order diagrams are still strongly suppressed. In this regime the
DP process can be approximated by neglecting higher-order loop diagrams. The resulting
propagator consists of a sum over n concatenated one-loop diagrams with n running from
zero to inﬁnity.
In momentum space the corresponding expression can be written as
a simple geometric series.
By carrying out the integration one obtains an ultraviolet-
divergent expression.
Hence, in order to regularize the propagator, an upper cutoﬀΩ
has to be introduced in momentum space. Physically this upper cutoﬀcorresponds to the
lattice spacing of the DP model. In other words, DP needs a lattice; there is no continuum
theory of DP.
In order to approximate the critical exponents, we use Wilson’s renormalization group
scheme [189] which consists of two steps (see Fig. 25). At ﬁrst the theory is coarse-grained
60

<!-- page 61 -->

Ω
step 2
step 1
Ω
Ω/Λ
Fig. 25:
Wilson’s renormalization group. Step 1: Scaling transformation in momentum space. Step 2:
Integration in a momentum shell (shaded region).
by a scaling transformation x →Λx with Λ < 1, leading to a change of the coeﬃcients in
the eﬀective action and a dilatation of momentum space (including the cutoﬀΩ). In the
second step the short-range ﬂuctuations are integrated out in a momentum shell. This can
be done by evaluating the Feynman diagrams in the range Ω≤k ≤Ω/Λ and absorbing the
resulting contributions in the coeﬃcients. The total change of the coeﬃcients determines
the RG ﬂow and therefore the critical exponents.
Let us ﬁrst consider the scaling transformation. Because of the duality symmetry of
DP (see Sec. 3.3) the action S is invariant under the replacement
φ(x, t) →−˜φ(x, −t) ,
˜φ(x, t) →−φ(x, −t) ,
(151)
implying that φ and ˜φ have exactly the same scaling behavior:
x →Λ x ,
t →Λz t ,
φ(x, t) →Λχ φ(Λx, Λzt) ,
˜φ(x, t) →Λχ ˜φ(Λx, Λzt) .
(152)
Under this scaling transformation the eﬀective action (142)-(143) turns into
S0[φ, ˜φ]
=
Z
ddx dt ˜φ(x, t)

τΛ2χ+d
| {z }
τ ′
∂t −DΛ2χ+d+z−2
|
{z
}
D′
∇2 −κΛ2χ+d+z
|
{z
}
κ′

φ(x, t) ,
Sint[φ, ˜φ]
=
1
2
Z
ddx dt Γ Λ3χ+d+z
|
{z
}
Γ′
φ(x, t)˜φ(x, t)

φ(, t) −˜φ(, t)

.
(153)
Thus, for an inﬁnitesimal dilatation Λ = 1 + l, the four coeﬃcients rescale as
τ ′
=
[1 + l (2χ + d)] τ ,
D′
=
[1 + l (2χ + d + z −2)] D ,
(154)
κ′
=
[1 + l (2χ + d + z)] κ ,
Γ′
=
[1 + l (3χ + d + z)] Γ .
In the second step of Wilson’s RG procedure the one-loop diagrams are integrated in a
momentum shell. The propagator is renormalized by diagram (d) in Fig. 24
G−1
0 (k, ω)′′ = G−1
0 (k, ω) −Γ2
2
Z
>
dk′ω′ G0(k
2 + k′, ω
2 + ω′) G0(k
2 −k′, ω
2 −ω′) ,
(155)
where ‘>’ denotes integration in the momentum shell Ω≤k ≤Ω/Λ. This equation can
be rewritten as
κ′′ −D′′k2 + iτ ′′ω = κ′ −D′k2 + iτ ′ω −Γ′2
2 JP ,
(156)
61

<!-- page 62 -->

where JP denotes the integral in Eq. (155). Integrating JP and expanding the result to
the lowest order in k and ω yields the series (see Appendix D)
JP = lKdΩd
2τ

1
Ω2D −κ −
Ω2D
4(Ω2D −κ)2 k2 +
iτ
2(Ω2D −κ)2 ω + . . .

.
(157)
Therefore, the coeﬃcients in Eq. (156) are renormalized to one-loop order by
τ ′′ = τ ′ −
Γ2lKd
8(Ω2D −κ)2 ,
D′′ = D′ −
Γ2lKdΩ2D
16τ(Ω2D −κ)2 ,
κ′′ = κ′ −
Γ2lKd
4τ(Ω2D −κ) .
(158)
Finally, we have to renormalize the coupling constant Γ. Because of the duality symme-
try (119) the cubic vertices renormalize identically [see diagrams (e) and (f) in Fig. 24].
For the cubic vertices it is suﬃcient to carry out the integration at k = ω = 0:
Γ′′ = Γ′ −2Γ3
Z
>
dkω G2
0(k, ω) G0(−k, −ω) = Γ′ −
lΓ3Kd
2τ(Ω2D −κ)2 .
(159)
Adding the changes of the coeﬃcients under rescaling (154) and the subsequent shell
integration (158)-(159) we obtain the RG ﬂow equations
∂lτ
=
τ

2χ + d −
Γ2KdΩd
8τ(DΩ2 −κ)2

,
∂lD
=
D

2χ + d + z −2 −
Γ2KdΩd
16τ(DΩ2 −κ)2

,
(160)
∂lκ
=
κ

2χ + d + z −
Γ2KdΩd
4κτ(DΩ2 −κ)

,
∂lΓ
=
Γ

3χ + d + z −
Γ2KdΩd
2τ(DΩ2 −κ)2

.
Two scaling combinations appear in these equations, namely
S1 =
Γ2KdΩd
16τ(DΩ2 −κ)2 ,
S2 =
Γ2KdΩd
4κτ(DΩ2 −κ) .
(161)
Two of the four parameters τ, D, κ, Γ can be chosen freely11. Here we ﬁx the coeﬃcients
of the spatial and temporal derivatives, i.e., we require τ and D to be invariant under RG
transformations. Thus the ﬁrst two equations read
4 −ǫ + 2χ −2S1 = 0 ,
2 −ǫ + 2χ + z −S1 = 0 ,
(162)
where d = 4 −ǫ. The RG ﬂow is then described by two diﬀerential equations:
∂lκ = κ

4 −ǫ + 2χ + z −S2

= κ

2 + S1 −S2

,
∂lΓ = Γ

4 −ǫ + 3χ + z −8S1

= Γ

ǫ/2 −6S1

.
(163)
11On the level of lattice models such as the DK model or the contact process, this freedom corresponds
to choosing the time scale and a point on the phase transition line.
62

<!-- page 63 -->

Γ
κ
Fig. 26:
Linearized RG ﬂow near the ﬁxed point of directed percolation. To reach the ﬁxed point, the
system has to be on the dashed line, i.e., there is one parameter in the model which has to be tuned to
criticality.
At the ﬁxed point (S∗
1, S∗
2), κ and Γ are invariant under RG transformations, i.e., 2+S∗
1 −
S∗
2 = 0 and ǫ/2 −6S∗
1 = 0. Therefore, the ﬁxed point is located at
S∗
1 = ǫ/12 ,
S∗
2 = 2 + ǫ/12 .
(164)
Inserting this solution into Eq. (162) we obtain two of the three critical exponents, namely
χ = −2+ 7ǫ/12 and z = 2−ǫ/12. The third exponent can be determined by investigating
the RG ﬂow in the vicinity of the ﬁxed point. Because of Eq. (161), the ﬁxed point values
for κ and Γ are given by
κ∗= 4DΩ2ǫ
24 + 5ǫ = DΩ2
6
ǫ + 0(ǫ2) ,
(165)
Γ2∗=
2D(24 + ǫ)
24 + 5ǫ
r ǫτ
3Kd
2
= 4D2τ
3Kd
ǫ + 0(ǫ2) ,
(166)
where we assumed that Ωd ≃Ω4. Close to the ﬁxed point, the RG ﬂow in Eq. (163)
can be linearized. As shown in Fig. 26, the ﬂow is attractive along the dashed line and
repulsive elsewhere. To ﬁrst order in ǫ the corresponding Jacobian is triangular. Hence
its eigenvalues are given by the diagonal elements
∂κκ(2 + S1 −S2)|κ=κ∗,Γ=Γ∗= 2 −ǫ/4 + 0(ǫ2) ,
∂ΓΓ(ǫ/2 −6S1)|κ=κ∗,Γ=Γ∗= −ǫ + 0(ǫ2) .
(167)
The positive eigenvalue corresponds to the repulsive eigenvector (dotted line in Fig. 26).
Since the parameter κ plays the role of the reduced percolation probability p −pc, this
eigenvalue is equal to ν−1
⊥, rendering the third critical exponent. Because of χ = −β/ν⊥
and z = ν||/ν⊥we thus obtain the critical exponents
β = 1 −ǫ/6 + 0(ǫ2) ,
ν⊥= 1/2 + ǫ/16 + 0(ǫ2) ,
ν|| = 1 + ǫ/12 + 0(ǫ2) .
(168)
63

<!-- page 64 -->

                                        







































                                                                                                                                                                                                                                                                                                                                                                                            





































(a)
(c)
(b)
t
Fig. 27:
(a) Schematic drawing of a DP cluster. (b) The same cluster with an absorbing boundary. (c)
The same cluster in a parabola-shaped geometry.
A two-loop approximation of these exponents (see Eq. (138)) can be found in Ref. [171].
Although the two-loop result is quite accurate in 3+1 dimensions, it cannot compete with
numerical methods in lower dimensions. For example, in 1+1 dimensions the approxi-
mation for the density exponent β diﬀers from the known numerical value by more than
40%. Even one-dimensional fermionic ﬁeld theories, which have been introduced recently
in Ref. [190], turn out to be inaccurate. Therefore, regarding quantitative results, ﬁeld-
theoretic methods are only of limited interest. However, in many cases they are extremely
useful to understand essential universal properties of the system. For example, various
scaling relations can only be proven by means of ﬁeld-theoretic considerations. In fact, the
ﬁeld-theoretic renormalization group is one of the most powerful tools of nonequilibrium
statistical mechanics.
3.6
Surface critical behavior
As in equilibrium statistical mechanics, nonequilibrium critical phenomena depend cru-
cially on the boundary conditions of the system. Because of long-range correlations, the
choice of the boundary conditions may aﬀect the physical properties of the entire system.
The critical behavior at surfaces of equilibrium models has been studied extensively (for
a review see Igl´oi et al. [191]). As suggested by Cardy [192], surface critical phenomena may
be described by introducing an additional surface exponent for the order parameter ﬁeld
which is generally independent of the other bulk exponents. A similar picture emerges
in nonequilibrium statistical physics. However, since in this case there is no symmetry
between space and time, we have to distinguish between spatial, temporal and mixed
surfaces. The simplest example of a spatial surface is a semi-inﬁnite system with a wall.
Close to the wall the scaling behavior of the order parameter is characterized by a surface
critical exponent βs whose value depends on the type of boundary condition. The most
important example of a temporal surface is the initial state of a nonequilibrium system.
As shown below, correlations in the initial state may in fact change the entire evolution
of a stochastic process. Finally, we will consider systems with mixed boundary conditions
such as DP in a parabola-shaped space-time geometry. Mixed boundary conditions may
be viewed as moving boundaries, i.e., the system size varies with time.
64

<!-- page 65 -->

DP with an absorbing wall
In DP an absorbing wall may be introduced by cutting all bonds crossing a given (d-1)-
dimensional hyperplane in space (see Fig. 27). Hence for p > pc the stationary density of
active sites close to the wall ρstat
s
is expected to be smaller than the density in the bulk.
In fact, the density at the wall is found to scale as
ρstat
s
∼(p −pc)βs
(169)
with a surface critical exponent βs > β. The problem of an absorbing wall was ﬁrst studied
in the simpler case of CDP where a surface exponent βCDP
s
= 2 was found [193,194]. In a
series of papers this scaling theory was later applied to DP with an absorbing wall [195–198]
(for a review see [199]). By means of series expansions and numerical simulations it was
observed that the mean survival time T of a cluster in the inactive phase next to the wall
scales as T ∼∆−τs, where ∆denotes the distance from criticality. In 1+1 dimensions
the exponent τs was estimated by 1.0002(3), leading to the remarkable conjecture βs =
ν∥−1 [196].
However, very recent series expansions favor the value 1.00014(2)̸ = 1,
disproving the conjecture [200]. In fact, in view of dimensional analysis it seems to be
unlikely that βs and ν∥are related by a simple linear scaling relation.
Moreover, in
2+1 dimensions the numerical value τs = 0.26(2) cannot be simply related to the other
exponents. Similarly, the ﬁeld-theoretic one-loop result [198]
τs = −1/2 + 11ǫ/48 + 0(ǫ2) ,
βs = 3/2 −7ǫ/48 + 0(ǫ2)
(170)
indicates that the surface exponent is generally independent of the other exponents (168).
The ﬁeld-theoretic analysis was also extended to DP with an absorbing edge [198]. A
closely related application is the study of spreading processes in narrow channels [201]. It
is also interesting to study DP with an active wall. This case is related to the problem of
local persistence and will be discussed below.
DP clusters in a parabola
Several years before the problem of an absorbing wall was investigated, Kaiser and Turban
considered the much more complicated problem of DP in a parabola-shaped geometry [202,
203]. Assuming an absorbing boundary of the form x = ±ctσ they proposed a general
scaling theory. It is based on the observation that the width c of the parabola c scales
as c →Λzσ−1c under rescaling (92). Therefore, the boundary is relevant for σ > 1/z and
irrelevant otherwise. To implement this scaling theory, the scaling forms (101)-(102) have
to be extended by an invariant argument of the form tσ−d/z/c. The survival probability
of a cluster (102), for example, has to be generalized by
P(t) ∼t−δ g(∆t1/ν∥, td/z/N, tσ−d/z/c) .
(171)
This scaling form is supported by numerical results and a mean ﬁeld approximation [203].
The authors also derived a conjecture for the fractal dimensions
d∥(σ) = 1 −zσ(d∥−1) ,
d⊥(σ) = d∥(σ)/σ ,
(172)
where d∥= (β + γ)/ν∥.
65

<!-- page 66 -->

10
0
10
1
10
2
10
3
10
4
t
0.1
1
ρ(t)
10
−2
10
−1
10
0
10
1
10
2
ρ0t
δ+θ
0.1
1
ρ(t)t
δ
ρ0=1
ρ0=0.3
ρ0=0.1
ρ0=0.03
ρ0=0.01
θ
−δ
(a)
(b)
Fig. 28:
“Critical initial slip” of the particle density measured in a (1+1)-dimensional directed bond
percolation process at criticality on a lattice with 104 sites. (a): Particle density ρ(t) for various initial
densities ρ0 as a function of time. The dashed lines indicate the slopes +θ and −δ. (b): Data collapse of
the same data according to Eq. (175).
Early-time behavior and critical initial slip
In Sec. (3.4) we reviewed two Monte Carlo techniques for systems with phase transitions
into absorbing states which diﬀer in their initial state. In simulations starting with a fully
occupied lattice the particle density at criticality decreases as ρ(t) ∼t−β/ν∥. On the other
hand, in dynamic simulations starting from a single particle (active seed), we observe an
increase of the average number of particles as N(t) ∼tθ. In general the exponent θ is
independent from the bulk exponents β, ν⊥, ν∥. In the case of DP, however, the duality
symmetry under time reversal (see Eq. (119)) implies the additional hyperscaling relation
θ = (dν⊥−2β)/ν∥.
(173)
An interesting crossover phenomenon between initial increase and asymptotic decay of the
number of particles emerges when a critical spreading process starts with a low-density
distribution of active sites. Fig. 28a shows the temporal behavior of the density of active
sites ρ(t) for various initial densities ρ0. The density ﬁrst increases as ρ(t) ∼tθ until it
reaches a maximum value at time tc when it crosses over to the usual asymptotic decay
ρ(t) ∼t−β/ν∥. This phenomenon is sometimes referred to as the critical initial slip of
nonequilibrium systems. As can be seen in Fig. 28a, the curves converge to a single one
after suﬃciently long time when the memory of the initial condition is lost. The crossover
time tc depends on the initial density ρ0 and scales as
tc ∼ρ
−1/(β/ν∥+θ)
0
.
(174)
In ﬁnite-size systems near criticality the critical initial slip may be described by adding
the scale-invariant argument ρ0tβ/ν∥+θ to the scaling form (101), i.e.
ρ(t) ∼t−β/ν∥f

∆t1/ν∥, td/z/N, ρ0tβ/ν∥+θ
.
(175)
The scaling function f behaves asymptotically as f(0, 0, u) ∼u for u →0 and f(0, 0, u) =
const for u →∞. To verify this scaling form at criticality, we have plotted ρ(t)tδ versus
ρ0tδ+θ in Fig. 28b. As can be seen, we obtain a convincing data collapse.
66

<!-- page 67 -->

The critical initial slip in a DP process can be interpreted as follows. In a low-density
initial state the active sites are separated by empty intervals of a certain average size
ξ0. As time evolves, they generate individual clusters of connected sites (see Sec. 3.1).
Initially these clusters are spatially separated; they do not interact and the particle number
therefore increases as tθ.
Only a fraction t−δ of these clusters survives, each of them
spanning a volume of ξd
⊥. These surviving clusters start touching each other when ξd
⊥∼
ρ−1
0 tδ
c. Therefore, we expect the crossover to take place at tc ∼ρ1/(δ−d/z)
0
. Insertion of
the DP hyperscaling relation (173) leads to Eq. (174). It is worth being mentioned that
dynamic simulations starting from a single particle represent the limit ρ0 →0. In this
case tc diverges and the critical initial slip extends to the entire temporal evolution of the
system.
Correlated initial conditions
The previously discussed early-time behavior shows that initial states with short-range
correlations may aﬀect the temporal evolution of a DP process for a limited time until the
system crosses over to the usual decay of the particle density. Let us now turn to initial
states with long-range correlations of the form
⟨sisi+r⟩∼rσ−d ,
(176)
where 0 ≤σ ≤d controls the power-law decay of the correlations on large scales. In
one dimension such states can be generated by creating uncorrelated empty intervals12
of length ℓwhich are algebraically distributed as P(ℓ) ∼ℓ−1−σ.
There are, however,
many possibilities to create such states because higher order correlations can be chosen
freely. Apart from cutoﬀs, the resulting particle conﬁgurations do not exhibit a speciﬁc
length scale ξ0, instead they are characterized by a fractal dimension df = σ. It turns out
that long-range correlations may change the entire temporal evolution of a DP process
(similar phenomena can be observed in other nonequilibrium critical systems such as in
the annihilation model [205]).
For σ = d the particles are homogeneously distributed, leading to the usual long-time
behavior ρ(t) ∼t−β/ν∥. For σ →0 the fractal dimension tends to zero, corresponding
to isolated particles where we expect the density to increase a ρ(t) ∼tθ.
In between
numerical simulations suggest that the decay exponent changes continuously by
ρ(t) ∼ρ0tα(σ) ,
α(σ) =
(
θ
for σ ≤σc
1
z(d −σ −β/ν⊥)
for σ > σc
,
(177)
where σc = β/ν⊥plays a role of a critical threshold above which correlations in the initial
state become relevant (see Fig. 29). Below σc the initial distribution of particles is so
sparse that interactions between growing clusters turn out to be irrelevant.
The numerical result can be proven by a simple ﬁeld-theoretic calculation [204]. In
order to take the initial state into account, the ﬁeld-theoretic action (142)-(143) has to be
extended by the term
Sic = µ
Z
ddx ˜φ(x, 0) φ0(x) .
(178)
12The generation of fractal distribution is crucial since for σ < d the particle density is zero. Therefore,
appropriate cutoﬀs have to be introduced, as described in detail in Ref. [204].
67

<!-- page 68 -->

0.0
0.2
0.4
0.6
0.8
1.0
σ
−0.2
−0.1
0.0
0.1
0.2
0.3
α
Fig. 29:
Correlated initial conditions in a DP process. Numerical estimates for the exponents α(σ) in one
spatial dimension. The theoretical prediction (177) is shown by the solid line. The dashed line indicates
the ‘natural’ correlations of DP.
Here the initial particle distribution is represented by a ﬁeld φ0(x) that is coupled to
the ‘creation operator’ ˜φ at time t = 0 via a coupling constant µ. The scaling behavior
of φ0(x) depends on the fractal dimension.
Obviously, homogeneous initial conditions
φ0(x) = const are invariant under rescaling, whereas a fully localized initial condition
φ0(x) = δd(x) has the scaling dimension −d. Therefore, the initial density should scale as
φ0(x) →Λdf −dφ0(x) = Λσ−dφ0(x) .
(179)
It is easy to verify that the contribution Sic does not lead to additional loop corrections in
the ﬁeld theory. Therefore, the coupling between the system and the initial condition will
not be renormalized. Moreover, it can be shown that higher-order contributions of the form
˜φk(x, 0) φ0(x) with k > 1 are irrelevant under renormalization [130, 204]. Consequently,
the coupling constant µ scales as
µ →Λσ+χµ ,
(180)
where χ = −β/ν⊥. Hence the correlations in the initial state are relevant if σ > σc = β/ν⊥.
Scaling invariance of the expression ρ(t) ∼ρ0tα(σ) implies that αz = d−σ +χ, completing
the proof of Eq. (177).
Interestingly, the above calculation does not depend on the speciﬁc form of correlations
in the initial state but only on the scaling dimension of the distribution. That is, no matter
how the particles are distributed – as long as the distribution scales as in Eq. (179), the
particle density decays according to the scaling form (177). Moreover, it is interesting to
note that a critical DP process itself generates two-point correlations ⟨sisi+r⟩∼r−β/ν⊥,
corresponding to the ‘natural’ fractal dimension df,⊥= d −β/ν⊥. Therefore, choosing
‘natural’ correlations σ = d −β/ν⊥the number of particles remains almost constant (see
dashed lines in Fig. 29).
Similar phenomena have been observed in the Glauber-Ising
model with correlated initial conditions [206].
68

<!-- page 69 -->

Persistence probability in a DP process
In the past few years it has been realized that certain ﬁrst passage quantities of critical
nonequilibrium processes exhibit a power law decay with non-trivial exponents. One of
these quantities is the local persistence probability Pl(t), deﬁned as the probability that
a local variable si(t) at a given site i does not change its state until time t during the
temporal evolution. In various systems it was found that
Pl(t) ∼t−θl ,
(181)
where θl is the so-called local persistence exponent [207–213]. A similar quantity, the global
persistence probability Pg(t), which is deﬁned as the probability that the global order
parameter does not change its sign up to time t, is also found to decay as a power law with
a global persistence exponent θg [214–216]. In general the exponents θl and θg are diﬀerent
and independent of the usual scaling exponents. Since the persistence probabilities depend
on the history of evolution as a whole13, it is generally hard to determine these exponents
analytically. In fact, only a few exact results have been obtained so far [208, 210, 216].
Persistence exponents are known to exhibit certain universal properties. For example, the
local persistence exponent of the two-dimensional Glauber model in the ordered phase
T < Tc does not depend on T [217–219], whereas it is non-universal with respect to the
initial magnetization [210].
Most researchers believe that persistence exponents are to
some extent ‘less universal’ than ordinary bulk exponents.
In a DP process the local persistence probability Pl(t) may be deﬁned as the probability
that a particular site never becomes active up to time t. Numerical simulations suggest
that the local persistence exponent is given by
θl = 1.50(1)
(182)
independent of the initial density of active sites [220] . Moreover, the numerical data are
in agreement with the scaling form
Pl(t, N, ∆) ∼t−θl f(∆t1/ν∥, N −z/dt) ,
(183)
where ∆= p−pc denotes the distance from criticality. The local persistence probability can
also be related to certain return probabilities in a DP process with an absorbing boundary
or an active source [220]. Similar measurements of the global persistence probability Pg(t)
suggest that θg > θl in agreement with recent ﬁeld-theoretic results [221].
3.7
The inﬂuence of quenched disorder
In DP models it is usually assumed that the percolation probability does not vary in space
and time. However, in realistic spreading processes the rate for oﬀspring production is
not homogeneous, rather it ﬂuctuates in space and time. For example, the local density of
open channels in a porous rock will vary because of inhomogeneities of the material. Sim-
ilarly, most spreading processes take place in inhomogeneous environments. It is therefore
important to investigate how quenched disorder aﬀects the critical properties of a spread-
ing process. It turns out that even weak disorder may aﬀect or even destroy the critical
behavior of DP.
13Persistence probabilities may be regarded as autocorrelation functions involing inﬁnitely many points.
69

<!-- page 70 -->

In the DP Langevin equation (139) the parameter κ plays the role of the percolation
probability. Quenched disorder may be introduced by random variations of κ, i.e., by
adding another Gaussian noise ﬁeld η:
κ →κ + η(x, t) .
(184)
Thus, the resulting Langevin equation reads
∂tρ(x, t) = κρ(x, t) −Λρ2(x, t) + D∇2ρ(x, t) + ζ(x, t) + ρ(x, t)η(x, t) .
(185)
The noise η is quenched in the sense that quantities like the particle density are averaged
over many independent realizations of the intrinsic noise ζ while the disorder ﬁeld η is
kept ﬁxed.
In the following we distinguish three diﬀerent types of quenched disorder,
namely spatially, temporally quenched, and spatio-temporally quenched disorder.
The
three variants of quenched disorder diﬀer in how far they aﬀect the critical behavior of
DP.
Spatially quenched disorder
For spatially quenched disorder, the disorder ﬁeld η is deﬁned through the correlations
η(x)η(x′) = γ δd(x −x′) ,
(186)
where the bar denotes the average over many independent realizations of the disorder ﬁeld
(in contrast to the ensemble average ⟨. . . ⟩over realizations of the intrinsic noise ζ). The
parameter γ is an amplitude controlling the intensity of disorder. In order to ﬁnd out
whether this type of noise aﬀects the critical behavior of DP, let us again consider the
properties of the Langevin equation under rescaling. At the critical dimension dc = 4 the
additional term ρη scales as
ρη →Λ−dc/2−χρη ,
(187)
i.e., spatially quenched disorder is a marginal perturbation. Therefore, it may seriously
aﬀect the critical behavior at the transition. The same result is obtained by considering the
ﬁeld-theoretic action. Without quenched noise, DP is described by the action of Reggeon
ﬁeld theory (see Sec. 3.5)
S =
Z
ddx dt ˜φ
h
∂t −κ −D∇2 + Γ
2 (φ −˜φ)
i
φ
(188)
where φ(x, t) represents the local particle density while ˜φ(x, t) denotes the Martin-Siggia-
Rosen response ﬁeld. As shown by Janssen [222], spatially quenched noise can be taken
into account by adding the term
S →S + γ
Z
ddx
Z
dt ˜φφ
2
.
(189)
By simple power counting one can prove that this additional term is indeed a marginal
perturbation. Janssen showed by a ﬁeld-theoretic analysis that the stable ﬁxed point is
shifted to an unphysical region, leading to runaway solutions of the ﬂow equations in the
physical region of interest. Therefore, spatially quenched disorder is expected to crucially
disturb the critical behavior of DP. The ﬁndings are in agreement with earlier numerical
70

<!-- page 71 -->

t
Fig. 30:
(1+1)-dimensional DP with spatially quenched disorder. Left: In the glassy phase quenched
disorder forces active sites to percolate in narrow ‘channels’ where the local percolation probability is high.
Right: Supercritical disordered DP process starting from a single seed, leading to avalanches (marked by
the arrows) where the spreading agent overcomes a local barrier.
results by Moreira and Dickman [223] who reported non-universal logarithmic behavior
instead of power laws. Later Caﬁero et al. [224] showed that DP with spatially quenched
randomness can be mapped onto a non-Markovian spreading process with memory, in
agreement with previous results.
From a more physical point of view, spatially quenched disorder in 1+1 dimensional
DP was studied by Webman et al. [225]. It turns out that even weak disorder drastically
modiﬁes the phase diagram. Instead of a single critical point one obtains an intermediate
phase of very slow glassy-like dynamics. The glassy phase is characterized by non-universal
exponents which depend on the percolation probability and the disorder amplitude. For
example, in a supercritical 1+1 dimensional DP process without quenched disorder the
boundaries of a cluster propagate at constant average velocity v. However, in the glassy
phase v decays algebraically with time. The corresponding exponent turns out to vary
continuously with the mean percolation probability. The power-law behavior is due to
‘blockades’ at certain sites where the local percolation probability is small (see Fig. 30).
Similarly, in the subcritical edge of the glassy phase, the spreading agent becomes localized
at sites with high percolation probability.
For d > 1, however, numerical simulations
indicate that a glassy phase does not exist.
Temporally quenched disorder
Temporally quenched disorder is deﬁned by the correlations
η(t)η(t′) = γ δ(t −t′) .
(190)
In this case the additional term is a relevant perturbation which scales as ρη →Λ−z/2−χρη.
Therefore, we expect the critical behavior and the associated critical exponents to change
entirely. In the ﬁeld-theoretic formulation this corresponds to adding a term of the form
S →S + γ
Z
dt
Z
ddx ˜φφ
2
.
(191)
The inﬂuence of temporally quenched disorder was investigated in detail in Ref. [116]. Em-
ploying series expansion techniques it was demonstrated that the three exponents β, ν⊥, ν∥
71

<!-- page 72 -->

x
t
Fig. 31:
Spreading process in a coupled map lattice for r = 3 and D = 0.57. Chaotic and laminar sites are
represented by black and white pixels, respectively. Deterministic synchronous updates lead to symmetric
clusters (left) with nonuniversal behavior. Asynchronous updates (right) destroy these correlations, leading
to typical DP clusters.
vary continuously with the disorder strength. Thus the transition no longer belongs to
the DP universality class.
Spatio-temporally quenched disorder
For spatio-temporally quenched disorder, the noise ﬁeld η is uncorrelated in both space
and time:
η(x, t)η(x′, t′) = γ δd(x −x′)δ(t −t′) .
(192)
In Reggeon ﬁeld theory, this would correspond to the addition of the term
S →S + γ
Z
ddx dt
h
˜φφ
i2
,
(193)
being an irrelevant perturbation. In fact, this noise has essentially the same properties as
the intrinsic noise and can be considered as being annealed. Spatio-temporally quenched
disorder is expected in systems where each time step takes place in a new spatial environ-
ment of the system. Examples include water in porous media subjected to a gravitational
ﬁeld as well as systems of ﬂowing sand on an inclined plane (see Sec. 3.9). In these cases
the critical behavior of DP should remain valid on large scales.
3.8
Related models
Directed percolation plays a role in various other contexts such as in coupled map lattices,
the problem of friendly walkers, real-valued spreading processes, models with particle
conservation, and even in systems with inﬁnitely absorbing states. In the following we
discuss some of these related models. Moreover, we investigate the special case of compact
directed percolation in more detail.
Spreading transitions in deterministic systems
Spreading transitions can also be observed in certain deterministic lattice models. Instead
of using random numbers, these models employ chaotic maps in order to generate random
72

<!-- page 73 -->

behavior. A simple example of such a chaotic map is given by
u(t + 1) = f
 u(t)

,
f(x) =





rx
if 0 ≤x < 1/2 ,
r(1 −x)
if 1/2 < x ≤1 ,
x
if 1 < x ≤r/2 ,
(194)
where r is a free parameter. The chaotic motion of f for x ≤1 is governed by a tent map
of slope r. However, if r exceeds the value 2, the map eventually reaches an absorbing
state with x > 1, the so-called ‘laminar’ state of the model. In a coupled map lattice [226]
many of these local maps ui(t) are coupled by a diﬀusive interaction of the form
ui(t + 1) = f
 ui(t)

+ D
2

f
 ui−1(t)

−2f
 ui(t)

+ f
 ui+1(t)

,
(195)
where D plays the role of a diﬀusion constant. The coupled map lattice evolves deter-
ministically by synchronous updates. By varying D it exhibits a nonequilibrium phase
transition from a ‘chaotic’ phase into a ‘laminar’ state. The existence of absorbing states
led Pomeau to the conjecture that the transition should belong to the DP universality
class [227], hoping that the apparent randomness of the chaotic maps would eﬀectively
lead to stochastic spreading of activity on large scales. However, subsequent numerical
simulations did not agree with this conjecture [228], in particular the exponents were found
to depend on r. The non-universal behavior of spreading transitions in deterministic sys-
tems is caused by subtle correlations emerging as artifacts of the deterministic update
rule. For example, a ‘cluster’ of chaotic sites starting from a single active seed remains
symmetric throughout the whole temporal evolution, leading to a qualitatively diﬀerent
spreading behavior (see Fig. 31). The consequences of these correlations are not yet fully
understood. However, replacing the synchronous dynamics of Eq. (195) by asynchronous
updates, the deterministic correlations are destroyed and the resulting phase transition
is indeed characterized by DP exponents [143].
Similar transitions of two-dimensional
coupled map lattices have been investigated in Ref. [229].
DP and the problem of ‘friendly walkers’
The so-called problem of ‘friendly walkers’ is deﬁned as follows. Consider the paths of m
random walkers on a diagonal square lattice. All walks originate in site (0, 0) and end in
site (x, t). While traveling the walkers may share the same bonds but they are not allowed
to cross each other (see Fig. 32). In the partition sum Zm(x, t) each possible conﬁguration
of random walks is weighted by a factor pk, where p > 0 is a parameter and k denotes the
number of bonds used by at least one of the walkers. For p < 1 it is therefore advantageous
for the walkers to be ‘friendly’ to each other, i.e., to share the same bonds.
Some time ago, Arrowsmith and Essam [230] suggested a close relationship between
DP and the problem of friendly walkers. More precisely, they showed that the partition
function Zm(x, t) is related to the pair-connectedness function c(x, t) of a directed bond
percolation process by [231]
c(x, t) = lim
m→0 Zm(x, t) ,
(196)
where p is the usual percolation probability. Here the limit m →0 has to be performed as
a suitable continuation of polynomial expressions. For example, let us consider m friendly
random walkers traveling from the origin (0, 0) to the point (x, t) = (0, 2) (see right part of
73

<!-- page 74 -->

m
n
m-n
(x,t)=(0,4)
m
(0,0)
...
...
Fig. 32:
The problem of friendly walkers on a diagonal square lattice. Left: A particular realization of
three random walkers traveling from the origin (0, 0) to the site (0, 4). Right: All m + 1 realizations of m
random walkers traveling from (0, 0) to (0, 2).
Fig. 32). There are m+1 possible conﬁgurations; two of them use only two bonds while the
others use four bonds. Hence the partition function is given by Zm(0, 2) = (m−1)p4+2p2.
Inserting m = 0 we obtain Z0(0, 2) = 2p2 −p4. In fact, this expression is exactly equal
to the pair connectedness function c(0, 2) in a directed bond percolation process. This
equivalence holds for any (x, t) and also in higher dimensions. Recently this result has
been generalized to friendly walkers with arbitrary interactions [232].
The problem of friendly walkers may also be interpreted as a ﬂow of integer numbers
on a diagonal square lattice. At the origin there is a source creating an integer number m.
While traveling on the directed lattice, this integer number may split up into several parts.
Finally there is a sink where all integers merge into a single one and disappear. Clearly,
the integers represent just the number of friendly walkers sharing the same bond.
Even more remarkably, it has been shown that DP is related to the partition sum
of a chiral Potts model [230, 231, 233], generalizing the well-known result of Fortuin and
Kasteleyn for isotropic percolation [234]. However, since the deﬁnition of the chiral Potts
model is rather cumbersome, this relation is not of immediate practical beneﬁt.
DP with real-valued degrees of freedom
DP models are usually deﬁned in terms of discrete local variables si = 0, 1 representing
inactive and active sites. An interesting variant of DP is ‘self-organized directed percola-
tion’, where real-valued local degrees of freedom are used [235–238]. To understand the
basic mechanism, let us consider directed bond percolation. Clearly, a given path between
two sites is conducting if all bonds along the path are open, i.e., all random numbers
generated along the path have to be larger than p. Thus, in order to ﬁnd out whether
a path is conducting, it is only necessary to keep track of the smallest random number
generated along this path. This number may be considered as the weight of the path,
being a measure of its weakest link. However, a pair of sites can be connected by many
diﬀerent paths. For the target site to become active, at least one of these paths has to
be conducting. Therefore, two sites are connected if the maximum of all weights is larger
than p.
Interestingly, the maximal weight can be computed by a local update rule which is
deﬁned in terms of real-valued degrees of freedom xi(t) ∈[0, 1]. Starting with the initial
74

<!-- page 75 -->

0
200
400
600
800
1000
i
0.6
0.7
0.8
0.9
1
xi(t=10
4)
1
10
100
1000
10000
t
0
0.2
0.4
0.6
0.8
min {xi(t)}
Fig. 33:
Self-organized directed percolation in 1+1 dimensions. Left: Spatial conﬁguration of the chain
after 104 updates. The arrow indicates the percolation threshold of directed bond percolation. Right:
Minimum value of all sites as a function of time, approaching the percolation threshold of bond DP
(dashed line).
condition xi(0) = 0 the system evolves according to
xi(t + 1) = min

max(z−
i , xi−1(t)), max(z+
i , xi−1(t))

,
(197)
where we used the notation of Eq. (80).
A typical spatial conﬁguration of a (1+1)-
dimensional chain after 104 updates is shown in Fig 33. Using this update rule, the binary
state si(t) of the corresponding directed bond percolation process can be retrieved by the
projection
si(t) = Θ(p −xi(t)) ,
(198)
where Θ denotes the heaviside step function. Remarkably, the update rule (197) does
not involve the percolation probability p.
Instead, it processes all values of p at once
until a particular value of p is selected by application of the projection rule (198). Thus,
‘self-organized directed percolation’ can be used as a tool for very eﬃcient oﬀ-critical
simulations.
Spreading process with particle conservation
Recently Broeker and Grassberger [239] introduced another interesting ‘self-organized’
variant of DP which is motivated as follows.
A gardener takes care of N plants in a
ﬂowerbed.
The ﬂowers are seized with a parasite.
Once a plant is struck, it perishes
irreversibly. Moreover, the parasite may spread to neighboring plants. However, if the
number of befallen plants exceeds a certain number, the gardener replaces one of them,
keeping the number of infected plants constant. In more technical terms, the number of
active sites is conserved by means of a global update rule. The update consists of two
steps. At ﬁrst one of the active sites activates a randomly chosen neighbor, modeling the
spreading of the parasite. If this move was successful, another randomly chosen active
site is deactivated, representing the global control of the gardener. Clearly, the number of
active sites M is conserved, i.e., the model has no absorbing states.
The number of active sites M is speciﬁed by the initial state. For example, we may
start with a compact domain of M active sites on an inﬁnite lattice. Initially, spreading
75

<!-- page 76 -->

occurs only at the edges of the domain. As time proceeds, the distribution of active sites
becomes more and more sparse, forming a diﬀusing cloud. Nevertheless, the cloud keeps
its integrity and reaches a typical size after some time. Amazingly, the dynamic processes
in the interior of the cloud are those of an almost critical DP process. In fact, as shown
in Ref. [239], most properties of the cloud can be explained in terms of DP scaling laws.
Considering a small region in the interior of the cloud, the relation to DP is quite obvious:
The two processes for oﬀspring production A →2A and self-destruction A →Ø occur
randomly in space, just as in a contact process with random-sequential updates. However,
the global control adjusts the ratio of the eﬀective rates and drives the system to criticality.
Branching Potts interfaces
Recently Cardy [240] studied a ﬁeld-theory for branching interfaces between ordered
domains of a q-state Potts model.
In two spatial dimensions these interfaces are one-
dimensional objects. For q < qc they become fractal with a vanishing interfacial tension
at the critical point, while for q > qc the interfacial width diverges at a ﬁnite value of
the tension, indicating a ﬁrst-order transition. In a certain limit, namely q →−∞, the
model becomes equivalent to a DP process. Therefore, the model provides a ﬁeld theory
of directed percolation that diﬀers from the standard ﬁeld theory discussed in Sec. 3.5.
Although both ﬁeld theories ‘intersect’ in one dimension, they are completely diﬀerent.
In particular, the loop expansion starts out from diﬀerent critical dimensions, namely
dc = 4 for Reggeon ﬁeld theory and dc = 2 for branching Potts interfaces. Consequently,
in the latter case the one-loop estimates for the critical exponents in d = 1 are much more
accurate.
DP models with inﬁnitely many absorbing states
According to the DP conjecture, phase transitions into a single absorbing state belong
generically to the DP universality class.
However, DP behavior may also be observed
in models with several or even inﬁnitely many absorbing states. An interesting example
is the dimer-trimer model for heterogeneous catalysis introduced by K¨ohler and ben-
Avraham [241]. This model generalizes the ZGB model and is deﬁned by the reaction
scheme
ØØ →AA
at rate p ,
ØØØ →BBB
at rate 1 −p ,
(199)
AB →ØØ
at rate ∞.
On an inﬁnite lattice this model has inﬁnitely many absorbing states. For example, con-
ﬁgurations of dimers and trimers separated by single vacant sites are absorbing.
The
dimer-trimer model displays a phase transition in 2+1 dimensions. Initially, the values
of the critical exponents were found to be diﬀerent from those of DP. Later reﬁned sim-
ulations conﬁrmed, however, that the dimer model still belongs to the DP universality
class [242].
The same result was found in a similar model for catalysis of dimers and
monomers [243,244]. Another important example is the pair contact process without dif-
fusion [245] which is deﬁned by the reaction scheme
2A →3A ,
2A →Ø .
(200)
In this model solitary particles neither react nor diﬀuse. Starting from random initial
conditions, the critical pair contact process evolves into certain frozen conﬁgurations, as
76

<!-- page 77 -->

t
Fig. 34:
The pair contact process 2A →3A, 2A →Ø at criticality.
demonstrated in Fig. 34. As can be seen, it is important that single particles are not
allowed to diﬀuse. In fact, by adding diﬀusion of individual particles the critical behavior
of the model changes entirely (see Sec. 4.5).
In all models with inﬁnitely many absorbing states and non-conserved order parameter
the critical exponents β, ν∥, ν⊥coincide with those of DP. This observation suggests an
extension of the DP conjecture to systems with several absorbing states which are char-
acterized by a non-conserved single-component order parameter [246]. However, it was
realized that the dynamic exponents δ and θ depend on the initial condition and even
violate the usual DP hyperscaling relation (120). Mendes et al. [129] resolved this prob-
lem by introducing the generalized hyperscaling relation (100). However, the sum δ + θ is
believed to be independent of the initial condition.
Recently Mu˜noz et al. proposed a Langevin equation for systems with inﬁnitely many
absorbing states [147]. It diﬀers from the usual DP Langevin equation (139) by an addi-
tional term:
∂tφ(x, t) =κφ(x, t) −λφ2(x, t) + D∇2φ(x, t) + ζ(x, t) +
+ αφ(x, t) exp
h
−w
Z t
0
dt′ φ(x, t′)
i
.
(201)
Here α and w are certain constants (the noise correlations are assumed to be the same as
in Eq. (140)). This Langevin equation is non-Markovian, i.e., it has a temporal memory.
The memory is local since the integral correlates ﬂuctuations at the same position in
space.
From the physical point of view, this memory encodes the local realization of
the absorbing state. As can be seen in Fig. 34, the emerging inactive domains have a
highly inhomogeneous structure which can be regarded as a ﬁngerprint of the history
of the spreading process. A detailed numerical analysis of the Langevin equation (201)
conﬁrmed that the exponenets β, ν∥, and ν⊥do belong to the DP class while δ and θ vary
with the density of the initial state [247]. In order to explain the apparent nonuniversality
of the spreading exponents, Grassberger developed a simple toy model that grasps the
main properties of such spreading processes [248]. In this toy model the spreading rate at
a given site changes irreversibly at the ﬁrst encounter with the spreading agent. Although
the model does not involve multiple absorbing states, it displays similar ‘nonuniversal’
properties.
77

<!-- page 78 -->

Fig. 35:
Snapshot of a critical epidemic process with immunisation grown from a single seed after 200
time steps. Active and immune sites are represented by bold and thin dots, respectively. An zone of
activity propagates outwards, leaving a cluster of immune sites behind. The resulting cluster belongs to
the universality class of isotropic percolation.
Another important example for systems with inﬁnitely many absorbing states is dam-
age spreading where two copies of a stochastic system evolve under the same realization
of thermal noise. The concept of damage spreading will be discussed in detail in Sec. 5.
Epidemic processes with immunization
As we have seen in Sec. 3.1, epidemic models without immunization belong generically
to the DP universality class.
In most cases, however, an infected individual becomes
increasingly immune after recovery, i.e., the susceptibility for a new infection decreases.
Cardy and Grassberger showed that epidemic models with immunization are in the same
universality class as dynamic percolation [249,250]. It is important to note that dynamic
percolation diﬀers signiﬁcantly from directed percolation. For example, let us consider
a spreading process with immunization in 2+1 dimensions starting from an initial state
where all sites are equally susceptible for infections. If a single site in the center is infected,
there is a ﬁnite probability that the disease will spread. However, since infected sites
become increasingly immune, a more or less irregular front of activity moves away from
the origin, leaving behind a certain cluster of immune sites (see Fig. 35). The morphology
of this cluster depends on the percolation parameter. In the supercritical case there is a
ﬁnite probability that the front moves to inﬁnity, whereas in the subcritical regime the
process stops after some time. At criticality it turns out that the generated cluster has
the same asymptotic properties as critical clusters of isotropic percolation [109] (cf. left
part of Fig. 10). Thus, dynamic percolation can be used as a tool to generate isotropic
percolation clusters and should not be confused with directed percolation. In particular,
the critical exponents turn out to be diﬀerent in both cases. Interestingly, even a small
degree of immunization suﬃces for a (2+1)-dimensional epidemic process to cross over
from directed to dynamic percolation (together with a shift of the critical point).
A
renormalized ﬁeld theory of dynamic percolation was studied in [251].
Compact directed percolation
Let us ﬁnally come back to compact directed percolation (CDP) [117] which characterizes
the critical behavior of the DK model at the upper terminal point of the phase transition
78

<!-- page 79 -->

line p1 = 1/2, p2 = 1 (see Fig. 14). The case p2 = 1 is special because there are two
symmetric absorbing states, namely the dry state s1 = . . . = sN = 0 and the entirely wet
state s1 = . . . = sN = 1. In contrast to DP, CDP has a global Z2 symmetry
si →1 −si ,
p1 →1 −p1 .
(202)
Since wet sites cannot spontaneously become dry, compact islands of active sites are
formed. In 1+1 dimensions CDP is fully equivalent to a zero temperature Ising model with
Glauber dynamics or the voter model [252]. Expressing the dynamic processes in terms
of kinks X between wet and dry domains, the kinks perform an annihilating random
walk X + X →Ø.
Therefore, (1+1)-dimensional CDP is exactly solvable [113].
The
corresponding critical exponents are given by
β = 0 ,
β′ = 1 ,
ν∥= 2 ,
ν⊥= 1 ,
δ = 1/2 ,
θ = 0 ,
˜z = 1 .
(203)
It should be noted that these exponents do not comply with the usual DP hyperscaling
relation (120).
However, as pointed out in Ref. [253], they satisfy the generalized hy-
perscaling relation (100). In fact, as can be veriﬁed easily, for CDP the backbone of a
two-point function (see Sec. 3.5) is no longer statistically invariant under time reversal.
Because of the vanishing exponent β, the CDP transition is discontinuous. In fact, for
p1 < 1/2, p2 = 1 the empty lattice is a stable stationary state while for p1 > 1/2 the fully
occupied lattice is stable. Various spreading models display a crossover from CDP to DP.
In these models, the rate for the reaction A →Ø is very small. Therefore, clusters appear
to be compact on small scales. On larger scales, however, clusters break up into several
active branches, leading to DP behavior in the asymptotic limit. This type of crossover
has been studied in detail in Refs. [254–257] and may also play a role in experiments of
ﬂowing granular matter (see Sec. 3.9).
3.9
Experimental realizations of directed percolation
So far we have seen that directed percolation is the generic universality class for nonequilib-
rium phase transitions into absorbing states. In fact, DP seems to be of similar importance
as the Ising model in equilibrium statistical mechanics. Despite this success in theoretical
statistical physics, the critical behavior of DP, especially the values of the critical expo-
nents, have not yet been conﬁrmed experimentally. The lack of experimental evidence is
indeed surprising, especially since a large number of possible experimental realizations have
been suggested in the past. As Grassberger emphasizes in a summary on open problems
in DP [106]:
”...there is still no experiment where the critical behavior of DP was seen. This
is a very strange situation in view of the vast and successive theoretical eﬀorts
made to understand it. Designing and performing such an experiment has thus
top priority in my list of open problems.”.
What might be the reason for the apparent lack of experimental evidence? It seems that
the basic features of DP, which can easily be implemented on a computer, are quite diﬃcult
to realize in nature. One of these idealized assumptions is the existence of an absorbing
79

<!-- page 80 -->

Fig. 36:
Left: The reaction rate on Pt(210) as a function of PCO (ﬁgure reprinted from Ref. [259]). Right:
STM image of a catalytic reaction on a Pt(111) crystal (ﬁgure reprinted from Ref. [260]).
state. In real systems, however, a perfect non-ﬂuctuating state cannot be realized. For
example, a poisoned catalytic surface is not completely frozen, instead it will always be
aﬀected by small ﬂuctuations. Although these ﬂuctuations are strongly suppressed, they
could still be strong enough to ‘soften’ the transition, making it impossible to quantify
the critical exponents.
Another reason might be the inﬂuence of quenched disorder due to spatial or temporal
inhomogeneities. In most experiments frozen randomness is expected to play a signiﬁcant
role. For example, a real catalytic surface is not fully homogeneous but characterized by
certain defects leading to spatially quenched disorder. As has been shown in Sec. 3.7, this
type of disorder may aﬀect or even destroy the critical behavior of DP.
In the following we summarize some of the most important experimental applications
which have been proposed so far [258].
Other experimental applications in systems of
growing interfaces will be discussed in Sec. 6.
Catalytic reactions
It is well known that under speciﬁc conditions certain catalytic reactions mimic the mi-
croscopic rules of DP models. For example, as shown in Fig. 17, the ZGB model for the
catalytic reaction CO + O →CO2 on a platinum surface displays a continuous transition
at y = y1 belonging to DP. In real catalytic reactions, however, only the discontinuous
transition at y = y2 can be observed. Fig. 36 shows the reaction rates as functions of
the CO pressure measured in a catalytic reaction on a Pt(210) surface [259]. Although
this experiment was designed in order to investigate the technologically interesting regime
of high activity close to the ﬁrst-order phase transition, the results clearly indicate that
poisoning with oxygen does not occur. Instead the reactivity increases almost linearly
with the CO pressure. Similar results were obtained for Pt(111) and for other catalytic
materials. Thus, so far there is no experimental evidence for DP transitions in catalytic
reactions.
One may speculate why the DP transition is obscured or even destroyed under exper-
imental conditions. On the one hand, the reaction chain in the experiment is much more
80

<!-- page 81 -->

φ
Fig. 37:
Simpliﬁed drawing of the Douady-Daerr experiment.
Pouring sand on top of a plane with
inclination φ, a thin layer settles and remains immobile. Perturbing the layer locally with a stick we can
release an avalanche of ﬂowing.
complicated than in the ZGB model [261]. Moreover, the O-poisoned system might not
be a perfect absorbing state, i.e., the surface can still adsorb CO molecules although it is
already saturated. Another possibility is thermal (nonreactive) desorption of oxygen, act-
ing as an external ﬁeld which drives the system away from criticality [262] [cf. Eq. (116)].
Finally, defects and inhomogeneities of the catalytic material could lead to an eﬀective
(spatially quenched) disorder.
For a long time the microscopic dynamics processes were diﬃcult to study experimen-
tally. However, in recent years novel techniques such as scanning tunneling microscopy
(STM) led to an enormous progress in the understanding of catalytic reactions, pointing
at various unexpected subtleties. For example, it was observed that reactions preferably
take place at the perimeter of oxygen islands [260].
Furthermore, it was realized that
adsorbed CO molecules on Pt(111) may form three diﬀerent rotational patterns represent-
ing the c(4×2) structure of CO on platinum, i.e., there are several competing absorbing
states [263]. Moreover, the STM technique allows one to trace individual molecular re-
actions and to determine the corresponding reaction rates. In addition, the inﬂuence of
defects such as terraces on catalytic reactions can be quantiﬁed experimentally [264]. We
may therefore expect a considerable progress in the understanding of catalytic reactions
in near future.
Flowing granular matter
Recently it has been shown that simple systems of ﬂowing sand on an inclined plane, such
as the experiments performed by Douady and Daerr [265,266], could serve as experimental
realizations of DP [267]. In the Douady-Daerr experiment glass beads with a diameter of
250-425 µm are poured uniformly at the top of an inclined plane covered by a rough velvet
cloth (see Fig. 37). As the beads ﬂow down, a thin layer settles and remains immobile.
Increasing the angle of inclination φ by ∆φ the layer becomes dynamically unstable, i.e.,
by locally perturbing the system at the top of the plane an avalanche of ﬂowing granular
matter will be released. In the experiment these avalanches have the shape of a fairly
regular triangle with an opening angle θ. As the increment ∆ϕ decreases, the value of θ
decreases, vanishing as
tan θ ∼(∆ϕ)x
(204)
81

<!-- page 82 -->

150 000 time steps
100 sites
800 sites
CDP:
3000 time steps
DP:
Fig. 38:
Typical clusters generated at criticality on small and large scales, illustrating the crossover from
CDP to DP.
with a certain critical exponent x. The experimental results suggest the value x = 1 [266].
In order to explain the experimentally observed triangular form of the avalanches,
Bouchaud et al. proposed a mean-ﬁeld theory based on deterministic equations, taking
the actual local thickness of the ﬂowing avalanche into account [268]. This theory predicts
the exponent x = 1/2. Another explanation assumes that ﬂowing sand may be interpreted
as a nearest-neighbor spreading process [267]. Here the avalanche is considered as a cluster
of active sites. Identifying the vertical coordinate of the plane with time and the increment
of inclination ∆ϕ with p −pc, the opening angle is expected to scale as
tan θ ∼ξ⊥/ξ∥∼(∆ϕ)ν∥−ν⊥,
(205)
where ν∥and ν⊥are the scaling exponents of the spreading process under consideration.
To support this scaling argument, a simple lattice model was introduced which mimics
the physics of ﬂowing sand [267]. The model exhibits a transition from an inactive to an
active phase with avalanches whose compact shapes reproduce the experimental observa-
tions. On laboratory scales the model predicts a transition belonging to the universality
class of compact directed percolation [see Eq. (203)], implying that
x = ν∥−ν⊥= 1.
(206)
The CDP behavior, however, is only transient and crosses over to DP after a very long
time. Thus the Douady-Daerr experiment – performed on suﬃciently large scales – may
serve as a physical realization of DP. Irregularities of the layers thickness can be considered
as spatio-temporally quenched disorder which is irrelevant on large scales (see Sec. 3.7).
Thus, in contrast to catalytic reactions, the problem of quenched disorder does not play a
major role in this type of experiments.
The crossover from CDP to DP is very slow and presently not accessible in experiments.
To illustrate the crossover, two avalanches are plotted on diﬀerent scales in Fig. 38. The
left one represents a typical avalanche within the ﬁrst few thousand time steps. As can be
seen, the cluster appears to be compact. However, as shown in the right panel of the ﬁgure,
the cluster breaks up into several branches after a very long time. Recent experimental
studies [269] conﬁrm that for high angles of inclination critical avalanches do split up
into several branches (see Fig. 39). Yet here the avalanches have no well deﬁned front, the
82

<!-- page 83 -->

Fig. 39:
Left: Triangular avalanches in the Douady-Daerr experiment. Right: Avalanches splitting up
into several branches for high angles of inclination (reprinted with kind permission from A. Daerr).
propagation velocity of separate branches rather depends on their thickness. It is therefore
no longer possible to interpret the vertical axis as a time coordinate. Another problem
is the kinetic energy of the grains.
According to arguments by Dickman et al. [238],
continuous phase transitions into absorbing states can only be observed if the inertia of
particles can be neglected.
Finally, it is not yet known how the spreading process depends on correlations in the
initial state. As shown in Sec. 3.6, such long-range correlations may change the values of
certain dynamic critical exponents. However, recent studies of a single rolling grain on an
inclined rough plane [270] support that there are presumably no long-range correlations
due to a ‘memory’ of rolling grains. By means of molecular dynamics simulations it was
shown that the motion of a rolling grain consists of many small bounces on each grain
of the supporting layer. Therefore, the rolling grain quickly dissipates almost all of the
energy gain from the previous step and thus forgets its history very fast. For this reason
it seems to be unlikely that quenched disorder of the prepared layer involves long-range
correlations. Therefore, ﬂowing granular matter seems to be a promising candidate for an
experimental realization of DP.
Porous media
DP is often motivated as a model for water percolating through a porous medium in a
gravitational ﬁeld. Due to an external driving force, the ﬂow in the medium is assumed
to be strictly unidirectional, i.e., the water can only ﬂow downwards (in contrast to the
depinning models of Sec. 6.2 where the water can ﬂow forth and back). Although this
application seems to be quite natural, it is diﬃcult to realize experimentally. As shown in
Ref. [271], porous media in nature are highly irregular. By cutting sandstone into slices
and digitizing the section images, the porosity distribution and the local connectivity were
measured and averaged over 99 samples. As expected, the pores have diﬀerent sizes and
are distributed irregularly. In addition, the percolation probability is found to depend on
the local porosity and the direction in space, i.e., sandstone is an anisotropic material.
But there are even more fundamental problems. On the one hand, water is a conserved
quantity, leading to unpredictable long-range correlations in the bulk. On the other hand,
water can always ﬂow against the gravity ﬁeld by means of capillary forces. Therefore,
it is quite diﬃcult or even impossible to verify scaling laws in such experiments and is
not yet clear whether the relation to DP is meaningful or simply a commonly accepted
83

<!-- page 84 -->

misconcept.
Epidemics
Another frequently quoted application of DP is the spreading of epidemics without immu-
nization [101, 111]. In an epidemic process infection and recovery resemble the reaction-
diﬀusion scheme of DP (81). If the rate of infection is very low, the infectious disease will
disappear after some time. If infections occur more frequently, the disease may spread
and survive for a very long time.
However, spreading processes in nature are usually
not homogeneous enough to reproduce the critical behavior of DP. Moreover, in many
realistic spreading processes short-range interactions are no longer appropriate. This sit-
uation emerges, for example, when an infectious disease is transported by insects. Such
long-range interactions may be described by L´evy ﬂights, leading to continuously varying
critical exponents (see Sec. 4.1).
Forest ﬁres
A closely related problem is the spreading of forest ﬁres [102,134]. Tephany et al. studied
the propagation of ﬂame fronts on a random lattice both under quiescent conditions and in
a wind tunnel [272]. The experimental estimates of the critical exponents at the spreading
transition are in rough agreement with the predictions of isotropic and directed percolation,
respectively. However, the accuracy of these experiments remains limited.
Calcium dynamics
DP transitions may also occur in certain kinetic models for the dynamics of Calcium ions in
living cells. Ca2+ ions play an important physiological role as second messenger for various
purposes ranging from hormonal release to the activation of egg cells by fertilization [273,
274]. The cell uses nonlinear propagation of increasing intracellular Ca2+ concentration,
so-called calcium waves, as a tool to transmit signals over distances that are much longer
than the diﬀusion length.
For example, propagating Ca2+ waves can be observed in
the immature Xenopus laevis oocyte [275].
So far theoretical work focused mainly on
deterministic reaction-diﬀusion equations in the continuum, explaining various phenomena
such as solitary and spiral waves [276]. Recently improved models have been introduced
which take also the stochastic nature of Calcium release into account [277]. As expected,
the transition in one of these models belongs to the DP universality class [278]. However,
from the experimental point of view it seems to be impossible to conﬁrm or disprove this
conjecture. On the one hand, the size of a living cell is only a few order of magnitude larger
than the diﬀusion length, leading to strong ﬁnite-size eﬀects in the experiment. On the
other hand, inhomogeneities as well as internal structures of the cell lead to a completely
unpredictable form of quenched noise. Therefore, it seems to be impossible to identify the
universality class of the transition in such experiments. It would be rather an achievement
to ﬁnd clear evidence for the very existence of a phase transition between survival and
extinction of propagating calcium waves.
Directed polymers
DP is also related to the problem of directed polymers [279]. In contrast to DP, which is
deﬁned as a local process, the directed polymer problem selects directed paths in a random
medium by global optimization. Under certain conditions, namely a bimodal distribution of
random numbers, both problems were shown to be closely related [280]. More speciﬁcally,
84

<!-- page 85 -->

the roughness exponent of the optimal path in a directed polymer problem is predicted to
cross over from the KPZ value 2/3 to the DP value ν∥/ν⊥≃0.63 at the transition point.
Directed polymers were used to describe the propagation of cracks [281]. However, in such
experiments it is usually impossible to verify the tiny crossover from KPZ to DP.
Turbulence
Finally, DP has also been considered as a toy model for turbulence.
As suggested in
Ref. [227], the front between turbulent and laminar ﬂow should exhibit the critical be-
havior of DP. For example, the velocity of the front should scale algebraically with a
combination of DP exponents. However, these predictions are based rather on heuristic
arguments than on rigorous results. In fact, in many respects turbulent phenomena show
a much richer behavior than DP. Nevertheless there are certain similarities between DP
and turbulence. Therefore, the study of DP could be helpful for a better understanding
of turbulent phenomena.
Summary and outlook
Directed percolation is keeping theoretical physicists fascinated since more than four
decades. Several reasons make directed percolation so appealing. First of all, DP is a
very simple model in terms of its dynamic rules. Nevertheless, the DP phase transition
turns out to be highly nontrivial. In fact, DP belongs to the very few critical phenomena
which have not yet been solved exactly in one spatial dimension. Therefore, the critical
exponents are not yet known analytically.
High-precision estimates indicate that they
might be given rather by irrational than by simple fractional values. Moreover, DP is
extremely robust.
It stands for a whole universality class of phase transitions from a
ﬂuctuating phase into absorbing states. In fact, a large variety of models displays phase
transitions belonging to the DP universality class. Thus, on the theoretical level, DP plays
the role of a standard universality class similar to the Ising model in equilibrium statistical
mechanics.
In spite of its simplicity, no experiment is known which conﬁrms the values of the
critical exponents quantitatively. An exception may be the wetting experiment performed
by Buldyrev et al. (see below in Sec. 6.2), where the value of the roughness exponent α
coincides with ν⊥/ν∥within less than 10%. However, since the results of similar experi-
ments are scattered over a wide range, further experimental eﬀort in this direction would
be needed in order to conﬁrm the existence of DP in this type of systems.
Apart from diﬃculties to realize a non-ﬂuctuating absorbing state, a fundamental prob-
lem of DP experiments is the emergence of quenched disorder caused by inhomogeneities
of the system. Depending on the type of disorder, even weak inhomogeneities might ob-
scure or even destroy the DP transition. Therefore, the most promising experiments are
those where quenched disorder is irrelevant on large scales. This is the case, for example,
in wetting experiments and systems of ﬂowing granular matter.
Although there is certainly a lack of experimental evidence, there is no reason to
believe that DP is a purely artiﬁcial model. To be optimistic, it is helpful to recall the
history of the Ising model, which has been introduced almost one century ago. Although
the Ising model is probably the best studied system in equilibrium statistical mechanics,
there are only few experiments in which the critical exponents have been reproduced (for
a review see Ref. [10]). For this reason, many physicists believe that DP should have a
85

<!-- page 86 -->

counterpart in reality as well, mostly because of its simplicity and robustness. In this
respect, Grassbergers message remains valid: The experimental realization of DP is an
outstanding problem of top priority.
86

<!-- page 87 -->

4
Other classes of spreading transitions
This Section discusses various other types of nonequilibrium phase transitions into ab-
sorbing states which do not belong to the universality class of directed percolation. In
particular we will address spreading processes with long-range interactions and additional
symmetries.
4.1
Long-range spreading processes
According to the DP conjecture (see Section 3.3) phase transitions in spreading models
with short range interactions generically belong to the DP universality class. In many
realistic spreading processes, however, short-range interactions do not appropriately de-
scribe the underlying transport mechanism. This situation emerges, for example, if an
infectious disease is transported by insects. Typically the motion of the insects is not a
random walk, one rather observes occasional ﬂights over long distances before the next
infection occurs. Similar phenomena are expected when the spreading agent is subjected
to a turbulent ﬂow. Intuitively it is clear that occasional spreading over long distances will
signiﬁcantly alter the spreading properties. On a theoretical level such a super-diﬀusive
transport may well be described by L´evy ﬂights [104], i.e., by uncorrelated random moves
over long distances r which are algebraically distributed as
P(r) ∼1/rd+σ ,
(σ > 0) .
(207)
The exponent σ is a free parameter that controls the characteristic shape of the distribu-
tion. The algebraic tale leads to occasional long-distance ﬂights, as shown in Fig. 40.
Anomalous directed percolation, as originally proposed by Mollison [101] in the context
of epidemic spreading, is a generalization of DP in which the spreading agent is transported
by L´evy ﬂights. As in the case of ordinary DP, we expect anomalous DP to be characterized
by certain universal critical exponents β, ν⊥, and ν||. The question is how these exponents
depend on σ, whether they are independent from one another, and how they cross over to
the exponents of ordinary DP. Based on ﬁeld-theoretic considerations, Grassberger [282]
claimed that the critical exponents of anomalous DP should depend continuously on the
control exponent σ. Very recently this work has been considerably clariﬁed and extended
by Janssen et al. [283], who presented a comprehensive ﬁeld-theoretic analysis of anomalous
spreading processes with and without immunization.
Anomalous directed percolation: Field-theoretic predictions
In order to include long-range spreading in the Langevin equation (139), the Laplacian
has to be replaced by a non-local expression. This term can be written as an integral that
describes L´evy ﬂights over the distance r according to the probability distribution P(r):
∂tφ(x, t) = κ φ(x, t) −λ φ2(x, t) + ζ(x, t)
+ D
Z
ddx′ P(|x −x′|)
h
φ(x′, t) −φ(x, t)
i
.
(208)
The two contributions in the integrand describe gain and loss processes, respectively.
Keeping the most relevant terms in a small momentum expansion [283], this equation may
87

<!-- page 88 -->

Fig. 40:
Ordinary random walk (σ = ∞) in comparison with L´evy ﬂights (σ < ∞).
be written as
∂tφ(x, t) =

DN∇2 + DA∇σ
A + κ

φ(x, t) −λφ2(x, t) + ζ(x, t) ,
(209)
where the noise correlations are assumed to be the same as in Eq. (140). DN and DA
are the rates for normal and anomalous diﬀusion, respectively. The anomalous diﬀusion
operator ∇σ
A describes moves over long distances and is deﬁned through its action in
momentum space
∇σ
A eik·x = −kσ eik·x ,
(210)
where k = |k|. The standard diﬀusive term DN∇2 takes the short-range component of
the L´evy distribution into account. Note that even if this term were not initially included,
it would still be generated under renormalization of the theory. The mean-ﬁeld theory
of anomalous DP is completely analogous to that of ordinary DP. For σ < 2 a scaling
analysis yields the mean ﬁeld results
dc = 2σ ,
βMF = 1 ,
νMF
⊥
= 1/σ ,
νMF
∥
= 1 .
(211)
For σ ≥2 these exponents cross over smoothly to the ordinary DP mean-ﬁeld expo-
nents (124).
The mean ﬁeld approximation is expected to be quantitatively accurate
above the upper critical dimension dc, while for d ≤dc ﬂuctuation eﬀects have to be taken
into account. By using standard techniques one can derive the eﬀective action
S[φ, ˜φ] =
Z
ddx dt
h
˜φ(∂t −κ −DN∇2 −DA∇σ
A)φ + Γ
2 (˜φφ2 −˜φ2φ)
i
.
This expression diﬀers from the usual action of Reggeon ﬁeld theory (142)-(143) by the
addition of a term representing anomalous diﬀusion. Simple power counting on this action
conﬁrms that the upper critical dimension is dc = 2σ, below which ﬂuctuation eﬀects
become important.
The ﬁeld-theoretic RG calculation basically follows the same lines
as in the case of DP. The resulting critical exponents to one-loop order in d = 2σ −ǫ
dimensions are given by
β = 1 −2ǫ/7σ + O
 ǫ2
,
ν⊥= 1/σ + 2ǫ/7σ2 + O
 ǫ2
,
(212)
ν|| = 1 + ǫ/7σ + O
 ǫ2
.
88

<!-- page 89 -->

t
x
=0.5
σ
=1.0
σ
8
=2.0
σ
=
σ
Fig. 41:
Critical anomalous directed percolation in 1+1 dimensions for diﬀerent values of σ. The ﬁgure
shows typical clusters starting in the center of the lattice. The case σ = ∞corresponds to ordinary DP. As
σ decreases, spatial structures become more and more smeared out until in the mean ﬁeld limit σ = 1/2
the particles appear to be randomly distributed over the whole system. For small values of σ the system
quickly reaches the absorbing state due to extremely strong ﬁnite size eﬀects.
Moreover, it can be shown that the hyperscaling relation (120) holds for arbitrary values
of σ. Thus, to one-loop order, θ and δ are given by
θ = ǫ/7σ + O
 ǫ2
,
δ = 1 −3ǫ/7σ + O
 ǫ2
.
(213)
Finally, since DA will not be renormalized, one can prove the additional exact scaling
relation
ν|| −ν⊥(σ −d) −2β = 0 .
(214)
This equation implies that anomalous DP is described by two rather than three indepen-
dent critical exponents. Moreover, it has another surprising consequence. Assuming that
β, ν⊥and ν|| change continuously with σ and cross over smoothly, it predicts for ﬁxed d
the value σc where the system should cross over to ordinary DP. In order to compute σc we
simply have to insert the numerically known values of the DP exponents into Eq. (214).
Surprisingly one obtains σc = 2.0766(2) in one, σc ≃2.2 in two, and σc = 2 + ǫ/12 in
d = 4−ǫ spatial dimensions. Thus, the crossover takes place at σc > 2 which collides with
the intuitive argument that the anomalous diﬀusion operator ∇σ
A should only be relevant
if σ < 2. But, as pointed out in Ref. [283], this naive argument may be wrong in an
interacting theory where the critical behavior is determined by a nontrivial ﬁxed point of
a RG transformation. The ﬁeld-theoretic calculation rather predicts anomalous diﬀusion
to be relevant in the range 2 ≤σ < σc(d) for d < 4. The surprising conclusion would be
that
∇2
A̸ = ∇2
(215)
in certain interacting theories. Loosely speaking, the tendency to correlate particles in
local spots of activity makes a DP process more sensitive to long-range ﬂights. Therefore,
the relevance of L´evy ﬂights sets in earlier than in the case of simple diﬀusion.
A lattice model for anomalous directed percolation
The ﬁeld-theoretic predictions can be veriﬁed numerically by studying a lattice model for
anomalous DP that generalizes directed bond percolation [284]. The model is deﬁned on
a tilted square lattice and evolves by parallel updates. As usual, a binary variable si(t)
89

<!-- page 90 -->

0
0.5
1
1.5
2
2.5
σ
0
0.5
1
1.5
2
ν||
ν
T
β
Fig. 42:
Numerical estimates for the dynamic critical exponents β, and the derived exponents ν⊥and ν||
in comparison with the ﬁeld-theoretic predictions (solid lines) and the DP exponents (dot-dashed lines).
is attached to each lattice site i. si = 1 means that the site is active (infected) whereas
si = 0 denotes an inactive (healthy) site. The dynamic rules depend on two parameters,
namely the control exponent σ > 0 and the bond probability 0 ≤p ≤1. For a given
conﬁguration {si(t)} at time t, the next conﬁguration {si(t+1)} is constructed as follows.
First the new conﬁguration is initialized by setting si(t + 1) := 0. Then a loop over all
active sites i in the previous conﬁguration is executed. In the (1+1)-dimensional case this
loop consists of the following steps:
1. Generate two random numbers zL and zR from a ﬂat distribution between 0 and 1.
2. Deﬁne two real-valued spreading distances rL = z−1/σ
L
and rR = z−1/σ
R
, for spreading
to the left (L) and to the right (R). The corresponding integer spreading distances
dL and dR are deﬁned as the largest integer numbers that are smaller than rL and
rR, respectively.
3. Generate two further random numbers yL and yR drawn from a ﬂat distribution
between 0 and 1, and assign si+1−2dL(t + 1) := 1 if yL < p and si−1+2dR(t + 1) := 1
if yR < p. In ﬁnite systems the arithmetic operations in the indices are carried out
modulo N by assuming periodic boundary conditions, i.e. si ≡si±N.
This model includes two special cases. For σ →∞it reduces to ordinary directed bond
percolation (see Sec. 3.1). On the other hand, for σ →0 the interaction becomes totally
random. In this case the mean-ﬁeld approximation becomes exact with a transition taking
place at pc = 1/2. In between, the spreading properties of the model change drastically, as
demonstrated in Fig. 41. As can be veriﬁed easily, the assignment r = z−1/σ reproduces
the normalized probability distribution
P(r) =
 σ/r1+σ
if
r > 1 ,
0
if
r ≤1 .
(216)
90

<!-- page 91 -->

0.5
1
1.5
2
2.5
3
σ
0.4
0.5
0.6
0.7
0.8
0.9
1
1.1
α(σ)
direct measurement
extrapolated exponents
field−theoretical result
Fig. 43:
Anomalous annihilation process. The graph shows direct and extrapolated estimates for the
decay exponent α, as a function of σ. The solid line represents the exact result (neglecting logarithmic
corrections at σ = 1).
As usual, the distribution has a lower cutoﬀat rmin = 1, representing the lattice spacing.
Yet, in contrast to other models [285,286], no upper cutoﬀis introduced. In order to reduce
ﬁnite-size eﬀects, the target site is determined by assuming periodic boundary conditions,
i.e., the particle may ‘revolve’ several times around the system.
An interesting aspect of anomalous DP is the possibility to choose σ in such a way
that the critical dimension dc = 2σ approaches the actual physical dimension where the
simulations are performed. Even in one spatial dimension this allows the one-loop re-
sults (212) to be veriﬁed. For example, if σ = 1/2+µ, the critical dimension of the system
is dc = 1 + 2µ. Hence the exponents in a (1+1)-dimensional system change to ﬁrst order
in µ as
β
=
1 −8µ/7 + O
 µ2
,
ν⊥
=
2 −12µ/7 + O
 µ2
,
ν||
=
1 + 4µ/7 + O
 µ2
,
z
=
1/2 + 5µ/7 + O
 µ2
,
δ
=
1 −12µ/7 + O
 µ2
,
θ
=
4µ/7 + O
 µ2
.
(217)
In Fig. 42, the predicted initial slopes are indicated by solid lines. Clearly they are in fair
agreement with the numerical estimates, conﬁrming the ﬁeld-theoretic results of Eq. (212).
This is one of the rare cases where we can directly ‘see’ the ﬁeld-theoretic results in the
simulation data.
Anomalous annihilation process
The more simple case of anomalous pair annihilation A + A →∅with long-range hop-
ping [287] can be solved exactly by a similar ﬁeld-theoretic analysis. In the ordinary an-
nihilation process [44] with short-range interactions, the average particle density is known
to decay as in Eq. (37). The L´evy-ﬂight case may be described theoretically by inserting
an additional operator ∇σ
A into the ﬁeld-theoretic action for pair annihilation [44]. The
resulting action, which can also be derived systematically [284], reads
S[φ, ˜φ] =
Z
ddx dt
n
˜φ(∂t −DN∇2 −DA∇σ
A)φ + 2λ˜φφ2 + λ˜φ2φ2 −φ0 ˜φδ(t)
o
,
(218)
91

<!-- page 92 -->

where φ0 represents the initial (homogeneous) density at t = 0. An analysis of this action
follows very closely that of Ref. [44]. For σ < 2, power counting reveals the upper critical
dimension of the model to be dc = σ < 2. For d > dc mean-ﬁeld theory is expected to
be quantitatively accurate, with an asymptotic density decay ∼t−1. Below dc, however,
the renormalized reaction rate ﬂows to an order ǫ = σ −d ﬁxed point. The decay of the
density can therefore be predicted via dimensional arguments (see Fig. 43):
ρ(t) ∼



t−d/σ
for d < σ ,
t−1 ln t
for d = dc = σ ,
t−1
for d > σ .
(219)
4.2
Absorbing-state transitions in systems with additional symmetries
As outlined previously, directed percolation is the canonical universality class for phase
transitions into a single absorbing states. According to the DP-conjecture, we may there-
fore expect non-DP behavior to occur in systems with several absorbing states. However,
it is important to note that the existence of several absorbing states alone does not auto-
matically lead to non-DP behavior at the transition point. As we have seen in Sec. 3.8,
even models with inﬁnitely many absorbing states may still belong to the DP universal-
ity class. Non-DP critical behavior emerges only if there is a symmetry among diﬀerent
absorbing states.
The ﬁrst examples of non-DP transitions into absorbing states were discovered by
Grassberger et al. [288,289] who observed “a new type of kinetic critical phenomenon” in
certain one-dimensional stochastic cellular automata. The density exponent β ≃0.6(2) in
1+1 dimensions was found to diﬀer signiﬁcantly from the usual DP exponent β ≃0.277.
Partially because of the complicated dynamic rules of these models it took almost ten
years until the mechanism behind this type of non-DP behavior was clearly identiﬁed.
Up to now two universality classes of spreading transitions with non-DP behavior
have been found, namely the so-called parity-conserving class (PC) and Z2-symmetric
directed percolation (DP2). The PC class is represented most prominently by branching-
annihilating random walks with even number of oﬀspring (BAWE) [290], where the number
of particles is preserved modulo 2. The DP2 class, on the other hand, which is also referred
to as the directed Ising class, introduces two symmetric absorbing states.
As we will
see below, both universality classes coincide in one spatial dimension wherefore they are
usually considered to be identical. However, it is important to note that they diﬀer from
each other in higher dimensions.
One may speculate whether systems with a symmetry among several spreading agents
will also be able to display novel critical properties. However, such multi-color spreading
processes were found to belong to the DP class as well [291]. A ﬁeld-theoretic analysis
conﬁrms this observation [292] and predicts that the symmetry among the spreading agents
may be spontaneously broken.
The parity-conserving universality class
The parity-conserving universality class is represented most prominently by branching
annihilating walks with an even number of oﬀspring [290,293–297]. These non-conserved
92

<!-- page 93 -->

t
c
λ=λ
λ<λ
c
c
λ>λ
Fig. 44:
The parity-conserving universality class. Typical space-time trajectories of a branching annihi-
lating random walk with two oﬀspring below, at and above the critical point.
random walks are deﬁned by the reaction-diﬀusion scheme
A →
λ (n + 1)A ,
2A →
α Ø ,
(220)
where n = 2, 4, 6, . . . denotes the number of oﬀspring. As an essential feature, this process
conserves the number of particles modulo 2. Early numerical studies in 1+1 dimensions
assumed instantaneous on-site annihilation α = ∞. In this case the model displays a
continuous phase transition only for n ≥4 [293–295], while there is no such transition for
n = 2. Later Zhong and ben-Avraham demonstrated that a phase transition also emerges
in the case of two oﬀspring, provided that the annihilation rate α is ﬁnite [290]. Fig. 44
shows a typical cluster grown from a single seed. In contrast to DP, a PC process starting
with an odd number of particles cannot reach the empty state since at least one particle is
left. Therefore, initial states with even and odd number of particles are expected to lead
to diﬀerent cluster morphologies.
The relaxational properties of PC models in the subcritical phase diﬀer signiﬁcantly
from the standard DP behavior. While the particle density in DP models decays exponen-
tially as ρ(t) ∼e−t/ξ∥, in PC models it decays algebraically in the long-time limit. More
precisely, the temporal evolution of PC processes in the inactive phase is governed by the
annihilation process 2A →Ø. Under RG transformations we therefore expect the system
to ﬂow towards the ﬁxed point of particle annihilation. Consequently, in the subcritical
phase the particle density decays algebraically as in Eq. (37).
A systematic ﬁeld theory for PC models has been presented in Refs. [296, 297], con-
ﬁrming the existence of such an annihilation ﬁxed point.
However, the ﬁeld-theoretic
treatment of the (1+1)-dimensional case poses considerable diﬃculties. They stem from
the presence of two critical dimensions: dc = 2, above which mean-ﬁeld theory applies,
and d′
c ≈4/3, where for d > d′
c (d < d′
c) the branching process is relevant (irrelevant) at
the annihilation ﬁxed point. Therefore, the physically interesting spatial dimension d = 1
cannot be accessed by a controlled ǫ-expansion down from upper critical dimension dc = 2.
Nevertheless the usual scaling theory is still valid below d′
c. Currently the best numerical
estimates of the critical exponents in 1+1 dimensions are
β = 0.92(2) ,
δ + θ = 0.286(2) ,
ν∥= 3.22(6) ,
2/z = 1.15(1) ,
(221)
ν⊥= 1.83(3) .
The actual values of δ and θ in dynamic simulations depend on the initial condition. If
the process starts with a single particle, it will never stop, hence δ = 0. On the other
93

<!-- page 94 -->

Fig. 45:
Temporal evolution of the (1+1)-dimensional interacting monomer-dimer model of Ref. [299].
The two absorbing states emerge as alternating conﬁgurations with A-particles at even and odd sites,
respectively.
hand, if the initial seed consists of two particles, one observers that the roles of δ and θ
are exchanged, i.e. θ = 0.
It has been customary to investigate whether the numerical estimates of the critical
exponents can be ﬁtted by simple rational numbers. In fact, the estimates of Eq. (221)
are in good agreement with the rational values [295] β = 12/13, ν∥= 42/13, ν⊥= 24/13,
δ + θ = 2/7, and ˜z = 8/7. In particular, β/ν⊥, the exponent for the decay of spatial
correlations at criticality, should be equal to 1/2. In the past, however, rational values
were also proposed for the DP exponents and later disproved by more accurate numerical
estimates [298].
The DP2 universality class
The DP2 universality class describes phase transitions in spreading models with two sym-
metric absorbing states. Since the two absorbing states compete one another, the resulting
critical behavior is diﬀerent from ordinary directed percolation14. In various models, for
example in certain cellular automata [288, 289] as well as in interacting monomer-dimer
models [299–302], the two absorbing states emerge as checkerboard-like conﬁgurations of
particles at even or odd sites, respectively (see Fig. 45). Other DP2 models explicitly intro-
duce two symmetric inactive states, as, for example, nonequilibrium Ising models [303,304],
Z2-symmetric generalizations of the DK model and the contact process [263], monomer-
monomer surface reaction models with two absorbing states [305], and certain cellular
automata [306]. Even in monomer-monomer models with three absorbing states a DP2
transition emerges at certain points in the parameter space [307,308].
On a phenomenological level, spreading transitions with several absorbing states may
be deﬁned by introducing a single active state A and n symmetric absorbing states
I1, . . . , In which can be regarded as having diﬀerent colors. The system evolves according
to the following descriptive rules:
1. Spreading of activity: Active sites turn their inactive nearest neighbors into the
active state.
14Without symmetry, one of the two absorbing states will dominate so that the critical behavior crosses
over to DP after suﬃciently long time.
94

<!-- page 95 -->

t
λ>λc
λ=λc
λ<λc
Fig. 46:
Simulation of the generalized contact process with n = 2 absorbing states starting from a random
initial condition. The two diﬀerent types of inactive domains are shown in black and grey. The active sites
between the domains are represented by white pixels.
2. Spontaneous recovery: Active sites turn spontaneously into an inactive state of a
randomly chosen color.
3. Boundaries between inactive domains of diﬀerent colors are free to separate again,
leaving active sites in between.
Rules 1 and 2 resemble the usual infection and recovery processes of DP. Rule 3 is new and
distinguishes diﬀerent colors. Roughly speaking, this rule ensures that inactive domains
of diﬀerent colors cannot stick together irreversibly, rather they will always be separated
by ﬂuctuating active ‘interfaces’.
The symmetry under global permutation of the col-
ors ensures that absorbing domains of diﬀerent colors compete one another, leading to
interesting critical behavior.
Following these descriptive rules, we can introduce a generalized version of the Domany-
Kinzel cellular automaton (see Sec. 3.2) with n absorbing states [263]. It uses the same
type of lattice and is deﬁned by the conditional probabilities
P(Ik | Ik, Ik) = 1 ,
P(A|A, A) = 1 −n P(Ik|A, A) = p2 ,
P(A|Ik, A) = P(A|A, Ik) = p1 ,
P(Ik|Ik, A) = P(Ik|A, Ik) = 1 −p1 ,
P(A|Ik, Il) = 1 ,
(222)
where k, l = 1, . . . , n and k̸ = l. Notice that rule 3 is implemented by transition IkIl →A,
ensuring that active sites are created between two inactive domains of diﬀerent colors.
For n = 1 the above model reduces to the original Domany-Kinzel model. For n = 2 it
displays a DP2 phase transition at the critical threshold p1,c = p2,c = 0.5673(5). Similarly
it is possible to deﬁne a generalized contact process by the rates
wAA→AIk = wAA→IkA = 1/2n ,
wAIk→IkIk = wIkA→IkIk = 1/2 ,
wAIk→AA = wIkA→AA = λ/2 ,
wIkIl→IkA = wIkIl→AIl = λ/2 ,
(223)
95

<!-- page 96 -->

Here the last equation implements rule 3. For n = 1 this model reduces to the usual
contact process introduced in Sec. 3.2. For n = 2 the model undergoes a DP2 transition
at the critical point λc = 1.592(5). A typical evolution of the generalized contact process
in 1+1 dimensions is shown in Fig. 46. In the active phase λ > λc small inactive islands of
random color are generated which survive only for a short time. Approaching the phase
transition their size and lifetime grows while the density of active sites decreases. Notice
that according to rule 3 a thin ﬁlm of active sites separates diﬀerent inactive domains.
As expected, the numerical estimates of the critical exponents are in agreement with
the PC exponents (221). Both models can easily be generalized to higher dimensions.
However, in higher dimensions the phase transition is presumably characterized by mean
ﬁeld behavior. Similarly, increasing the number of absorbing state does not lead to new
universality classes. Simulations with n ≥3 symmetric absorbing states in 1+1 dimensions
indicate that the system is again described by mean ﬁeld exponents.
If the Z2 symmetry of DP2 models is broken by an external ﬁeld, the critical behavior
at the transition crosses over to ordinary DP [137]. Roughly speaking the external ﬁeld
favors one of the absorbing states so that pairs of kinks between oppositely oriented inactive
domains form bound ‘dipoles’ of a certain size. Interpreting these dipoles as composite
particles, they recombine and produce a single oﬀspring at certain rates, resembling an
ordinary DP process on large scales.
The diﬀerence between the PC and DP2 universality classes
In the DP2 class the ‘kinks’ between diﬀerently colored domains may be interpreted as
branching-annihilating particles with even number of oﬀspring. Although the number of
active sites is generally not conserved modulo 2, we may associate with each active island
between diﬀerently colored inactive domains a particle X. Obviously, these particles per-
form an eﬀective branching-annihilating random walk X →3X, 2X →Ø. Therefore, the
DP2 class and the PC class coincide in 1+1 dimensions. However, it is important to note
that they are diﬀerent in higher dimensions. Active sites of PC models in d ≥2 dimensions
can be considered as branching-annihilating walkers, whereas DP2 models describe the dy-
namics of branching-annihilating interfaces between oppositely oriented inactive domains
(see Fig. 47). Therefore, the corresponding ﬁeld theories are expected to be diﬀerent. A
ﬁeld theory for the PC class was presented in Ref. [296,297], whereas the development of
ﬁeld theories for branching-annihilating interfaces started only recently [240].
In order to understand the diﬀerence between PC and DP2, it is helpful to consider two
other universality classes which also coincide in 1+1 dimensions, namely the annihilation
process A + A →Ø and the voter model [252]. The voter model is a two-state model
with spins si = ±1. It evolves by random-sequential updates +−→+ + / −−with equal
rates. Interpreting +−kinks as particles A, the voter model and the annihilation process
coincide in 1+1 dimensions. However, in higher dimensions they are diﬀerent. In fact,
even their Langevin equations turn out to be diﬀerent. As shown in Sec. 2.6, the Langevin
equation of the annihilation process reads
∂tρ(x, t) = −λρ2(x, t) + D∇2ρ(x, t) + ζ(x, t) ,
(224)
⟨ζ(x, t)ζ(x′, t′)⟩= Γ ρ2(x, t) δd(x −x′) δ(t −t′) ,
where ρ(x, t) represents the coarse-grained density of A-particles. On the other hand, the
96

<!-- page 97 -->

Fig. 47:
Conﬁgurations of an almost critical PC (left) and a DP2 process (right) in d = 2 dimensions.
Langevin equation of the voter model [253] is given by
∂tρ(x, t) = D∇2ρ(x, t) + ζ(x, t)
(225)
⟨ζ(x, t)ζ(x′, t′)⟩= Γ ρ(x, t) [1 −ρ(x, t)] δd(x −x′)δ(t −t′) ,
where ρ(x, t) ∈[0, 1] represents the local orientation of the domain. Obviously, the two
equations stand for diﬀerent universality classes.
DP2 surface critical behavior
The inﬂuence of an absorbing wall in systems with a DP2 transition was studied in
Ref. [309].
Analyzing the generalized Domany-Kinzel model (222) with two absorbing
states in a semi-inﬁnite geometry, it turned out that absorbing and reﬂective boundary
conditions play complementary roles. Moreover, since β and β′ may be diﬀerent in the
DP2 class, one has to introduce two diﬀerent surface exponents βs and β′
s. For absorbing
boundary conditions the numerical estimates in a (1+1)-dimensional DP2 process are
βs = 1.34(2) ,
β′
s = 2.06(2) .
(226)
For reﬂecting boundary conditions the two values are simply exchanged. This property
is related to a duality transformation in parity-conserving processes [310].
It is quite
remarkable that the two surface exponents seem to obey the scaling relation
1
2(βs + β′
s) = ν∥−1 ,
(d = 1)
(227)
generalizing the conjecture βs = ν∥−1 in the case of DP.
4.3
Activated random walks
So far we considered nonequilibrium phase transitions where a parameter (e.g. the perco-
lation probability p) has to be tuned to criticality. Other systems with conserved dynamics
can be tuned to criticality by varying the particle density in the initial state. An inter-
esting example is the activated random walk of n particles [238]. In this model each site
can be occupied with arbitrarily many particles. Sites with at least two particles are ac-
tive, i.e., their particles may move independently to randomly selected neighbors. Sites
with only one particle are frozen. On an inﬁnite lattice, this model has inﬁnitely many
97

<!-- page 98 -->

absorbing states. The control parameter is the density of particles. For a low density, the
model quickly evolves into one of the absorbing states, whereas it remains active for high
particle densities. Near the transition, the stationary density of active sites ρstat scales as
ρstat ∼(ζ −ζc)β ,
ζ = n/N ,
(228)
with the critical point ζc ≃0.9486 in 1+1 dimensions [238]. However, unlike other models
with inﬁnitely many absorbing state (see Sec. 3.8) the transition does not belong to the
DP universality class. For example, in 1+1 dimensions the measured exponent β = 0.43(1)
diﬀers signiﬁcantly from the expected DP value βDP ≃0.277. Obviously, this deviation
is due to the conservation law. Hence, activated random walks provide a new universality
class of phase transitions into absorbing states, caused by an additional symmetry, namely
particle conservation.
4.4
Absorbing phase transitions and self-organized criticality
In contrast to ordinary transitions into absorbing states, self-organized critical systems
exhibit long-range correlations and power laws without being tuned to a certain critical
point.
In some sense the annihilation process A + A →Ø (see Sec. 2.4) can be con-
sidered as a simple example of self-organized criticality. Without tuning of a parameter
the annihilation process generates long-range correlations with power-law characteristics.
Like many other coarsening processes and growth phenomena the driving mechanism is a
stationary state which is approached but never reached.
The term ‘self-organized criticality’ (SOC) refers to a diﬀerent type of models that
are attracted to a stationary critical state without being tuned to a critical point. The
chief examples are the sandpile model [36] and the Bak-Sneppen model [311]. The concept
of SOC has been used to explain the large variety of power laws observed in nature, for
example 1/f noise, the distribution of earthquakes, and the dynamics of ﬁnancial markets,
to name only a few [312]. For more than one decade SOC was considered as a quite separate
ﬁeld of theoretical statistical physics, being more or less unrelated to conventional phase
transitions with a tuning parameter. Recently, it was pointed out [313–315] that SOC
is in fact closely related to ordinary phase transitions into (inﬁnitely many) absorbing
states (for a survey see [238]).
More precisely, two classes of SOC models have to be
distinguished. The ﬁrst class of SOC models, exempliﬁed by the Bak-Sneppen model,
employs extremal dynamics. In this case the site with an extremal value is selected for
the next update, i.e., the dynamic rules provide a mechanism of global supervision. In the
second class, which is represented by the sandpile model, the bulk dynamics is conserved.
Here a slow driving force competes with the loss of particles at the systems boundaries
and drives the system to criticality. Hence the process of self-organization is characterized
by a separation of time scales for avalanches and driving.
In the latter case, SOC is related to a conventional absorbing-state transition as follows.
As explained in Ref. [238], any system with conserved local dynamics and a continuous
absorbing-states transition can be converted into a SOC model by (1) adding a process for
increasing the density in inﬁnitesimal steps, and (2) implementing a process for decreasing
the density inﬁnitesimally while the system is active. For example, let us consider the
activated random walk. Adding a process for random deposition of particles and a process
for loss of particles during avalanches at the boundaries, we obtain the so-called Manna
98

<!-- page 99 -->

sandpile model [316] which is known to exhibit SOC. Using a deterministic variant of the
same model, one obtains the famous Bak-Tang-Wiesenfeld model [36]. However, in order
to create a SOC counterpart for directed percolation, a slightly diﬀerent recipe has to be
used, as demonstrated in Ref. [314,317].
4.5
The annilation/ﬁssion process
As shown in Sec. 2.6, the inﬂuence of ﬂuctuations in (1+1)-dimensional systems can be
very diﬀerent. In the annihilation process, for example, the particles become anticorre-
lated, i.e., they try to be far away from each other. As shown in Ref. [43], this type of
ﬂuctuations corresponds to ‘imaginary’ noise in the Langevin equation. In a DP process,
however, the particles are highly correlated and the noise turns out to be real. Three
years ago Howard and T¨auber posed the question whether it is possible to interpolate
between real and imaginary noise. As a prototype of such a transition, they considered
the annihilation/ﬁssion process
2A →3A ,
2A →Ø
(229)
with diﬀusion of single particles. Note that this is a binary process, i.e., at least two par-
ticles are required to meet at the same (or at neighboring) sites in order to self-destruct
or to create oﬀspring. Moreover, there are no exceptional symmetries such as parity con-
servation on the microscopic level. The model exhibits a nonequilibrium phase transition
between an active phase and two non-symmetric absorbing states, namely the empty lat-
tice and the state with only one diﬀusing particle. Performing a ﬁeld-theoretic Howard
and T¨auber argued that this transition should belong to an independent yet unknown
universality class of phase transitions.
The annihilation/ﬁssion process may be interpreted as a pair contact process plus
diﬀusion of single particles. As a consequence, the static background shown in Fig. 34
begins to ﬂuctuate.
Since the decay of the particle density in the subcritical phase is
governed by the annihilation process 2A →Ø, it is natural to expect that the transition
does not belong to the DP class. Similarly, the possibility of DP2 critical behavior seems
to be unlikely since there is no parity conservation or Z2-symmetry in the system [163,318].
Thus, the annihilation/ﬁssion process may represent a new type of nonequilibrium critical
phenomenon which has not yet been studied before.
4.6
The Lipowski model
A particularly interesting model with a non-DP transition in d = 2 spatial dimensions has
been introduced recently by Lipowski [319]. It has inﬁnitely many absorbing states and is
deﬁned by extremely simple dynamic rules.
The model is deﬁned as follows. The state of the system is speciﬁed by real numbers
yij ∈(0, 1) which reside on the bonds of a d-dimensional square lattice (see Fig. 48).
Initially all these variables are randomly distributed between 0 and 1. The model evolves
by random sequential updates depending on a control parameter p which plays the role of
a percolation probability. For each update attempt a site i is selected at random. Its local
‘ﬁeld’ hi is deﬁned by the sum of the y-variables of the four connecting bonds. If hi > p
the four variables are replaced by random numbers drawn from a ﬂat distribution between
99

<!-- page 100 -->

Fig. 48:
The Lipowski model. Each bond of a square lattice carries a weight between zero and one. A
plaquette of four bonds (hollow circles) is randomly selected. If the sum of the four weights exceeds a
certain threshold, they are replaced by four independent random numbers between zero and one.
0 and 1. Otherwise the update attempt is abandoned. As usual in models with random-
sequential dynamics, each update attempt corresponds to a time increment of ∆t = 1/N,
where N is the total number of sites. A site is called ‘active’ if it is susceptible for a
replacement, i.e. hi > p. Measuring the density of active sites in a numerical simulation,
the model displays a continuous phase transition between a ﬂuctuating active and a frozen
phase. Clearly, the model has inﬁnitely many absorbing states.
The Lipowski model poses a puzzle: In one dimension the static exponent β coincides
with the DP exponent β ≃0.276. In two dimensions, however, the static exponent β is
found to be much smaller than the expected value 0.58. Surprisingly, the measured value
seems to coincide with the DP value in one dimension. Lipowski argued that the model
should provide a mechanism for dimensional reduction, placing the two-dimensional critical
phenomenon into a one-dimensional universality class. However, the observed coincidence
may well be accidental. Moreover, it is not obvious how such a mechanism should work. It
is therefore an interesting open question whether dimensional reduction can be observed
in stochastic lattice models.
100

<!-- page 101 -->

5
Damage spreading
One of the central problems of dynamic system theory is the dependence of the system’s
temporal evolution on the initial conditions. It is well known that nonlinear systems with
deterministic dynamics may be extremely sensitive to small perturbations of the initial
state. Even in simple systems such as in a periodically driven pendulum, a small variation
of the initial parameters can change the entire trajectory completely.
If the distance
between two inﬁnitesimally diﬀerent trajectories diverges during the temporal evolution,
a dynamic system is said to exhibit chaotic behavior.
The notion of chaos has been introduced in the context of deterministic systems where
a trajectory is uniquely determined by the initial condition. In random processes, however,
trajectories are not uniquely determined; the time evolution of such a system is not repro-
ducible and therefore the usual deﬁnition of chaotic behavior does not apply. Nevertheless,
one may pose the question how the system responds to changes in the initial condition.
It would be interesting to know, for example, how the biological evolution on earth would
have been aﬀected if the initial conditions were slightly diﬀerent.
In order to address
this question, Kauﬀman introduced the concept of damage spreading (DS) [320, 321]. In
damage spreading simulations [322] two copies (replicas) S, S′ of a stochastic model evolve
simultaneously. Initially the two copies diﬀer only on a small number of sites. This diﬀer-
ence is considered as a small perturbation (damage) in one of the two systems. Moreover,
it is assumed that the replicas evolve under identical realizations of thermal noise, i.e.,
both copies use the same sequence of random numbers in the simulation. If the number
of sites in diﬀerent states, the so-called Hamming distance
∆(t) =
X
i
∆i(t) =
X
i
1 −δsi(t),s′
i(t)
(230)
does not go to zero in the long-time limit, damage is said to spread, indicating high
sensitivity with respect to the initial condition. Otherwise, if ∆(t) vanishes, damage is
said to heal, indicating a weak inﬂuence of the initial condition. In order to get statistically
meaningful results, the Hamming distance has to be averaged over many realizations of
randomness.
Damage spreading ﬁrst appeared in physics literature in the mid eighties [323–327]
and attracted considerable interest and attention. The main reason behind this initial
enthusiasm was the hope that damage may spread in some regions of a system’s parame-
ter space and disappear elsewhere, indicating the existence of chaotic and regular phases
in stochastic systems. The initial enthusiasm abated during subsequent years, the main
reason being an apparent lack of an objective measure whether damage does or does not
spread in a given system. More precisely, it was realized that the location of the phase
boundaries may depend on details of the algorithmic implementation. However, if spread-
ing or healing of damage indicated some intrinsic property of the system, one would not
expect the result to depend on details of the algorithm used to generate its dynamics.
Meanwhile DS has been applied to a large variety of models (see Table 3). In view of
the vast literature on DS simulations, it is therefore necessary to carefully analyze the
conceptual problems of this technique.
In the following we discuss several aspects of DS. First we present a simple example in
order to explain how a DS simulation depends on the algorithmic implementation. From a
more mathematical point of view, this phenomenon can also be understood by analyzing
101

<!-- page 102 -->

Equilibrium models:
Ref.
Ising ferromagnet
[325,326,328–331]
Heisenberg model
[332]
XY model
[333]
Potts and Ashkin-Teller models
[334–336]
discrete N-vector ferromagnets
[337]
spin glasses
[338–345]
Ising model with microcanonical constraints
[346]
diluted Ising model
[347,348]
Ising model in quenched random ﬁelds
[349]
Frustrated Ising systems
[350]
Layered Ising systems
[351]
Nonequilibrium models:
Kauﬀman model
[320,321,352,353]
The game of life
[354,355]
ZGB model
[356]
Ohta-Jasnov-Kawasaki and Ginzburg-Landau model
[357]
Irreversible reaction-diﬀusion processes
[358,359]
Models for surface reactions
[360]
Restricted solid on solid growth models
[361]
Models of the PC class
[362]
SOC models
[363]
Traveling salesman problem
[364]
Boolean random networks
[365]
Table 3:
Some applications of damage spreading.
the joint master equation of the two replicas. Because of their algorithmic dependence,
DS simulations are ambiguous and cannot be used as a criterion for ‘chaotic’ and regular
phases. However, to some extent the ambiguity of DS can be overcome by an algorithm-
independent deﬁnition of DS phases.
In this approach the entire family of physically
legitimate algorithms for a given dynamic system is considered as a whole. Furthermore,
we summarize what is known about the universal properties of DS transitions. Finally we
discuss several applications of DS simulations.
5.1
Damage spreading phases
Algorithmic dependence of damage spreading simulations
The conceptual problem of DS was ﬁrst discovered in the Domany-Kinzel cellular automa-
ton (cf. Fig. 14). Martins et al. [366] observed that in a certain region of the active phase
damage spreads and heals elsewhere. Subsequently several other authors determined the
boundary of this region with increasing accuracy [114, 367–369]. Independently, mean-
ﬁeld type approximations of varying complexity conﬁrmed the existence of this ‘chaotic
phase’ [368,370–372]. Its boundary, however, was shown to depend on the manner in which
the dynamic procedure of the DK model is implemented on a computer [370,371], while the
evolution of a single replica is completely insensitive to the algorithmic implementation.
This prompted Grassberger [369] to observe that
“it is misleading to speak of diﬀerent phases in the DK automaton, ...instead
these are diﬀerent phases for very speciﬁc algorithms for simulating pairs of
such automata”.
102

<!-- page 103 -->

equilibrium
dynamic
rule
algorithmic
implementation
model
...
Ising model
Potts model
Ashkin-Teller
N-vector
heat bath
Metropolis
Kawasaki
Q2R
Wolf
...
...
spin-orienting (heat bath)
spin-flipping (Glauber)
mixtures of Glauber / heat bath
Fig. 49:
The thermodynamic ensemble of an equilibrium model (left) may be generated by diﬀerent
dynamic procedures (center) that are characterized by speciﬁc transition probabilities.
Each of these
dynamic procedures in turn can be realized by diﬀerent algorithms (right). These algorithms are fully
equivalent and cannot be distinguished by an observer of a single system. However, their damage spreading
properties turn out to be diﬀerent (see text).
To understand the algorithmic dependence of DS simulations, let us consider a much
simpler system, namely the Ising model.
In this context it is important to note that
there are diﬀerent levels of variety in stochastic lattice models (see Fig. 49). As explained
in Sec. 2.10, the equilibrium ensemble of the Ising model can be generated by various
diﬀerent dynamic rules. For example, heat bath and Metropolis dynamics represent two
diﬀerent dynamic systems which have the same stationary state. Initially it was hoped
that DS would not depend on the intrinsic dynamics, allowing regular and chaotic phases
to be identiﬁed as equilibrium properties [325–327]. However, later it was realized that
diﬀerent dynamic procedures (such as Glauber versus Metropolis [328,331], Q2R [329], or
Kawasaki dynamics [330]) exhibit diﬀerent DS properties. Yet, this is not surprising since
the dynamic rules, although generating the same equilibrium ensemble, represent diﬀerent
dynamic systems. Since DS is a dynamic phenomenon it is quite natural that it depends
on the dynamics under consideration. Similarly, it is not surprising that DS depends on
the type of updates [373,374] and the interaction range [375].
The real conceptual problem of DS occurs at the second level in Fig. 49: On a com-
puter each dynamic rule can be realized by several algorithms. Regarding a single system,
these algorithms are fully equivalent and cannot be distinguished. However, in DS simula-
tions they lead to diﬀerent results. To understand this apparent paradox, let us consider
the Ising model with spin-orienting (standard heat bath) and spin-ﬂipping (Glauber) dy-
namics. Both algorithmic prescriptions represent the same dynamic rule which mimics
the contact of the Ising model with a thermal reservoir by means of local spin dynamics.
More precisely, an observer of a single system who analyzes the trajectories would be
unable to decide whether its temporal evolution was generated by standard heat bath or
Glauber dynamics. This can already be seen in the example of a single-spin Ising model
at inﬁnite temperature. Since in this case the spin σ(t) = ±1 changes randomly in time,
the transition rates are simply given by w−1→+1 = w+1→−1 = 1. These transition proba-
bilities deﬁne our dynamic system. However, this system can be realized by two diﬀerent
algorithms, namely by spin-orienting updates (standard heat bath dynamics)
z=rnd(0,1);
if (z<0.5) σ(t + dt)=1; else σ(t + dt)=-1;
and by spin-ﬂipping updates (Glauber dynamics)
103

<!-- page 104 -->

z=rnd(0,1);
if (z<0.5) σ(t + dt)=−σ(t); else σ(t + dt)=σ(t)
.
It is obvious that both procedures are fully equivalent on a single replica, i.e., an observer
would be unable to decide which of the two algorithms has been used. The diﬀerence
between the two algorithms may become evident only if we observe the evolution of two
replicas S and S′ of the system in a DS simulation: For spin-orienting updates an initial
‘damage’ σ(0) = −σ(0)′ heals immediately while it is preserved when spin-ﬂipping updates
are used.
A similar algorithmic dependence can be observed in the full Ising model (58) at ﬁnite
temperature. Deﬁning the transition probability (cf. Eq. (59))
pi(t) =
ehi(t)
ehi(t) + e−hi(t)
(231)
we may express the update rules of heat bath dynamics by
σi(t + 1) = sign[pi(t) −zi(t)] ,
(232)
where zi(t) ∈(0, 1) are equally distributed random numbers. The corresponding update
rule for Glauber dynamics is given by
σi(t + 1) =
 +sign[pi(t) −z]
if σi(t) = +1
−sign[1 −pi(t) −z]
if σi(t) = −1 .
(233)
It is easy to verify that for given {σi−1(t), σi(t), σi+1(t)}, the probability to get σi(t+1) =
+1 is the same in both cases. Hence, by observing the temporal evolution of a single Ising
system, one cannot tell which of the two methods was used to generate the trajectory in
conﬁguration space. The two algorithms can only be distinguished when two copies are
simulated in parallel, i.e., by studying damage spreading.
Investigating the two-dimensional Ising model with Glauber dynamics, Stanley et
al. [326] and Mariz et al. [328] found that damage spreads in the disordered phase T > Tc
and heals elsewhere. Performing more precise simulations, Grassberger [376] realized that
the DS transition occurs slightly below Tc. This observation was also be supported by a
mean ﬁeld theory [377]. Moreover, the critical point of the DS transition was found to vary
continuously when mixtures of Glauber and heat bath dynamics are used [378]. Very simi-
lar properties are observed in the three-dimensional Glauber model [376,379–381]. Turning
to the Ising model with spin orienting (standard heat bath) dynamics, it is possible to
prove that damage does not spread at any temperature in any dimension [327,328].
The master equation of damage spreading simulations
In order to understand the algorithmic dependence of DS simulations from a more funda-
mental point of view, let us consider two copies S1 and S2 of an arbitrary nonequilibrium
system with asynchronous dynamics (random-sequential updates). Each of the two sys-
tems evolves according to a master equation with a given Liouville operator L. The total
system S = (S1, S2) constitutes a new nonequilibrium system that evolves according to
a joint master equation with a certain Liouville operator M. If both systems were using
independent random numbers, M would be given by the tensor product L ⊗L. However,
in DS simulations the use of the same sequence of random numbers leads to nontrivial
correlations between S1 and S2.
104

<!-- page 105 -->

According to the deﬁnition of damage spreading, the Liouville operator M of the total
system (S1, S2) is restricted by certain physical constraints. On the one hand, each replica
is required to evolve according to its own natural dynamics. Hence, by integrating out
the degrees of freedom of one of the replicas we obtain the Liouville operator of the other
system:
X
s1
⟨s1, s2|M|s′
1, s′
2⟩= ⟨s2|L|s′
2⟩
for all s2, s′
1, s′
2 ,
(234)
X
s2
⟨s1, s2|M|s′
1, s′
2⟩= ⟨s1|L|s′
1⟩
for all s1, s′
1, s′
2 .
(235)
These restrictions already imply probability conservation for the total system. On the
other hand, the trajectories of the two replicas, once they have reached the same state (no
damage), have to be identical:
⟨s1, s2|M|s′, s′⟩= ⟨s1|L|s′⟩δs1,s2 .
(236)
The ambiguity of damage spreading simulations is due to the fact that the restrictions
(234)-(236) do not fully determine the Liouville operator M of the total system. This can
easily be veriﬁed by counting the degrees of freedom. For a system with n conﬁgurations
the restrictions determine less than 3n3 of the n4 matrix elements of M so that damage
spreading is ambiguous for n ≥3. But even in the case of a single-spin Ising model with
n = 2 states one can show that only 14 of 24 equations are linearly independent so that two
out of 16 matrix elements of M can be chosen freely. The remaining degrees of freedom are
the origin of the algorithmic dependence of DS simulations. A similar ambiguity occurs
in the joint transfer matrix of models with parallel dynamics [382].
An algorithm-independent deﬁnition of damage spreading
In order to overcome these conceptual diﬃculties, let us consider the entire family of
dynamic procedures consistent with certain physically dictated constraints [382]. Then
for any particular system one of the three possibilities may hold:
1. Damage spreads for any member of the family of dynamic procedures.
2. Damage heals for any member of this family.
3. Damage spreads for a subset of the possible dynamic procedures, and heals for the
complementing subset.
Obviously this allows us to classify damage spreading properties in an algorithm-indepen-
dent manner. To this end we must, however, consider simultaneously the entire set of
possible algorithms (dynamic procedures) which are consistent with the physics of the
model under consideration. This set of physically ‘legitimate’ algorithms may be deﬁned
by certain restrictions that are dictated by the dynamics of the single evolving system:
105

<!-- page 106 -->

0.5
0.6
0.7
0.8
0.9
1
0
0.2
0.4
0.6
0.8
1
Inactive
2
3
1
phase
Fig. 50:
Algorithm-independent damage spreading phases in the Domany-Kinzel model. 1: Damage never
spreads. 2: Damage may spread, depending on the algorithm. 3: Damage always spreads.
a) Deﬁnition of DS: The dynamic rules for the evolution of the pair of replicas are such
that a single replica evolves according to its ‘natural’ dynamics. Once both replicas
have reached the same conﬁguration, their temporal evolution will be identical.
b) Interaction range: The transition probabilities for the combined system at site i may
depend only on those sites that aﬀect the evolution of site i under the dynamic rules
of a single system.
c) Symmetry: The rules that govern the evolution for the pair of systems do not break
any of the symmetries of the single-replica dynamics.
The ﬁrst restriction is simply a verbal formulation of Eqs. (234)-(236). The second con-
dition tells us that the interaction range in the combined system of two replicas must
not exceed the interaction range of a single system. The third rule implies, for example,
that if there is a left-right symmetry in the evolution of a single system, the same must
hold for the pair of replicas. It can be shown that these conditions suﬃce to unambigu-
ously determine a set of physically ‘legitimate’ algorithms. This was demonstrated for the
one-dimensional Domany-Kinzel model [382] for which the phase diagram in Fig. 50 was
found.
Clearly, the subjectivity in deﬁning DS phases, as described before, has now been
shifted to selecting the restrictions deﬁning which DS procedure is ‘legitimate’. However,
the speciﬁcation of such a family by means of physically motivated criteria appears to be
less arbitrary than choosing, at random, one out of a continuum of physically equivalent
update procedures. It should also be emphasized that the algorithm-independent deﬁnition
of DS phases does not mean that DS is reﬂected in the dynamic behavior of a single system,
so that Grassberger’s observation still holds: DS is a property of a pair of stochastic
systems.
106

<!-- page 107 -->

damage
replica 2
t
replica 1
Fig. 51:
Damage spreading in the (1+1)-dimensional Domany-Kinzel cellular automaton evolving in the
stationary active state at p1 = 0, p2 = 0.313. An initial damage is introduced at the lattice site in the
center.
5.2
Universality of damage spreading transitions
The DP conjecture for damage spreading
As can be seen in Fig. 51, spreading of damage is in many respects similar to spreading
of activity in a DP process. As in DP damage spreads to nearest neighbors and heals
spontaneously. Once both copies have reached the same conﬁguration (no damage), their
evolution will be identical, i.e., they are conﬁned to some ‘absorbing’ subspace from where
they cannot escape. In contrast to DP the spreading process depends crucially on the
actual evolution of the two replicas, providing a ﬂuctuating background in which the
spreading process takes place. Nevertheless DS follows the same spirit as the spreading of
activity in a DP model. This observation led Grassberger to the conjecture that damage
spreading transitions belong generically to the directed percolation universality class [369].
So far analytical support for this conjecture came from approximate mean-ﬁeld argu-
ments [371] and an exact statement by Kohring and Schreckenberg [370], who noted that
on the p2 = 0 line the dynamics of damage spreading in the DK automaton is precisely
identical to the evolution of the DK automaton itself, hence on this line DS is trivially in
the DP universality class. To prove this statement, consider two replicas of S and S′ of a
DK automaton evolving according to Eq. (83) with equal random numbers zi(t) = z′
i(t).
The probabilities PD(∆i = 1|si−1, si+1; s′
i−1, s′
i+1) to generate damage at site i are listed
in Table 4. For p2 = 0 they may be expressed as
PD(∆i = 1|si−1, si+1; s′
i−1, s′
i+1) =



p1 if si−1 = si+1 and s′
i−1̸ = s′
i+1 ,
p1 if si−1̸ = si+1 and s′
i−1 = s′
i+1 ,
0
otherwise
(237)
or equivalently
PD(∆i = 1|si−1, si+1; s′
i−1, s′
i+1) = p1(1 −δ∆i−1,∆i+1) .
(238)
Therefore, the damage variables ∆i(t) evolves precisely according to the probabilistic
rules of a single DK automaton. Hence for p2 = 0 the DS transition belongs to the DP
universality class. This mapping of DS to DP was later extended to other regions in the
phase diagram of the DK model [382]. Although such an exact mapping is usually not
available for other models, various numerical simulations show that most DS transitions
107

<!-- page 108 -->

σ′
i−1, σ′
i+1
σi−1, σi+1
00
01
10
11
00
0
p1
p1
p2
01
p1
0
0
max(p1, p2)
10
p1
0
0
max(p1, p2)
11
p2
max(p1, p2)
max(p1, p2)
0
Table 4:
Probabilities PD(∆i = 1|si−1, si+1; s′
i−1, s′
i+1) for the generation of damage in the DK model.
are indeed characterized by DP exponents, supporting Grassberger’s conjecture.
The
same applies to deterministic cellular automata with chaotic behavior if a small noise is
added [383].
Diﬀerent DS properties may be expected in models with cluster dynamics [384]. How-
ever, as pointed out in Refs. [385, 386], it is diﬃcult to extend the deﬁnition of ‘using
the same random numbers’ to cluster algorithms since the number of clusters may diﬀer
on both replicas. This diﬃculty can be overcome by introducing a random background
ﬁeld [387]. Assigning a local random number to each lattice site the new orientation of a
cluster depends on whether the sum of its random numbers positive or negative. Although
cluster algorithms involve long-range correlations, DS transitions still belong to the DP
universality class unless they coincide with the thermodynamic phase transition of the
Ising system.
DS transitions with non-DP behavior
The critical properties of a DS transition are expected to change if one of the four con-
ditions of the DP conjecture is violated.
For example, DS transitions in models with
frozen randomness (such as the Kauﬀman model [320,321]) do not belong to the DP class.
Non-DP behavior is also expected when the DS order parameter exhibits an additional
Z2 symmetry (cf. Sec. 4.2). In the context of damage spreading it is important to note
that such a symmetry should be a property of the DS order parameter, i.e., the Hamming
distance. For example, the Z2 symmetry of Ising systems is not suﬃcient – inverting all
spins in both replicas does not change the Hamming distance between the two conﬁgura-
tions. Similarly, models with a non-DP transition do not automatically exhibit non-DP
damage spreading. For example, Grassbergers cellular automaton A, which has a DP2
phase transition, displays an ordinary DS transition belonging to DP [362].
Non-DP behavior at the DS transition can be observed in systems with a symmetry
between two ‘absorbing states’ of ∆i(t), one without damage ∆= 0 and the other with
full damage ∆= 1. The simplest example of such a symmetry is the Ising model with
Glauber dynamics. In fact, if both replicas are in opposite states (full damage), they will
always evolve through opposite conﬁgurations. Unfortunately, there is no DS transition
in the one-dimensional Glauber model. However, by exploiting the algorithmic freedom,
it is possible to construct a modiﬁed dynamic rule which exhibits a DP transition [388].
This rule depends on a parameter λ and is deﬁned by
σi(t + 1) =
 +sign(pi(t) −z)
if y = 1
−sign(1 −pi(t) −z)
if y = −1
,
(239)
y = 1
2σi
h
1 + σi−1(t)σi+1(t)

+

1 −σi−1(t)σi+1(t)

sign(λ −˜z)
i
,
108

<!-- page 109 -->

where z, ˜z ∈(0, 1) are two independent random numbers. On a single replica this update
rule is fully equivalent to Glauber and heat bath dynamics for all 0 ≤λ ≤1. However, the
eﬀective rate for spreading of damage depends on λ. For ﬁxed temperature J/kBT = 0.25
a DS transition with non-DP exponents occurs at the critical value λ∗= 0.82(1). This
example demonstrates that additional symmetries of the DS order parameter may lead to
non-DP behavior at the DS transition. A similar situation emerges in the nonequilibrium
Ising model introduced by Menyh´ard [362].
5.3
Applications of damage spreading
Measurement of critical exponents in equilibrium models
Damage spreading simulations can also be used to determine certain static and dynamic
properties of systems at thermal equilibrium. This application of DS was ﬁrst demon-
strated by Coniglio et al. [389] who showed that there exists an exact relation between the
Hamming distance ∆and certain correlation functions. The essential idea is to consider
two copies S, S′ of an Ising model with spins σi, σ′
i = ±1 and introducing a small damage
by keeping a single spin in one of the systems ﬁxed during the whole temporal evolution,
say σ′
0 = −1. Since this perturbation breaks the symmetry between the two copies, it is
important to distinguish two diﬀerent types of damage at site i, namely σi = 1, σ′
i = −1
and σi = −1, σ′
i = 1. The probabilities of ﬁnding these diﬀerent types of damage in the
stationary state is given by
d+−= ⟨(1 −si)s′
i⟩,
d−+ = ⟨si(1 −si)′⟩,
(240)
where si = 1
2(σi −1) (see also Ref. [322]). Notice that these quantities can be understood
as two-point correlation functions between the two copies. However, their diﬀerence
Γi = d+−
i
−d−+
i
= ⟨si⟩−⟨s′
i⟩
(241)
is a combination of one-point functions, i.e., by taking a certain combination of the dam-
age probabilities one obtains quantities that describe the properties of a single system.
Therefore, such quantities do not depend on the algorithmic implementation.
In fact,
using detailed balance and ergodicity one can prove that
Γi =
C0,i
2(1 −m) ,
(242)
where C0,i is the two-point correlation function and m the magnetization of the Ising
model at thermal equilibrium. This relation is exact and does not depend on the speciﬁc
algorithmic implementation used in the simulation. Moreover, it can be shown by mono-
tonicity arguments that for standard heat bath dynamics the total damage ∆is related
to the static susceptibility χ by
χ = 2(1 −m)
X
i
Γi .
(243)
Thus, damage spreading simulations can be used to determine the static exponents β and
ν and the critical temperature of the Ising model (see Refs. [390–392]). Similar relations
between Hamming distance and correlation functions were found in certain lattice models
with absorbing states [393].
109

<!-- page 110 -->

Fig. 52:
Identiﬁcation of domain walls. Left: The two-dimensional Ising model evolving by heat bath
dynamics. Right: Corresponding domain walls detected by the observable ∆i deﬁned in Eq. (244). Figure
reprinted from Ref. [219].
DS is also used as a tool for accurate measurements of the dynamic exponent z at the
phase transition of equilibrium systems, provided that the thermodynamic transition and
the DS transition coincide [394,395]. In this case the spreading exponent z is not given by
the DP exponent ν∥/ν⊥, instead it is expected to coincide with the dynamic exponent z of
the model under consideration. The knowledge of z is important in order to estimate the
critical slowing down of given dynamic system. In a series of papers the dynamic exponent
of the Ising model with heat bath dynamics has been measured in two [396–400] and three
dimensions [397, 399, 401, 402]. After an initial controversy, Grassberger [403] compared
several simulation methods and found the estimates z = 2.172(6) in two and z = 2.032(4)
in three dimensions. Later reﬁned simulations conﬁrmed these results [404].
Identiﬁcation of domain walls in coarsening systems
Damage spreading techniques can also be used to identify domain walls in coarsening
systems. Coarsening phenomena [405] occur in various dynamic systems as, for example,
in the Ising model with Glauber dynamics. In the ordering phase, starting with random
initial conditions, patterns of ferromagnetic domains are formed whose typical size grows
with time as t1/2. For zero temperature, these domains are completely ordered and the
domain walls can be identiﬁed as bonds between oppositely oriented spins. For nonzero
temperature, however, it is diﬃcult to deﬁne domain walls as one has to distinguish be-
tween ‘true’ domains and islands of the minority phase generated by thermal ﬂuctuations.
In order to identify coarsening domains for 0 < T < Tc, Derrida [217] proposed to compare
two identical copies evolving under the same thermal noise. One copy starts with random
initial conditions and begins to coarsen, whereas the other copy starts from a fully magne-
tized state and remains ordered as time evolves. It is assumed that all spin ﬂips occurring
in ordered replica can be regarded as thermal ﬂuctuations. Therefore, when a spin ﬂip
occurs simultaneously in both replicas, it can be considered as a thermal ﬂuctuation, while
it is a signature the coarsening process otherwise. This method was used to determine the
fraction of persistent spins of the Ising model [207] as a function of time, conﬁrming that
the persistence exponent does not change for 0 ≤T < Tc.
Since Derrida’s method allows only one type of domains to be detected, it is impossible
110

<!-- page 111 -->

to identify the precise location of domain walls. This can be overcome by comparing three
replicas instead of two [219]. As before, the ﬁrst copy starts with random initial conditions
and serves as the master copy in which the coarsening process takes place. The other
copies start from fully ordered initial conditions with positive and negative magnetization,
respectively. Fluctuations in the ﬁrst copy, which do not occur in the other two replicas,
indicate the presence of a domain wall. More precisely, domain walls may be detected by
the observable
∆i(t) =

1 −
Y
j
1 + σ(1)
j (t)σ(2)
j (t)
2

1 −
Y
j
1 + σ(1)
j (t)σ(3)
j (t)
2

,
(244)
where the superscripts denote the replica. As shown in Fig. 52, this method works sur-
prisingly well.
However, it turns out that the appearance depends on the algorithmic
implementation, i.e., Glauber and heat bath dynamics leads to diﬀerent results. This is
not surprising as the method relies on the same ideas as damage spreading.
5.4
Damage spreading and experiments
In this Section we have seen that the concept of damage spreading depends on the algorith-
mic implementation and therefore lacks a well-deﬁned meaning. The suggested algorithm-
independent deﬁnition of DS resolves these diﬃculties only partly, since it is based on
certain (physically motivated) ad hoc assumptions.
Therefore, DS does not provide a
strict deﬁnition of ‘chaotic’ and ‘regular’ phases in stochastic systems. Nevertheless DS
has been useful to estimate critical exponents of certain nonequilibrium phase transitions
and to identify domain walls in coarsening system at nonzero temperature. Concerning
the question whether DS is relevant for experiments, there is a clear answer: DS is an
artiﬁcial concept which does not exist in nature. In particular, there is no meaning of
‘using the same sequence of random numbers’ in an experimental system. DS is rather a
simulation technique for a pair of systems, taking advantage of our ability to use determin-
istic pseudo random number generators. It is indeed not surprising that such a concept,
for which there is no correspondence in nature, expresses its incompleteness by certain
inconsistencies.
111

<!-- page 112 -->

6
Interface growth
A further important ﬁeld of statistical physics is the study of crystal growth and transitions
between diﬀerent morphologies of moving interfaces [406]. During the last two decades
there has been an enormous progress in the understanding of growth processes (for reviews
see e.g. Refs. [33,407–409]). In this Section we will focus on certain classes of depinning
transitions which are closely related to nonequilibrium phase transitions into absorbing
states.
6.1
Roughening transitions – a brief introduction
Models of growing interfaces may be realized either on a discrete lattice or by continuum
equations. Discrete solid-on-solid (SOS) models are usually deﬁned on a d-dimensional
square lattice of N = Ld sites associated with integer variables hi marking the actual
height of the interface at site i. Clearly, this description does not allow for overhangs
of the interface.
The set of all heights {hi} determines the state of the system which
evolves according to certain stochastic rules for adsorption and desorption. In restricted
solid-on-solid models (RSOS) the dynamic rules satisfy the additional constraint
|hi −hi+1| ≤1 ,
(245)
i.e., adjacent sites may diﬀer by at most one height step. To characterize the evolution of
the interface, it is useful to introduce the mean height ¯h(t) and the width w(t):
¯h(t) = 1
N
X
i
hi(t) ,
w(t) =
h 1
N
X
i
h2
i (t) −
 1
N
X
i
hi(t)
2i1/2
.
(246)
A surface is called ‘smooth’ if the heights hi are correlated over arbitrarily large distances.
Otherwise, if distant heights become uncorrelated, the surface is said to be rough. A rough
surface is typically characterized by a diverging width when the system size is taken to
inﬁnity.
In many cases growing interfaces exhibit simple scaling laws. In a ﬁnite system, starting
from a ﬂat interface, the width ﬁrst increases algebraically as w ∼tγ, where γ > 0 is the
growth exponent of the system15.
In the long-time limit, when the correlation length
reaches the system size, the width saturates at some constant value wsat ∼Lα, where
α is the roughness exponent. This type of scaling behavior is known as Family-Vicsek
scaling [410], as described by the scaling form
w(N, t) ∼N αf(t/N z) ,
(247)
where z = α/γ is the dynamic exponent of the model. f(u) is a universal scaling function
with the asymptotic behavior f(u) ∼uγ for u →0 and f(u) = const for u →∞. Notice
that this scaling form is invariant under rescaling
x →Λx ,
t →Λzt ,
h(x, t) →Λαh(Λx, Λzt) ,
(248)
15In most textbooks the growth exponent is denoted by β. Here we use a diﬀerent symbol in order to
distinguish this exponent from the density exponent of DP.
112

<!-- page 113 -->

where Λ is a dilatation parameter. Notice that in contrast to the order parameter of a
spreading process (92), the ﬂuctuations in the heights increase under rescaling.
The critical exponents γ and z are used to categorize various universality classes of
roughening interfaces.
Typically each of these universality classes is characterized by
certain symmetry properties of the system and may be associated with a speciﬁc stochastic
diﬀerential equation. This equation describes the growth of a continuous height ﬁeld h(x, t)
and consists of the most relevant operators under rescaling (248) that are consistent with
the symmetries of the system. For example, postulating the symmetries
1. translational invariance in space x →x + ∆x,
2. translational invariance in time t →t + ∆t,
3. translational invariance in height direction h →h + ∆h,
4. reﬂection invariance in space x →−x,
5. up/down symmetry h →−h,
one is led to the Edwards-Wilkinson (EW) universality class [411], as described by the
equation
∂h(x, t)
∂t
= v + σ∇2h(x, t) + ζ(x, t) ,
(249)
where v denotes the mean velocity and σ the surface tension. ζ(x, t) is a zero-average
Gaussian noise ﬁeld with variance
⟨ζ(x, t)ζ(x′, t′)⟩= 2Dδd(x −x′)δ(t −t′)
(250)
taking the stochastic nature of deposition into account. This equation is linear and thus
exactly solvable. Scaling invariance (248) implies that the critical exponents are given by
α = 1 −d/2 ,
γ = 1/2 −d/4 ,
z = 2 .
(EW)
(251)
Edwards-Wilkinson growth processes are invariant under up/down reﬂection of the in-
terface h →−h. However, if atoms are adsorbed from a gas phase above the interface
there is no particular reason for the system to be up/down symmetric. In that case the
above equation has to be extended by the most relevant term that breaks the up/down
symmetry, leading to the Kardar-Parisi-Zhang (KPZ) equation [186,412]
∂h(x, t)
∂t
= v + σ∇2h(x, t) + λ(∇h(x, t))2 + ζ(x, t) .
(252)
For a one-dimensional interface the critical exponents of the KPZ universality class are
given by
α = 1/2 ,
γ = 1/3 ,
z = 3/2 .
(KPZ in 1d)
(253)
whereas in d ≥2 dimensions only numerical estimates are known (see Ref. [33]).
It is particularly interesting to study roughening transitions between a smooth and
a rough phase. A roughening transition is usually accompanied by a diverging spatial
correlation length ξ⊥. In the smooth phase this correlation length provides a typical scale
below which the interface appears to be rough. On larger scales, however, the interface
113

<!-- page 114 -->

turns out to be smooth. Approaching the roughening transition the correlation length ξ⊥
diverges whereby the entire interface becomes rough. One of the simplest models displaying
a roughening transition is the two-dimensional discrete Gaussian SOS model [33]. Another
important example is the KPZ equation which exhibits a roughening transition in d > 2
spatial dimensions.
In the following we will focus on certain growth models with depinning transitions.
In particular we will discuss depinning transitions in random media, polynuclear growth
processes, and solid-on-solid growth processes with evaporation at the edges of plateaus. In
the pinned phase of these models the interface is smooth and does not propagate. Varying
a control parameter the interface undergoes a depinning transition; it starts moving and
evolves into a rough state. As we will see below, various depinning transitions are closely
related to phase transitions into absorbing states.
6.2
Depinning transitions of driven interfaces
An interesting class of depinning transitions can be observed in experiments of driven
interfaces in random media [33].
In these experiments a liquid is pumped through a
porous medium. If the driving force F is suﬃciently low the liquid cannot move through
the medium since the air/liquid interface is pinned at certain pores.
Above a critical
threshold, however, the interface starts moving through the medium with an average
velocity v. Close to the transition, v is expected to scale as
v ∼(F −Fc)θ ,
(254)
where θ is the velocity exponent. One of the ﬁrst experiments in 1+1 dimensions was
performed by Buldyrev et al., who studied the wetting of paper in a basin ﬁlled with
suspensions of ink or coﬀee [413]. Here the driving force F is a result of capillary forces
competing with the total weight of the absorbed suspension. Consequently, the interface
becomes pinned at a certain height where F ≃Fc. Once the interface has stopped, the
width should scale as
w(ℓ, t) ∼ℓαf(t/ℓ˜z) ,
(255)
where ℓis a box size and α the roughening exponent.
Measuring the interface width
Buldyrev et al. found the roughness exponent α = 0.63(4). In various other experiments
the values are scattered between 0.6 and 1.25. This is surprising since the Kardar-Parisi-
Zhang (KPZ) class [186] predicts the much smaller value α = 1/2.
It is believed that the large values of α are due to inhomogeneities of the porous
medium. Because of these inhomogeneities, the interface does not propagate uniformly by
local ﬂuctuations as in the KPZ equation, it rather propagates by avalanches. In the lit-
erature two universality classes for this type of interfacial growth have been proposed. In
case of linear growth the interface should be described by a random ﬁeld Ising model [414],
leading to the exponents α = 1, ˜z = 4/3, and θ = 1/3 in 1+1 dimensions. In the presence
of a KPZ-type nonlinearity, however, the roughening process should exhibit a depinning
transition which is related to DP [415,416]. The underlying DP mechanism diﬀers signif-
icantly from an ordinary directed percolation process in a porous medium subjected to a
gravitational ﬁeld (cf. Sec. 3.9). In ordinary DP the spreading agent represents active sites
and percolates along the given direction. In the present case, however, water may ﬂow not
114

<!-- page 115 -->

1
q
q
1
Fig. 53:
Simple model for interface depinning in random media. The pores of the medium are represented
by cells on a diagonal square lattice.
The permeability across the edges of the cells depends on the
direction of ﬂow: In downward direction all edges are permeable, whereas in upward direction they are
permeable with probability q and impermeable otherwise. The right panel of the ﬁgure shows a particular
conﬁguration of open (dashed) and closed (solid) edges.
Pumping in water from below, the interface
becomes pinned along a directed path of solid lines, leading to a ﬁnite cluster of wet cells (shaded region).
The dashed arrow represents an open path in order to illustrate a possible ﬂow.
only in the direction of the pumping force but also in opposite direction. Moreover, the
DP process runs perpendicular to the direction of growth, as will be explained below.
A simple model exhibiting a depinning transition is shown in Fig. 53. In this model
the pores of the material are represented by cells of a diagonal square lattice. The liquid
can ﬂow to neighboring cells by crossing the edges of a cell. Depending on the direction of
ﬂow these edges can either be permeable or impermeable. For simplicity we assume that
all edges are permeable in downwards direction, whereas in upwards direction they can
only be crossed with a certain probability q. Thus, by starting with a horizontal row of
wet cells at the bottom, we obtain a compact cluster of wet cells, as illustrated in Fig. 53.
The unrestricted ﬂow downwards ensures that the cluster has no overhangs. Clearly, the
size of the cluster (and therefore the penetration depth of the liquid) depends on q. If q is
large enough, the cluster is inﬁnite, corresponding to a moving interface. If q is suﬃciently
small, the cluster is bound from above, i.e., the interface becomes pinned.
The depinning transition is related to DP as follows.
As can seen in the ﬁgure, a
pinned interface may be interpreted as a directed path along impermeable edges running
from one boundary of the system to the other. Obviously, the interface becomes pinned
only if there exists a directed path of impermeable bonds connecting the boundaries of
the system.
Hence the depinning transition is related to an underlying directed bond
percolation process with probability p = 1 −q running perpendicular to the direction of
growth. The pinning mechanism is illustrated in Fig. 54, where a supercritical DP cluster
propagates from left to right. The cluster’s backbone, consisting of bonds connecting the
two boundaries, is indicated by bold dots. The shaded region denotes the resulting cluster
of wet cells. As can be seen, the interface becomes pinned at the lowest lying branch
of the DP backbone. Therefore, the roughening exponent coincides with the meandering
exponent of the backbone
α = ν⊥/ν∥.
(256)
Moreover, by analyzing the dynamics of the moving interface, it can be shown that the
dynamic critical exponents are given by θ = α and ˜z = 1. Thus, depinning transitions in
inhomogeneous porous media may serve as possible experimental realizations of the DP
universality class.
The prediction (256) matches surprisingly well with the experimental result α = 0.63(4)
115

<!-- page 116 -->

Fig. 54:
Pinned interface and the underlying DP process. The ﬁgure is explained in the text.
obtained by Buldyrev et al. [413]. Therefore, it is near at hand to regard this experiment
as a ﬁrst quantitative experimental evidence of DP exponents. However, only one ex-
ponent has been veriﬁed, and it is not fully clear how accurate and reproducible these
exponents are. Moreover, pressure diﬀerences may cause long-range correlations in the
bulk, leading to a ﬂat interface on large scales. This means that gravity could destroy
the asymptotic critical behavior. In fact, in subsequent experiments the estimates for the
critical exponents are scattered over a wide range. For example, Dougherty and Carle
measured the dynamical avalanche distribution of an air/water interface moving through
a porous medium made of glass beads [417]. According to the DP hypothesis, the distribu-
tion P(s) of avalanche sizes s should decrease algebraically. In the experiment, however, a
stretched exponential behavior P(s) ∼s−be−s/L is observed even for small ﬂowing rates.
The estimates for the exponent b are inconclusive; they depend on the time window of
the measurement and vary between −0.5 and 0.85. Even more recently Albert et al. pro-
posed to identify the universality class by measuring the propagation velocity of locally
tilted parts of the interface [418]. Their results suggest that interfaces propagating in glass
beads are not described by a DP depinning process, instead they seem to be related to
the random-ﬁeld Ising model. Altogether the emerging picture is not yet fully transparent
and further experimental eﬀort in this direction would be desirable.
Depinning experiments were also carried out in 2+1 dimensions with a spongy-like
material used by ﬂorists, as well as ﬁne-grained paper rolls [419]. In this case, however, the
exponent α is not related to 2+1-dimensional DP, instead it corresponds to the dynamic
exponent of percolating directed interfaces in 2+1 dimensions. In experiments as well as
in numerical simulations a roughness exponent α = 0.50(5) was obtained.
6.3
Polynuclear growth
A completely diﬀerent type of depinning transition takes place in models for polynuclear
growth (PNG) [141,420–422]. A key feature of PNG models is the use of parallel updates,
leading to a maximal propagation velocity of one monolayer per time step. For a high
adsorption rate the interface of PNG models is smooth and propagates at maximal velocity
v = 1. Roughly speaking, this means that the interface is ‘pinned’ at the light cone of the
dynamics. Decreasing the adsorption rate below a certain critical threshold, PNG models
exhibit a roughening transition to a rough phase with v < 1. In contrast to equilibrium
roughening transitions, which only exist in d ≥2 dimensions, PNG models display a
roughening transition even in one spatial dimension.
116

<!-- page 117 -->

t
t+1/2
t+1
Fig. 55:
Polynuclear growth model. In the ﬁrst half time step atoms are deposited with probability p. In
the second half time step islands grow deterministically by one step and coalesce.
One of the simplest PNG model investigated so far is deﬁned by the following dynamic
rules [141]. In the ﬁrst half time step atoms ‘nucleate’ stochastically at the surface by
hi(t + 1/2) =
 hi(t) + 1
with prob. p ,
hi(t)
with prob. 1 −p .
(257)
In the second half time step the islands grow deterministically in lateral direction by one
step. This type of growth may be expressed by the update rule
hi(t + 1) = max
j∈<i>

hi(t + 1/2), hj(t + 1/2)

,
(258)
where j runs over the nearest neighbors of site i.
The relation to DP can be established as follows. Starting from a ﬂat interface hi(0) =
0, let us interpret sites at maximal height hi(t) = t as active sites of a DP process.
Obviously, the adsorption process (257) turns active sites into the inactive state with
probability 1 −p, while the process (258) resembles oﬀspring production. Therefore, if
p is large enough, the interface is smooth and propagates with maximal velocity v = 1.
This situation corresponds to the active phase of DP. Approaching the phase transition,
we expect the density of sites at maximal height to scale as
nmax = 1
N
X
i
δhi−t ∼(p −pc)β,
(259)
where N denotes the system size and β ≃0.277 is the density exponent of DP. Below a
critical threshold, however, the density of active sites at the maximum height hi(t) = t
vanishes after some time, the growth velocity is smaller than 1, and the interface evolves
into a rough state.
Although this mapping to DP is not exact, numerical simulations
strongly support the validity of Eq. (259).
As will be shown below, PNG models are
actually a realization of unidirectionally coupled DP processes [423,424].
Polynuclear growth models have also been used to describe the growth of colonial or-
ganisms such as fungi and bacteria [425]. This study was motivated by recent experiments
with the yeast Pichia membranaefaciens on solidiﬁed agarose ﬁlm [426]. By varying the
concentration of polluting metabolites, diﬀerent front morphologies were observed. The
model proposed in [425] aims to explain these morphological transitions on a qualitative
level. It is easy to verify that this model follows the same spirit as the PNG model de-
ﬁned above. They both employ parallel dynamics and exhibit a DP-related roughening
transition.
Concerning experimental realizations of PNG models, one major problem – apart from
quenched disorder – is the use of parallel updates. The type of updates in these models
is crucial; by using random-sequential updates the transition is lost since in this case
117

<!-- page 118 -->

                         




            


                    




            


                         




                    




                    




                    




1-q
q
q
q
q
q=0.15
q=0.23
q=0.189
(1-q)/2
(1-q)/2
Fig. 56:
Restricted solid-on-solid growth model exhibiting a roughening transition from a non-moving
smooth to a moving rough phase. Monomers are randomly deposited whereas desorption takes place only
at the edges of plateaus.
there is no maximum velocity.
However, in realistic experiments atoms do not move
synchronously, rather the adsorption events are randomly distributed in time. Therefore,
random sequential updates might be more appropriate to describe such experiments. It
thus remains an open question to what extent PNG processes can be realized in nature.
6.4
Growth with evaporation at the edges of plateaus
DP-related roughening transitions can also be observed in certain solid-on-solid growth
processes with random-sequential updates [427, 428]. As a key feature of these models,
atoms may desorb exclusively at the edges of existing layers, i.e., at sites which have at
least one neighbor at a lower height. This corresponds to a physical situation where the
binding energy of atoms at the edges is much smaller than the binding energy of atoms in
completed layers. By varying the growth rate, such growth processes display a roughening
transition from a non-moving smooth phase to a moving rough phase.
A simple solid-on-solid model for this type of growth is deﬁned by the following dynamic
rules [427]: For each update a site i is chosen at random and an atom is adsorbed
hi →hi + 1 with probability q
(260)
or desorbed at the edge of a plateau
hi →min(hi, hi+1) with probability (1 −q)/2 ,
hi →min(hi, hi−1) with probability (1 −q)/2 .
(261)
Moreover, the growth process is assumed to be restricted, i.e., updates are only carried
out if the resulting conﬁguration obeys the constraint (245).
The qualitative behavior of this model is illustrated in Fig. 56. For small q the desorp-
tion processes (261) dominate. If all heights are initially set to the same value, this level
will remain the bottom layer of the interface. Small islands will grow on top of the bottom
layer but will be quickly eliminated by desorption at the island edges. Thus, the interface
is eﬀectively anchored to its bottom layer and a smooth phase is maintained. The growth
velocity v is therefore zero in the thermodynamic limit. As q is increased, more islands on
top of the bottom layer are produced until above qc ≃0.189, the critical value of q, they
merge forming new layers at a ﬁnite rate, giving rise to a ﬁnite growth velocity.
118

<!-- page 119 -->

Interestingly, this model can be interpreted as a driven diﬀusion process of two oppo-
sitely charged types of particles. The charges
ci,i+1 = hi+1 −hi ∈{−1, 0, +1}
(262)
are bond variables and represent a change of height between two adjacent interface sites.
In this representation the dynamic rules (260)-(261) can be implemented by randomly
selecting two neighboring bonds and performing the following processes with probabilities
indicated on the arrows:
Ø+
q
⇄
(1−q)/2
+Ø ,
ØØ
q
⇄
1−q
+−,
−Ø
q
⇄
(1−q)/2
Ø−,
−+
q→ØØ .
(263)
In the smooth phase q < qc, the charges are arranged as closely bound +−dipoles.
For q > qc, the dipoles become unbound wherefore the ﬂuctuations in the total charge,
measured over a distance of order N, diverge with N. Thus the transition can be described
in terms of correlations between charged particles.
Relation to directed percolation
At the transition the dynamics of the model is related to DP as follows. Starting with a ﬂat
interface at zero height, let us consider all sites with hi = 0 as particles A of a DP process.
Growth according to Eq. (260) corresponds to spontaneous annihilation A →Ø. Con-
versely, desorption may be regarded as a particle creation process. However, since atoms
may only desorb at the edges of plateaus, particle creation requires a neighboring active
site to be present. This process therefore corresponds to oﬀspring production A →2A.
Clearly, these processes resemble the dynamic rules of a DP process. In contrast to PNG
models, the DP process takes place at the bottom layer of the interface. Moreover, the
roughening transition does not depend on the use of either parallel or random-sequential
updates.
In the model without the restriction (245), the relation to DP can be proven exactly. In
this case the processes at the bottom layer decouple from the evolution of the interface at
higher levels. Introducing local variables si = δhi,0 it is possible to map Eqs. (260)–(261)
exactly onto the dynamic rules
if si = 1
si →0 with prob. q ,
if {si−1, si, si+1} = {0, 0, 1}
si →1 with prob. (1 −q)/2 ,
if {si−1, si, si+1} = {1, 0, 0}
si →1 with prob. (1 −q)/2 ,
if {si−1, si, si+1} = {1, 0, 1}
si →1 with prob. 1 −q .
These rules deﬁne a contact process on a square lattice (cf. Fig. 16) in which particles
self-annihilate at unit rate and create oﬀspring with rate λ/2 = (1 −q)/2q. Since the
percolation threshold of the contact process λc = 3.29785(8) is known, we can predict the
critical point of the unrestricted growth to be given by qc = 0.232675(5). model. For the
restricted model there is no exact mapping to a contact process. Nevertheless numerical
simulations strongly suggest that the bottom layer still exhibits DP behavior.
The underlying DP transition is the origin of the roughening transition. In the active
phase of DP the interface ﬂuctuates close to the bottom layer so that the interface is
smooth. Approaching criticality the bottom layer occupation n0 vanishes as
n0 ∼(qc −q)β ,
(264)
119

<!-- page 120 -->

10
100
1000
0.1
1
n=1
n=2
n=3
n=4
n=5
t
n
t
<M  >
Fig. 57:
Spontaneous symmetry breaking. Left: Simulation in the smooth phase starting from a random
interface conﬁguration. Black (white) pixels indicate even (odd) heights. As time evolves the systems
’coarsens’, providing a robust mechanism for the elimination of minority islands. Right: Order parameters
Mn versus time at criticality starting from a ﬂat interface.
where β is the density exponent of DP. In the inactive phase of DP the interface detaches
from the bottom layer and evolves into a rough state. Therefore, the front velocity v for
q > qc is proportional to the characteristic survival time of a DP process in the inactive
phase
v ∼1/ξ∥∼(q −qc)1/ν∥.
(265)
With respect to the microscopic rules the roughening transition in these models seems to
be as robust as a DP transition. Including the dynamics at higher levels of the interface,
the models turn out to be described by unidirectionally coupled DP processes (see below).
Scaling of the interface width
Turning to the scaling properties of the interface width (246) the emerging picture is less
clear and may indicate the presence of many length scales. With only a single length scale
one would expect the saturation width in ﬁnite systems to scale as wsat(N) ∼N α in the
growing phase and wsat(N) ∼N α′ at criticality, where α and α′ are generally diﬀerent
critical exponents. Thus the expected scaling form would read
wsat(N, q) = N α′ g
 N 1/ν⊥(q −qc)

,
(266)
where g(u) is a universal scaling function with the asymptotic behavior
g(u) =





|u|−α′ν⊥
for u →−∞,
const
for u →±0 ,
|u|(α−α′)ν⊥
for u →+∞.
(267)
This scaling form implies that the width in a ﬁnite size system saturates at
wsat(N, q) =





(qc −q)−α′ν⊥
if q < qc ,
N α′
if q = qc ,
N α(qc −q)(α−α′)ν⊥
if q < qc .
(268)
120

<!-- page 121 -->

However, the width at criticality actually increases as w(t) ∼(ln t)γ until it saturates at
wsat ∼(ln N)κ ,
(269)
where κ ≃0.43, suggesting that α′ = 0.
Moreover, as the growth rate is increased
and the interface starts to move, the width diverges as wsat(q) −wsat(qc) ∼(q −qc)0.95,
suggesting that α ≃0.9. This value diﬀers from the expected value α = 1/2 for growing
interfaces in one spatial dimension. Therefore, the scaling behavior of the interface width is
presumably characterized by many diﬀerent length scales. This point of view is conﬁrmed
by a numerical analysis of correlation functions at diﬀerent levels [428].
Spontaneous symmetry breaking
The roughening transition in growth models with evaporation at edges of islands is accom-
panied by spontaneous symmetry breaking of translational invariance in height direction.
This symmetry is associated with a family of non-conserved magnetization-like order pa-
rameters
Mn =
 1
N
N
X
j=1
exp
2πihj
n + 1
 .
(270)
The order parameter M1 = | 1
N
P
j(−1)hj|, for example, measures the diﬀerence between
the densities of sites at even and odd heights, respectively (see Fig. 57). In the smooth
phase the order parameters Mn tend towards stationary positive values. Approaching the
phase transition these values vanish algebraically as
⟨Mn⟩∼(qc −q)θn ,
(271)
where θ1 ≈0.65, θ2 ≈0.40, and θ3 ≈0.23. At criticality, starting from a ﬂat interface, the
order parameters decay as Mn(t) ∼t−θn/ν∥, as shown in the right hand graph of Fig. 57.
In the rough phase q > qc, where the heights hi become uncorrelated over large distances,
Mn = 0 in the thermodynamic limit.
Experimental realizations
With respect to experimental realizations of the dynamic rules (260)-(261) it is important
to note that atoms are not allowed to diﬀuse on the surface. This assumption is rather
unnatural since in most experiments the rate for surface diﬀusion is much higher than
the rate for desorption back into the gas phase. Therefore, it will be diﬃcult to realize
this type of homoepitaxial growth experimentally.
However, in a diﬀerent setup, the
above model could well be relevant [429]. As illustrated in Fig. 58, a laterally growing
monolayer could resemble the dynamic rules (260) and (261) by identifying the edge of the
monolayer with the interface contour of the growth model. In this case ‘surface diﬀusion’,
i.e. diﬀusion of atoms along the edge of the monolayer, is highly suppressed. Moreover,
in single-step systems (such as fcc(100) surfaces) it would also be possible to implement
the restriction (245).
6.5
Unidirectionally coupled directed percolation processes
So far we have seen that the scaling properties of quantities involving only the bottom
level of the interface may be adequately described by using DP exponents. In particular it
121

<!-- page 122 -->

Fig. 58:
Possible experimental realization of the growth process deﬁned in Eqs. (260) and (261) by lateral
growth of a monolayer on a substrate.
10
−3
10
−2
10
−1
0.1
1
10
−3
10
−2
10
−1
0.1
1
10
−3
10
−2
10
−1
0.1
1
k=0
k=1
k=2
k=3
mk
mk
ε
(R)
ε
(U)
k
nk
(a)
k
p-p
n
c
(b)
q -q
c
(c)
q -q
c
n
Fig. 59:
Coupled directed percolation. Stationary densities nk as functions of the distance from criticality
in (a) the unrestricted growth model, (b) the restricted growth model, and (c) a unidirectionally coupled
sequence of directed bond percolation processes.
was shown that the density of exposed sites at the bottom layer decreases as n0 ∼(qc−q)β
until it vanishes at the transition. Let us now turn to the scaling properties of the ﬁrst
few layers k = 1, 2, . . . above the bottom layer. Since the scaling properties at the bottom
layer k = 0 in the smooth phase are completely characterized by the three DP exponents
β, ν⊥, and ν||, it is natural to assume that the next layers obey similar scaling laws with
analogous exponents, β(k), ν(k)
⊥, and ν(k)
|| . For example, the densities
nk = 1
N
k
X
j=0
X
i
δhi,j ,
k = 0, 1, 2, . . .
(272)
of sites at height ≤k above the bottom layer are expected to scale as
nk ∼(qc −q)β(k) .
k = 1, 2, 3, . . .
(273)
In principle all these exponents could be diﬀerent and independent from each other. How-
ever, extensive numerical simulations and ﬁeld-theoretic considerations [428] suggest that
the scaling exponents ν(k)
⊥
and ν(k)
||
are identical on all levels and equal to the DP ex-
ponents ν⊥and ν||. This remarkable property implies that the growth process at criti-
cality is characterized by a single dynamic exponent z = ν||/ν⊥. The density exponents
β(1), β(2), β(3), . . . , however, turn out to be diﬀerent and considerably reduced compared
β(0) = β. These exponents appear to be non-trivial in the sense that they are not simply
related to DP exponents.
In order to explain the reduced values of β(k) it is useful to introduce the following
particle interpretation. Assuming that the minimal height of the interface is zero, lattice
122

<!-- page 123 -->

sites at height hi ≤0, 1, 2, . . . may be associated with particles A, B, C, . . . , respectively.
By deﬁnition, particles of diﬀerent species are allowed to occupy the same site simulta-
neously. For example, the presence of an A-particle induces simultaneous presence of all
other particle species at the same site. As shown before, the A-particles of the unrestricted
model evolve independently according to the dynamic rules of a contact process. Similarly,
in absence of A-particles, the B-particles will evolve according to the rules of a contact
process. In presence of an A-particle, however, B-particles are instantaneously created
at the same site, giving rise to an eﬀective reaction A →A + B at inﬁnite rate. As this
reaction does not modify the conﬁguration of the A-particles, it couples the two subsys-
tems in one direction without feedback. Similarly, the B-particles induce the creation of
C-particles by an eﬀective reaction B →B +C. Thus we obtain a unidirectionally coupled
sequence of contact processes corresponding to the reaction-diﬀusion scheme
A ↔2A
B ↔2B
C ↔2C
A →A + B
B →B + C
C →C + D
etc.
Notice that this sequence can be truncated at any level without changing the dynamics of
the lower levels. In fact, numerical simulations of the growth model and a unidirectionally
coupled sequence of three (1+1)-dimensional directed bond percolation processes yield
compatible estimates for the critical exponents (cf. Fig. 59):
β(0) = 0.28(1) ,
β(1) = 0.13(2) ,
β(2) = 0.05(1) .
(274)
Thus, ‘coupled DP’ is a universal phenomenon and should play a role even in more gen-
eral contexts, namely whenever DP-like processes are coupled unidirectionally without
feedback.
Mean ﬁeld approximation
The simplest set of Langevin equations for unidirectionally coupled DP reads [428]
∂tφk(x, t) = κkφk(x, t) −λφ2
k(x, t) + D∇2φk(x, t) + µφk−1(x, t) + ζk(x, t) ,
(275)
where ζk are independent multiplicative noise ﬁelds with correlations
⟨ζk(x, t)⟩= 0 ,
⟨ζk(x, t)ζl(x′, t′)⟩= 2Γ φk(x, t) δk,l δd(x −x′) δ(t −t′) .
(276)
Assuming that φ−1 ≡0, the lowest equation for k = 0 reduces to the ordinary Langevin
equation of DP (139). The parameters κk control the rates for oﬀspring production at
level k. In principle all these parameters could be diﬀerent, leading to interesting mul-
ticritical behavior [423]. Here we will restrict to the simplest case where the parameters
κk = κ coincide.
The reduced values of β(1), β(2), . . . can already be understood on the level of a sim-
ple mean ﬁeld approximation.
Determining the stationary solutions of the mean ﬁeld
equations
∂tφk = κφk −λφ2
k + µφk−1
(277)
and expanding the result for small values of κ, the ﬁelds φk are found to scale asymptoti-
cally with the mean ﬁeld exponents
β(k)
MF = 2−k .
(278)
123

<!-- page 124 -->

In addition, dimensional analysis of Eqs. (275)-(276) reveals that the mean ﬁeld scaling
exponents ν(k)
⊥
and ν(k)
∥
coincide with the DP exponents νMF
⊥
= 1/2 and νMF
∥
= 1 on all
levels, i.e., they do not depend on k.
Field-theoretic treatment
The above mean ﬁeld approximation is expected to hold above the critical dimension
dc = 4. In less than four dimensions ﬂuctuation corrections to β(k) have to be taken into
account. These corrections can be approximated by ﬁeld-theoretic renormalization group
techniques [423,424]. To this end the master equation is mapped onto a second-quantized
bosonic operator representation, leading to the eﬀective action (cf. Sect 3.5)
S[φ0, ˜φ0, φ1, ˜φ1, . . . ] =
K−1
X
k=0
Z
ddx dt
(
˜φk(x, t)

τ∂t −D∇2 −κ

φk(x, t) −µ˜φkφk−1
+ Γ
2
˜φk(x, t)

φk(x, t) −˜φk(x, t)

φk(x, t)
)
,
(279)
where K is the total number of levels in the hierarchy.
Because of the unidirectional
structure the truncation of the hierarchy at ﬁnite K does not aﬀect the temporal evolution
of the subsystems k ≤K. The ﬁeld-theoretic treatment turns out to be rather complex,
even a one-loop calculation for a two-level system involves 34 diﬀerent diagrams. Details
of these calculations are given in Ref. [424], whose main results we summarize here. A
one-loop calculation in the inactive phase reveals that the RG ﬂow is characterized by
two ﬁxed lines. One of them is unstable and corresponds to a situation where the two
systems are decoupled. The other one is stable and describes the interacting case. One
can show that strongly ultraviolet-divergent contributions cancel along the stable ﬁxed
line. To determine the exponent β(1), similar calculations have to be performed in the
active phase. Choosing a particular point along the ﬁxed line it is possible to derive the
result
β(0) = 1 −ǫ/6 ,
β(1) = 1/2 −ǫ/8 ,
(280)
where ǫ = 4 −d. This result explains the downward correction of the critical exponent
β(1) in less than four dimensions.
The ﬁeld-theoretic analysis involves various technical and conceptual problems. On the
one hand, in the active phase several infrared-divergent diagrams are encountered [430],
without being clear to what extent they will aﬀect the physical properties of coupled DP.
On the other hand, the coupling constant µ between diﬀerent levels is shown to be a
relevant quantity, i.e., it grows and ﬁnally diverges under RG transformations. This may
be the cause for numerically observed violations of scaling. In fact, the curves in Fig. 59
for k ≥1 are not perfectly straight but slightly bent. This curvature is neither related
to transients nor to ﬁnite size eﬀects. A careful analysis shows that the ﬁeld-theoretic
prediction seems to apply in an intermediate scaling regime, whose size depends on the
coupling constant µ. The breakdown of scaling may be an artifact of the lattice realization
which, in contrast to the ﬁeld-theoretic prediction, seems to limit the value of µ.
124

<!-- page 125 -->

                         




                         




                         




                         




                    



                         




                         




               


q
q=0.2
q=0.317
q=0.7
1-q
1-q
q
q
q
1-q
Fig. 60:
Roughening transition in a model for dimer adsorption and desorption. Left: Dynamic rules for
adsorption and desorption at the edges of terraces. Right: Typical interface conﬁgurations below, at, and
above criticality.
6.6
Parity-conserving roughening transitions
So far we discussed two examples of depinning transitions related to DP, namely interface
depinning in driven random media, polynuclear growth, and growth without evaporation
at the edges of terraces. It is therefore interesting to ask whether it is possible to ﬁnd
examples of roughening transitions with non-DP behavior in one spatial dimension.
As discussed in Sec. 4.2 non-DP phase transitions into absorbing states can be observed
in systems with additional symmetries.
For example, non-DP behavior is observed in
parity-conserving branching processes. Therefore, the question arises whether it is possible
to replace the underlying DP transition in the previously discussed growth models by a
parity-conserving mechanism. As shown in Ref. [431], this can be done by considering
the growth of dimers with evaporation at edges of plateaus. As dimers consist of two
atoms, the number of particles at each height level is preserved modulo 2. This deﬁnition
of the dynamic rules mimics a BAWE at the bottom layer of the interface. It turns out
that some of the critical properties of the roughening transition in this model are indeed
characterized by PC exponents (cf. Eq. (221)).
As shown in Fig. 60, the model generalizes the growth model introduced in Sec. 6.4,
replacing monomers by dimers. In each attempted update two adjacent sites i and i + 1
are selected randomly. If the heights hi and hi+1 are equal, one of the following moves is
carried out. Either a dimer is adsorbed with probability p
hi →hi + 1 ,
hi+1 →hi+1 + 1
(281)
or desorbed with probability 1 −p:
hi, hi+1 →min(hi−1, hi, hi+1, hi+2) .
(282)
An update will be rejected if it leads to a violation of the RSOS restriction (245).
The transition of the model is illustrated in Fig. 60: If q is very small, only a few
dimers are adsorbed for a short time so that the interface is smooth and pinned at the
(spontaneously selected) bottom layer. As q increases, a growing number of dimers covers
the surface and large islands of several layers stacked on top of each other are formed.
When q exceeds the critical value qc ≃0.317(1), the mean size of the islands diverges and
the interface evolves into a rough state.
125

<!-- page 126 -->

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
t
0.1
1
n0
n1
n2
10
1
10
2
10
3
10
4
t
0.1
nA
nB
nC
10
−6
10
−4
10
−2
10
0
10
2
10
4
t / L
z
1
0.3
t
−α exp(W
2)
(a)
(c)
(b)
Fig. 61:
Dimer model: (a) Finite size scaling of the width. (b) Decay of the densities n0, n1, n2 as a
function of time. (c) Decay of the corresponding densities in unidirectionally coupled BAWE’s.
Naively one may expect the interface to detach from the bottom layer in the rough
phase, resulting in a ﬁnite propagation velocity. In the present case it turns out, how-
ever, that the interface remains pinned to the initial height, i.e., it does not propagate
at constant velocity. This is due to the fact that a stochastic deposition process cannot
create a dense packing of dimers. Instead the emerging conﬁgurations are characterized
by a certain density of defects (solitary sites at the bottom layer) where dimers cannot
be adsorbed. Because of the RSOS condition (245) these defects act as ‘pinning centers’
which prevent the interface from growing, leading to triangular conﬁgurations as shown
in Fig. 60. The pinning centers cannot disappear spontaneously, they can only diﬀuse by
interface ﬂuctuations and recombine in pairs so that their number is expected to decrease
very slowly. Therefore, the interface of an inﬁnite system in the rough phase does not
propagate at constant velocity. Instead the squared width and the average height grow
logarithmically with time as ¯h(t) ∼w(t) ∼log t. At criticality, this behavior crosses over
to w2(t) ∼¯h(t) ∼log t. Thus the expected scaling form reads
w2(L, t) ≃log
h
t−γ F(t/Lz)
i
,
(283)
where F is a universal scaling function and γ plays the role of a growth exponent. For
γ = 0.172(10) and z = 1.76(5) we obtain a fairly accurate data collapse (see Fig. 61a),
supporting the supposition that z is the dynamic exponent of the PC universality class.
Similarly the densities nk (see Eq. (272)) are found to decay at criticality as nk(t) ∼t−δk
with
δ0 = 0.280(10) , δ1 = 0.200(15) , δ2 = 0.120(15) ,
(284)
where δ0 = β/ν∥is the usual cluster survival exponent of the PC class.
Using the particle interpretation of Sec. 6.5, the temporal evolution of the A-particles
resembles a BAWE. Similarly, the B-particles perform an eﬀective BAWE on top of inac-
tive islands of the A-system. Since A-particles instantaneously create B-particles, the two
subsystems are coupled by the eﬀective reaction A →A + B at inﬁnite rate. As this reac-
tion does not modify the conﬁguration of the A-particles, it couples the two subsystems
only in one direction without feedback. On the other hand, the RSOS condition (245)
introduces an eﬀective feedback so that the A-particles are not completely decoupled from
the B-particles. However, it seems that the inhibiting inﬂuence of the B-particles does not
aﬀect the critical behavior of the A-particles. Similarly, the C-particles are coupled to the
126

<!-- page 127 -->

                                                                



substrate
bulk phase
interface
surface phase
Fig. 62:
Wetting of a surface – schematic illustration.
B-particles by the eﬀective reaction B →B + C. Therefore, the dimer model resembles a
unidirectionally coupled sequence of BAWE’s, corresponding to the reaction scheme
A →3A
B →3B
C →3C
2A →Ø
2B →Ø
2C →Ø
A →A + B
B →B + C
C →C + D
etc.
In fact, as shown in Fig. 61 b-c, a coupled hierarchy of three BAWE’s shows the same
critical behavior as the dimer model at the ﬁrst few layers.
6.7
Nonequilibrium wetting transitions
Wetting phenomena are observed in various physical systems where a bulk phase (e.g. a
gas or liquid) is brought into contact with a wall or a substrate. Because of interactions
between bulk phase and surface, a thin layer of another phase may be formed which is
attracted to the substrate. The thickness of the layer ﬂuctuates and may depend on various
parameters such as temperature or chemical potential. As some of these parameters are
varied, the thickness of the layer may diverge, leading to a wetting transition. Theoretical
models for wetting phenomena ignore the details of molecular interactions between surface,
layer and bulk phase.
Instead, they characterize the system by an interface without
overhangs that separates the two phases. The conﬁguration of the interface is given in
terms of the height h(x) ≥0 of the interface at point x on the surface.
Within this
approach, wetting transitions may be viewed as the unbinding of the interface from the
wall.
Wetting transitions at thermal equilibrium have been theoretically studied and exper-
imentally observed in a large variety of systems (for a review, see Ref. [432]). Models for
equilibrium wetting are usually deﬁned by an eﬀective Hamiltonian of the form
H =
Z
dd−1x
σ
2 (∇h)2 + V [h(x)]

,
(285)
where σ denotes the surface tension of the interface and d −1 the dimension of the in-
terface [433–435]. The potential V [h(x)] yields the eﬀective interaction between wall and
interface. Usually it contains an attractive component binding the interface to the wall.
However, as temperature or other parameters are varied, the attractive component of the
potential may become so weak that the potential is longer able to bind the interface,
127

<!-- page 128 -->

0.0
0.5
1.0
1.5
2.0
p
0.0
0.5
1.0
1.5
q
special case q=p=1
λ=0
DP
exactly soluble case
q=p=1  (EW class)
p
growing phase
bound phase
KPZ
KPZ
q
Fig. 63:
Phase diagram of the nonequilibrium wetting model. The wetting transition takes place along
the solid line. Along the dashed line the growth velocity in the model without a wall does not depend on
a global tilt of the interface, indicating that the eﬀective KPZ nonlinearity vanishes. Both lines intersect
in the point q = p = 1, where the microscopic processes obey detailed balance.
leading to a wetting transition. In d = 2 dimensions one usually distinguishes between
critical and complete wetting. Critical wetting refers to the divergence of the interface
width when the wetting transition is approached by moving along the coexistence curve of
the bulk and surface phases. On the other hand, complete wetting refers to the divergence
of the interface width when the chemical potential diﬀerence between the two phases is
varied, moving towards the coexistence curve. Critical and complete wetting transitions
are associated with generally diﬀerent critical exponents.
While equilibrium wetting transitions are well understood, the investigation of wetting
transitions under nonequilibrium conditions has started only recently. Here the surface
layer is adsorbed to the wall by a growth process whose dynamics, unlike in equilibrium
processes, does not obey detailed balance. As expected, the resulting critical exponents
diﬀer from those observed at thermal equilibrium.
A lattice model for nonequilibrium wetting
A simple lattice model for nonequilibrium wetting can be deﬁned as follows [436]. As in
Sec. 6.4, the interface is given by height variables hi = 0, 1, . . . , ∞. For each update a site
i is selected randomly and one of the following moves is carried out:
hi →hi + 1
with prob. q/(q + p + 1) ,
hi →min(hi−1, hi, hi+1)
with prob. 1/(q + p + 1) ,
(286)
hi →hi −1
if hi−1 = hi = hi+1
with prob. p/(q + p + 1) .
The selected move will be rejected if it would result in a violation of the RSOS constraint
|hi −hi+1| ≤1. In addition, a hard-core wall at zero height h = 0 is introduced, i.e., a
process is only carried out if the resulting interface heights are non-negative. Generally
the processes deﬁned above do not satisfy detailed balance.
128

<!-- page 129 -->

By varying the relative rates of these processes, a transition from a binding to a non-
binding phase is found. This wetting transition can be understood as follows. Without the
wall for ﬁxed p > 0, the parameter q controls the mean growth velocity of the interface.
This velocity may be positive or negative and vanishes at some critical value q = qc. On
large time scales a lower wall will only aﬀect the interface dynamics if the interface moves
downwards, i.e.
q ≤qc, leading to a smooth interface.
In the moving phase q > qc,
however, the interface does not feel the wall and evolves into a rough state (see Fig. 63).
It is obvious that in the moving phase the interface velocity scales as v ∼(q −qc)y with
y = 1. In the smooth phase q < qc, the expected scaling for bottom layer occupation and
width is given by
ρ0 ∼(qc −q)x0 ,
w ∼(qc −q)−γ ,
(287)
where x0 and γ are certain critical exponents.
The above model includes two special
cases, namely p = 0, where the interface cannot move below its actual minimum height,
and p = 1, where the dynamic rules satisfy detailed balance.
The exactly solvable case p = 1
For p = 1 the dynamic rules (286) can be mapped onto an exactly solvable equilibrium
model which exhibits a transition to complete wetting. For q < 1 the probability of ﬁnding
the interface in the conﬁguration {h1, . . . , hN} is given by the distribution
P(h1, . . . , hN) = P(H) = Z−1 qH(h1,... ,hN) ,
(288)
where H = H(h1, . . . , hN) = PN
i=1 hi is the sum of all heights and Z = P
h1,... ,hN qH
denotes the partition sum running over all possible interface conﬁgurations. Eq. (288) can
be proven by verifying the detailed-balance condition
w(σH →σH+1)/w(σH+1 →σH) = q .
(p = 1)
(289)
which is consistent with the hard wall constraint hi ≥0.
The steady state distribu-
tion (288) does not exist for q > 1 where the interface propagates at constant velocity.
The critical exponents x0 and γ can be computed by analyzing the transfer matrix acting
in spatial direction [437,438]
Th,h′ =
 qh
if |h −h′| ≤1
0
otherwise
,
(290)
where h, h′ ≥0.
Steady state properties can be derived from the eigenvector φ that
corresponds to the largest eigenvalue µ of the transfer matrix P∞
h′=0 Th,h′φh′ = µ φh.
From the squares of the eigenvector components various steady state quantities can be
derived. For example, the probability ρh of ﬁnding the interface at height h is given by
ρh = φ2
h/ P
h′ φ2
h′. Close to criticality, where ǫ = 1 −q is small, one can carry out the
continuum limit φh →φ(˜h), replacing the discrete heights h by real-valued heights ˜h.
Then, the above eigenvalue problem turns into a diﬀerential equation [437] which, to
leading order in ǫ, is given by
 ∂2
∂˜h2 + (3 −µ) −3ǫ˜h

φ(˜h) = 0 .
(291)
129

<!-- page 130 -->

Simple dimensional analysis indicates that the height variables scale as h ∼ǫ−1/3 wherefore
the width diverges as w2 ∼ǫ−2/3. Similarly one can show by elementary calculations that
ρ0 ∼ǫ. The critical exponents for p = 1 are thus given by
x0 = 1 ,
γ = 1/3 .
(292)
The generic case p̸ = 1
In the nonequilibrium case p̸ = 1 the critical exponents can only be approximated numer-
ically. For p = 0, where the interface cannot move below its actual minimum height, the
hard-core wall becomes irrelevant and the model reduces to a growth process similar to
the one discussed in Sec. 6.4, belonging to the universality class of unidirectionally coupled
directed percolation. For p = 1 the numerical estimates are consistent with the previously
derived exact results. For 0 < p < 1 a very slow crossover to a diﬀerent behavior with
critical exponents y = 1.00(3), x0 = 1.5(1), and γ = 0.41(3) is observed. Similar results
are obtained for p > 1 except for x0 which is close to 1 in this case.
It is believed that this type of nonequilibrium wetting can be modeled by the KPZ
equation (252) with an additional term for the eﬀective interaction between the wall and
the interface [439,440]
∂h(r, t)
∂t
= v + σ∇2h(r, t) + λ(∇h(r, t))2 + ζ(r, t) −∂V [h(r, t)]
∂h(r, t)
.
(293)
Dimensional analysis suggests that the width exponent described by this equation should
be given by γ = (2 −z)/(2z −2), where z = 3/2 is the dynamic exponent of the KPZ
universality class, yielding γ = 1/2. The numerical estimates of γ ≃0.41 are smaller,
presumably because of the very slow crossover to the exactly solvable case γ = 1/3. In
addition, Eq. (293) suggests that the bottom layer occupation ρ0 should be proportional
to the inverse correlation length, hence x0 = ν⊥= 1/(2z −2) = 1. However, this scaling
argument seems to hold only for p > 1, whereas for 0 < p < 1 much larger values for x0 are
observed. As shown in Ref. [436], the changing sign of the coeﬃcient λ in Eq. (293) leads
to diﬀerent universality classes in both cases, corresponding to the distinction between an
‘upper’ and a ‘lower’ wall in Ref. [440]. In fact, comparing the growth velocities of a ﬂat
and a tilted interface in absence of a wall it is found that the nonlinear term (∇h(r, t))2
is indeed non-vanishing in the (p, q) plane, except for the dashed line shown in Fig. 63.
Therefore, moving along the transition line, the sign of λ changes precisely at the integrable
point, leading to a diﬀerent exponent x0.
Nonequilibrium wetting of an attractive substrate
In the above wetting model the substrate is introduced as a hard wall at zero height.
Therefore, the model neglects interactions between the substrate in the surface layer.
Loosely speaking, the free energies of the exposed and the wetted substrate are assumed
to be equal. In order to describe more realistic situations, the model has to be generalized
by taking interactions between the substrate and the surface layer into account.
This
can be done by introducing a modiﬁed growth rate q0 at zero height [441]. Obviously the
attractive short-range interaction at the bottom layer is a surface eﬀect. Therefore, the
critical point qc remains unchanged. However, if the interaction is strong enough, the
transition becomes discontinuous. In the equilibrium case p = 1 it can be proven by using
transfer matrix methods that for q0 < 2/3 there is a ﬁrst-order wetting transition.
130

<!-- page 131 -->

0
0.5
1
1.5
0
0.2
0.4
0.6
0.8
1
1.2
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
























moving phase
tricritical point
phase coexistence
pinned phase
                                                                              


                                                                                                        



                                                                              


                                                                              


a)
d)
c)
b)
Fig. 64:
Phase coexistence in a nonequilibrium wetting model. Left: Phase diagram of the wetting model
with an attractive substrate q0 = 0.2. Right: Elimination of an island in the coexistence phase. The
ﬂat island (a) grows quickly until it reaches a triangular shape (b). Because of the KPZ nonlinearity,
the deposition rate decreases, leading to a laterally shrinking island (c). Finally the island is eliminated,
ensuring phase coexistence.
In the nonequilibrium case the morphology of the phase transition depends on the
sign of the KPZ nonlinearity.
For p > 1 the emerging picture is essentially the same
as for p = 1, although with diﬀerent critical exponents. For p < 1, however, there is a
whole region in the phase diagram where the pinned and the moving phase coexist. As
illustrated in Fig. 64, islands generated by ﬂuctuations quickly grow until they reach an
almost triangular shape. Since the KPZ nonlinearity is negative, adsorption processes at
the inclined edges of the islands are strongly suppressed, allowing the attractive interaction
to reduce the size of the island in lateral direction. This ensures the stability of the pinned
phase in parts of the phase diagram where a free interface would move away from the wall.
This model demonstrates that nonequilibrium eﬀects do not only lead to diﬀerent critical
exponents but may also change the whole phase structure of a model.
Acknowledgments
This review is based on my Habilitation thesis written at the Max-Planck-Institut f¨ur
Physik komplexer Systeme and submitted to the Freie Universit¨at Berlin in June 1999. I
would like to thank I. Peschel, who made this work possible, as well as E. Domany and
D. Mukamel, who introduced me to the ﬁeld of nonequilibrium phase transitions. I am also
grateful to V. Rittenberg, whom I owe my experience in the ﬁeld of integrable equilibrium
systems. I am indebted to U. Alon, M. Antoni, M. R. Evans, M. Falcke, Y. Y. Gold-
schmidt, M. Henkel, M. J. Howard, H. K. Janssen, A. Jim´enez-Dalmaroni, H. M. Kodu-
vely, K. Krebs, R. Livi, G. ´Odor, M. Pfannm¨uller, A. Politi, Y. Rozov, S. Ruﬀo, S. Sandow,
G. Schliecker, H. Simon, D. Stauﬀer, U. C. T¨auber, B. Wehefritz, and J. S. Weitz, who es-
sentially contributed as my collaborators to the presented results. I also express my thanks
to many other colleagues for numerous stimulating discussions. Particularly I would like
to thank P. Fulde and the Max-Planck-Institut f¨ur Physik komplexer Systeme for generous
support.
131

<!-- page 132 -->

Appendix A: Vector space notation and tensor products
A one-dimensional lattice model, whose sites i = 1, . . . , N are either occupied (si = 1)
or vacant (si = 0), can be in 2N diﬀerent states s = {s1, s2, . . . , sN}. In order to represent
the probability distribution Pt(s) as a vector in a 2N-dimensional vector space let us deﬁne
an orthogonal set of basis vectors |s⟩corresponding to the conﬁgurations of the system.
Using the representation
|0⟩=
1
0

,
|1⟩=
0
1

(A.1)
the basis vectors are given by
|s⟩= |s1⟩⊗|s2⟩⊗. . . ⊗|sN⟩,
(A.2)
where ’⊗’ denotes the tensor product of two vectors:
a1
a2

⊗
b1
b2

=




a1b1
a1b2
a2b1
a2b2



.
(A.3)
Row vectors ⟨s| are the transposed vectors of |s⟩. In this vector space the probability
distribution Pt(s) can be represented by the vector
|Pt⟩=
X
s
Pt(s) |s⟩.
(A.4)
Deﬁning the sum vector over all states
⟨1| =
X
s
⟨s| = (1, 1)⊗N = (1, 1, . . . , 1)
(A.5)
the normalization of the probability distribution can be simply expressed as ⟨1|Pt⟩= 1.
Similarly the ensemble average ⟨A(t)⟩of any observable A can be expressed as
⟨A(t)⟩= ⟨1|A|Pt⟩.
(A.6)
The empty lattice is represented by the vector |vac⟩= |0⟩⊗N. Local operators act only on
a ﬁnite number of adjacent sites. For example, a single-site operator Ai can be written as
Ai = 1 ⊗1 ⊗. . .
⊗A⊗
| {z }
i-th position
. . . ⊗1 ,
(A.7)
where 1 and A are 2 × 2 matrices. Here the tensor product of two matrices is deﬁned by
a1
a2
a3
a4

⊗
b1
b2
b3
b4

=




a1b1
a1b2
a2b1
a2b2
a1b3
a1b4
a2b3
a2b4
a3b1
a3b2
a4b1
a4b2
a3b3
a3b4
a4b3
a4b4



.
(A.8)
Similarly one can deﬁne two-site operators Bi,i+1, where B is a 4 × 4 matrix. The above
notations can be easily generalized to systems with n > 1 particle species by introducing
local vectors with n components in Eq. (A.1).
132

<!-- page 133 -->

Appendix B: Derivation of the eﬀective action
In order to derive the eﬀective action of Reggeon ﬁeld theory by integration of the
noise, we ﬁrst introduce a response ﬁeld ˜φ(x, t). This allows the δ-function to be expressed
as an oscillating integral
Z ∼
Z
Dζ P[ζ]
Z
Dφ D ˜φ I[φ, ˜φ] exp
h
i
Z
ddx dt ˜φ
 ∂tφ −D∇2φ −κφ + λφ2 −ζ
i
,
(B.1)
where I[φ, ˜φ] denotes the Jacobian. After a Wick rotation in the complex plane the φ-
dependent noise contribution can be separated by
Z ∼
Z
Dφ D ˜φ I[φ, ˜φ] exp
h
−
Z
ddx dt ˜φ
 ∂tφ −D∇2φ −κφ + λφ2i
×
Z
Dζ P[ζ] exp
Z
ddx dt ˜φζ

.
(B.2)
Because of the correlations (140) the probability distribution P[ζ] is given by
P[ζ] = f[φ] exp

−
Z
ddx dt
ζ2(x, t)
2Γφ(x, t)

,
(B.3)
where f[φ] is a (ﬁeld-dependent) normalization factor. This allows the noise to be inte-
grated
Z
Dζ P[ζ] exp
Z
ddx dt ˜φζ

= f[φ]
Z
Dζ exp
hZ
ddx dt

˜φζ −ζ2
2Γφ
i
= f[φ]
Z
Dζ exp
hZ
ddx dt
Γ
2
˜φ2φ −ζ2
2Γφ
i
= ¯f[φ] exp
Γ
2
Z
ddx dt φ˜φ2
,
(B.4)
where we used Gaussian integration of the form
Z +∞
−∞
dη
1
√2πζφ exp

˜φη −η2
2ζφ

= exp
1
2ζφ˜φ2
.
(B.5)
The resulting partition function reads
Z ∼
Z
Dφ D ˜φ I′[φ, ˜φ] exp
h
−
Z
ddx dt

˜φ [∂tφ −D∇2φ −κφ] + λ˜φφ2 −Γ
2
˜φ2φ
i
.
(B.6)
It is convenient to symmetrize the cubic terms in the partition function. To this end we
rescale the ﬁelds by
φ′ =
p
2λ/Γφ ,
˜φ′ =
p
Γ/2λ˜φ ,
Γ′ =
√
2Γλ .
(B.7)
In order to keep the action symmetrized during the RG procedure, one has to introduce
an additional coeﬃcient τ in front of the time derivative. Dropping the primes the eﬀec-
tive action S = S0 + Sint is given by the expressions (144) and (145). Alternatively the
action may directly be derived from the master equation of a contact process by introduc-
ing bosonic creation and annihilation operators (for this standard procedure we refer to
Refs. [24,25,442]).
133

<!-- page 134 -->

d
1
2
3
4
5
6
Sd
2
2π
4π
2π2
8π2/3
π3
Kd
1
π
1
2π
1
2π2
1
8π2
1
12π3
1
64π3
Vd
2
π
4π/3
π2/2
8π2/15
π3/6
Table 5:
Surface area Sd, Kd = Sd/(2π)d, and the volume Vd of a d-dimensional sphere.
Appendix C: Shell integration
The one-loop integrals in Wilson’s renormalization group approach take the form
I(k) =
1
(2π)d
Z
>
ddk′ f(k2, k·k′, k′2)
(C.1)
where ’>’ denotes integration in the momentum shell Ω(1 −l) < |k′| ≤Ω. This integral
can be written as
I(k) =
lΩd
(2π)d Sd−1
Z π
0
dθ sind−2 θ f(k2, Ω|k| cos θ, Ω2) ,
(C.2)
where Sd and Vd denote surface area and volume of a d-dimensional sphere:
Sd =
2 πd/2
Γ(d/2) ,
Vd = Sd
d .
(C.3)
We also use the notation Kd = Sd/(2π)d. For easy reference we listed some of the values
for Sd, Kd, and Vd in Table 5.
To evaluate Eq. (C.2) it is often helpful to use the formulas
Z π
0
dθ sind−2 θ
=
Sd
Sd−1
,
(C.4)
Z π
0
dθ sind−2 θ cos2 θ
=
Sd
d Sd−1
.
(C.5)
In particular, if the function f does not depend on θ the integral I(k) reduces to
I(k) =
1
(2π)d
Z
>
ddk′ f(k2, k′2) = lKdΩd f(k2, Ω2).
(C.6)
134

<!-- page 135 -->

Appendix D: One-loop integrals for directed percolation
The integral for propagator renormalization reads
JP
=
Z
>
Dk′ω′ G0(k
2 + k′, ω
2 + ω′) G0(k
2 −k′, ω
2 −ω′)
(D.1)
=
Z
>
Dk′ω′
1

D(k
2 + k′)2 −κ −iτ(ω
2 + ω′)
 
D(k
2 −k′)2 −κ −iτ(ω
2 −ω′)
 .
Denoting A± = D(k
2 ± k′)2 −κ −iτ ω
2 and integrating with respect to the pole iτω′ = A+
we obtain
JP
=
1
(2π)d+1
Z
>
ddk′
Z ∞
−∞
dω′
1
(A+ −iτω′) (A−+ iτω′)
=
2πi
(2π)d+1iτ
Z
>
ddk′
1
A+ + A−
=
1
2 (2π)d τ
Z
>
ddk′
1
1
4Dk2 + Dk′2 −κ −i
2τω
(D.2)
=
lKdΩ4
2τ

1
4Dk2 + Ω2D −κ −i
2τω
 .
The integral for vertex renormalization is given by
JV =
Z
>
DkωG2
0(k, ω) G0(−k, −ω)
=
1
(2π)d+1
Z
>
ddk
Z ∞
−∞
dω
1
(Dk2 −κ −iτω)2(Dk2 −κ + iτω) .
(D.3)
Integration with respect to the poles iτω = Dk2 −κ yields the result
JV =
2πi
(2π)d+1iτ
Z
>
ddk
1
4(Dk2 −κ)2 =
lKdΩd
4τ(Ω2D −κ)2 .
(D.4)
135

<!-- page 136 -->

Appendix E: Notations
Frequently used symbols
A, B, C, . . .
particle species
Ø
vacant site
α
roughness exponent
β
density exponent
β(k)
density exponents at diﬀerent hierarchy levels
δ
exponent for temporal decay
D
diﬀusion constant
Γ
noise amplitude
d
spatial dimension of the system
D, E, ¯D, ¯E
matrices for matrix product states
∆
distance from criticality
∆(t)
Hamming distance (damage)
γ
growth exponent
Γ
noise amplitude
hi(t)
interface height at site i
H
Hamiltonian, internal energy
< i >
set of nearest neighbors of site i
θ
critical initial slip exponent
θn
critical exponents of the order parameters Mn
Iℓ
probability for empty interval of length ℓ
J
coupling constant in equilibrium systems
κ
rate for oﬀspring production in Langevin equations
L
length of the system
λ
coeﬃcient for nonlinear term in Langevin equations
ρ(x, t)
particle density
Λ
dilatation parameter for scaling transformations
µ
coupling constant between diﬀerent levels
L
Liouville operator
Mn
order parameters for spontaneous symmetry breaking
nk(t)
density of sites at height k above the bottom layer
N
system size (total number of sites N = Ld)
ν∥
temporal scaling exponent
ν⊥
spatial scaling exponent
p
percolation probability
Pt(s)
probability to ﬁnd the system in state s at time t
Pl(t)
local persistence probability
Pg(t)
global persistence probability
q
probability for interface growth
s
state of the model
si
local state si = 0, 1, . . . of lattice site i
136

<!-- page 137 -->

S
ﬁeld-theoretic action
σ
control exponent for anomalous DP
σi
local Ising spin σi = ±1 at site i
T
temperature
v
interface velocity
ws→s′
transition rate from state s to state s′
w(N, t)
interface width
Ω
cutoﬀin momentum space
ζ(x, t)
noise ﬁeld in Langevin equations
ξ⊥
spatial correlation length
ξ∥
temporal correlation length
z
dynamic critical exponent
Z
partition sum
Abbreviations
BAWE
branching-annihilating random walk with even number of oﬀspring
CDP
compact directed percolation
DK
Domany-Kinzel (model)
DMRG
density matrix renormalization group
DP
directed percolation
DP2
directed percolation with two absorbing states
DS
damage spreading
EW
Edwards-Wilkinson
IMF
improved mean ﬁeld
IPDF
interparticle distribution function
KPZ
Kardar-Parisi-Zhang
MC
Monte Carlo
MF
mean ﬁeld
MPS
matrix product state
PNG
polynuclear growth
RG
renormalization group
RSOS
restricted solid on solid
SOS
solid on solid
137

<!-- page 138 -->

References
[1] J. W. Gibbs, The elementary principles in statistical mechanics, Scribner, New York, 1902.
[2] R. J. Baxter, Exactly solved models in statistical mechanics, Academic Press, New York, 1982.
[3] A. J. Guttmann, Asymptotic analysis of power-series expansions, in Phase transitions and critical
phenomena, edited by C. Domb and J. L. Lebowitz, volume 13, Academic Press, London, 1989.
[4] D. J. Amit,
Field theory, the renormalization group, and critical phenomena,
World Scientiﬁc,
Singapore, 1984.
[5] T. M. Liggett, Interacting particle systems, Springer, Berlin, 1985.
[6] H. Spohn, Large scale dynamics of interacting particles, Springer, Berlin, 1991.
[7] B. M. McCoy and T. T. Wu, The two-dimensional Ising model, Harvard University Press, Cambridge,
1973.
[8] A. M. Polyakov, Conformal symmetry of critical ﬂuctuations, Sov. Phys. JETP Lett. 12, 381 (1970).
[9] J. L. Cardy, Conformal Invariance, in Phase transitions and critical phenomena, edited by C. Domb
and J. L. Lebowitz, volume 11, Academic Press, New York, 1987.
[10] M. Henkel, Conformal Invariance and Critical Phenomena, Springer Verlag, Berlin, 1999.
[11] B. Schmittmann and R. K. P. Zia,
Statistical mechanics of driven diﬀusive systems,
in Phase
transitions and critical phenomena, edited by C. Domb and J. L. Lebowitz, volume 17, Academic
Press, New Yorkls, 1995.
[12] H. Risken, The Fokker-Planck equation, Springer Verlag, Berlin, 1988.
[13] S. Alexander, J. Bernasconi, W. R. Schneider, and R. Orbach,
Excitation dynamics in random
one-dimensional systems, Rev. Mod. Phys. 53, 175–198 (1981).
[14] R. J. Glauber, Time-dependent statistics of the Ising model, J. Math. Phys. 4, 294 (1963).
[15] B. U. Felderhof, Spin relaxation of the Ising chain, Rep. Math. Phys. 1, 215–234 (1971).
[16] E. Lieb, T. Schultz, and D. Mattis, Ann. Phys. 16, 407 (1961).
[17] M. B. Salamon, Phase transitions in ionic conductors, Springer-Verlag, Berlin, 1979.
[18] W. Dieterich, P. Fulde, and I. Peschel, Theoretical models for superionic conductors, Adv. Phys.
29, 527–605 (1980).
[19] J. K¨arger, K. Hahn, V. Kukla, and C. R¨odenbeck, Molek¨ule im G¨ansemarsch, Phys. Bl¨att. 54,
811–815 (1998).
[20] N. G. van Kampen,
Stochastic processes in physics and chemistry,
North-Holland, Amsterdam,
1981.
[21] M. Schreckenberg, A. Schadschneider, K. Nagel, and M. Ito, Discrete stochastic models for traﬃc
ﬂow, Phys. Rev. E 51, 2939–2949 (1995).
[22] R. Kroon, H. Fleurent, and R. Sprik, Diﬀusion-limited exciton fusion reaction in one-dimensional
tetramethylammonium manganese trichloride (TMMC), Phys. Rev. E 47, 2462–2472 (1993).
[23] L. P. Kadanoﬀand J. Swift,
Transport coeﬃcients near the critical point: a master equation
approach, Phys. Rev. 165, 310–322 (1968).
[24] M. Doi, Stochastic theory of diﬀusion-controlled reactions, J. Phys. A 9, 1479–1495 (1976).
[25] P. Grassberger and M. Scheunert, Fock-space methods for identical classical objects, Fortschr. Phys.
28, 547–578 (1980).
[26] S. Alexander and T. Holstein, Lattice diﬀusion and the Heisenberg ferromagnet, Phys. Rev. B 18,
301–302 (1978).
[27] F. C. Alcaraz, M. Droz, M. Henkel, and V. Rittenberg, Reaction-diﬀusion processes, critical dynam-
ics, and quantum chains, Ann. Phys. (N.Y.) 230, 250–302 (1994).
[28] G. M. Sch¨utz, Integrable stochastic many-body systems, Habilitationsschrift, Universit¨at Bonn, 1998.
[29] B. Derrida and M. Evans, The asymmetric exclusion model: exact results through a matrix ap-
proach, in Nonequilibrium statistical mechanics in one dimension, edited by V. Privman, Cambridge
University Press, U.K., 1997.
[30] A. S. Mikhailov, Foundations of synergetics I, Springer Verlag, Berlin, 1994.
[31] M. v. Smoluchowski,
Versuch einer mathematischen Theorie der Koagulationskinetik kolloider
L¨osungen, Z. Phys. Chem. 92, 129 (1917).
[32] V. Privman, editor, Nonequilibrium statistical mechanics in one dimensions, Cambridge University
Press, U.K., 1997.
138

<!-- page 139 -->

[33] A. L. Barab´asi and H. Stanley, Fractal concepts in surface growth, Cambridge University Press,
U.K., 1995.
[34] D. E. Wolf, M. Schreckenberg, and A. Bachem, Traﬃc and granular ﬂow, World Scientiﬁc, Singapore,
1996.
[35] W. Kinzel, Percolation structures and processes, in Ann. Isr. Phys. Soc., edited by G. Deutscher,
R. Zallen, and J. Adler, volume 5, Adam Hilger, Bristol, 1983.
[36] P. Bak, C. Tang, and K. Wiesenfeld, Self-organized criticality, Phys. Rev. Lett. 59, 381–384 (1987).
[37] G. T. Barkema, M. J. Howard, and J. L. Cardy,
Reaction-diﬀusion front for A+B to 0 in one
dimension, Phys. Rev. E 53, R2017–R2020 (1996).
[38] M. R. Evans, Y. Kafri, H. M. Koduvely, and D. Mukamel,
Phase separation and coarsening in
one-dimensional driven diﬀusive systems: local dynamics leading to long-range hamiltonians, Phys.
Rev. E 58, 2764–2778 (1998).
[39] S. K. Scott, Oscillations, waves, and chaos in chemical kinetics, Oxford University Press, U.K.,
1994.
[40] M. R. Evans, N. Rajewsky, and E. R. Speer, Exact solution of a cellular automaton for traﬃc, J.
Stat. Phys. 95, 45–96 (1999).
[41] J. L. Cardy, Scaling and renormalization in statistical physics, Cambridge University Press, Cam-
bridge, U.K., 1996.
[42] D. C. Mattis and M. L. Glasser, The uses of quantum ﬁeld theory in diﬀusion-limited reactions,
Rev. Mod. Phys. 70, 979–1001 (1998).
[43] M. J. Howard and U. C. T¨auber, ‘Real’ versus ‘imaginary’ noise in diﬀusion-limited reactions, J.
Phys. A 30, 7721–7731 (1997).
[44] B. P. Lee, Renormalization group calculation for the reaction kA→0, J. Phys. A 27, 2633–2652
(1994).
[45] G. ´Odor, A. Krikelis, G. Vesztergombi, and F. Rohrbach,
Eﬀective Monte Carlo simulation on
System-V massively parallel associative string processing architecture, in Proceedings of the 7-th
Euromicro workshop on parallel and distributed processing, Funchal (Portugal), edited by B. Werner,
IEEE Computer society press, Los Alamitos, 1999.
[46] H. Simon, Concentration for one and two species one-dimensional reaction-diﬀusion systems, J.
Phys. A 28, 6585–6603 (1995).
[47] J. E. Santos, The duality relation between Glauber dynamics and the diﬀusion-annihilation model
as a similarity transformation, J. Phys. A 30, 3249–3259 (1997).
[48] M. Henkel, E. Orlandini, and J. Santos,
Reaction-diﬀusion processes from equivalent integrable
quantum chains, Ann. Phys. (N.Y.) 259, 163–231 (1997).
[49] P. Martin, String-like lattice models and Hecke algebras, J. Phys. A 22, 3103–3112 (1989).
[50] P. Martin and V. Rittenberg, A template for quantum spin chain spectra, Int. J. Mod. Phys. A 7,
707–730 (1992).
[51] M. Henkel, E. Orlandini, and G. Sch¨utz, Equivalences between stochastic systems, J. Phys. A 28,
6335–6344 (1995).
[52] L. H. Gwa and H. Spohn, Bethe solution for the dynamical-scaling exponent of the noisy Burgers
equation, Phys. Rev. A 46, 844–854 (1992).
[53] E. Siggia, Pseudospin formulation of kinetic Ising models, Phys. Rev. B 16, 2319–2320 (1977).
[54] K. Krebs, M. Pfannm¨uller, B. Wehefritz, and H. Hinrichsen,
Finite size scaling studies of one-
dimensional Reaction-Diﬀusion Systems, J. Stat. Phys. 7, 1429–1470 (1995).
[55] H. Hinrichsen, K. Krebs, and I. Peschel, Solution of a one-dimensional Diﬀusion-Reaction Model
with spatial asymmetry, Z. Phys. B 100, 105–114 (1996).
[56] F. C. Alcaraz and V. Rittenberg,
Reaction-diﬀusion processes as physical realizations of Hecke
algebras, Phys. Lett. B 314, 377–380 (93).
[57] S. Sandow, Partially asymmetric exclusion process with open boundaries, Phys. Rev. E 50, 2660–
2667 (1994).
[58] V. Pasquier and H. Saleur, Common structures between ﬁnite systems and conformal ﬁeld theories
through quantum groups, Nucl. Phys. B 330, 523–556 (1990).
[59] H. N. V. Temperley and E. H. Lieb, Relations between the ‘percolation’ and ‘colouring’ problem
and other graph-theoretical problems associated with regular planar lattices: some exact results for
the ’percolation’ problem, Proc. Roy. Soc. A 322, 251–280 (1971).
[60] F. C. Alcaraz, Exact steady states of asymmetric diﬀusion and two-species annihilation with back
reaction from the ground state of quantum spin models, Int. J. Mod. Phys. B 8, 3349–3361 (1994).
139

<!-- page 140 -->

[61] D. ben Avraham, The method of inter-particle distribution functions for diﬀusion-reaction systems
in one dimension, Mod. Phys. Lett. B 9, 895–919 (1995).
[62] P. A. Alemany and D. ben Avraham, Inter-particle distribution functions for one-species diﬀusion-
limited annihilation, Phys. Lett. A 206, 18–25 (1995).
[63] D. ben Avraham, The coalescence process, A + A →A, and the method of interparticle distribu-
tion functions, in Nonequilibrium Statistical Mechanics in One Dimension, edited by V. Privman,
Cambridge University Press, U.K., 1997.
[64] C. R. Doering and D. ben Avraham,
Interparticle distribution functions and rate equations for
diﬀusion-limited reactions, Phys. Rev. A 38, 3035–3042 (1988).
[65] C. R. Doering and D. ben Avraham, Diﬀusion-limited coagulation in the presence of particle input:
exact results in one dimension, Phys. Rev. Lett. 62, 2563–2566 (1989).
[66] M. A. Burschka, C. R. Doering, and D. ben Avraham, Transition in the relaxation dynamics of a
reversible diﬀusion-limited reaction, Phys. Rev. Lett. 63, 700–703 (1989).
[67] D. ben Avraham, M. A. Burschka, and C. R. Doering, Statics and dynamics of a diﬀusion-limited
reaction: anomalous kinetics, nonequilibrium self-ordering, and a dynamic transition, J. Stat. Phys.
60, 695–728 (1990).
[68] C. R. Doering, M. A. Burschka, and W. Horsthemke, Fluctuations and correlations in a diﬀusion-
reaction system – exact hydrodynamics, J. Stat. Phys. 65, 953–970 (1991).
[69] C. R. Doering and M. A. Burschka, Long crossover times in a ﬁnite system, Phys. Rev. Lett. 64,
245–248 (1990).
[70] H. Hinrichsen, V. Rittenberg, and H. Simon,
Universality properties of the stationary states in
the one-dimensional coagulation-diﬀusion model with external particle input,
J. Stat. Phys. 86,
1203–1236 (1997).
[71] R. Kroon and R. Sprik, Diﬀusion-limited exciton kinetics in one-dimensional systems, in Nonequi-
librium statistical mechanics in one dimensions, edited by V. Privman, Cambridge University Press,
U.K., 1997.
[72] J. Prasad and R. Kopelman, Molecular reaction kinetics inside channel pores – delayed ﬂuorescence
of naphthalene in methanol, Chem. Phys. Lett. 157, 535–538 (1989).
[73] R. Kopelman, C. S. Li, and Z. Y. Shi, One-dimensional exciton fusion kinetics in dilute polymer
blends, J. Luminescence 45, 40–42 (1990).
[74] B. Derrida, E. Domany, and D. Mukamel, An exact solution of the asymmetric exclusion model with
open boundaries, J. Stat. Phys. 69, 667–687 (1992).
[75] G. Sch¨utz and E. Domany, Phase transition in an exactly soluble one-dimensional exclusion process,
J. Stat. Phys. 72, 277–296 (1993).
[76] B. Derrida, M. R. Evans, V. Hakim, and V. Pasquier, Exact solution of a one-dimensional asymmetric
exclusion model using a matrix formulation, J. Phys. A 26, 1493–1517 (1993).
[77] F. H. L. Essler and V. Rittenberg, Representation of the quadratic algebra and partially asymmetric
diﬀusion with open boundaries, J. Phys. A 29, 3375–3407 (1996).
[78] K. Mallick and S. Sandow, Finite-dimensional representations of the quadratic algebra: Applications
to the exclusion process, J. Phys. A 30, 4513–4526 (1997).
[79] A. Schadschneider and M. Schreckenberg, Cellular automaton models and traﬃc ﬂow, J. Phys. A
26, L679–L683 (1993).
[80] G. M. Sch¨utz,
Time-dependent correlation functions in a one-dimensional asymmetric exclusion
process, Phys. Rev. E 47, 4265–4277 (1993).
[81] S. Yukawa, M. Kikuchi, and S. Tadaki, Dynamical phase transition in one-dimensional traﬃc ﬂow
model with blockage, J. Phys. Soc. Jap. 63, 3609–3618 (1994).
[82] H. Hinrichsen, Matrix product ground states for exclusion processes with parallel dynamics, J. Phys.
A 29, 3659–3567 (96).
[83] N. Rajewsky, A. Schadschneider, and M. Schreckenberg,
The asymmetric exclusion model with
sequential update, J. Phys. A 29, L305–L309 (1996).
[84] A. Honecker and I. Peschel, Matrix product states for a one-dimensional lattice gas with parallel
dynamics, J. Stat. Phys. 88, 319–345 (1997).
[85] N. Rajewsky, L. Santen, A. Schadschneider, and M. Schreckenberg,
The asymmetric exclusion
process: comparison of update procedures, J. Stat. Phys. 92, 151–194 (1998).
[86] H. Hinrichsen, S. Sandow, and I. Peschel, On matrix product ground states for reaction-diﬀusion
models, J. Phys. A 29, 2643–2649 (1996).
140

<!-- page 141 -->

[87] I. Peschel, V. Rittenberg, and U. Schultze, Spectra of quantum chains without the Yang-Baxter
equation, Nucl. Phys. B 430, 633–655 (1994).
[88] A. Kl¨umper, A. Schadschneider, and J. Zittartz, Groundstate properties of a generalized VBS-model,
Z. Phys. B 87, 281–287 (1992).
[89] M. Fannes, B. Nachteraele, and R. F. Werner, Exact antiferromagnetic ground states of quantum
spin chains, Europhys. Lett. 10, 633–637 (1989).
[90] A. Kl¨umper, A. Schadschneider, and J. Zittartz, Matrix product ground states for one-dimensional
spin-1 quantum antiferromagnets, Europhys. Lett. 24, 293–297 (1993).
[91] M. R. Evans, D. P. Foster, C. Godr`eche, and D. Mukamel, Asymmetric exclusion model with two
species: spontaneous symmetry breaking, J. Stat. Phys. 80, 69–102 (1995).
[92] H. Hinrichsen and S. Sandow,
Deterministic exclusion process with a stochastic defect: matrix-
product ground states, J. Phys. A 30, 2745–2756 (1997).
[93] K. Mallick, Shocks in the asymmetric exclusion model with an impurity, J. Phys. A 29, 5375–5386
(1996).
[94] B. Derrida and M. R. Evans, Bethe ansatz solution for a defect particle in the asymmetric exclusion
process, J. Phys. A 32, 4833–4850 (1999).
[95] F. C. Alcaraz, S. Dasmahapatra, and V. Rittenberg, N-species stochastic models with boundaries
and quadratic algebras, J. Phys. A 31, 845–878 (1998).
[96] V. Karimipour, Multispecies asymmetric simple exclusion process and its relation to traﬃc ﬂow,
Phys. Rev. E 59, 205–212 (1999).
[97] R. Stinchcombe and G. M. Sch¨utz, Application of operator algebras to stochastic dynamics and the
Heisenberg chain, Phys. Rev. Lett. 75, 140–143 (1995).
[98] G. M. Sch¨utz, Dynamic matrix ansatz for integrable reaction-diﬀusion processes, Euro. Phys. J. B
5, 589–597 (1998).
[99] S. Sandow and K. Krebs, Matrix product eigenstates for one-dimensional stochastic models and
quantum spin chains, J. Phys. A 30, 3165–3173 (1997).
[100] M. Kaulke and I. Peschel, A DMRG study of the q-symmetric Heisenberg chain, Euro. Phys. J. B
5, 727–734 (1998).
[101] D. Mollison, Spatial contact models for ecological and epidemic spread, J. Roy. Stat. Soc. B 39,
283 (1977).
[102] E. V. Albano, Critical behaviour of a forest ﬁre model with immune trees, J. Phys. A 27, L881–L886
(1994).
[103] S. Havlin and D. ben Avraham, Diﬀusion in disordered media, Adv. Phys. 36, 695–798 (1987).
[104] J.-P. Bouchaud and A. Georges, Anomalous diﬀusion in disordered media: statistical mechanics,
models and physical applications, Phys. Rep. 195, 127–293 (1990).
[105] W. Kinzel, Phase transitions of cellular automata, Z. Phys. B 58, 229–244 (1985).
[106] P. Grassberger,
Directed percolation: results and open problems,
in Nonlinearities in complex
systems, proceedings of the 1995 Shimla conference on complex systems, edited by S. Puri et al.,
Narosa Publishing, New Dehli, 1997.
[107] J. Marro and R. Dickman, Nonequilibrium phase transitions in lattice models, Cambridge University
Press, Cambridge, 1999.
[108] J. W. Essam, Percolation theory, Rep. Prog. Phys. 43, 833–912 (1980).
[109] D. Stauﬀer and A. Aharony, Introduction to Percolation Theory, Taylor & Francis, London, 1992.
[110] S. R. Broadbent and J. M. Hammersley, Proc. Camb. Phil. Soc. 53, 629 (1957).
[111] T. E. Harris, Contact interactions on a lattice, Ann. Prob. 2, 969 (1974).
[112] S. Wolfram, Statistical mechanics of cellular automata, Rev. Mod. Phys. 55, 601–644 (1983).
[113] E. Domany and W. Kinzel, Equivalence of cellular automata to Ising models and directed percolation,
Phys. Rev. Lett. 53, 311–314 (1984).
[114] G. F. Zebende and T. J. P. Penna, The Domany-Kinzel cellular automaton phase diagram, J. Stat.
Phys. 74, 1273–1279 (1994).
[115] A. Y. Tretyakov and N. Inui, Critical behaviour for mixed site-bond directed percolation, J. Phys.
A 28, 3985–3990 (1995).
[116] I. Jensen, Temporally disordered bond percolation on the directed square lattice, Phys. Rev. Lett.
77, 4988–4991 (1996).
141

<!-- page 142 -->

[117] J. W. Essam, Directed compact percolation: cluster size and hyperscaling, J. Phys. A 22, 4927–4937
(1989).
[118] R. Dickman and I. Jensen, Time-dependent perturbation theory for nonequilibrium lattice models,
Phys. Rev. Lett. 67, 2391–2394 (1991).
[119] I. Jensen and R. Dickman, Time-dependent perturbation theory for nonequilibrium lattice models,
J. Stat. Phys. 71, 89–127 (1993).
[120] R. Dickman and J. K. da Silva, Moment ratios for absorbing-state phase transitions, Phys. Rev. E
58, 4266–4270 (1998).
[121] R. M. Ziﬀ, E. Gulari, and Y. Barshad, Kinetic phase transitions in an irreversible surface-reaction
model, Phys. Rev. Lett. 56, 2553–2556 (1986).
[122] G. Grinstein, Z. W. Lai, and D. A. Browne,
Critical phenomena in a nonequilibrium model of
heterogeneous catalysis, Phys. Rev. A 40, 4820–4823 (1989).
[123] P. Meakin, Simple models for heterogeneous catalysis with a poisoning transition, J. Chem. Phys.
93, 2903–2910 (1990).
[124] I. Jensen, H. C. Fogedby, and R. Dickman, Critical exponents for an irreversible surface reaction
model, Phys. Rev. A 41, 3411–3414 (1990).
[125] C. A. Voigt and R. M. Ziﬀ,
Epidemic analysis of the second-order transition in the Ziﬀ-Gulari-
Barshad surface-reaction model, Phys. Rev. E 56, R6241–R6244 (1997).
[126] T. Aukrust, D. A. Browne, and I. Webman, Critical behavior of an autocatalytic reaction model,
Phys. Rev. E 41, 5294–5301 (1990).
[127] P. Grassberger and A. de la Torre, Reggeon ﬁeld theory (Schl¨ogl’s ﬁrst model) on a lattice; Monte
Carlo calculations of critical behaviour, Ann. Phys. (N.Y.) 122, 373–396 (1979).
[128] M. A. Mu˜noz, G. Grinstein, and Y. Tu,
Survival probability and ﬁeld theory in systems with
absorbing states, Phys. Rev. E 56, 5101–5105 (1997).
[129] J. F. F. Mendes, R. Dickman, M. Henkel, and M. C. Marques, Generalized scaling for models with
multiple absorbing states, J. Phys. A 27, 3019–3028 (1994).
[130] F. van Wijland, K. Oerding, and H. J. Hilhorst,
Wilson renormalization of a reaction-diﬀusion
process, Physica A 251, 179–201 (1998).
[131] R. Dickman and M. Burschka, Nonequilibrium critical poisoning in a single-species model, Phys.
Lett. A 127, 132–137 (1988).
[132] P. Grassberger, On the critical behavior of the general epidemic process and dynamical percolation,
Math. Biosci. 62, 157–172 (1982).
[133] J. Nahmias, H. T´ephany, and E. Guyon,
Propagation of combustion on a heterogeneous two-
dimensional network, Revue de Physique Appliqu´e 24, 773–777 (1989).
[134] E. V. Albano, Spreading analysis and ﬁnite-size scaling study of the critical behavior of a forest ﬁre
model with immune trees, Physica A 216, 213–216 (1995).
[135] F. Schl¨ogl, Chemical reaction models for non-equilibrium phase transitions, Z. Phys. A 253, 147–161
(1972).
[136] D. A. Browne and P. Kleban, Critical phenomena in nonequilibrium phase transitions, Appl. Phys.
A 51, 194–202 (1990).
[137] H. Park and H. Park, Critical behavior of an absorbing phase transition in an interacting monomer-
dimer model, Physica A 221, 97–103 (1995).
[138] H. Takayasu and A. Y. Tretyakov, Extinction, survival, and dynamical phase transition of branching
annihilating random walk, Phys. Rev. Lett. 68, 3060–3063 (1992).
[139] I. Jensen, Critical behavior of branching annihilating random walks with an odd number of oﬀspring,
Phys. Rev. E 47, R1–R4 (1993).
[140] E. V. Albano, M. Hoyuelos, and H. O. Martin, One-dimensional generalized branching annihilating
random walker process with stochastic generation of oﬀsprings, Physica A 239, 531–541 (1997).
[141] J. Kert´esz and D. E. Wolf,
Anomalous roughening in growth processes,
Phys. Rev. Lett. 62,
2571–2574 (89).
[142] T. Nagatani,
Fractal structure and pattern in 2-species growth-model – a generalized directed
percolation, J. Phys. Soc. Jap. 62, 1085–1088 (1992).
[143] J. Rolf, T. Bohr, and M. H. Jensen, Directed percolation universality in asynchronous evolution of
spatio-temporal intermittency, Phys. Rev. E 57, R2503–R2506 (1998).
142

<!-- page 143 -->

[144] H. K. Janssen, On the nonequilibrium phase transition in reaction-diﬀusion systems with an absorb-
ing stationary state, Z. Phys. B 42, 151–154 (1981).
[145] P. Grassberger, On phase transitions in Schl¨ogl’s second model, Z. Phys. B 47, 365–374 (1982).
[146] E. V. Albano, Irreversible phase transitions into non-unique absorbing states in a multicomponent
reaction system, Physica A 214, 426–434 (1995).
[147] M. A. Mu˜noz, G. Grinstein, R. Dickman, and R. Livi,
Critical behavior of systems with many
absorbing states, Phys. Rev. Lett. 76, 451–454 (1996).
[148] M. A. Mu˜noz, G. Grinstein, R. Dickman, and R. Livi, Inﬁnite numbers of absorbing states: critical
behavior, Physica D 103, 485–490 (1997).
[149] M. A. Mu˜noz, G. Grinstein, and R. Dickman, Phase structure of systems with inﬁnite numbers of
absorbing states, J. Stat. Phys. 91, 541–569 (1998).
[150] P. Meakin and D. J. Scalapino, Simple models for heterogeneous catalysis: phase transition-like
behavior in nonequilibrium systems, J. Chem. Phys. 87, 731–741 (1987).
[151] H. Park, J. K¨ohler, I. M. Kim, D. ben Avraham, and S. Redner,
Excluded volume eﬀects in
heterogeneous catalysis – reactions between dollars and dimes, J. Phys. A 26, 2071–2079 (1993).
[152] M. Hoyuelos, E. V. Albano, and H. O. M´artin, A multilayer contact process, J. Phys. A 30, 431–443
(1997).
[153] J. Zhuo, S. Redner, and H. Park, Critical behavior of an interacting surface-reaction model, J. Phys.
A 26, 4197–4213 (1993).
[154] M. Henkel, Local scale invariance and strongly anisotropic equilibrium critical systems, Phys. Rev.
Lett. 78, 1940–1943 (1997).
[155] E. Ben-Naim and P. L. Krapivsky, Cluster approximation for the contact process, J. Phys. A 27,
L481–L487 (1994).
[156] G. ´Odor, Estimation of the order parameter exponent of critical cellular automata using the enhanced
coherent anomaly method, Phys. Rev. E 51, 6261–6262 (1995).
[157] M. Henkel and G. Sch¨utz, Finite-lattice extrapolation algorithms, J. Phys. A 21, 2617–2633 (1988).
[158] P. Grassberger, Directed percolation in 2+1 dimensions, J. Phys. A 22, 3673–3679 (1989).
[159] D. ben Avraham, R. Bidaux, and L. S. Schulman, Transfer-matrix methods and results for directed
percolation, Phys. Rev. A 43, 7093–7096 (1991).
[160] S. R. White, Density matrix formulation for quantum renormalization groups, Phys. Rev. Lett. 69,
2863–2866 (1992).
[161] J. Hieida, Application of the density matrix renormalization group method to a non-equilibrium
problem, J. Phys. Soc. Jap. 67, 369–372 (1998).
[162] E. Carlon, M. Henkel, and U. Schollw¨ock,
Density matrix renormalization group and reaction-
diﬀusion processes, Euro. Phys. J. B 12, 99–106 (1999).
[163] E. Carlon, M. Henkel, and U. Schollw¨ock, Critical properties of the reaction-diﬀusion model 2A →
3A, 2A →Ø, eprint cond-mat/9912347, to appear in Phys. Rev. E.
[164] J. W. Essam, K. de’Bell, J. Adler, and F. M. Bhatti, Analysis of extended series for bond percolation
on the directed square lattice, Phys. Rev. B 33, 1982–1986 (1986).
[165] J. W. Essam, A. J. Guttmann, and K. De’Bell, On two-dimensional directed percolation, J. Phys.
A 21, 3815–3832 (1988).
[166] I. Jensen and R. Dickman, Time-dependent perturbation theory for diﬀusive non-equilibrium lattice
models, J. Phys. A 26, L151–L157 (1993).
[167] I. Jensen and A. J. Guttmann, Series expansions of the percolation probability for directed square
and honeycomb lattices, J. Phys. A 28, 4813–4833 (1995).
[168] I. Jensen, Low-density series expansions for directed percolation, I. A new eﬃcient algorithm with
applications to the square lattice, J. Phys. A 32, 5233–5249 (1999).
[169] I. Jensen and A. J. Guttmann, Series expansions for two-dimensional directed percolation, Nucl.
Phys. B 47, 835–837 (1996).
[170] I. Jensen, Critical behavior of the three-dimensional contact process, Phys. Rev. E 45, R563–R566
(1992).
[171] J. B. Bronzan and J. W. Dash, Higher order epsilon terms in the renormalization group approach
to Reggeon ﬁeld theory, Phys. Lett. B 51, 496–498 (1974).
[172] H. D. I. Abarbanel and J. B. Bronzan, Structure of the Pomeranchuk singularity in Reggeon ﬁeld
theory, Phys. Rev. D 9, 2397–2410 (1974).
143

<!-- page 144 -->

[173] H. D. I. Abarbanel, J. Bartels, J. B. Bronzan, and D. Sidhu, Inclusive cross sections and angular
distributions in Reggeon ﬁeld theory, Phys. Rev. D 12, 2798–2814 (1975).
[174] H. D. I. Abarbanel, J. B. Bronzan, R. L. Sugar, and A. R. White, Reggeon ﬁeld theory: formulation
and use, Phys. Rep. C 21, 121–182 (1975).
[175] H. D. I. Abarbanel, J. B. Bronzan, A. Schwimmer, and R. L. Sugar, Intercept of the Pomeron, Phys.
Rev. D 14, 632–646 (1976).
[176] W. R. Frazer, H. Hoﬀman, J. R. Fulco, and R. L. Sugar, Representation of the Pomeron propagator,
Phys. Rev. D 14, 2387–2401 (1976).
[177] R. C. Brower, M. A. Furman, and M. Moshe, Critical exponents for the Reggeon quantum spin
model, Phys. Lett. B 76, 213–219 (1978).
[178] M. Moshe, Recent developments in Reggeon ﬁeld theory, Phys. Rep. C 37, 257–341 (1978).
[179] P. Grassberger and K. Sundermeyer, Reggeon ﬁeld theory and Markov processes, Phys. Lett. B 77,
220–222 (1978).
[180] S. P. Obukhov, The problem of directed percolation, Physica A 101, 145–155 (1980).
[181] J. L. Cardy and R. L. Sugar,
Directed percolation and Reggeon ﬁeld theory,
J. Phys. A 13,
L423–L427 (1980).
[182] V. N. Gribov, Sov. Phys. JETP 26, 414 (1968).
[183] V. N. Gribov and A. A. Migdal, Sov. Phys. JETP 28, 784 (1969).
[184] G. Grinstein, M. A. Mu˜noz, and Y. H. Tu, Phase structure of systems with multiplicative noise,
Phys. Rev. Lett. 76, 4376–4379 (1996).
[185] M. A. Mu˜noz, Nature of diﬀerent types of absorbing states, Phys. Rev. E 57, 1377–1383 (1998).
[186] G. P. M. Kardar and Y.-C. Zhang, Dynamic scaling of growing interfaces, Phys. Rev. Lett. 56,
889–892 (1986).
[187] R. Dickman, Numerical study of a ﬁeld theory for directed percolation, Phys. Rev. E 50, 4404–4409
(1994).
[188] G. Huber, M. H. Jensen, and K. Sneppen, Distributions of self-interactions and voids in (1+1)-
dimensional directed percolation, Phys. Rev. E 52, R2133–R2136 (1995).
[189] K. G. Wilson and J. Kogut, The renormalization group and the ǫ expansion, Phys. Rep. 12 C,
75–200 (1974).
[190] V. Brunel, K. Oerding, and F. van Wijland, Fermionic ﬁeld theory for directed percolation in (1+1)
dimensions, J. Phys. A 33, 1085–1097 (2000).
[191] F. Igl´oi, I. Peschel, and L. Turban, Inhomogeneous systems with unusual critical behaviour, Adv.
Phys. 42, 683–740 (1993).
[192] J. L. Cardy, Critical behaviour at an edge, J. Phys. A 16, 3617–3628 (1983).
[193] J. W. Essam and D. TanlaKishani, Directed compact percolation near a wall: I. Biased growth, J.
Phys. A 27, 3743–3750 (1994).
[194] J. W. Essam and A. J. Guttmann, Directed compact percolation near a wall: II. Cluster length and
size, J. Phys. A 28, 3591–3598 (1995).
[195] H. K. Janssen, B. Schaub, and B. Schmittmann, The eﬀects of surfaces on directed percolation and
related stochastic evolution processes, Z. Phys. B 72, 111–114 (1988).
[196] J. W. Essam, A. J. Guttmann, I. Jensen, and D. TanlaKishani, Directed percolation near a wall, J.
Phys. A 29, 1619–1628 (1996).
[197] K. B. Lauritsen, K. Sneppen, M. Markosov´a, and M. H. Jensen,
Directed percolation with an
absorbing boundary, Physica A 247, 1–9 (1997).
[198] P. Fr¨ojdh, M. Howard, and K. B. Lauritsen, Directed percolation with a wall or edge, J. Phys. A
31, 2311–2320 (1998).
[199] M. Howard, P. Frojdh, and K. B. Lauritsen,
Surface Critical Behavior in Systems with Non-
Equilibrium Phase Transitions, Phys. Rev. E 61, 167–183 (2000).
[200] I. Jensen, Low-density series expansions for directed percolation II: The square lattice with a wall,
J. Phys. A 32, 6055–6062 (1999).
[201] E. V. Albano, Displacement of inactive phases by the reactive regime in a lattice gas model for a
dimer-monomer irreversible surface reaction, Phys. Rev. E 55, 7144–7152 (1997).
[202] C. Kaiser and L. Turban, Fractal dimensions of conﬁned clusters in two-dimensional directed per-
colation, J. Phys. A 27, L579–L583 (1994).
144

<!-- page 145 -->

[203] C. Kaiser and L. Turban, Surface shape and local critical behaviour in two-dimensional directed
percolation, J. Phys. A 28, 351–363 (1995).
[204] H. Hinrichsen and G. ´Odor, Correlated initial conditions in directed percolation, Phys. Rev. E 58,
311–317 (1998).
[205] P. A. Alemany and D. H. Zanette, Eﬀect of an initial fractal distribution of particles in coagulation
and annihilation reaction-diﬀusion models, Chaos, Solitons & Fractals 6, 11–15 (1995).
[206] A. J. Bray, K. Humayun, and T. J. Newman, Kinetics of ordering for correlated initial conditions,
Phys. Rev. B 43, 3699–3702 (1991).
[207] B. Derrida, A. J. Bray, and C. Godr`eche, Non-trivial exponents in the zero temperature dynamics
of the 1d Ising and Potts models, J. Phys. A 27, L357–L361 (1994).
[208] A. J. Bray, B. Derrida, and C. Godr`eche, Non-trivial algebraic decay in a soluble model of coarsening,
Europhys. Lett. 27, 175–180 (1994).
[209] J. L. Cardy, Proportion of unaﬀected sites in a reaction-diﬀusion process, J. Phys. A 28, L19–L24
(1995).
[210] B. Derrida, V. Hakim, and V. Pasquier, Exact exponent for the number of persistent spins in the
zero-temperature dynamics of the one-dimensional Potts model, J. Stat. Phys. 85, 763–797 (1996).
[211] S. N. Majumdar, C. Sire, A. J. Bray, and S. J. Cornell, Nontrivial exponent for simple diﬀusion,
Phys. Rev. Lett. 77, 2867–2870 (1996).
[212] B. Derrida, V. Hakim, and R. Zeitak, Persistent spins in the linear diﬀusion approximation of phase
ordering and zeros of stationary Gaussian processes, Phys. Rev. Lett. 77, 2871–2874 (1996).
[213] N. Menyh´ard and G. ´Odor, Non-Markovian persistence at the parity conserving point of a one-
dimensional nonequilibrium kinetic Ising model, J. Phys. A 30, 8515–8521 (1997).
[214] S. N. Majumdar, A. J. Bray, S. J. Cornell, and C. Sire, Global persistence exponent for nonequilib-
rium critical dynamics, Phys. Rev. Lett. 77, 3704–3407 (1996).
[215] L. Sch¨ulke and B. Zheng, Monte Carlo measurement of the global persistence exponent, Phys. Lett.
A 233, 93–98 (1997).
[216] B. P. Lee and A. D. Rutenberg, Persistence, poisoning, and autocorrelations in dilute coarsening,
Phys. Rev. Lett. 79, 4842–4845 (1997).
[217] B. Derrida, How to extract information from simulations of coarsening at ﬁnite temperature, Phys.
Rev. E 55, 3705–3707 (1997).
[218] D. Stauﬀer, Universality of Derrida coarsening in Ising models, Int. J. Mod. Phys. C 8, 361–364
(1997).
[219] H. Hinrichsen and M. Antoni, Identiﬁcation of domain walls in coarsening systems at ﬁnite temper-
ature, Phys. Rev. E 57, 2650–2655 (1998).
[220] H. Hinrichsen and H. M. Koduvely, Numerical study of local and global persistence in directed
percolation, Euro. Phys. J. B 5, 257–264 (1998).
[221] K. Oerding and F. van Wijland, Global persistence in directed percolation, J. Phys. A 31, 7011–7021
(1998).
[222] H. K. Janssen, Renormalized ﬁeld theory of the Gribov process with quenched disorder, Phys. Rev.
E 55, 6253–6256 (1997).
[223] A. G. Moreira and R. Dickman, Critical dynamics of the contact process with quenched disorder,
Phys. Rev. E 54, R3090–R3093 (1996).
[224] R. Caﬁero, A. Gabrielli, and M. A. Mu˜noz, Disordered one-dimensional contact process, Phys. Rev.
E 57, 5060–5068 (1998).
[225] I. Webman, D. ben Avraham, A. Cohen, and S. Havlin, Dynamical phase transitions in a random
environment, Phil. Mag. B 77, 1401–1412 (1998).
[226] H. Chat´e and P. Manneville, Spatio-temporal intermittency in coupled map lattices, Physica D 32,
409–422 (1988).
[227] Y. Pomeau, Front motion, metastability and subcritical bifurcations in hydrodynamics, Physica D
23, 3–11 (1986).
[228] M. Houlrik, I. Webman, and M. H. Jensen, Mean-ﬁeld theory and critical behavior of coupled map
lattices, Phys. Rev. E 41, 4210–4222 (1990).
[229] Y. Cuche, R. Livi, and A. Politi, Phase transitions in 2D linearly stable coupled map lattices, Physica
D 103, 369–380 (1997).
145

<!-- page 146 -->

[230] D. K. Arrowsmith and J. W. Essam,
Extension of the Kasteleyn-Fortuin formulas to directed
percolation, Phys. Rev. Lett. 65, 3068–3071 (1990).
[231] D. K. Arrowsmith, P. Mason, and J. W. Essam, Vicious walkers, ﬂows and directed percolation,
Physica A 77, 267–272 (1991).
[232] J. Cardy and F. Colaiori, Directed percolation and generalized friendly random walkers, Phys. Rev.
Lett. 82, 2232–2235 (1999).
[233] T. Tsuchiya and M. Katori, Chiral Potts models, friendly walkers and directed percolation problem,
J. Stat. Phys. 67, 1655–1666 (1998).
[234] C. M. Fortuin and P. W. Kasteleyn, On the random-cluster model. I. Introduction and relation to
other models, Physica A 57, 536–564 (1972).
[235] A. Hansen and S. Roux, Application of ‘logical transport’ to determine the directed and isotropic
percolation thresholds, J. Phys. A 20, L873–L878 (1987).
[236] P. Grassberger and Y.-C. Zhang, ‘Self-organized’ formulation of standard percolation phenomena,
Physica A 224, 169–179 (1996).
[237] S. Maslov and Y.-C. Zhang, Self-organized critical directed percolation, Physica A 223, 1–6 (1996).
[238] R. Dickman, M. A. Mu˜noz, A. Vespinani, and S. Zapperi, Paths to self-organized criticality, Braz.
J. Phys. 30, 27–41 (2000).
[239] H. M. Br¨oker and P. Grassberger, SOC in a population model with global control, Physica A 267,
453–470 (1999).
[240] J. Cardy, Renormalisation group theory of branching potts interfaces, 1998, unpublished, eprint
cond-mat/9806098.
[241] J. K¨ohler and D. ben Avraham, The dimer-trimer model for heterogeneous catalysis: a new univer-
sality class, J. Phys. A 24, L621–L627 (1991).
[242] I. Jensen, Critical behavior of nonequilibrium models with inﬁnitely many absorbing state, Int. J.
Mod. Phys. B 8, 3299–3331 (1994).
[243] K. Yaldram, K. M. Khan, N. Ahmed, and M. A. Khan,
A new class of critical exponents in a
dimer-monomer catalytic surface reaction, J. Phys. A 26, L801–L805 (1993).
[244] I. Jensen, Critical behaviour of a surface reaction model with inﬁnitely many absorbing states, J.
Phys. A 27, L61–L68 (1994).
[245] I. Jensen, Critical behavior of the pair contact process, Phys. Rev. Lett. 70, 1465–1468 (1993).
[246] I. Jensen and R. Dickman, Nonequilibrium phase transitions in systems with inﬁnitely many ab-
sorbing states, Phys. Rev. E 48, 1710–1725 (1993).
[247] C. L´opez and M. A. Mu˜noz, Numerical analysis of a Langevin equation for systems with inﬁnite
absorbing states, Phys. Rev. E 56, 4864–4867 (1997).
[248] P. Grassberger, H. Chat´e, and G. Rousseau, Spreading in media with long-time memory, Phys. Rev.
E 55, 2488–2495 (1997).
[249] P. Grassberger, On the critical behavior of the general epidemic process and dynamical percolation,
Math. Biosc. 63, 157–172 (1983).
[250] J. L. Cardy and P. Grassberger,
Epidemic models and percolation,
J. Phys. A 18, L267–L271
(1985).
[251] H. K. Janssen, Renormalized ﬁeld theory of dynamical percolation, Z. Phys. B 58, 311–317 (1985).
[252] R. Durrett, Lecture Notes on Particle Systems and Percolation, Wadsworth, Paciﬁc Grove, CA,
1988.
[253] R. Dickman and A. Y. Tretyakov, Hyperscaling in the Domany-Kinzel cellular automaton, Phys.
Rev. E 52, 3218–3220 (1995).
[254] E. Frey, U. C. T¨auber, and F. Schwabl, Crossover from isotropic to directed percolation, Phys. Rev.
E 49, 5058–5072 (1994).
[255] V. L. Nguyen and A. Rubio, Percolation in 2-D random systems: from isotropic to directed perco-
lation, Solid State Comm. 95, 833–840 (1995).
[256] J. F. F. Mendes, R. Dickman, and H. Herrmann, Crossover from directed percolation to compact
directed percolation, Phys. Rev. E 54, R3071–R3074 (1996).
[257] P. Fr¨ojdh and M. den Nijs, Crossover from isotropic to directed percolation, Phys. Rev. Lett. 78,
1850–1853 (1997).
[258] H. Hinrichsen, On possible experimental realizations of directed percolation, Braz. J. Phys. 30,
69–82 (2000).
146

<!-- page 147 -->

[259] M. Ehsasi, M. Matloch, O. Frank, J. H. Block, K. Christmann, F. S. Rys, and W. Hirschwald,
Steady and nonsteady rates of reaction in a heterogeneously catalyzed reaction: oxydation of CO on
platinum, experiments and simulations, J. Chem. Phys. 91, 4949–4960 (1989).
[260] J. Wintterlin, S. V¨olkening, T. V. W. Janssens, T. Zambelli, and G. Ertl, Atomic and macroscopic
reaction rates of a surface-catalyzed reaction, Science 278, 1931–1933 (1997).
[261] R. M. Eiswirth, K. Krischer, and G. Ertl, Nonlinear dynamics in the CO-oxidation on Pt single
crystal surface, Appl. Phys. A 51, 79–90 (1990).
[262] suggested by R. Dickman, private communication (1999).
[263] H. Hinrichsen, Stochastic lattice models with several absorbing states, Phys. Rev. E 55, 219–226
(1997).
[264] T. Zambelli, J. Wintterlin, J. Trost, and G. Ertl, Identiﬁcation of the “active sites” of a surface-
catalyzed reaction, Science 273, 1688–1690 (1996).
[265] S. Douady and A. Daerr, Avalanches and inclined chutes, in Physics of Dry Granular Media, edited
by H. J. Herrmann et al., page 339, Kluwer Academic Publishers, New York, 1998.
[266] A. Daerr and S. Douady, Two types of avalanche behaviour in granular media, Nature 399, 241–243
(1999).
[267] H. Hinrichsen, A. Jim´enez-Dalmaroni, Y. Rozov, and E. Domany, Flowing sand – a physical real-
ization of directed percoaltion, Phys. Rev. Lett. 83, 4999–5002 (1999).
[268] J. P. Bouchaud and M. E. Cates, Triangular and uphill avalanches of a tilted sandpile, Granular
Matter 1, 101–103 (1998).
[269] A. Daerr, private communication (1999).
[270] S. Dippel, G. G. Batrouni, and D. E. Wolf,
Collision-induced friction in the motion of a single
particle on a bumpy inclined line, Phys. Rev. E 54, 6845–6856 (1996).
[271] R. Hilfer, T. Rage, and B. Virgin, Local percolation probabilities for a natural sandstone, Physica
A 241, 105–110 (1997).
[272] H. Tephany, J. Nahmias, and J. A. M. S. Duarte, Combustion on heterogeneous media. A critical
phenomenon, Physica 242, 57–69 (1997).
[273] M. J. Berridge, Inositol trisphosphate and calcium signaling, Nature 361, 315–325 (1993).
[274] J. Keizer, Y.-X. Li, J. Rinzel, and S. Stojilkovic, InsP3-induced Ca2+ excitability of the endoplasmic
reticulum, Mol. Biol. Cell 6, 945–951 (1995).
[275] J. Lechleiter and D. Clapham, Molecular mechanisms of intracellular Ca2+ excitability in Xenopus
laevis oocytes, Cell 69, 283 (1992).
[276] L. F. Jaﬀe, Classes and mechanisms of calcium waves, Cell Calcium 14, 736 (1993).
[277] J. Keizer and G. D. Smith,
Spark-to-wave transition: saltatory transmission of Ca2+ waves in
cardiac myocytes, Bioph. Chem. 72, 87–100 (1998).
[278] M. B¨ar, M. Falcke, H. Levine, and L. Tsimring, Discrete stochastic modeling of calcium channel
dynamics, 2000, to appear in Phys. Rev. Lett.
[279] D. A. Huse and C. Henley, Pinning and roughening of domain walls in Ising systems due to random
impurities, Phys. Rev. Lett. 54, 2708–2711 (1985).
[280] E. Perlsman and S. Havlin, The directed polymer – directed percolation transition, Europhys. Lett.
46, 13–17 (1999).
[281] E. B. J. P. Bouchaud and G. Lapasset, Models of fractal cracks, Phys. Rev. Lett. 71, 2240–2243
(1993).
[282] P. Grassberger, in Fractals in physics, edited by L. Pietronero and E. Tosatti, Elsevier, Amsterdam,
1986.
[283] H. K. Janssen, K. Oerding, F. van Wijland, and H. J. Hilhorst, L´evy-ﬂight spreading of epidemic
processes leading to percolating clusters, Euro. Phys. J. B 7, 137–145 (1999).
[284] H. Hinrichsen and M. Howard, A model for anomalous directed percolation, Euro. Phys. J. B 7,
635–643 (1999).
[285] E. V. Albano, Branching annihilating L´evy ﬂights: Irreversible phase transitions with long-range
interactions, Europhys. Lett. 34, 97–102 (1996).
[286] S. A. Cannas, Phase diagram of a stochastic cellular automaton with long-range interactions, Physica
A 258, 32–44 (1998).
[287] G. Zumofen and J. Klafter, Reaction dynamics controlled by enhanced diﬀusion, Phys. Rev. E 50,
5119–5122 (1994).
147

<!-- page 148 -->

[288] P. Grassberger, F. Krause, and T. von der Twer, A new type of kinetic critical phenomenon, J.
Phys. A 17, L105–L109 (1984).
[289] P. Grassberger, Some further results on a kinetic critical phenomenon, J. Phys. A 22, L1103–L1107
(1989).
[290] D. Zhong and D. ben Avraham, Universality class of two-oﬀspring branching annihilating random
walks, Phys. Lett. A 209, 333–337 (1995).
[291] A. Lipowski, Monte Carlo study of a 2-species contact process, J. Phys. A 29, L355–L361 (1996).
[292] H. K. Janssen, Spontaneous symmetry breaking in directed percolation with many colors: diﬀeren-
tiation of species in the Gribov process, Phys. Rev. Lett. 78, 2890–2893 (1997).
[293] D. ben Avraham, F. Leyvraz, and S. Redner, Propagation and extinction in branching annihilating
random walks, Phys. Rev. E 50, 1843–1850 (1994).
[294] I. Jensen, Conservation laws and universality in branching annihilating random walks, J. Phys. A
26, 3921–3930 (1993).
[295] I. Jensen,
Critical exponents for branching annihilating random walks with an even number of
oﬀspring, Phys. Rev. E 50, 3623–3633 (1994).
[296] J. Cardy and U. C. T¨auber, Theory of branching and annihilating random walks, Phys. Rev. Lett.
77, 4870 (1996).
[297] J. Cardy and U. C. T¨auber, Field theory of branching and annihilating random walks, J. Stat. Phys.
90, 1–56 (1998).
[298] I. Jensen, Low-density series expansions for directed percolation on square and triangular lattices,
J. Phys. A 29, 7013–7040 (1996).
[299] M. H. Kim and H. Park, Critical behavior of an interacting monomer-dimer model, Phys. Rev. Lett.
73, 2579–2582 (1994).
[300] H. Park, M. H. Kim, and H. Park,
Dynamic scaling behavior of an interacting monomer-dimer
model, Phys. Rev. E 52, 5664–5666 (1995).
[301] W. Hwang, S. Kwon, H. Park, and H. Park, Critical phenomena of nonequilibrium dynamical systems
with two absorbing states, Phys. Rev. E 57, 6438–6540 (1998).
[302] P. Bhattacharyya, Dynamical critical properties of a one-dimensional probabilistic cellular automa-
ton, Euro. Phys. J. B 3, 247–252 (1998).
[303] N. Menyh´ard, One-dimensional non-equilibrium kinetic Ising models with branching annihilating
random walk, J. Phys. A 27, 6139–6146 (1994).
[304] N. Menyh´ard and G. ´Odor,
Non-equilibrium phase transitions in one-dimensional kinetic Ising
models, J. Phys. A 29, 7739–7755 (1996).
[305] K. S. Brown, K. E. Bassler, and D. A. Browne, Mean-ﬁeld analysis and Monte Carlo study of an
interacting two-species reaction model, Phys. Rev. E 56, 3953–3958 (1997).
[306] F. Bagnoli, N. Boccara, and P. Palmerini, Phase transitions in a probabilistic cellular automaton
with two absorbing states, in Dynamical modeling in biotechnologies, edited by F. Bagnoli, P. Lio,
and S. Ruﬀo, World Scientiﬁc, Singapore, 1998.
[307] K. E. Bassler and D. A. Browne, Nonequilibrium critical dynamics of a three species monomer-
monomer model, Phys. Rev. Lett. 77, 4094–4097 (1996).
[308] K. E. Bassler and D. A. Browne, Three-species monomer-monomer model: A mean-ﬁeld analysis
and Monte Carlo study, Phys. Rev. E 55, 5225–5233 (1997).
[309] K. B. Lauritsen, P. Fr¨ojdh, and M. Howard, Surface critical behavior in systems with absorbing
states, Phys. Rev. Lett. 81, 2104–2107 (1998).
[310] K. Mussawisade, J. E. Santos, and G. M. Sch¨utz,
Branching-annihilating random walks in one
dimension: some exact results, J. Phys. A 31, 4381–4394 (1998).
[311] P. Bak and K. Sneppen, Punctuated equilibrium and criticality in a simple model of evolution, Phys.
Rev. Lett. 71, 4083–4086 (1993).
[312] P. Bak, How Nature works, Copernicus, New York, 1996.
[313] O. Narayan and A. A. Middleton, Avalanches and the renormalization group for pinned charge-
density waves, Phys. Rev. B 49, 244–256 (1994).
[314] D. Sornette, A. Johansen, and I. Dornic, Mapping self-organized criticality onto criticality, Journal
de Physique 5, 325–335 (1995).
[315] R. Dickman, A. Vespignani, and S. Zapperi, Self-organized criticality as an absorbing-state transition,
Phys. Rev. E 57, 5095–5105 (1995).
148

<!-- page 149 -->

[316] S. S. Manna, Two-state model of self-organized criticality, J. Phys. A 24, L363–L369 (1991).
[317] D. Sornette and I. Dornic, Parallel Bak-Sneppen model and directed percolation, Phys. Rev. E 54,
3334–3338 (1996).
[318] H. Hinrichsen,
The pair contact process with diﬀusion – a new type of nonequilibrium critical
behavior?, eprint cond-mat/0001177, to appear in Phys. Rev. E.
[319] A. Lipowski, Dimensional reduction in a model with inﬁnitely many absorbing states, Phys. Rev.
E 60, R6255–R6258 (1999).
[320] S. A. Kauﬀman, J. Theor. Biol. 22, 437 (1969).
[321] S. A. Kauﬀman, Emergent properties in random complex automata, Physica D 10, 145–156 (1984).
[322] H. J. Herrmann, Damage spreading, Physica A 168, 516–528 (1990).
[323] O. Martin, Lyapunov exponents of stochastic dynamical systems, J. Stat. Phys. 41, 249–261 (1985).
[324] B. Derrida and D. Stauﬀer,
Phase transitions in two-dimensional Kauﬀman cellular automata,
Europhys. Lett. 2, 739–745 (1986).
[325] M. Creutz, Deterministic Ising dynamics, Ann. Phys. 167, 62–72 (1986).
[326] H. Stanley, D. Stauﬀer, J. Kert´esz, and H. Herrmann, Dynamics of spreading phenomena in two-
dimensional Ising model, Phys. Rev. Lett. 59, 2326–2328 (1987).
[327] B. Derrida and G. Weisbuch, Dynamical phase transitions in 3-dimensional spin glasses, Europhys.
Lett. 4, 657–662 (1987).
[328] A. M. Mariz, H. J. Herrmann, and L. de Arcangelis, Comparative study of damage spreading in
the Ising model using heat-bath, Glauber and Metropolis dynamics, J. Stat. Phys. 59, 1043–1050
(1990).
[329] S. C. Glotzer, D. Stauﬀer, and S. Sastry, Damage spreading in the Q2R Ising model, Physica A
164, 1–11 (1990).
[330] S. C. Glotzer and N. Jan, Damage spreading in the Kawasaki Ising model, Physica A 173, 325–331
(1991).
[331] N. Jan and L. de Arcangelis, Computational Aspects of Damage Spreading, in Annual Review of
Computational Physics, edited by D. Stauﬀer, World Scientiﬁc, Singapore, 1995.
[332] E. N. Miranda and N. Parga, Dynamical phase transitions in the classical Heisenberg model, J.
Phys. A 24, 1059–1064 (1991).
[333] J. Chiu and S. Teitel, The eﬀect of dynamics on damage spreading in the two-dimensional classical
XY model, J. Phys. A 23, L891–L894 (1990).
[334] A. M. Mariz, Damage spreading in the Potts and Ashkin-Teller models: exact results, J. Phys. A
23, 979–986 (1990).
[335] L. da Silva, F. A. Tamarit, and A. C. N. Magalh˜aes, Damage spreading in a two-dimensional Potts
ferromagnet in a magnetic ﬁeld, J. Phys. A 30, 2329–2344 (1997).
[336] M. F. de A. Bibiano, F. G. B. Moreira, and A. M. Mariz, Damage spreading in the q-state Potts
model, Phys. Rev. E 55, 1448–1450 (1997).
[337] A. M. Mariz, A. M. C. Souza, and C. Tsallis, Spreading of damage in the discrete N-vector ferro-
magnet: exact results, J. Phys. A 26, L1007–L1010 (1993).
[338] L. de Arcangelis, A. Coniglio, and H. J. Herrmann, Damage spreading in spin glasses, Europhys.
Lett. 9, 749–754 (1989).
[339] H. R. da Cruz, U. M. S. Costa, and E. M. F. Curado, Spreading of damage in a 3D Ising spin glass,
J. Phys. A 22, L651–L655 (1989).
[340] I. A. Campbell and L. de Arcangelis, On the damage spreading in Ising spin glasses, Physica A
178, 29–43 (1991).
[341] I. A. Campbell,
Damage spreading in ±J and Gaussian Ising spin glasses,
Europhys. Lett. 21,
959–964 (1993).
[342] I. A. Campbell and L. Bernardi, Dynamics of Ising spin glasses above the freezing temperature,
Phys. Rev. B 50, 12643–12646 (1994).
[343] F. Wang, N. Kawashima, and M. Suzuki, The dynamical critical behaviour of the ±J spin glass at
the damage-spreading threshold, Europhys. Lett. 33, 165–169 (1996).
[344] R. M. C. de Almeida, L. Bernardi, and I. A. Campbell, Damage spreading in ±J asymmetric Ising
spin glasses, Journal de Physique 5, 355 (1995).
[345] T. Wappler, T. Vojta, and M. Schreiber, Monte Carlo simulations of the dynamical behavior of the
Coulomb glass, Phys. Rev. B 55, 6272–6277 (1997).
149

<!-- page 150 -->

[346] F. G. B. Moreira, A. J. F. Souza, and A. M. Mariz, Damage spreading in the Ising model with a
microcanonical constraint, Phys. Rev. E 53, 332–334 (1996).
[347] P. Grassberger, Magnetization decay in the diluted Ising model, Int. J. Mod. Phys. 7, 89–97 (1996).
[348] E. S. de Sousa, A. M. Mariz, F. D. Nobre, and U. M. S. Costa, The two-dimensional site-diluted
Ising ferromagnet. A damage-spreading analysis, Physica A 241, 469–480 (1997).
[349] T. Vojta, Damage spreading in a random ﬁeld system, J. Phys. A 30, L643–L649 (1997).
[350] A. A. J´unior and F. D. Nobre, Damage spreading in the fully frustrated square-lattice Ising model,
Physica A 243, 58–66 (1997).
[351] A. V. Lima, M. L. Lyra, and U. M. S. Costa, Damage spreading scaling and temporal evolution near
the surface-bulk transition of Ising layered ﬁlms, J. Magn. Mag. Mat. 171, 329–335 (1997).
[352] A. Coniglio, D. Stauﬀer, and N. Jan, Search for multifractality in damage spreading for Kauﬀman
cellular automata, J. Phys. A 20, L1103–L1106 (1987).
[353] D. Stauﬀer, Hunting for the fractal dimension of the Kauﬀman model, Physica D 38, 341–344
(1989).
[354] F. Bagnoli, R. Rechtman, and S. Ruﬀo, Some facts of life (hiring cell densities), Physica A 171,
249–264 (1991).
[355] R. A. Monetti and E. V. Albano, Critical edge between frozen extinction and chaotic life, Phys.
Rev. E 52, 5825–5831 (1995).
[356] E. V. Albano, Damage spreading in the Ziﬀ-Gulariu-Barshad model, Phys. Rev. E 50, 1129–1134
(1994).
[357] I. Graham, E. Hern´andez-Garc´ia, and M. Grant, Damage spreading during domain growth, Phys.
Rev. E 49, R4763–R4766 (1994).
[358] E. V. Albano, Dynamics of damage spreading in irreversible reaction processes, Phys. Rev. Lett.
72, 108–111 (1994).
[359] E. V. Albano, Damage spreading in a single-component irreversible reaction process: Dependence
of the system’s immunity on the Euclidean dimension, J. Stat. Phys. 78, 1147–1155 (1995).
[360] E. V. Albano, Study of damage spreading in dimer-dimer surface reaction models, Physica A 215,
451–460 (1995).
[361] J. M. Kim, Y. Lee, and I.-M. Kim,
Spreading of damage in the ballistic deposition and larger
curvature models, Phys. Rev. E 54, 4603–4606 (1996).
[362] G. ´Odor and N. Menyh´ard,
Damage spreading for one-dimensional nonequilibrium models with
parity conserving phase transitions, Phys. Rev. E 57, 5168–5177 (1998).
[363] A. Bhowal, Damage spreading in the ‘sandpile’ model of SOC, Physica A 247, 327–330 (1997).
[364] E. N. Miranda and N. Parga,
Damage spreading and multifractality in the traveling salesman
problem, J. Phys. A 22, L907L912 (1989).
[365] B. Luque and R. V. Sole, Phase transitions in random networks: simple analytic determination of
critical points, Phys. Rev. E 55, 257–260 (1997).
[366] M. L. Martins, H. F. V. de Resende, C. Tsallis, and A. C. N. de Magalh˜aes, Evidence for a new
phase in the Domany-Kinzel cellular automaton, Phys. Rev. Lett. 66, 2045–2048 (1991).
[367] M. L. Martins, G. F. Zebende, T. J. P. Penna, and C. Tsallis, The Domany-Kinzel cellular automaton:
relaxation time, susceptibility and constrained dynamics, J. Phys. A 27, 1–8 (1994).
[368] A. S. H. Rieger and M. Schreckenberg, Re-entrant behaviour in the Domany-Kinzel cellular automa-
ton, J. Phys. A 27, L423–L429 (1994).
[369] P. Grassberger, Are damage spreading transitions generically in the universality class of directed
percolation?, J. Stat. Phys. 79, 13–23 (1995).
[370] G. A. Kohring and M. Schreckenberg, The Domany-Kinzel cellular automaton revisited, Journal de
Physique 2, 2033–2047 (1992).
[371] F. Bagnoli, On damage spreading transitions, J. Stat. Phys. 85, 151–164 (1996).
[372] T. Tom´e, Spreading of damage in the Domany-Kinzel cellular automaton: a mean-ﬁeld approach,
Physica A 212, 99–109 (1994).
[373] F. D. Nobre, A. M. Mariz, and E. S. Sousa, Spreading of damage: An unexpected disagreement
between the sequential and parallel updating in Monte Carlo simulations, Phys. Rev. Lett. 69, 13–16
(1992).
150

<!-- page 151 -->

[374] S. C. Glotzer, P. H. Poole, and N. Jan, Comment on ”Spreading of damage: an unexpected dis-
agreement between the sequential and parallel updatings in Monte Carlo simulations”, Phys. Rev.
Lett. 70, 2046–2047 (1993).
[375] S. S. Manna,
Damage spreading in Ising-model with increased range of interaction,
Journal de
Physique 51, 1261–1265 (1990).
[376] P. Grassberger,
Damage spreading in the Ising model with Glauber dynamics,
J. Phys. A 28,
L67–L72 (1995).
[377] T. Vojta, Damage spreading and dynamic stability of kinetic Ising models, J. Phys. A 30, L7–L13
(1997).
[378] A. M. Mariz and H. J. Herrmann, Mixing heat-bath and Glauber dynamics: damage spreading in
the Ising model, J. Phys. A 22, L1081–L1084 (1989).
[379] G. Le Ca¨er, Damage spreading in a three-dimensional Ising model in a magnetic ﬁeld, J. Phys. A
22, L647–L650 (1989).
[380] G. Le Ca¨er, Damage spreading in a ferromagnetic cubic Ising system as a function of temperature
and magnetic ﬁeld, Physica A 159, 329–342 (1989).
[381] U. Gropengiesser, Cluster distribution of damage in the 3D Ising model, Physica A 207, 492–429
(1994).
[382] H. Hinrichsen, J. S. Weitz, and E. Domany, An algorithm-independent deﬁnition of damage spreading
- application to directed percolation, J. Stat. Phys. 88, 617–636 (1997).
[383] F. Bagnoli, R. Rechtman, and S. Ruﬀo,
Damage spreading and Lyapunov exponents in cellular
automata, Phys. Lett. A 172, 34–38 (1992).
[384] S. Wang and R. Swendsen, Cluster Monte Carlo algorithms, Physica A 167, 565–579 (1990).
[385] D. Stauﬀer, Damage spreading in 3d Ising model with Swendsen-Wang dynamics, Physica A 162,
27–30 (1989).
[386] D. Stauﬀer,
Dynamics of Swendsen-Wang clusters in 2d Ising model,
Physica A 171, 471–474
(1991).
[387] H. Hinrichsen, E. Domany, and D. Stauﬀer, Damage spreading in a 2d Ising model with Swendsen-
Wang dynamics, J. Stat. Phys. 91, 807–814 (1998).
[388] H. Hinrichsen and E. Domany, Damage spreading in the Ising model,
Phys. Rev. E 56, 94–98
(1997).
[389] A. Coniglio, L. de Arcangelis, H. J. Herrmann, and N. Jan, Exact relations between damage spreading
and thermodynamical properties, Europhys. Lett. 8, 315–320 (1989).
[390] L. R. da Silva and H. J. Herrmann, Damage spreading in a gradient, J. Stat. Phys. 52, 463–470
(1988).
[391] N. Boissin and H. J. Herrman, Gradient method for thermal phase transitions, J. Phys. A 24,
L43–L47 (1991).
[392] G. G. Batrouni and A. Hansen, Static critical exponents from the dynamics of damage spreading
and overhangs in the Ising model with temperature gradient, J. Phys. A 25, L1059–L1064 (1992).
[393] M. Droz and L. Frachebourg,
Damage spreading in cellular automata of nonequilibrium phase
transitions, Phys. Lett. A 148, 447–451 (1990).
[394] S. C. Glotzer, P. H. Poole, and N. Jan,
Time-dependent thermodynamic properties of the Ising
model from damage spreading, J. Stat. Phys. 68, 895–910 (1992).
[395] F. Wang and M. Suzuki, Time evolution of damage in the Ising model, Physica A 223, 34–49
(1996).
[396] K. McIsaac and N. Jan, On the dynamic exponent of the two-dimensional Ising model, J. Phys. A
25, 2139–2145 (1992).
[397] D. L. Hunter, L. de Arcangelis, R. Matz, P. H. Poole, and N. Jan, Time-dependent critical properties
of Ising models by damage spreading, Physica 196, 188–208 (1993).
[398] D. Stauﬀer, Eﬀective analysis of damage spreading in Ising models, J. Phys. A 26, L599–L604
(1993).
[399] F. Wang, N. Hatano, and M. Suzuki, Study on dynamical critical exponents of the Ising model using
the damage spreading method, J. Phys. A 28, 4543–4552 (1995).
[400] F. Montani and E. V. Albano, Dynamics of damage spreading in the two-dimensional Ising magnet
at criticality, Phys. Lett. A 202, 253–257 (1995).
151

<!-- page 152 -->

[401] P. H. Poole and N. Jan, Dynamics properties of the two- and three-dimensional Ising models by
‘damage spreading’, J. Phys. A 23, L453–L459 (1990).
[402] R. Matz, D. L. Hunter, and N. Jan, The dynamic critical exponent of the 3D Ising model, J. Stat.
Phys. 74, 903–908 (1994).
[403] P. Grassberger, Damage spreading and critical exponents for ‘model A’ Ising dynamics, Physica A
214, 547–559 (1995).
[404] U. Gropengiesser, Damage spreading and critical exponents for ‘model A’ Ising dynamics, Physica
A 215, 308–310 (1995).
[405] P. C. Hohenberg and B. I. Halperin, Theory of dynamic critical phenomena, Rev. Mod. Phys. 49,
435–479 (1977).
[406] P. Meakin, The growth of fractal aggregates and their fractal measures, in Phase transitions and
critical phenomena, edited by C. Domb and J. L. Lebowitz, volume 12, Academic Press, New York,
1988.
[407] J. Krug and H. Spohn, Kinetic roughening of growing surfaces, in Solids far from equilibrium, edited
by C. Godr`eche, Cambridge University Press, U.K., 1991.
[408] A. J. McKane, M. Droz, J. Vannimenus, and D. Wolf, editors,
Scale Invariance, Interfaces and
Non-Equilibrium Dynamics, volume B 344, Plenum press, New York, 1995.
[409] J. Krug, Origins of scale invariance in growth processes, Adv. Phys. 46, 139–282 (1997).
[410] F. Family and T. Vicsek, Scaling of the active zone in the Eden process on percolation networks and
the ballistic deposition model, J. Phys. A 18, L75–L81 (1985).
[411] S. F. Edwards and D. R. Wilkinson, The surface statistics of a granular aggregate, Proc. R. Soc.
Lond. A 381, 17–31 (1982).
[412] T. Halpin-Healy and Y.-C. Zhang,
Kinetic roughening phenomena, stochastic growth, directed
polymers and all that. Aspects of multidisciplinary statistical mechanics, Phys. Rep. 254, 215–414
(1995).
[413] S. V. Buldyrev, A. L. Barab´asi, F. Caserta, S. Havlin, H. E. Stanley, and T. Vicsek, Anomalous
interface roughening in porous media: Experiment and model,
Phys. Rev. A 45, R8313–R8316
(1992).
[414] T. Nattermann, S. Stepanow, and L. H. Tang,
Dynamics of interface depinning in a disordered
medium, Journal de Physique 2, 1483–1488 (1992).
[415] L. H. Tang and H. Leschhorn, Pinning by directed percolation, Phys. Rev. A 45, R8309–R8312
(1992).
[416] H. A. Makse, S. Buldyrev, H. Leschhorn, and H. E. Stanley, The pinning paths of an elastic interface,
Europhys. Lett. 41, 251–256 (1998).
[417] A. Dougherty and N. Carle, Distribution of avalanches in interfacial motion in a porous medium,
Phys. Rev. E 58, 2889–2893 (1998).
[418] R. Albert, A. L. Barab´asi, N. Carle, and A. Dougherty,
Driven interfaces in disordered media:
determination of universality classes from experimental data, Phys. Rev. Lett. 81, 2926–2929 (1998).
[419] S. V. Buldyrev, A. L. Barab´asi, S. Havlin, J. Kert´esz, H. E. Stanley, and H. S. Xenias, Anomalous
interface roughening in porous media: experiment and model, Physica A 191, 220–226 (1992).
[420] C. Lehner, N. Rajewsky, D. Wolf, and J. Kert´esz,
Anomalous roughening in three-dimensional
polynucleation growth, Physica A 164, 81–90 (1990).
[421] A. Toom, On critical phenomena in interacting growth systems. I. General, J. Stat. Phys. 74, 91–109
(1994).
[422] A. Toom, On critical phenomena in interacting growth systems. II. Bounded growth, J. Stat. Phys.
74, 111–130 (1994).
[423] U. C. T¨auber, M. J. Howard, and H. Hinrichsen, Multicritical behavior in coupled directed percola-
tion processes, Phys. Rev. Lett. 80, 2165–2168 (1998).
[424] Y. Y. Goldschmidt, H. Hinrichsen, M. Howard, and U. C. T¨auber, Novel non-equilibrium critical
behavior in unidirectionally coupled stochastic processes, Phys. Rev. E 59, 6381–6408 (1999).
[425] J. L´opez and H. Jensen, Nonequilibrium roughening transition in a simple model of fungal growth
in 1+1 dimension, Phys. Rev. Lett. 81, 1734–1737 (1998).
[426] T. Sams, K. Sneppen, M. H. Jensen, C. Ellegaard, B. E. Christensen, and U. T. Thrane, Morpho-
logical instabilities in a growing yeast colony: experiment and theory, Phys. Rev. Lett. 79, 313–316
(1997).
[427] U. Alon, M. R. Evans, H. Hinrichsen, and D. Mukamel, Roughening transition in a one-dimensional
growth process, Phys. Rev. Lett. 76, 2746–2749 (1996).
152

<!-- page 153 -->

[428] U. Alon, M. Evans, H. Hinrichsen, and D. Mukamel, Smooth phases, roughening transitions, and
novel exponents in one-dimensional growth models, Phys. Rev. E 57, 4997–5012 (1998).
[429] suggested by O. Biham, private communication (1999).
[430] Y. Y. Goldschmidt, Comment on ”Multicritical behavior in coupled directed percolation processes”,
Phys. Rev. Lett. 81, 2178 (1998).
[431] H. Hinrichsen and G. ´Odor, Roughening transition in a model for dimer adsorption and desorption,
Phys. Rev. Lett. 82, 1205 (1999).
[432] S. Dietrich, Wetting phenomena, in Phase transitions and critical phenomena, edited by C. Domb
and J. L. Lebowitz, volume 12, Academic Press, London, 1988.
[433] D. M. Kroll and R. Lipowsky, Pinning transitions in d-dimensional Ising ferromagnets, Phys. Rev.
B 26, 5289–5292 (1982).
[434] E. Br´ezin, B. I. Halperin, and A. Leibler, Critical wetting in three dimensions, Phys. Rev. Lett. 50,
1387–1390 (1983).
[435] D. S. Fisher and D. A. Huse, Wetting transitions: a functional renormalization-group approach,
Phys. Rev. B 32, 247–256 (1985).
[436] H. Hinrichsen, R. Livi, D. Mukamel, and A. Politi, A model for nonequilibrium wetting transitions
in two dimensions, Phys. Rev. Lett. 79, 2710–2713 (1997).
[437] J. M. J. van Leeuwen and H. J. Hilhorst, Pinning of a rough interface by an external potential,
Physica A 107, 319–329 (1981).
[438] T. W. Burkhardt, Localisation-delocalisation transition in a solid-on-solid model with a pinning
potential, J. Phys. A 14, L63–L68 (1981).
[439] Y. Tu, G. Grinstein, and M. A. Mu˜noz, Systems with multiplicative noise: critical behavior from
KPZ equation and numerics, Phys. Rev. Lett. 78, 274–277 (1997).
[440] M. A. Mu˜noz and T. Hwa, On nonlinear diﬀusion with multiplicative noise, Europhys. Lett. 41,
147–152 (1998).
[441] H. Hinrichsen, R. Livi, D. Mukamel, and A. Politi, First order phase transition in a 1+1-dimensional
nonequilibrium wetting process, Phys. Rev. E 61, R1032–R1035 (2000).
[442] L. Peliti, Path integral approach to birth-death processes on a lattice, Journal de Physique 46,
1469–1483 (1984).
153
