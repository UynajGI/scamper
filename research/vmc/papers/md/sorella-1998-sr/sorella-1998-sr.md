# Green Function Monte Carlo with Stochastic Reconfiguration

Sandro Sorella

Istituto Nazionale di Fisica della Materia and International School for Advanced Studies Via

Beirut 4, 34013 Trieste, Italy

(November 26, 2024)

## Abstract

A new method for the stabilization of the sign problem in the Green Function Monte Carlo technique is proposed. The method is devised for real lattice Hamiltonians and is based on an iterative ”stochastic reconfiguration” scheme which introduces some bias but allows a stable simulation with constant sign. The systematic reduction of this bias is in principle possible. The method is applied to the frustrated $J _ { 1 } - J _ { 2 }$ Heisenberg model, and tested against exact diagonalization data. Evidence of a finite spin gap for $J _ { 2 } / J _ { 1 } > \sim 0 . 4$ is found in the thermodynamic limit.

02.70.Lq,75.10.Jm,75.40.Mg

As well known the Green Function Monte Carlo method (GFMC) allows to obtain the exact ground state properties of a many body hamiltonian with a statistical method. One of the most severe restriction is that only positive definite Green function GF can be sampled, otherwise the method is facing the well known ”sign problem”. Approximate techniques like the fixed node approximation (FN) have been developed to circumvent the sign problem but at the very least they cannot be systematically improved to achieve the exact answer within statistical errors. This property has severely limited the applications of GFMC to fermions and frustrated boson models. In this letter I propose a new approach to stabilize the sign problem, the GFMC with stochastic reconfiguration (GFMCSR), which will be shortly described below, revisiting also the basic steps of the standard GFMC on a lattice. [1,2]

In order to filter out the ground state of a given lattice hamiltonian H the standard power method may be applied iteratively :

$$
\psi _ { n + 1 } ( x ^ { \prime } ) = \sum _ { x } ( \Lambda \delta _ { x ^ { \prime } , x } - H _ { x ^ { \prime } , x } ) \psi _ { n } ( x )\tag{1}
$$

where x represents conventionally the index of a complete basis $| x > , H _ { x ^ { \prime } , x }$ being the corresponding matrix elements of the hamiltonian which in the following are assumed real, and Λ is a positive constant that allows the convergence of $\psi _ { n }$ to the ground state $\psi _ { 0 } ( x )$ , for large n. In numerical calculations of interesting lattice hamiltonians the dimension of the basis grows exponentially with the size and the particle number, though the matrix itself is very sparse and all its elements $H _ { x ^ { \prime } , x }$ , for given x, can be generally computed even for large system size. In this case an exact application of (1) is impossible unless for few steps. A way out is to use a stochastic approach , like GFMC ,which is particularly simple on a lattice.

In order to implement stochastically the iteration (1) the corresponding lattice GF

$$
\begin{array} { r } { G _ { x ^ { \prime } , x } = \Lambda \delta _ { x ^ { \prime } , x } - H _ { x ^ { \prime } , x } } \end{array}\tag{2}
$$

may be decomposed in the following way:

$$
G _ { x ^ { \prime } , x } = s _ { x ^ { \prime } , x } p _ { x ^ { \prime } , x } b _ { x }\tag{3}
$$

where $p _ { x ^ { \prime } , x }$ is a normalized stochastic matrix, $b _ { x } \geq 0$ is a normalization constant and the matrix s takes into account the sign of the GF. The typical choice is to take $p _ { x ^ { \prime } , x } = | G _ { x ^ { \prime } , x } | / b _ { x }$ $\begin{array} { r } { b _ { x } = \sum _ { x ^ { \prime } } \left| G _ { x ^ { \prime } , x } \right| } \end{array}$ and $s _ { x ^ { \prime } , x } = \operatorname { s g n } G _ { x ^ { \prime } , x }$ , which is identically one if there is no sign problem.

In the GFMC method the so called ”walker“ is defined by a weight w and a configuration x.. At a given iteration n the walker is assumed to sample statistically the state $\psi _ { n } ( x )$ in Eq.(1), in the sense that the probability $P _ { n } ( w , x )$ to have the walker with weight w (not restricted to be positive) in a given configuration x satisfies: R $\dot { { d w P } _ { n } } ( w , x ) w = \psi _ { n } ( x )$ . Then the matrix multiplication (1) can be implemented statistically , in the precise sense that $\begin{array} { r } { \int d w P _ { n + 1 } ( w , x ) w = \psi _ { n + 1 } ( x ) } \end{array}$ , by the following three steps:

1. scale the walker weight by $b _ { x } \colon w ^ { \prime } = b _ { x } w$

2. select randomly a new configuration x′ according to the stochastic matrix $p _ { x ^ { \prime } , x }$

3. finally multiply the weight of the walker by the sign factor $s _ { x ^ { \prime } , x } \colon w ^ { \prime } \to w ^ { \prime } s _ { x ^ { \prime } , x }$ (MI)

In principle the previous Markov process determines, for large n, the ground state of H even with a single walker. In practise it is convenient to use a large number M of walkers, which I indicate by $\left( w _ { j } , x _ { j } \right) j = 1 , \cdot \cdot \cdot M$ , shorthand in the following also by vector notations w, x.

If there is sign problem the average walker sign $< s > _ { n } = \frac { < \sum _ { j } w _ { j } > _ { n } } { < \sum _ { j } | w _ { j } | > _ { n } }$ decreases exponentially to zero as the Markov iteration MI is repeatedly applied and it is basically impossible to reach a reasonably large value of n.

Recently a remarkable progress in GFMC on a lattice was the extension of the FN to this case. The method is based on a definition of an effective GF $G _ { x ^ { \prime } , x } ^ { f }$ which is always positive definite but yields a good variational estimate of the energy. For later purposes we define this effective GF in a slightly different way, by introducing a parameter γ: which allows to sample also the negative elements of the GF :

$$
G _ { x ^ { \prime } , x } ^ { f } = \left\{ \begin{array} { c c } { - H _ { x ^ { \prime } , x } } & { \mathrm { i f } H _ { x ^ { \prime } , x } \leq 0 } \\ { \gamma H _ { x ^ { \prime } , x } } & { \mathrm { i f } H _ { x ^ { \prime } , x } > 0 } \\ { \Lambda - H _ { x , x } - ( 1 + \gamma ) \mathcal { V } _ { \mathrm { s f } } ( x ) } & { \mathrm { i f } x = x ^ { \prime } } \end{array} \right.\tag{4}
$$

where the diagonal sign-flip contribution is given by [3,4]:

$$
\mathcal { V } _ { \mathrm { s f } } ( x ) = \sum _ { H _ { \chi ^ { \prime } , \chi } > 0 , x ^ { \prime } \neq x } H _ { \chi ^ { \prime } , \chi }\tag{5}
$$

For $\gamma = 0$ the usual formulation [4] is recovered, whereas for $\gamma > 0 \ [ 5 ]$ the crossing to the negative sign region is allowed so that the exact GF can be written as $G _ { x ^ { \prime } , x } = s _ { x ^ { \prime } , x } G _ { x ^ { \prime } , x } ^ { f }$ where $s _ { x ^ { \prime } , x }$ is finite and non zero and is determined by the ratio $\frac { G _ { x ^ { \prime } , x } } { G _ { x ^ { \prime } , x } ^ { f } }$ with G and $G ^ { f }$ given by Eq.(2) and Eq.(4) respectively. The value of the constant $\gamma$ necessary to cross the ”nodal surface” was chosen to be $1 / 2$ in all forthcoming applications.

In the basic decomposition (3) the stochastic matrix $p _ { x ^ { \prime } , x } = G _ { x ^ { \prime } , x } ^ { f } / b _ { x }$ and the normalization coefficient $\begin{array} { r } { b _ { x } = \sum _ { x ^ { \prime } } G _ { x ^ { \prime } , x } ^ { f } } \end{array}$ are instead determined only by $G ^ { f }$

By omitting the last step $w ^ { \prime }  w s _ { x ^ { \prime } , x }$ in the Markov iteration process MI, the state $\psi _ { n }$ is indeed propagated trough the positive GF $G ^ { f }$ . The main property used in the following is that at any Markov iteration n we can have a statistic knowledge of both the state $\psi _ { n } ( x )$ obtained with the exact GF and of $\psi _ { n } ^ { f } ( x )$ obtained instead with the approximate but positive definite one $G ^ { f }$ . To this purpose the $j ^ { t h }$ walker is defined by two weights $\boldsymbol { w } _ { j } ^ { f }$ and $w _ { j }$ corresponding to the propagation of the walker by $G ^ { f }$ and $G$ respectively. These weights act on the same configuration $x _ { j }$ . Hereafter the vector w represents therefore a shorthand notation for the 2M components $w _ { j } , w _ { j } ^ { f }$ for j = 1, . . . M .

The walker vector w, x allows to determine statistically the state:

$$
\psi _ { n } ( x ) = \int d [ \underline { { w } } ] \sum _ { \underline { { x } } } P _ { n } ( \underline { { w } } , \underline { { x } } ) \sum _ { j } \delta _ { x , x _ { j } } w _ { j } / M\tag{6}
$$

and analogously $\psi _ { n } ^ { f } ( x )$ by replacing the weights $w _ { j }$ with the positive ones $\boldsymbol { w } _ { j } ^ { f }$ in the previous equation. In this way the configurations generated by the described Markov process MI, if weighted with the constants $\boldsymbol { w } _ { j } ^ { f }$ , are distributed for large n, according to the variational state corresponding to $G ^ { f }$ . This is a reasonable variational wavefunction (WF), which will be the initial approximation to which systematic corrections will be applied, as described later on.

Apart for the previous technical definitions, we can explain in few words the basic idea used for the stabilization of the sign problem, The iteration MI converges to the ground state, but due to the sign problem, only few iterations can be performed with a reasonable statistical accuracy. However, the representation of the state $\psi _ { n } ( x )$ in terms of the walker population $x _ { j } , w _ { j }$ is not unique. In fact it is perfectly possible to represent the same state $\psi _ { n } ( x )$ either with a walker population with very small average sign or with a one with a very large average sign. If such reconfigurations are possible each few $k _ { p }$ steps, the average sign may be stabilized to a large value during the iteration (1) and there will be no difficulty to sample the ground state for $n \longrightarrow \infty .$ with no sign problem.

I will show that this reconfiguration is well defined and indeed possible. The set of M walkers $( { \underline { { w } } } , { \underline { { x } } } )$ are defined via their probability function $P _ { n } ( { \underline { { w } } } , { \underline { { x } } } )$ which in turn defines the state $\psi _ { n } ( x )$ by Eq.(6). The task is to change $P _ { n }$ onto a new probability distribution $P _ { n } ^ { \prime }$ corresponding to a steadily high sign for the walker population. This without changing the information content, the state $\psi _ { n } ( x )$

Let us define the new state $\psi _ { n } ^ { \prime } ( x )$ , as the one obtained by averaging over $P _ { n } ^ { \prime }$ in Eq.(6), then the reconfiguration is exact if $P _ { n } ^ { \prime }$ is such that:

$$
\psi _ { n } ^ { \prime } ( x ) = \psi _ { n } ( x ) \mathrm { f o r ~ a l l } \ x\tag{7}
$$

In general it is difficult or impractical to realize all these conditions (7) as their number equals the dimension of the Hilbert space. I consider therefore a set of operators $O ^ { k } , k =$ $1 , \cdots p < < M$ and require only $p + 1$ stochastic reconfiguration conditions:

$$
\sum _ { x ^ { \prime } , x } O _ { x ^ { \prime } , x } ^ { k } \psi _ { n } ^ { ' } ( x ) = \sum _ { x ^ { \prime } , x } O _ { x ^ { \prime } , x } ^ { k } \psi _ { n } ( x )\tag{8}
$$

for $k = 1 , \cdots p _ { \mathrm { : } }$ , beyond the normalization one $\begin{array} { r } { \sum _ { x } \psi ^ { \prime } ( x ) = \sum _ { x } \psi _ { n } ( x ) } \end{array}$

The previous equations (8) mean that the so called ”mixed averages” of the operators $O ^ { k }$ coincide before and after the reconfiguration. [6]

The main idea of this work is that these $p + 1$ conditions can be fulfilled exactly (for chosen operators) by defining the reconfiguration in the following form:

$$
P _ { n } ^ { \prime } ( \underline { { w } } ^ { \prime } , \underline { { x } } ^ { \prime } ) = \int d [ \underline { { w } } ] \sum _ { \underline { { x } } } \prod _ { i = 1 } ^ { M } \left\{ \frac { \sum _ { j } | p _ { x _ { i } } | \delta _ { X _ { i } ^ { \prime } , X _ { j } } } { \sum _ { j } | p _ { x _ { j } } | } \delta ( w _ { i } ^ { \prime } - \frac { \sum _ { j } w _ { j } } { \beta M } \mathrm { s g n } p _ { x _ { i } ^ { \prime } } ) \delta ( w _ { i } ^ { f \prime } - | w _ { i } ^ { \prime } | ) \right\} P _ { n } ( \underline { { w } } , \underline { { x } } )\tag{9}
$$

where $\beta = \frac { \sum _ { j } p _ { x _ { j } } } { \sum _ { j } | p _ { x _ { j } } | }$ is the average sign after the reconfiguration which is supposed to be much higher to stabilize the process. The new configurations $x _ { i } ^ { \prime }$ are taken randomly among the old ones $\{ x _ { j } \}$ , according to the table $p _ { x _ { j } }$ , defined below.. The positive weights $\boldsymbol { w } _ { j } ^ { f }$ represent a good starting point for the definition of a reconfiguration with large $\beta .$ Though there is some arbitrariness in the definition of the coefficients $p _ { x _ { j } } ,$ a simple and convenient choice is:

$$
p _ { x _ { j } } = w _ { j } ^ { f } ( 1 + \sum _ { k } \alpha _ { k } ( O _ { j } ^ { k } - \hat { O } _ { f } ^ { k } ) )
$$

where $\begin{array} { r } { \bar { O } _ { f } ^ { k } = \frac { \sum _ { j } w _ { j } ^ { f } O _ { j } ^ { k } } { \sum _ { j } w _ { j } ^ { f } } } \end{array}$ are the averages over the positive weights $\boldsymbol { w } _ { j } ^ { f }$ of the mixed estimates $\begin{array} { r } { O _ { j } ^ { k } = \sum _ { x ^ { \prime } } O _ { x ^ { \prime } , x _ { j } } ^ { k } } \end{array}$ corresponding to the operator $O ^ { k }$ and the configuration $x _ { j }$

Then, in order to satisfy the WF conditions (8), by using the definition (9), it is sufficient that the coefficients $p _ { x _ { j } }$ satisfy the following Markovian conditions:

$$
\frac { \sum _ { j } p _ { x _ { j } } O _ { J } ^ { k } } { \sum _ { j } p _ { x _ { j } } } = \frac { \sum _ { j } w _ { j } O _ { j } ^ { k } } { \sum _ { j } w _ { j } }\tag{10}
$$

which in turn determine the unknown variables $\alpha _ { k } .$ for $k = 1 , \cdots p _ { ; }$ for given w, x.

For hamiltonian not affected by the sign problem $( G ^ { f } = G \alpha _ { k } = 0$ and $\beta = 1 )$ this reconfiguration was already used to control the walker population size without introducing any source of systematic error. [7] The present more general reconfiguration (9) can be easily and efficiently implemented in a similar way.

Obviously the reconfiguration conditions (8) are equivalent to the exact ones (7), when the number p of linearly independent operators considered in (8) is equal to the large dimension of the Hilbert space. An important applicative issue is whether GFMCSR converges, within a reasonable accuracy, even with a small number p of meaningful operators $O ^ { k }$

We consider the frustrated $J _ { 1 } - J _ { 2 }$ Heisenberg spin $1 / 2$ model on a finite square lattice with L sites and periodic boundary conditions (tilted by 45 degrees for the L = 32 size only). The model hamiltonian is determined by an antiferromagnetic coupling $J _ { 1 } > 0$ between nearest neighbor spins and a frustrating coupling $J _ { 2 } > 0$ between next neighbor ones. [8–10] In all forthcoming examples the stochastic reconfigurations were applied frequently enough to maintain the average sign before reconfiguration $\sim 0 . 8$ , condition that minimize the statistical fluctuations. Moreover in each simulation it is important to work with a fairly large number of walkers, since in the $M \to \infty$ limit, the GFMCSR results are practically independent of the frequency of reconfigurations, as well as the overall constant energy shift Λ.

The accuracy of GFMCSR for the ground state is displayed in Tab.I , and compared with other methods. The variational WF (used also for GFMC importance sampling [6]) contains a Jastrow like factor

$$
E x p ( \frac { \eta } { 2 } \sum _ { R , R ^ { \prime } } v ( R - R ^ { \prime } ) S _ { R } ^ { z } S _ { R ^ { \prime } } ^ { z } )
$$

to mimic the interaction between the spins $S _ { R } ^ { z } = \pm 1 / 2$ at sites R, R′, where η is a variational parameter and the two-spin interaction v can be derived by using the method described in [11], yielding an explicit Fourier transform for v:

$$
v _ { q } / 2 = 1 - \sqrt { \frac { 2 - \sigma ( 1 - \cos q _ { x } \cos q _ { y } ) + \cos q _ { x } + \cos q _ { y } } { 2 - \sigma ( 1 - \cos q _ { x } \cos q _ { y } ) - \cos q _ { x } - \cos q _ { y } } }
$$

with $\sigma = 2 J _ { 2 } / J _ { 1 }$ . This potential is not defined for $J _ { 2 } / J _ { 1 } = 1 / 2$ , and in such case I have chosen to work with $\sigma = 0 . 8$ . Restriction to any subspace of total spin projection $\begin{array} { r } { S _ { t o t } ^ { z } = \sum _ { R } S _ { R } ^ { z } } \end{array}$ allows to evaluate the spin gap by performing two simulations for $S _ { t o t } ^ { z } = 0$ and $S _ { t o t } ^ { z } = 1$ Henceforth I will use the the same potential v in both subspaces, by optimizing η for the $S _ { t o t } ^ { z } = 0$ energy.

As shown in the table the accuracy of the variational WF is rather poor, and is considerably improved by the FN, at least for small $J _ { 2 }$ . This kind of accuracy is however not enough to determine the rapid increase of the spin gap as $J _ { 2 } / J _ { 1 }$ approaches the value $1 / 2$ of the classical transition. Instead, as shown in Fig.(1) the GFMCSR allows to achieve a good accuracy also on this delicate quantity by considering in the reconfigurations only the energy and the spin structure factor $\begin{array} { r } { S _ { q } ^ { z } = \sum _ { R , R ^ { \prime } } e ^ { i q ( R - R ^ { \prime } ) } S _ { R } ^ { z } S _ { R ^ { \prime } } ^ { z } } \end{array}$ symmetrized over all directions and for all non equivalent wavevectors q. Remarkably also mixed averages of correlation functions that are not included in such reconfiguration conditions (8) are also significantly improved (see table).

The way GFMCSR reaches the large n limit (at fixed number of operators p) is displayed in Fig.(2) where the initial n = 0 distribution was obtained by the FN for $\gamma = 0$ . For fixed p the algorithm is Markovian and reaches an equilibrium distribution for $n \to \infty$ , independent of the initial one (see example in Fig.2 where p was changed at the iteration indicated by the arrow), this in turn will converge to the ground state distribution for large p. A comparison with the standard ”release nodes” estimate is also shown in the picture. It is clear that there is no hope to obtain meaningful results in this case by the direct sampling of the sign. On the contrary this method looks very stable and, though approximate, a convergence to a reasonable accuracy is obtained even with a very small number of operators, compared to the dimension of the Hilbert space.

The data shown in the table and in the picture indicate that the accuracy of GFMCSR may become rather size independent with a relatively small increase of the operator number p. The error to work at finite small p is systematic. Thus there is a considerable cancellation of this error for the determination of the spin gap displayed in Fig.(1).

The calculation was therefore extended to the large size system up to L = 100 where exact diagonalization is not possible. The spin gap as a function of the system size is displayed in Fig.(3). This figure is consistent with the opening of a finite spin gap for $J _ { 2 } / J _ { 1 } \ge \sim 0 . 4$ This gap is certainly not an artifact of the variational WF, which is obviously gapless, as also confirmed numerically in the same figure. The present numerical results confirm that the transition to a spin liquid state with a finite spin gap but no classical order parameter should be close to $J _ { 2 } / J _ { 1 } = 0 . 4$ . [10]

This work was supported in part by INFM (PRA HTSC) and CINECA grant.

## REFERENCES

[1] N. Trivedi, D. Ceperley, Phys. Rev. B 41, 4552 (1990)

[2] K. Runge, Phys. Rev. B 44, 122252 (1992)

[3] H. van Bemmel et al. Phys. Rev. Lett. 72, 2442 (1994).

[4] D. ten Haaf et al. Phys. Rev. B 51, 13039 (1995).

[5] It is possible to prove that the method is variational also for $\gamma > 0$ , with a simple extension of the proof in [4].

[6] All the analysis remains unchanged if a guiding function $\psi _ { G } ( x )$ is used for importance sampling. The matrix elements of all the operators $O _ { x ^ { \prime } , x }$ (including the GF) have to be accordingly changed : $O _ { x ^ { \prime } , x }  \psi _ { G } ( x ^ { \prime } ) O _ { x ^ { \prime } , x } / \psi _ { G } ( x )$

[7] M. Calandra and S. Sorella, to appear in Phys. Rev. B .

[8] E. Dagotto, A. Moreo Phys. Rev. Lett. 63, 2148 (1989).

[9] T. Nakamura et al. J. Phys. Soc. Japan 61, 3494 (1992).

[10] J. Schulz et al. J. de Phys. 6, 675 (1996)..

[11] F. Franjic, S. Sorella, Prog. Theor. Phys. 97, 399 (1997)

<!-- image-->  
FIG. 1. Estimate of the spin gap for several methods: variational (empty triangles), FN (empty squares), GFMCSR p = 1 (empty dots) , GFMCSR (full dots) as in the table for $L = 1 6$ (upper points) and L = 32 (lower ones). The exact results are connected by continuous lines.

<!-- image-->  
FIG. 2. Energy per site vs. n for GFMCSR with p = 1 (upper curve to the left of the arrow) and $p = 9$ (remaining curves). The triangles represent the standard method with sign problem, i.e. with large error bars already for $n > 1 5$

<!-- image-->  
FIG. 3. Size scaling of the spin gap. The dashed lines are linear fit of the GFMCSR data with $p = 9 , 1 4 , 2 0$ for L = 36, 64, 100 respectively. Lower curves are the variational estimates and continuous lines are guides to the eye.

TABLES
<table><tr><td rowspan=1 colspan=1> $J _ { 2 } / J _ { 1 }$ </td><td rowspan=1 colspan=1>L</td><td rowspan=1 colspan=1>η</td><td rowspan=1 colspan=1>% VMC</td><td rowspan=1 colspan=1>% FN</td><td rowspan=1 colspan=1>% SRe</td><td rowspan=1 colspan=1>% SR</td></tr><tr><td rowspan=1 colspan=1>0.1</td><td rowspan=1 colspan=1>16</td><td rowspan=1 colspan=1>1.2</td><td rowspan=1 colspan=1>2.84 (2.2)</td><td rowspan=1 colspan=1>0.17 (0.1)</td><td rowspan=1 colspan=1>-0.03 (0.0)</td><td rowspan=1 colspan=1>0.02 (0.0)</td></tr><tr><td rowspan=1 colspan=1>0.2</td><td rowspan=1 colspan=1>16</td><td rowspan=1 colspan=1>1.15</td><td rowspan=1 colspan=1>2.80 (2.5)</td><td rowspan=1 colspan=1>0.41 (0.4)</td><td rowspan=1 colspan=1>0.00 (0.2)</td><td rowspan=1 colspan=1>0.03 (0.0)</td></tr><tr><td rowspan=1 colspan=1>0.3</td><td rowspan=1 colspan=1>16</td><td rowspan=1 colspan=1>1.1</td><td rowspan=1 colspan=1>3.25 (2.5)</td><td rowspan=1 colspan=1>0.87 (0.7)</td><td rowspan=1 colspan=1>0.12 (0.8)</td><td rowspan=1 colspan=1>0.05 (0.1)</td></tr><tr><td rowspan=1 colspan=1>0.4</td><td rowspan=1 colspan=1>16</td><td rowspan=1 colspan=1>0.8</td><td rowspan=1 colspan=1>3.38 (2.4)</td><td rowspan=1 colspan=1>1.76 (3.2)</td><td rowspan=1 colspan=1>0.56 (4.5)</td><td rowspan=1 colspan=1>0.26 (0.2)</td></tr><tr><td rowspan=1 colspan=1>0.5</td><td rowspan=1 colspan=1>16</td><td rowspan=1 colspan=1>0.85</td><td rowspan=1 colspan=1>5.65 (10.9)</td><td rowspan=1 colspan=1>3.84 (8.9)</td><td rowspan=1 colspan=1>2.08 (8.9)</td><td rowspan=1 colspan=1>0.66 (1.1)</td></tr><tr><td rowspan=1 colspan=1>0.1</td><td rowspan=1 colspan=1>32</td><td rowspan=1 colspan=1>1</td><td rowspan=1 colspan=1>1.55 (2.5)</td><td rowspan=1 colspan=1>0.22 (0.3)</td><td rowspan=1 colspan=1>0.05 (0.1)</td><td rowspan=1 colspan=1>0.02 (0.0)</td></tr><tr><td rowspan=1 colspan=1>0.2</td><td rowspan=1 colspan=1>32</td><td rowspan=1 colspan=1>1</td><td rowspan=1 colspan=1>1.78 (2.5)</td><td rowspan=1 colspan=1>0.48 (0.6)</td><td rowspan=1 colspan=1>0.15 (0.6)</td><td rowspan=1 colspan=1>0.05 (0.1)</td></tr><tr><td rowspan=1 colspan=1>0.3</td><td rowspan=1 colspan=1>32</td><td rowspan=1 colspan=1>1</td><td rowspan=1 colspan=1>2.23 (2.1)</td><td rowspan=1 colspan=1>0.85 (0.91)</td><td rowspan=1 colspan=1>0.30 (1.4)</td><td rowspan=1 colspan=1>0.10 (0.0)</td></tr><tr><td rowspan=1 colspan=1>0.4</td><td rowspan=1 colspan=1>32</td><td rowspan=1 colspan=1>0.8</td><td rowspan=1 colspan=1>3.07 (4.0)</td><td rowspan=1 colspan=1>1.61 (3.1)</td><td rowspan=1 colspan=1>0.26 (5.6)</td><td rowspan=1 colspan=1>0.21 (0.1)</td></tr><tr><td rowspan=1 colspan=1>0.5</td><td rowspan=1 colspan=1>32</td><td rowspan=1 colspan=1>0.9</td><td rowspan=1 colspan=1>4.51 (10.0)</td><td rowspan=1 colspan=1>2.92 (7.2)</td><td rowspan=1 colspan=1>1.52 (7.7)</td><td rowspan=1 colspan=1>0.46 (0.9)</td></tr><tr><td rowspan=1 colspan=1>0.1</td><td rowspan=1 colspan=1>36</td><td rowspan=1 colspan=1>1.1</td><td rowspan=1 colspan=1>1.86 (2.8)</td><td rowspan=1 colspan=1>0.21 (0.2)</td><td rowspan=1 colspan=1>0.1 (0.12)</td><td rowspan=1 colspan=1>0.02 (0.1)</td></tr><tr><td rowspan=1 colspan=1>0.2</td><td rowspan=1 colspan=1>36</td><td rowspan=1 colspan=1>1.1</td><td rowspan=1 colspan=1>2.22 (2.8)</td><td rowspan=1 colspan=1>0.47 (0.5)</td><td rowspan=1 colspan=1>0.16 (0.5)</td><td rowspan=1 colspan=1>0.07 (0.1)</td></tr><tr><td rowspan=1 colspan=1>0.3</td><td rowspan=1 colspan=1>36</td><td rowspan=1 colspan=1>1</td><td rowspan=1 colspan=1>2.31 (2.8)</td><td rowspan=1 colspan=1>0.91 (1.4)</td><td rowspan=1 colspan=1>0.35 (2.0)</td><td rowspan=1 colspan=1>0.11 (0.1)</td></tr><tr><td rowspan=1 colspan=1>0.4</td><td rowspan=1 colspan=1>36</td><td rowspan=1 colspan=1>0.8</td><td rowspan=1 colspan=1>3.34 (5.5)</td><td rowspan=1 colspan=1>1.74 (4.5)</td><td rowspan=1 colspan=1>0.51 (6.8)</td><td rowspan=1 colspan=1>0.26 (0.3)</td></tr><tr><td rowspan=1 colspan=1>0.5</td><td rowspan=1 colspan=1>36</td><td rowspan=1 colspan=1>0.9</td><td rowspan=1 colspan=1>5.09 (14.4)</td><td rowspan=1 colspan=1>3.34 (11.1)]</td><td rowspan=1 colspan=1>1.83 (11.8)</td><td rowspan=1 colspan=1>0.62 (2.1)</td></tr></table>

TABLE I. Percentage error of the energy (square antiferromagnetic order parameter $\vec { m } ^ { 2 }$ as in [7]) for the various methods: variational (VMC), fixed node (FN) , p = 1 GFMCSR (SRe) with the energy alone and $p = 5 , 8 , 9$ GFMCSR estimate (SR) with the energy and $S _ { q } ^ { z }$ for $L = 1 6 , 3 2 , 3 6$ The statistical errors are about one place in the last digit.