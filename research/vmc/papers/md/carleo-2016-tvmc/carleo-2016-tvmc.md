# Unitary dynamics of strongly-interacting Bose gases with time-dependent variational Monte Carlo in continuous space

Giuseppe Carleo

Institute for Theoretical Physics, ETH Zurich, Wolfgang-Pauli-Str. 27, 8093 Zurich, Switzerland

Lorenzo Cevolani

Laboratoire Charles Fabry, Institut d’Optique, CNRS, Univ. Paris Sud 11, 2 avenue Augustin Fresnel, F-91127 Palaiseau cedex, France

Laurent Sanchez-Palencia

Laboratoire Charles Fabry, Institut d’Optique, CNRS, Univ. Paris Sud 11, 2 avenue Augustin Fresnel, F-91127 Palaiseau cedex, France and Centre de Physique Théorique, Ecole Polytechnique, CNRS, Univ Paris-Saclay, F-91128 Palaiseau, France

Markus Holzmann

LPMMC, UMR 5493 of CNRS, Université Grenoble Alpes, F-38042 Grenoble, France and Institut Laue-Langevin, BP 156, F-38042 Grenoble Cedex 9, France

We introduce time-dependent variational Monte Carlo for continuous-space Bose gases. Our approach is based on the systematic expansion of the many-body wave-function in terms of multibody correlations and is essentially exact up to adaptive truncation. The method is benchmarked by comparison to exact Bethe-ansatz or existing numerical results for the integrable Lieb-Liniger model. We first show that the many-body wave-function achieves high precision for ground-state properties, including energy and first-order as well as second-order correlation functions. Then, we study the out-of-equilibrium, unitary dynamics induced by a quantum quench in the interaction strength. Our time-dependent variational Monte Carlo results are benchmarked by comparison to exact Bethe ansatz results available for a small number of particles, and also compared to quench action results available for non-interacting initial states. Moreover, our approach allows us to study large particle numbers and general quench protocols, previously inaccessible beyond the mean-field level. Our results suggest that it is possible to find correlated initial states for which the long-term dynamics of local density fluctuations is close to the predictions of a simple Boltzmann ensemble.

## I. INTRODUCTION

The study of equilibration and thermalization properties of complex many-body systems is of fundamental interest for many areas of physics and natural sciences [1]. For systems governed by classical physics, an exact solution of Newton’s equations of motion is often numerically feasible, using for instance moleculardynamics simulations. For quantum systems, the mathematical structure of the time-dependent Schrödinger equation is instead fundamentally more involved. Quantum Monte Carlo algorithms, the de facto tool for simulating quantum many-body systems at thermal equilibrium [2–4], cannot be directly used to study timedependent unitary dynamics. Out-of-equilibrium properties are then often treated on the basis of approximations drastically simplifying the microscopic physics. Irreversibility is either enforced with an explicit breaking of unitarity, e.g. within the quantum Boltzmann approach, or the dynamics is reduced to mean-field description using time-dependent Hartree-Fock and Gross-Pitaevskii approaches. Although these approaches may qualitatively describe thermalization [5, 6], their range of validity cannot be assessed because genuine quantum correlations and entanglement are ignored.

For specific systems, exact dynamical results can be derived. This is the case for integrable 1D models, for which Bethe ansatz (BA) solutions exist [7]. However, also in this case many open questions still persist. For example, the exact evaluation of correlation functions for out-ofequilibrium dynamics is at present an unsolved problem. As a result, despite important theoretical and experimental progress [8–15], a complete picture of thermalization (or its absence), e.g. based on general quench protocols [16], is still missing.

Numerical methods for strongly-interacting systems face important challenges as well. Numerical Renormalization Group (NRG) and Density-Matrix Renormalization Group (DMRG) approaches provide an essentially exact description of arbitrary 1D lattice systems in- and out-of-equilibrium [17–20], but they have less predictive power when applied to continuos-space systems. On one hand, multi-scale extensions of the DMRG optimization scheme to the limit of continuous-space lattices [21] are so-far limited to relatively small system sizes [22]. On the other hand, efficient ground-state optimization schemes for continuous quantum field matrix product states (c-MPS) [23] have been introduced only very recently [24], and applications to quantum dynamics are still to be realized. A further formidable challenge is the efficient extension of these approaches to higher dimensions, which is a fundamentally hard problem.

Another class of methods for strongly-interacting systems is based on variational Monte Carlo (VMC), combining highly-entangled variational states with robust stochastic optimization schemes [25]. Such approaches have been successfully applied to the description of continuous quantum systems, in any dimension and not only in 1D [26, 27]. More recently, out-of-equilibrium dynamics has become accessible with the extension of these methods to real-time unitary dynamics, within timedependent variational Monte Carlo (t-VMC) [28, 29]. So far, the t-VMC approach has been developed for lattice systems with bosonic [28, 29], spin [30–32], and fermionic [33] statistics, yielding a description of dynamical properties with an accuracy often comparable with MPS-based approaches.

In this Paper, we extend t-VMC to access dynamical properties of interacting quantum systems in continuous space. Our approach is based on a systematic expansion of the wave function in terms of few-body Jastrow correlation functions. Using the 1D Lieb-Liniger model as a test case, we first show that the inclusion of highorder correlations allows us to systematically approach the exact BA ground state energy. Our results improve by orders of magnitude on previously published VMC and c-MPS results, and are in line with latest state-ofthe-art developments in the field. We further compute single-body and pair correlation functions, hardly accessible by current BA methods. We then calculate the time evolution of the contact pair correlation function following a quench in the interaction strength. For the non-interacting initial state, we benchmark our results to exact BA calculations available for a small number of bosons and further compare to the quench action approach for large systems approaching the thermodynamic limit. We finally apply our method to the study of general quenches from arbitrary initial states, for which no exact results in the thermodynamic limit are currently available.

## II. METHOD

## A. Expansion of the Many-Body Wave-Function

Consider a non-relativistic quantum system of N identical bosons in d dimensions, and governed by the firstquantization Hamiltonian

$$
\mathcal { H } = - \frac { 1 } { 2 } \sum _ { i = 1 } ^ { N } \nabla _ { i } ^ { 2 } + \sum _ { i = 1 } ^ { N } v _ { 1 } ( \vec { x } _ { i } ) + \frac { 1 } { 2 } \sum _ { i \neq j } v _ { 2 } ( \vec { x } _ { i } , \vec { x } _ { j } ) ,\tag{1}
$$

where $v _ { 1 } ( \vec { x } )$ and $v _ { 2 } ( \vec { x } , \vec { y } )$ are, respectively, a one-body external potential and a pair-wise inter-particle interaction [34]. Without loss of generality, a time-dependent N-body state can be written as $\Phi ( { \dot { \bf X } } , t ) = \exp \left[ { \bar { U } } ( { \bf X } ; t ) \right]$ , where $\mathbf { X } = { \vec { x } } _ { 1 } , { \vec { x } } _ { 2 } \dots { \vec { x } } _ { N }$ is the ensemble of particle positions and U is a complex-valued function of the $N _ { - }$ particle coordinates, $\mathbb { R } ^ { N \times d } \to \mathbb { C } .$ . Since the Hamiltonian (1) contains only two-body interactions, it is expected that an expansion of U in terms of few-body Jastrow functions containing at most m-body terms, rapidly converges towards the exact solution. Truncating this expansion up to a certain order $M \leq N$ , leads to the Bijl-Dingle-Jastrow-Feenberg expansion [35–38]

$$
\begin{array} { c } { { \displaystyle { \cal U } ^ { ( M ) } ( { \bf X } ; t ) = \sum _ { i = 1 } u _ { 1 } ( \vec { x } _ { i } ; t ) + \frac { 1 } { 2 ! } \sum _ { i \neq j } u _ { 2 } ( \vec { x } _ { i } , \vec { x } _ { j } ; t ) + } } \\ { { + \ldots \cdot \displaystyle { \frac { 1 } { M ! } \sum _ { i _ { 1 } \neq i _ { 2 } \neq \ldots i _ { M } } u _ { M } ( \vec { x } _ { i _ { 1 } } , \vec { x } _ { i _ { 2 } } \ldots \vec { x } _ { i _ { M } } ; t ) } , } } \end{array}\tag{2}
$$

where $u _ { m } ( \mathbf { r } ; t )$ are functions of m particle coordinates, $\mathbf { r } = \Vec { x } _ { i _ { 1 } } , \Vec { x } _ { i _ { 2 } } \dots \Vec { x } _ { i _ { m } } ,$ , and of the time t. A global constraint on the function $\ddot { U } ( { \mathbf { X } } , t )$ is given by particle statistics. In the bosonic case, we demand that $U ( { \vec { x } } _ { 1 } , { \vec { x } } _ { 2 } \dots { \vec { x } } _ { N } ) =$ $U ( \vec { x } _ { \sigma ( 1 ) } , \vec { x } _ { \sigma ( 2 ) } \ldots \vec { x } _ { \sigma ( N ) } )$ , for all particle permutations σ. In general, the functions $u _ { m } ( \mathbf { r } ; t )$ can have an arbitrarily complex dependence on the m particle coordinates, which can prove problematic for practical applications. Nonetheless, a simplified functional dependence can often be imposed, resulting from the two-body character of the interactions in the original Hamiltonian. For $m \geq 3 .$ $u _ { m } ( \mathbf { r } ; t )$ can be conveniently factorized in terms of general two-particle vector and tensor functions following Ref. [26]. Details of this approach and the present implementation are presented in Appendix A and F.

An appealing property of the many-body expansion (2) is that it is able to describe intrinsically non-local correlations in space. For instance the two-body $u _ { 2 } ( \vec { x } _ { i } , \vec { x } _ { j } ; t )$ , as well as any m-body function $u _ { m } ( \mathbf { r } ; t )$ , can be long range in the particle separation $\lvert \vec { x } _ { i } - \vec { x } _ { j } \rvert$ . This non-local spatial structure allows for a correct description of gapless phases, where a two-body expansion may already capture all the universal features, in the sense of the renormalization group approach [39]. This is in contrast with the MPS decomposition of the wave-function, which is intrinsically local in space.

## B. Time-Dependent Variational Monte Carlo

The time evolution of the variational state (2) is entirely determined by the time-dependence of the Jastrow functions $u _ { m } ( \mathbf { r } ; t )$ . In order to establish optimal equations of motion for the variational parameters, we start noticing that the functional derivative of $U ( \mathbf { X } ; t )$ with respect to the variational, complex-valued, Jastrow functions $u _ { m } ( \mathbf { r } ; t )$ ,

$$
\frac { \delta U ( \mathbf { X } ; t ) } { \delta u _ { m } ( \mathbf { r } ; t ) } \equiv \rho _ { m } ( \mathbf { r } ) ,\tag{3}
$$

yields the m-body density operators

$$
\rho _ { m } ( { \bf r } ) = \frac { 1 } { m ! } \sum _ { i _ { 1 } \neq i _ { 2 } \neq . . . i _ { m } } \prod _ { j } \delta ( x _ { i _ { j } } - r _ { j } ) .\tag{4}
$$

The expectation values of the operators $\rho _ { m }$ over the state $| \Phi ( t ) \rangle$ give the instantaneous m-body correlations. For instance, $\begin{array} { r } { \langle \rho _ { 2 } ( r _ { 1 } , r _ { 2 } ) \rangle _ { t } ~ = ~ \sum _ { i < j } \langle \delta ( x _ { i } - r _ { 1 } ) \delta ( x _ { j } - r _ { 2 } ) \rangle _ { t } } \end{array}$ where $\langle . . . \rangle _ { t } = \langle \Phi ( t ) | . . . | \Phi ( t ) \bar { \rangle } / \langle \Phi ( t ) | \Phi ( t ) \rangle$ , is proportional to the two-point density-density correlation function.

We can then express the time derivative of the truncated variational state ${ \cal U } ^ { ( M ) } ( { \bf X } ; t )$ using the functional derivatives, Eq. (3), as a sum of the few-body density operators up to the truncation M , i.e.

$$
\partial _ { t } U ^ { ( M ) } ( { \bf { X } } , t ) = \sum _ { m = 1 } ^ { M } \int d { \bf { r } } \rho _ { m } ( { \bf { r } } ) \partial _ { t } u _ { m } ( { \bf { r } } ; t ) .\tag{5}
$$

The exact wave-function satisfies the Schrödinger equation $\begin{array} { r l r } { i \partial _ { t } U ( \mathbf { X } ; t ) } & { { } = } & { E _ { \mathrm { l o c } } ( \mathbf { X } , t ) } \end{array}$ , where $\begin{array} { r l } { { \cal E } _ { \mathrm { l o c } } ( { \bf \bar { X } } ; t ) } & { { } = } \end{array}$ $\frac { \langle \mathbf { X } | \mathcal { H } | \Phi ( t ) \rangle } { \langle \mathbf { X } | \Phi ( t ) \rangle }$ is the so-called local energy. The optimal time evolution of the truncated Bijl-Dingle-Jastrow-Feenberg expansion (2) can be derived imposing the Dirac-Frenkel time-dependent variational principle [40, 41]. In geometrical terms, this amounts to minimizing the Hilbert-space norm of the residuals $R ^ { ( M ) } ( { \bf { X } } ; t ) \equiv $ $\left\| i \partial _ { t } U ^ { ( M ) } ( { \mathbf { X } } ; t ) - E _ { \mathrm { { l o c } } } ^ { ( M ) } ( { \mathbf { X } } ; t ) \right\|$ , thus yielding a variational many-body state as close as possible to the exact one [42]. The minimization can be performed explicitly and yields a closed set of integro-differential equation for the Jastrow functions $u _ { m } ( \mathbf { r } ; t )$

$$
\sum _ { p = 1 } ^ { M } \int d \mathbf { r } ^ { \prime } \frac { \delta \left. \rho _ { m } ( \mathbf { r } ) \right. _ { t } } { \delta u _ { p } ( \mathbf { r } ^ { \prime } ; t ) } \partial _ { t } u _ { p } ( \mathbf { r } ^ { \prime } ; t ) = - i \frac { \delta \left. \mathcal { H } \right. _ { t } } { \delta u _ { m } ( \mathbf { r } ; t ) } .\tag{6}
$$

In practice, these equations are numerically solved for the time derivatives $\bar { \partial } _ { t } u _ { p } ( \mathbf { r } ^ { \prime } ; t )$ at each time step t. The expectation values taken over the time-dependent state, $\langle \cdot \rangle _ { t }$ , which enter Eq. (6), are found via a stochastic sampling of the probability distribution $\Pi ( { \bf X } , t ) = \left| \Phi ( { \bf X } , t ) \right| ^ { 2 }$ This is efficiently achieved by means of the Metropolis-Hastings algorithm, as per conventional Monte Carlo schemes (See Appendix B for details). It then yields the full time evolution of the truncated Bijl-Dingle-Jastrow-Feenberg state (2) after time integration.

The t-VMC approach as formulated here provides, in principle, an exact description of the real-time dynamics of the N-body system. The essential approximation lies however in the truncation of the Bijl-Dingle-Jastrow-Feenberg expansion to the M most relevant terms. In practical applications, the $M = 2 \mathrm { o r } M = 3$ truncation is often sufficient. Systematic improvement beyond $M = 3$ is possible [26], but may require a substantial computational effort.

## III. LIEB-LINIGER MODEL

As a first application of the continuous-space t-VMC approach, we consider the Lieb-Liniger model [44]. On one hand, some exact results and numerical data are available, allowing us to benchmark the Jastrow expansion and the t-VMC approach. On the other hand, several aspects of the out-of-equilibrium dynamics of this model are unknown, which we compute here for the first time using t-VMC.

The Lieb-Liniger model describes N interacting bosons in one dimension with contact interactions. It corresponds to the Hamiltonian (1) with $v _ { 1 } ( x ) ~ = ~ 0$ and $v _ { 2 } ( x , y ) = g \delta ( x - y )$ , where $g$ is the coupling constant. Here, we consider periodic boundary conditions over a ring of length $L ,$ and the particle density $n \ = \ N / L$ The density dependence is as usual expressed in terms of the dimensionless parameter $\gamma = m g / \hbar ^ { 2 } n$ . The Lieb-Liniger model is the prototypal model of continuous onedimensional strongly-correlated gas exactly solvable by Bethe ansatz [44]. This model is experimentally realized in ultracold atomic gases strongly confined in onedimensional optical traps, and several studies on out-ofequilibrium physics have been already realized [14, 45– 49].

As a result of translation invariance, we have $u _ { 1 } ( x ; t ) =$ const, and the first non-trivial term in the many-body expansion (2) is the two-body, translation invariant, function $u _ { 2 } ( | x - y | ; t )$ . To compute the functional derivatives of the many-body wave function, we proceed with a projection of the continuous Jastrow fields $u _ { m } ( \mathbf { r } ; t )$ onto a finite basis set. Here, we have found convenient to represent both the field Hamiltonian (1) and the Jastrow functions, onto a uniform mesh of spacing $^ { a , }$ leading to $\begin{array} { r } { n _ { \mathrm { v a r } } = \frac { L } { 2 a } } \end{array}$ variational parameters for the 2-body Jastrow term. In the following our results are extrapolated to the continuous limit, corresponding to $a  0 .$ The finitebasis projection as well as the numerical time-integration of Eq. (6) are detailed in the Appendix C and benchmarked against exact diagonalization results in the three particle case in Appendix E.

## A. Ground-state properties

To assess the quality of truncated Bijl-Dingle-Jastrow-Feenberg expansions for the Lieb-Liniger model, we start our analysis considering ground-state properties. An exact solution can be found from the Bethe ansatz (BA) and gives access to exact ground state energies and local properties [44]. Other non-local properties are substantially more difficult to extract from the BA solution, and unbiased results for ground state correlation functions have not been reported so far. In order to determine the best possible variational description of the ground-state within our many-body expansion, different strategies are possible. A first possibility is to consider the imaginarytime evolution $| \dot { \Psi } ( \tau ) \rangle = \dot { e } ^ { - \tau { \mathcal { H } } } | \Psi _ { 0 } \rangle$ which systematically converges to the exact ground-state in the limit $\tau \gg \Delta _ { 1 }$ where $\Delta _ { 1 } = E _ { 1 } - E _ { 0 }$ is the gap with the first excited state on a finite system and provided that the trial state $\Psi _ { 0 }$ is non-orthogonal to the exact ground state. Imaginarytime evolution in the variational subspace can be implemented considering the formal substitution $t \to - i \tau$ in the t-VMC equations (6). The resulting equations are equivalent to the stochastic reconfiguration approach [50]. However, direct minimization of the variational energy can be significantly more efficient, in particular for systems becoming gapless in the thermodynamic limit, where $\Delta _ { 1 } \sim \mathrm { p o l y } ( 1 / N )$ . Given the gapless nature of the Lieb-Liniger model, we have found computationally more efficient to adopt a Newton method to minimize the energy variance [51].

<!-- image-->

<!-- image-->

<!-- image-->  
Figure 1. Ground-state properties of the Lieb-Liniger model as obtained from different variational approaches : (a) relative accuracy of the ground-state energy , (b) one-body density matrix, (c) density-density pair correlation functions. $\dot { U } ^ { ( 2 ) }$ and $U ^ { ( 3 ) }$ denote results for the 2, and 3-body expansion (present work), $U _ { \mathrm { A G } } ^ { ( 2 ) }$ is the parametrized 2-body Jastrow state of Ref. [43], $\mathrm { c { - } M P S } .$ 1 results are from Ref. [23], and ${ \mathrm { c } } { \mathrm { - M P S } } _ { 2 }$ are those very recently reported in [24]. Distances r are in units of the inverse density $1 / n .$ Our variational results have been obtained for $N = 1 0 0$ particles. Finite-size corrections on local quantities are negligible, and very mildly affect the reported large-distances correlations. Overall statistical errors are of the order of symbol sizes, for (a), and line widths, for (b) and (c).

For the ground state, the many-body expansion truncated at $M = 2$ is exact not only in the non-interacting limit $\gamma = 0$ but also in the Tonks-Girardeau limit $\gamma \to \infty$ In this fermionized limit the wave-function can be written as the modulus of a Vandermonde determinant of plane waves, corresponding to the two-body Jastrow function $u _ { 2 } ( r ) = \log { ( \sin { r \pi } / L ) }$ [52]. To assess the overall quality of pair wave functions for ground state properties, we start comparing the variational ground-state energies $E$ obtained for $M \ = \ 2$ with the exact BA result. In Fig. 1(a) we show the relative error $\Delta E / E$ as a function of in the interaction strength γ. We find that the relative error is lower than $1 0 ^ { - 4 }$ for all values of $\gamma$ and the accuracy of our two-body Jastrow function is superior to previously published variational results based either on c-MPS [23] or VMC [43] [see Fig. 1(a) for a quantitative comparison]. Notice that the improvement with respect to previous VMC results is due to the larger variational freedom of our $u _ { 2 } ( r )$ function, which is not restricted to any specific functional form as done in Ref. [43].

Even though the accuracy reached by the two-body Jastrow function may already be sufficient for most practical purposes for all values of $\gamma _ { ; }$ , we have also considered higher order terms with $M = 3 .$ , shown in Fig. $\mathrm { ( a ) }$ . The introduction of the third-order term yields a sizable improvement such that the maximum error is about three orders of magnitude smaller than original c-MPS results [23], and feature a similar accuracy of recently reported c-MPS results [24]. Overall, our approach reaches a precision on a continuous-space system, which is comparable to state-of-the-art MPS/DMRG results for gapless systems on a lattice [53].

Finally, to further assess the quality of our ground state ansatz beyond the total energy, we have also studied non-local properties of the ground state wave function, which are not accessible by existing exact BA methods. In Figs. 1(b) and (c) we show, respectively, our results for the off-diagonal part of the one-body density matrix, $g _ { 1 } ( r ) \propto \langle \Psi ^ { \dagger } ( r ) \Psi ( 0 ) \rangle$ , and for the pair correlation function, $g _ { 2 } ( r ) \propto \langle \Psi ^ { \dagger } ( r ) \Psi ( r ) \Psi ^ { \dagger } ( 0 ) \Psi ( 0 ) \rangle$ , where $\Psi ( r )$ is the bosonic field operator. We find an overall excellent agreement with the results that have been obtained with c-MPS in Ref. [24], except for some small deviations at large values of r which we attribute to residual finite-sizeeffects in our approach. We found that the addition of the 3-body terms does not change significantly correlation functions. Already at the 2-body level, the present results are statistically indistinguishable from exact results obtained using our implementation [54] of the worm algorithm [55] and for the same system (not shown).

<!-- image-->

<!-- image-->  
Figure 2. Time-dependent expectation value of local two-body correlations after a quantum quench from a non-interacting state, $\gamma _ { i } = 0$ , to $\gamma _ { f } : \mathbf { \Gamma } ( \mathrm { a } )$ t-VMC results are compared with BA results obtained for a small number of particles [56, 57]. The correlation function is rescaled to have $g _ { 2 } ( 0 , 0 ) = 1 ;$ (b) t-VMC results for $N = 1 0 0$ particles and $\gamma _ { f } = 1 , 2 , 4 , 8$ (from top to bottom) compared to the quench action predictions from Ref. [58] (dashed lines), to the Boltzmann thermal averages at the effective temperature $T ^ { \star }$ (dotted-dashed lines), and to the GGE thermal averages prediction (rightmost dashed lines). Statistical error bars on t-VMC data are of the order of lines width.

## B. Quench dynamics

Having assessed the quality of the ansatz for local and non-local ground-state properties, we now turn to the study of the out-of-equilibrium properties of the Lieb-Liniger model. We focus on the description of the unitary dynamics induced by a global quantum quench of the interaction strength, from an initial value $\gamma _ { i }$ to a final value $\gamma _ { f }$ . Exact BA results are available only in the case of a non-interacting initial state $( \gamma _ { i } ~ = ~ 0 )$ . Even in this case, the dynamical BA equations can be exactly solved only for a modest number of particles with further truncation in the number of energy eigen-modes [56, 58], $N \lesssim 1 0$ , since the complexity of the BA solution increases exponentially with the number of particles. Simplifications in the thermodynamic limit are exploited by the quench action [59], and have been recently applied to quantum quenches starting from a non-interacting initial state [58]. In the following we first compare our t-VMC results to these existing results, and then present new results for quenches following a non-vanishing initial interaction strength.

To assess the quality of the time-dependent wave function we compare our results for the evolution of local density-density correlations, $g _ { 2 } ( 0 , t )$ , with the truncated BA results obtained in Refs. [56, 57] for a small number of particles, $N \simeq 6$ . Appendix E provides also further validation of our method accessing $g _ { 2 } ( r , t )$ at nonvanishing distances. The comparison shown in Fig. 2(a) shows an overall good agreement. The t-VMC and BA results are indistinguishable for weak interactions $( \gamma _ { f } = 1 )$ For larger interactions, we notice systematic but small differences between BA and t-VMC with $M \ = \ 2$ or $M = 3$ These differences amount to a small increase in the amplitude of the oscillations. This effect tends to increase with the interaction strength, being hardly visible for $\gamma = 1$ and more pronounced for $\gamma = 1 0$ . However, these oscillations result at large times from the discrete mode structure due to the very small number of particles. They vanish in the physical thermodynamic limit. In turn, the comparison at small particle numbers indicates an accuracy better than a few percent for timeaveraged quantities in the asymptotic large time limit up to $\gamma = 1 0$ , with results at $M = 3$ systematically improving on the $M = 2$ case. Concerning small and intermediate time scales, we do not observe systematic deviations between the t-VMC results and the BA solution. In particular, the relaxation times are remarkably well captured by the t-VMC approach. On the basis of this comparison and of the comparison for three particles with exact diagonalization results presented in Appendix E, we conclude that t-VMC allows accurate quantitative studies of both the relaxation and equilibration dynamics. This careful benchmarking now allows us to confidently apply the t-VMC approach regimes that are inaccessible to exact BA namely large but finite N and long times, as well as the case of non-vanishing initial interactions.

Let us consider relaxation of density correlations for a large number of particles, close to the thermodynamic limit (here we use N = 100). As shown in Fig. 2(b), we notice that the amplitudes of the large-time oscillations, attributed to the discrete mode spectrum, are now drastically suppressed compared to the quenches with $N = 6$ After an initial relaxation phase, the quantity $g _ { 2 } ( 0 , t )$ approaches a stationary value. Comparing our curves with those obtained with the quench action method, we find a qualitatively good agreement, albeit a general tendency to underestimate the quench action predictions is observed.

We now turn to quenches from interacting initial states $( \gamma _ { i } \neq 0 )$ to different interacting final states for which no results have been obtained by means of exact BA nor simplified quench action method so far. In Fig. 3 we show the asymptotic equilibrium values obtained with our t-VMC approach for quantum quenches from $\gamma _ { i } = 1$ (Left panel), and $\gamma _ { i } = 4$ (Right panel) to several values of $\gamma _ { f } .$ . Since, by the variational theorem, the ground state of $\mathcal { H } _ { i }$ gives an upper bound for the ground state energy of $\mathcal { H } _ { f }$ , the system is pushed into a linear combination of excited states of the final hamiltonian. For systems able to thermalize to the Boltzmann ensemble (BE), relaxation to a stationary state described by the density matrix $\rho _ { T ^ { \star } } = e ^ { - \mathcal { H } _ { f } / T ^ { \star } }$ , at an effective temperature $T ^ { \star }$ , would occur. Comparing the stationary value, $\bar { g _ { 2 } } ( 0 )$ , of our t-VMC calculations at long times to the thermal values of the pair correlation functions, $g _ { 2 } ^ { T ^ { \star } } ( 0 )$ , a necessary condition for simple Boltzmann thermalization is given by $\bar { g _ { 2 } } = g _ { 2 } ^ { T ^ { \star } }$ The effective temperature, $T ^ { \star }$ , is determined by imposing the energy expectation value of the final Hamiltonian, $\mathcal { H } _ { f }$ , in the ground state, $\Phi _ { 0 } ( \gamma _ { i } )$ of the initial Hamiltonian, $\langle \mathcal { H } _ { \mathrm { f } } \rangle _ { T ^ { \star } } = \langle \Phi _ { 0 } ( \gamma _ { i } ) | \mathcal { H } _ { \mathrm { f } } | \Phi _ { 0 } ( \gamma _ { i } ) \rangle$ . Here, the thermal expectation value, $\langle \mathcal { H } _ { \mathrm { f } } \rangle _ { T } ,$ ? at the equilibrium temperature $T ^ { \star }$ is computed from the Yang-Yang BA equations [60]. The quantity $\langle \mathcal { H } _ { \mathrm { f } } \rangle _ { T ^ { \star } }$ ? then depends on a single parameter $T ^ { \star }$ , that is fitted to match the value of $\left. \Phi _ { 0 } ( \gamma _ { i } ) \vert \mathcal { H } _ { \mathrm { f } } \vert \Phi _ { 0 } ( \gamma _ { i } ) \right.$

As shown in Fig.2-(b), Boltzmann thermalization certainly does not occur in the case for the Lieb-Liniger model when quenching from a non-interacting state, $\gamma _ { i } =$ $0 ,$ , where we find $\bar { g _ { 2 } } \ \bar { \neq } \ g _ { 2 } ^ { T ^ { \star } }$ . This can be understood in terms of the existence of dynamically conserved charges (beyond energy and density conservation) which can yield an equilibrium value substantially different from the BE prediction. In particular, it is widely believe that the Generalized Gibbs Ensemble (GGE) is the correct thermal distribution approached after the quench [61, 62]. Several constructive approaches for the GGE have been put forward in past years [8, 9, 58], and the quench action predictions reported in Fig.2-(b) converge to the GGE predictions for the thermal values. In Fig.2-(b) we also show the thermal GGE values $g _ { 2 } ^ { \mathrm { G G E } }$ (rightmost dashed lines), and notice that our results are much closer to the GGE predictions than the simple BE. Deviations from the asymptotic GGE results are observed at large $\gamma _ { f } ,$ a regime in which the accuracy of our approach is still sufficient to resolve the difference between the BE and the long-term equilibration value.

For correlated initial ground states, $\gamma _ { i } \neq 0 ,$ , GGE predictions are fundamentally harder to obtain than for the non-interacting initial states, and the BE is the only reference thermal distribution we can compare with at this stage. From our results we observe that the difference between $\bar { g _ { 2 } }$ and the simple BE prediction, $g _ { 2 } ^ { T ^ { \star } }$ , is quantitatively reduced, see Fig. 3. In particular, for $\gamma _ { i } = 4 .$ the stationary values $\bar { g _ { 2 } }$ are quantitatively close to the ones predicted by the Boltzmann thermal distribution at the effective temperature $T ^ { \star }$ . Even though this quantitative agreement is likely to be coincidental, the regimes of parameters quenches studied here provide guidance for future experimental studies. In particular, it will be of great interest to understand whether a cross-over from a strongly non-Boltzmann to a close-to-Boltzmann thermal behavior might occur as a function of the initial interaction strength also for other local observables.

<!-- image-->  
Figure 3. Time-dependent expectation value of local twobody correlations after a quantum quench from the interacting ground state at $\gamma _ { i } = 1 \ \mathrm { ( a ) }$ , and $\gamma _ { i } = 4 ~ \mathrm { ( b ) }$ , long-term dynamical averages (red continuous lines) are compared to thermal averages at the effective temperature set by energy conservation (black dashed lines). Uncertainties on the thermal averages of the order of lines width, and are larger for small values of $\gamma _ { f }$

## IV. CONCLUSIONS

In this paper we have introduced a novel approach to the dynamics of strongly-correlated quantum systems in continuous space. Our method is based on correlated many-body wave-function systematically expanded in terms of reduced m-body Jastrow functions. The unitary dynamics in the subspace of these correlated states was realized using time-dependent variational Monte Carlo. We have demonstrated the possibility or performing calculations up to the three-body level, $m \le 3$ , for the Lieb-Liniger model, for both static and dynamical properties. The improvement from m = 2 to m = 3 provides an internal criterium to judge the validity of our results whenever exact results are unavailable. Benchmarking t-VMC with exact or numerical approaches whenever available, we have found a very good agreement with existing results. For static properties, our approach is at the level of state-of-the-art MPS techniques in lattice systems and of latest c-MPS results for interacting gases. For dynamical properties, we have investigated for the first time general interaction quenches which are at the moment unaccessible to Bethe-ansatz approaches. Since the general structure of our t-VMC method does not depend on the dimensionality of the system, it can be directly applied to bosonic systems in higher dimensions with a polynomial increase in computational cost. The methods presented here therefore pave the way to accurate outof-equilibrium dynamics of two- and three-dimensional quantum gases and fluids beyond mean field approximations.

## ACKNOWLEDGMENTS

We acknowledge discussions with J. De Nardis, M. Dolfi, M. Fagotti, T. Osborne, and M. Troyer. We thank M. Ganahl for providing us the c-MPS results in Fig. 1, J. Zill for the Bethe ansatz results in Fig. 2 (a), and J. De Nardis for the quench actions results in Fig. 2 (b). This research was supported by the Marie Curie IEF program (FP7/2007-2013 - Grant Agreement No. 327143), the European Research Council Starting Grant "ALoGlaDis" (FP7/2007-2013 Grant Agreement No. 256294) and Advanced Grant "SIM-COFE" (FP7/2007-2013 Grant Agreement No. 290464), the European Commission FET-Proactive QUIC (H2020 grant No. 641122), the French ANR-16-CE30-0023-03 (THERMOLOC), and the Swiss National Science Foundation through NCCR QSIT. It was performed using HPC resources from GENCI-CCRT/CINES (Grant c2015056853).

[1] Nature Physics Insight on Non-equilibrium physics, Nature Physics 11, 103 (2015).

[2] D. Ceperley, Reviews of Modern Physics 67, 279 (1995).

[3] L. Pollet, Reports on Progress in Physics 75, 094501 (2012).

[4] T. Plisson, B. Allard, M. Holzmann, G. Salomon, A. Aspect, P. Bouyer, and T. Bourdel, Physical Review A 84, 061606 (2011).

[5] N. G. Berloff and B. V. Svistunov, Physical Review A 66, 013603 (2002).

[6] N. Navon, A. L. Gaunt, R. P. Smith, and Z. Hadzibabic, Nature 539, 72 (2016).

[7] M. Gaudin, La fonction d’onde de Bethe (Paris; New Tork: Masson, 1983).

[8] J.-S. Caux and R. M. Konik, Physical Review Letters 109, 175301 (2012).

[9] J.-S. Caux and F. H. L. Essler, Physical Review Letters 110, 257203 (2013).

[10] M. Collura, S. Sotiriadis, and P. Calabrese, Physical Review Letters 110, 245301 (2013).

[11] S. Trotzky, Y.-A. Chen, A. Flesch, I. P. McCulloch, U. Schollwoeck, J. Eisert, and I. Bloch, Nature Physics 8, 325 (2012).

[12] T. Langen, R. Geiger, M. Kuhnert, B. Rauer, and J. Schmiedmayer, Nature Physics 9, 640 (2013).

[13] D. Greif, G. Jotzu, M. Messer, R. Desbuquois, and T. Esslinger, Physical Review Letters 115, 260401 (2015).

[14] M. Gring, M. Kuhnert, T. Langen, T. Kitagawa, B. Rauer, M. Schreitl, I. Mazets, D. A. Smith, E. Demler, and J. Schmiedmayer, Science 337, 1318 (2012).

[15] M. Cominotti, D. Rossini, M. Rizzi, F. Hekking, and A. Minguzzi, Physical Review Letters 113, 025301 (2014).

[16] P. Calabrese and J. Cardy, Physical Review Letters 96, 136801 (2006).

[17] S. R. White, Physical Review Letters 69, 2863 (1992).

[18] S. R. White and A. E. Feiguin, Physical Review Letters 93, 076401 (2004).

[19] A. J. Daley, C. Kollath, U. Schollwock, and G. Vidal, Journal of Statistical Mechanics-Theory and Experiment , P04005 (2004).

[20] F. B. Anders and A. Schiller, Physical Review Letters 95, 196801 (2005).

[21] M. Dolfi, B. Bauer, M. Troyer, and Z. Ristivojevic, Physical Review Letters 109, 020604 (2012).

[22] M. Dolfi, A. Kantian, B. Bauer, and M. Troyer, Physical Review A 91, 033407 (2015).

[23] F. Verstraete and J. I. Cirac, Physical Review Letters 104, 190405 (2010).

[24] M. Ganahl, J. Rincon, and G. Vidal, arXiv:1611.03779 (2016).

[25] C. J. Umrigar, J. Toulouse, C. Filippi, S. Sorella, and R. G. Hennig, Physical Review Letters 98, 110201 (2007).

[26] M. Holzmann, B. Bernu, and D. M. Ceperley, Physical Review B 74, 104510 (2006).

[27] M. Taddei, M. Ruggeri, S. Moroni, and M. Holzmann, Physical Review B 91, 115106 (2015).

[28] G. Carleo, F. Becca, M. Schiro, and M. Fabrizio, Scientific Reports 2, 243 (2012).

[29] G. Carleo, F. Becca, L. Sanchez-Palencia, S. Sorella, and M. Fabrizio, Physical Review A 89, 031602 (2014).

[30] L. Cevolani, G. Carleo, and L. Sanchez-Palencia, Physical Review A 92, 041603 (2015).

[31] B. Blaß and H. Rieger, Scientific Reports 6, 38185 (2016).

[32] G. Carleo and M. Troyer, Science 355, 602 (2017).

[33] K. Ido, T. Ohgoe, and M. Imada, Physical Review B 92, 245106 (2015).

[34] We have conveniently set the particle mass and the reduced Planck constant to unity.

[35] A. Bijl, Physica 7, 869 (1940).

[36] R. B. Dingle, The London, Edinburgh, and Dublin Philosophical Magazine and Journal of Science 40, 573 (1949).

[37] R. Jastrow, Physical Review 98, 1479 (1955).

[38] E. Feenberg, Theory of quantum fluids, Pure and applied physics (Academic Press, 1969).

[39] C. L. Kane, S. Kivelson, D. H. Lee, and S. C. Zhang, Physical Review B 43, 3255 (1991).

[40] P. a. M. Dirac, Mathematical Proceedings of the Cambridge Philosophical Society 26, 376 (1930).

[41] I. Frenkel, Wave Mechanics: Advanced General Theory, The International series of monographs on nuclear energy: Reactor design physics No. v. 2 (The Clarendon Press, 1934).

[42] The natural norm induced by a quantum Hilbert space is the Fubini-Study norm, which is gauge invariant and therefore insensitive to the unknown normalizations of the quantum states we are dealing with here.

[43] G. E. Astrakharchik and S. Giorgini, Physical Review A 68, 031602 (2003).

[44] E. H. Lieb and W. Liniger, Physical Review 130, 1605 (1963).

[45] H. Moritz, T. Stoferle, M. Kohl, and T. Esslinger, Physical Review Letters 91, 250402 (2003).

[46] T. Kinoshita, T. Wenger, and D. S. Weiss, Nature 440, 900 (2006), wOS:000236736700033.

[47] J. P. Ronzheimer, M. Schreiber, S. Braun, S. S. Hodgman, S. Langer, I. P. McCulloch, F. Heidrich-Meisner, I. Bloch, and U. Schneider, Physical Review Letters 110, 205301 (2013).

[48] B. Fang, G. Carleo, A. Johnson, and I. Bouchoule, Physical Review Letters 113, 035301 (2014).

[49] G. Boéris, L. Gori, M. D. Hoogerland, A. Kumar, E. Lucioni, L. Tanzi, M. Inguscio, T. Giamarchi, C. D’Errico, G. Carleo, G. Modugno, and L. Sanchez-Palencia, Physical Review A 93, 011601 (2016).

[50] S. Sorella, Physical Review B 64, 024512 (2001).

[51] C. J. Umrigar and C. Filippi, Physical Review Letters 94, 150201 (2005).

[52] M. Girardeau, Journal of Mathematical Physics 1, 516 (1960).

[53] P. Pippan, S. R. White, and H. G. Evertz, Physical Review B 81, 081103 (2010).

[54] G. Carleo, G. Boéris, M. Holzmann, and L. Sanchez-Palencia, Physical Review Letters 111, 050406 (2013).

[55] M. Boninsegni, N. Prokof’ev, and B. Svistunov, Physical Review Letters 96, 070601 (2006).

[56] J. C. Zill, T. M. Wright, K. V. Kheruntsyan, T. Gasenzer, and M. J. Davis, Physical Review A 91, 023611 (2015).

[57] J. C. Zill, T. M. Wright, K. V. Kheruntsyan, T. Gasenzer, and M. J. Davis, New Journal of Physics 18, 045010 (2016).

[58] J. De Nardis, L. Piroli, and J.-S. Caux, Journal of Physics A: Mathematical and Theoretical 48, 43FT01 (2015).

[59] J.-S. Caux and F. H. L. Essler, Physical Review Letters 110, 257203 (2013).

[60] C. N. Yang and C. P. Yang, Journal of Mathematical Physics 10, 1115 (1969).

[61] M. Rigol, Physical Review Letters 116, 100601 (2016).

[62] M. Kollar and M. Eckstein, Physical Review A 78, 013626 (2008).

[63] M. Holzmann, D. M. Ceperley, C. Pierleoni, and K. Esler, Physical Review E 68, 046707 (2003).

## Appendix A: Functional Structure of Many-Body Terms

The local residuals ${ \cal { R } } ^ { ( M ) } ( { \bf { X } } ; t ) \ = \ i \partial _ { t } U ^ { ( M ) } ( { \bf { X } } ; t ) \ - \ $ $E _ { \mathrm { l o c } } ^ { ( M ) } ( { \bf X } ; t )$ are vanishing if the Schrödinger equation is exactly satisfied by the many-body wave-function truncated at some order M. The local energy $E _ { \mathrm { l o c } } ^ { ( M ) } ( { \bf X } , t )$ however, may contain effective interaction terms involving a number of bodies larger than M, which leads to a systematic error in the truncation. However, the structure of these additional terms stemming from the local energy can be systematically used to deduce the functional structure of the higher order terms. For example, the one-body truncated local energy reads, for one-

dimensional particles,

$$
\begin{array} { c } { { \displaystyle E _ { \mathrm { l o c } } ^ { ( 1 ) } ( { \bf X } ; t ) = - \frac { 1 } { 2 } \sum _ { i } \left\{ \left[ \partial _ { x _ { i } } u _ { 1 } ( x _ { i } ; t ) \right] ^ { 2 } + \partial _ { x _ { i } } ^ { 2 } u _ { 1 } ( x _ { i } ; t ) \right\} + } } \\ { { + \sum _ { i } v _ { 1 } ( x _ { i } ) + \displaystyle \frac { 1 } { 2 } \sum _ { i \neq j } ^ { N } v _ { 2 } ( x _ { i } , x _ { j } ) , \qquad ( \mathrm { A 1 } } } \end{array}
$$

and contains a 2-body term which cannot be accounted for exactly by $u _ { 1 }$ . Introduction of a symmetric two-body Jastrow factor $u _ { 2 } ( x _ { i } , x _ { j } ; t )$ , then leads to

$$
\begin{array} { l } { { \displaystyle E _ { \mathrm { l o c } } ^ { ( 2 ) } ( { \bf X } ; t ) = E _ { \mathrm { l o c } } ^ { ( 1 ) } ( { \bf X } ; t ) + } \ ~ } \\ { { \displaystyle ~ - \frac { 1 } { 2 } \sum _ { i \neq j } [ \partial _ { x _ { i } } u _ { 1 } ( x _ { i } ; t ) \partial _ { x _ { i } } u _ { 2 } ( x _ { i } , x _ { j } ; t ) ] + } \ ~ } \\ { { \displaystyle ~ - \frac { 1 } { 2 } \sum _ { i \neq j } \partial _ { x _ { i } } ^ { 2 } u _ { 2 } ( x _ { i } , x _ { j } ; t ) + } \ ~ } \\ { { \displaystyle ~ - \frac { 1 } { 2 } \sum _ { i \neq j } \sum _ { k \neq i } \partial _ { x _ { i } } u _ { 2 } ( x _ { i } , x _ { j } ; t ) \partial _ { x _ { i } } u _ { 2 } ( x _ { i } , x _ { k } ( \delta ) 2 ) } } \end{array}
$$

In the latter expression, one recognizes an effective twobody term which can be accounted for by $u _ { 2 }$ and an additional three-body term in the form of product of twobody functions. The functional form of the three-body Jastrow can be therefore deduced from this additional term and formed accordingly:

$$
u _ { 3 } ( x _ { i } , x _ { j } , x _ { k } ; t ) = \bar { u } _ { 3 } ( x _ { i } , x _ { j } ; t ) \bar { u } _ { 3 } ( x _ { j } , x _ { k } ; t ) ,\tag{A3}
$$

with two-body functions $\bar { u } _ { 3 } ( x _ { i } , x _ { j } ; t )$ containing new variational parameters to be determined. Upon pursuing this approach, the expansion can be systematically pushed to higher orders and the functional structure of the higher order functions inferred. The same constructive approach we have discussed here is also valid for the Schrödinger equation in imaginary-time $\partial _ { \tau } U ( { \bf { X } } ; \tau ) = - E _ { \mathrm { { l o c } } } ( { \bf { x } } ; \tau )$ 2 and has been successfully used to infer the functional structure for ground-state properties [63].

## Appendix B: Monte Carlo Sampling

In order to solve the t-VMC equations of motion, Eq. (6), expectation values of some given operator  need to be computed over the many-body wave-function $\Phi ( { \mathbf { X } } , t )$ This is achieved by means of Monte Carlo sampling of the probability distribution $\Pi ( { \bf X } ) = | \Phi ( { \bf X } ) | ^ { 2 }$ (In the following we omit explicit reference to the time t, assuming that all expectation values are taken over the wave-function at a given fixed time). An efficient way of sampling the given probability distribution is to devise a Markov chain of configurations $\mathbf { X } ( 1 ) , \mathbf { X } ( 2 ) , \ldots \mathbf { X } ( N _ { c } - 1 ) , \mathbf { X } ( N _ { c } )$ which are distributed according to $\Pi ( \mathbf { X } )$ . Quantum expectation values of a given operator can then be obtained as statistical expectation values over the Markov chain as

$$
\frac {  \Phi | \mathcal { O } | \Phi  } {  \Phi | \Phi  } \simeq \frac { 1 } { N _ { c } } \sum _ { i = 1 } ^ { N _ { c } } \mathcal { O } _ { \mathrm { l o c } } ( \mathbf { X } ( i ) ) ,\tag{B1}
$$

where $\begin{array} { r } { \mathcal { O } _ { \mathrm { l o c } } ( \mathbf { X } ) = \frac { \langle \mathbf { X } | \mathcal { O } | \Phi \rangle } { \langle \mathbf { X } | \Phi \rangle } } \end{array}$ , and the equivalence is achieved in the limit $N _ { c }  \dot { \infty }$

The Markov chain is realized by the Metropolis-Hastings algorithm. Given the current state of the Markov chain, $\mathbf { X } ( i )$ , a configuration $\mathbf { X } ^ { \prime }$ is generated according to a given transition probability $T ( \mathbf { X } ( i ) \ $ $\mathbf { X } ^ { \prime } )$ . The proposed configuration is then accepted (i.e. $\mathbf { X } ( i + 1 ) = \mathbf { \bar { X } } ^ { \prime } )$ with probability

$$
\begin{array} { r l r } {  { A ( \mathbf { X } ( i ) \to \mathbf { X } ^ { \prime } ) = } } \\ & { } & { = \operatorname* { m i n } [ 1 , \frac { \Pi ( \mathbf { X } ^ { \prime } ) } { \Pi ( \mathbf { X } ( i ) ) } \frac { T ( \mathbf { X } ^ { \prime } \to \mathbf { X } ( i ) ) } { T ( \mathbf { X } ( i ) \to \mathbf { X } ^ { \prime } ) } ] , } \end{array}\tag{B2}
$$

otherwise it is rejected and $\mathbf { X } ( i + 1 ) = \mathbf { X } ( i )$

In the present t-VMC calculations we use simple transition probabilities in which a single particle is displaced, while leaving all the other particles positions unchanged. In particular, a particle index $p$ is chosen with uniform probability $1 / N$ and the position of particle $p$ is then displaced according to $x _ { p } ^ { \prime } = x _ { p } + \eta _ { \Delta }$ , where $\eta _ { \Delta }$ is a random number uniformly distributed in $[ - \frac { \Delta } { 2 } , \frac { \Delta } { 2 } ]$ The amplitude $\Delta$ is an adjustable parameter and it can be typically chosen to be of the order of the average inter-particle distance. With this choice, the transition probability is simply

$$
T ( x _ { p }  x _ { p } ^ { \prime } ) = \frac { 1 } { N \Delta } ,\tag{B3}
$$

and the acceptance probability is therefore given by the mere ratio of the probability distributions, $\frac { \mathbf { \partial } _ { \mathbf { \overline { { X } } } } ( \mathbf { \partial } _ { i } ) } { \Pi ( \mathbf { \vec { X } } ( i ) ) }$ Π(X0)

## Appendix C: Finite Basis Projection

The numerical solution of the equations of motion (6) requires the projection of the Jastrow fields $u _ { m } ( \mathbf { r } ; t )$ onto a finite basis. The continuous variable r is reduced to a finite set of $P$ values for each order $m ,$ $( m , \mathbf { r } )  ( r _ { 1 , m } , r _ { 2 , m } \ldots r _ { P , m } )$ We introduced a superindex K spanning all possible values of the discrete variables $r _ { i , m } .$ The complete set of variational parameters resulting from the projection on the finite basis can then be written as $u _ { K } ( t )$ and the associated functional derivatives read $\rho _ { K } ( t )$

The integro-differential equations (6) are then brought to the algebraic form

$$
\sum _ { K ^ { \prime } } S _ { K , K ^ { \prime } } \dot { u } _ { K ^ { \prime } } ( t ) = - i \left. E _ { \mathrm { l o c } } ( t ) \rho _ { K } ( t ) \right. ,\tag{C1}
$$

where we have introduced the Hermitian correlation matrix

$$
\begin{array} { r l r } { S _ { K , K ^ { \prime } } ( t ) = \displaystyle \frac { \partial \left. \rho _ { K } ( t ) \right. _ { t } } { \partial u _ { K ^ { \prime } } } = } & { } & \\ { = \left. \rho _ { K } ( t ) \rho _ { K ^ { \prime } } ( t ) \right. _ { t } - \left. \rho _ { K } ( t ) \right. _ { t } \left. \rho _ { K ^ { \prime } } ( t ) \right. _ { t } , } \end{array}\tag{.(C2}
$$

At a given time, all the expectations values in Eq. (C1) can be explicitly computed with the stochastic approach described in Appendix B. We are therefore left with a linear system in the $n _ { \mathrm { v a r } }$ unknowns $\dot { u } _ { K } ( t )$ , which needs to be solved at each time t.

In the presence of a large number of variational parameters, $n _ { \mathrm { v a r } }$ , the solution of the linear system can be achieved using iterative solvers e.g., conjugate gradient methods, which do not need to explicitly form the matrix S. Calling $n _ { \mathrm { i t e r } }$ the number of iterations needed to obtain a solution for the linear system, the computational cost to solve $\mathrm { ( C 1 ) }$ is $\mathcal { O } ( M \times n _ { \mathrm { v a r } } \times n _ { \mathrm { i t e r } } )$ as opposed to the $\mathcal { O } ( M \times n _ { \mathrm { v a r } } ^ { 2 } )$ operations needed by a standard solver in which the matrix S is formed explicitly. In the present work we resort to the Minimal Residual (Min-Res) method, which is a variant of the Lanczos method, working in the Krylov subspace spanned by the repeated action of the matrix S onto an initial vector. In typical applications we obtain that $n _ { \mathrm { i t e r } } \ll n _ { \mathrm { v a r } }$ and several thousands of variational parameters can be efficiently treated. This is of fundamental importance when the continuous (infinite-basis) limit must be taken, for which $n _ { \mathrm { v a r } }  \infty$

Once the unknowns $\dot { u } _ { K } ( t )$ are determined, we can solve numerically the first-order differential equations given in Eq. (6) for given initial conditions $u _ { K } ( 0 )$ . In the present work we have adopted an adaptive 4th order Runge-Kutta scheme for the integration of the differential equations.

## Appendix D: Lattice Regularization For The Lieb-Liniger Model

We consider a general wave-function $\Psi ( x _ { 1 } , x _ { 2 } \dots x _ { N } )$ for N one-dimensional particles, governed by the Lieb-Liniger Hamiltonian. By means of Variational Monte Carlo, we want to sample $| \Phi ( { \mathbf { X } } ) | ^ { 2 }$ , this is achieved via a lattice regularization, i.e.

$$
\int d { \bf X } | \Phi ( { \bf X } ) | ^ { 2 } \simeq \sum _ { l _ { 1 } , l _ { 2 } . . . l _ { N } } | \Phi ( l _ { 1 } , l _ { 2 } . . . l _ { N } ) | ^ { 2 }
$$

where $l _ { i } = \{ 0 , a , \dots L - a \}$ are discrete particle positions, a the lattice spacing, L the box size and $N _ { s } = 1 + L / a$ the number of lattice sites. As a discretized Hamiltonian we take

$$
\begin{array} { l } { { \displaystyle H _ { a } \Phi ( l _ { 1 } \dots l _ { N } ) = - \frac { \hbar ^ { 2 } } { 2 m a ^ { 2 } } \sum _ { i } \Biggl \{ \frac { 4 } { 3 } \left[ \Phi ( l _ { 1 } \dots l _ { i } - a , \dots l _ { N } ) + \right. } } \\ { { \left. \Phi ( l _ { 1 } \dots l _ { i } + a , \dots l _ { N } ) \right] - } } \\ { { \displaystyle \frac { 1 } { 1 2 } \left[ \Phi ( l _ { 1 } \dots l _ { i } - 2 a , \dots l _ { N } ) + \Phi ( l _ { 1 } \dots l _ { i } + 2 a , \dots l _ { N } ) \right] + } } \\ { { \left. \qquad - \frac { 5 } { 2 } \Phi ( l _ { 1 } \dots l _ { i } , \dots l _ { N } ) \right\} + \Phi ( l _ { 1 } \dots l _ { N } ) \frac { g } { a } \sum _ { i < j } \delta ( l _ { i } , l _ { j } ) } } \end{array}
$$

The first terms constitute just the fourth-order approximation of the laplacian via central finite differences, whereas the last term corresponds the two-body delta interaction part.

With this discretization, a two-body Jastrow factor reads

$$
u _ { 2 } ( x _ { i j } ; t ) = u _ { 2 } ( l _ { i } , l _ { j } ; t ) ,
$$

where $u _ { 2 } ( a , b ; t )$ is a time-dependent matrix of size $N _ { s } \times$ $N _ { s }$ which, in 1D and in the presence of translational symmetry, depends only on dist $( a - b )$ , i.e. it has $N _ { s } / 2$ variational parameters.

## Appendix E: Benchmark Study for $N = 3$ on a Lattice

Here, we use exact diagonalization of a Hamiltonian within a given finite basis for a quantitative test of our method. Exact diagonalization is limited to very small systems on a finite basis, and we haven chosen a system containing $N = 3$ particles on $L / a = 4 0$ lattice sites as a simple, but highly non-trivial reference. In contrast to our comparison with BA methods, all observables can be accessed by exact diagonalization and we have used the off-diagonal one body density matrix $g _ { 1 } ( r , t )$ and the pair correlation function $g _ { 2 } ( r , t )$ at different distances $r =$ $\vert x _ { 1 } - x _ { 2 } \vert$ of two particles after time t where the system is quenched from the non-interacting initial state, $\gamma _ { i } = 0 .$ to a final interaction $\gamma _ { f } > 0$ , to provide a benchmark on a more general observable.

We first benchmark the influence of the time-step lattice size discretization $\Delta t$ error on $g _ { 2 }$ . From $\mathrm { F i g . 4 ( a ) }$ we see that the t-VMC dynamics is stable over a long time and the time step error can be brought to convergence. Further, we see that for final interaction $\gamma _ { f } = 4$ the truncation at the two-body level, $U ^ { ( 2 ) }$ , introduces only a small systematic error, mainly a dephasing effect, which is almost negligible at the scale of the figure. Due to the stochastic noise of the Monte Carlo integration, t-VMC introduces additional high frequency oscillations which are, however, well separable from the deterministic propagation. The amplitude of these high frequency oscillations also quantifies the purely statistical error of our data.

Whereas exact diagonalization is limited to rather small basis sets, we can access much a larger basis within t-VMC. In Fig.4 (b) we show results within the $U ^ { ( 2 ) }$ approximation with $\dot { L } / a = 8 0$ and $L / a = 1 6 0$ with time discretization $\Delta t \sim ( \dot { L ^ { \prime } } a ) ^ { 2 }$ We see that the basis set truncation in general introduces a dephasing at large enough time.

The systematic error of $U ^ { ( 2 ) }$ increases towards quenches to stronger interaction strength and becomes more visible for $\gamma _ { f } = 8$ shown in Fig.5(a). However, even in this case, the most important effect remains to be a simple dephasing, a small shift of averaged quantities is probable, but difficult to quantify precisely. Introducing a general three-body Jastrow fields, $U ^ { ( 3 ) }$ , described in detail in Appendix F, the systematic error for $N = 3$ can be fully eliminated.

In figure ${ \bf \Pi } ( \bf { b } )$ , we also benchmark the possibility of calculating the off-diagonal one-body density matrix after a quench. Here the systematic error of $U ^ { ( 2 ) }$ is more pronounced at smaller $\gamma _ { f }$ in the long range and time regime.

From our study of the three particle problem, we conclude that truncation of the many-body wave function at the level of $U _ { 2 }$ may provide an excellent approximation for $g _ { 2 } ( r , t )$ for quenches involving not too strong interaction strengths, $\gamma \lesssim 5 .$ . The systematic error due to the $U ^ { ( 2 ) }$ truncation is mainly a dephasing at large times involving small relative errors of time averaged quantities. Similar systematic dephasing errors will occur for too large time discretization or basis set truncation.

Since our method provides a parametrization of the full wave function for a given time, many different observables can be evaluated via usual Monte Carlo methods. However, the quality of different observables may vary and depend more sensitive on the inclusion of higher order correlations $U ^ { ( n ) }$ with $n > 2$ as in the case of the single body density matrix. Although these higher order terms are computationally expensive, the scaling is not exponential, and we have explicitly shown that calculations with $n = 3$ are feasible. We notice that the computational complexity may be further reduced by functional forms adapted to the problem [26].

## Appendix F: General Structure of $U ^ { ( 3 ) }$

For a general time dependent wave function, we have to go beyond the usual ground state structure of the threebody Jastrow given in Eq. (A3). Here, we provide details of our three-body term in a general form beyond the present application in one dimension.

Introducing M basis functions, $b ^ { a } ( r ) , \ a \ = \ 1 , \dots M .$ we can introduce many-body vectors [26], $\begin{array} { r l } { B _ { i \alpha } ^ { a } } & { { } = } \end{array}$ $\textstyle \sum _ { j } \mathbf { r } _ { i j } ^ { \alpha } b ^ { a } ( r _ { i j } )$ , where $\alpha = 1 , \ldots D$ indicates the summation over directions and $i = 1 , \dots N$ The variational parameters of a general three-body structure can then be written in terms of a matrix $w _ { a b }$ , such that

$$
\sum _ { i \neq j \neq k } u _ { 3 } ( \mathbf { r } _ { 1 } , \mathbf { r } _ { 2 } , \mathbf { r } _ { 3 } ) = \sum _ { a b } w _ { a b } W ^ { a b } , \quad W ^ { a b } = \sum _ { i \alpha } B _ { i \alpha } ^ { a } B _ { i \alpha } ^ { b }\tag{F1}
$$

In order to reduce the variational parameters $( \sim M ^ { 2 } )$ we may perform a singular value decomposition of the matrix $w _ { a b }$ to reduce the effective degrees of freedom.

<!-- image-->

<!-- image-->  
Figure 4. Time-dependent expectation value of the two-body correlations after a quantum quench from a non-interacting state, $\gamma _ { i } = 0 ,$ to $\gamma _ { f } = 4$ at three different distance $| x _ { 1 } - x _ { 2 } | = 0 , L / 1 0 , L / 4$ Here, the system is on a lattice with $L / a = 4 0$ lattice sites, the full line is obtained by exact diagonalization of the Hamiltonian, the other curves are from tVMC truncated at the level of $U ^ { ( 2 ) }$ . In the left figure (a), we show the convergence with different time step discretization. On the right figure (b), we show the approach to the continuum for t-VMC simulations using $U ^ { ( 2 ) }$ for discretizations $L / a = 4 0$ , 80 and 160.

<!-- image-->

<!-- image-->  
Figure 5. (a) Time-dependent expectation value of the two-body correlations, $g _ { 2 } ( r , t )$ , after a quantum quench from a noninteracting state, $\gamma _ { i } = 0 ,$ , to $\gamma _ { f } = 8$ at three different distance $| x _ { 1 } - x _ { 2 } | = 0 , L / 1 0 , L / 4$ Here, the system is on a lattice with $L / a = 4 0$ lattice sites, the full line is obtained by exact diagonalization of the Hamiltonian, the other curves are from t-VMC truncated at the level of $U ^ { ( 2 ) }$ or $U ^ { ( 3 ) }$ . In contrast to $\gamma _ { f } = 4$ shown in Fig.4, systematic differences of $U ^ { ( 2 ) }$ compared to the exact results are more visible here, the exact dynamics is recovered by inclusion of three-body terms, $U ^ { ( 3 ) }$ , into the tVMC wave function. (b) shows the off-diagonal single particle density matrix, $g _ { 1 } ( r , t )$ , at three different distances, $r = L / 1 0 , L / 4$ and $L / 2 .$