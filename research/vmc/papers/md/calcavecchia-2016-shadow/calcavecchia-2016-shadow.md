# Metal-Insulator Transition of Solid Hydrogen by the Antisymmetric Shadow Wave Function

Francesco Calcavecchia∗

LPMMC, UMR 5493 of CNRS, Universit´e Grenoble Alpes, 38042 Grenoble, France Institute of Physics, Johannes Gutenberg-University, Staudingerweg 7, D-55128 Mainz, Germany and Graduate School of Excellence Materials Science in Mainz, Staudingerweg 9, D-55128 Mainz, Germany

Thomas D. K¨uhne†

Dynamics of Condensed Matter, Department of Chemistry,

University of Paderborn, Warburger Str. 100, D-33098 Paderborn, Germany and

Paderborn Center for Parallel Computing and Institute for Lightweight Design,

Department of Chemistry, University of Paderborn,

Warburger Str. 100, D-33098 Paderborn, Germany

(Dated: November 9, 2018)

We revisit the pressure-induced metal-insulator-transition of solid hydrogen by means of variational quantum Monte Carlo simulations based on the antisymmetric shadow wave function. In order to facilitate studying the electronic structure of large-scale fermionic systems, the shadow wave function formalism is extended by a series of technical improvements, such as a revised optimization method for the employed shadow wave function and an enhanced treatment of periodic systems with long-range interactions. It is found that the superior accuracy of the antisymmetric shadow wave function results in a significantly increased transition pressure.

## I. INTRODUCTION

In 1935 Wigner and Huntington predicted that, at very high pressure, solid molecular hydrogen will dissociate and become an atomic metallic solid [1]. Because of its relevance to astrophysics [2], but in particular due to the possible high-Tc superconductivity [3] and the existence of a metallic liquid ground state [4] , the importance to grasp metallic hydrogen can hardly be overstated [5, 6]. Due to the fact that it is still not possible to reach the static compression (> 450 GPa) required to dissociate solid hydrogen, recently alternative routes to metallic hydrogen, though at lower pressure have been proposed [7]. On the one hand, the negative slope of the melting line [8] immediately suggests the possibility of producing liquid metallic hydrogen at reduced pressure, when exposed to finite temperature [9–12]. On the other hand, due to the persistence of the molecular phase, it has been predicted that metallization through bandgap closure may be possible even in the paired state [13, 14], which would be very consequential since it facilitates potential high-Tc superconductivity in molecular metallic hydrogen [15, 16]. However, computational studies recently demonstrated that even though the pairing structure is indeed persistent over the whole pressure range of Phase III, it is more importantly throughout insulating [17–20]. This is to say that metallization due to dissociation into atomic solid hydrogen may precede eventual bandgap closure.

Thus, in this paper, we investigate the molecularatomic metal-insulator transition in solid hydrogen. Due to the small energy differences between the various phases of high-pressure hydrogen, instead of the effective singleparticle density functional theory (DFT) [21, 22], the more accurate Quantum Monte Carlo (QMC) method is employed here [23–26].

The remainder of the paper is organized as follows. In section II we outline the variational Monte Carlo method and introduce the shadow wave function, as well as its antisymmetric variant. Section III contains the computational details, whereas in section IV we describe our implementation for extended systems. The eventual results are discussed in section V. The last section is devoted to the conclusions.

## II. VARIATIONAL MONTE CARLO

Variational Monte Carlo (VMC) [27], is a QMC method that permits to approximately solve the manybody Schr¨odinger equation. The main concepts underlying VMC are the application of the Rayleigh-Ritz variational principle and importance sampled Monte Carlo (MC) to efficiently evaluate high-dimensional integrals in order to compute the total energy [28, 29]. However, in contrast to quantum-chemical electronic structure methods [30], where the computational complexity grows rapidly with the number of electrons N, the formal scaling of VMC is similar to that of effective singleparticle theories such as Hartree-Fock (HF) or DFT [31]. Furthermore, as many-body correlation effects are explicitly taken into account by a prescribed trial wave function (WF), VMC is throughout more accurate than typical mean-field techniques and allows to treat even strongly correlated systems.

Nevertheless, since the exact WF of the electronic ground state is generally unknown, it is approximated by a trial WF $\psi _ { \mathrm { T } } ( R , \alpha )$ , where $R \equiv \left( \mathbf { r } _ { 1 } , \mathbf { r } _ { 2 } , \ldots , \mathbf { r } _ { N } \right)$ are the particle coordinates. The variational parameters $\alpha \equiv ( \alpha _ { i } ) _ { i = 1 , \dots n } .$ which corresponds to the lowest variational energy

$$
E = \frac { \int d R \psi _ { \mathrm { T } } ^ { * } ( R , \alpha ) H \psi _ { \mathrm { T } } ( R , \alpha ) } { \int d R \psi _ { \mathrm { T } } ^ { * } ( R , \alpha ) \psi _ { \mathrm { T } } ( R , \alpha ) } ,\tag{1}
$$

represents the best possible approximation of the electronic ground state within the given trial WF, while H is the system’s Hamiltonian. Therefore, the accuracy of a VMC simulation depends critically on how well the particular trial WF mimics the exact ground state WF.

For the purpose to efficiently evaluate the highdimensional integral of Eq. (1), it is convenient to rewrite it as

$$
E = \frac { \int d R | \psi _ { \mathrm { T } } ( R , \alpha ) | ^ { 2 } \frac { H \psi _ { \mathrm { T } } ( R , \alpha ) } { \psi _ { \mathrm { T } } ( R , \alpha ) } } { \int d R | \psi _ { \mathrm { T } } ( R , \alpha ) | ^ { 2 } } .\tag{2}
$$

This facilitates to compute E using the MC method by sampling M points from the probability density function

$$
\rho ( R ) = \frac { | \psi _ { \mathrm { T } } ( R , \alpha ) | ^ { 2 } } { \int d R | \psi _ { \mathrm { T } } ( R , \alpha ) | ^ { 2 } } .\tag{3}
$$

Employing the $\mathrm { { M } ( R T ) ^ { 2 } }$ algorithm (also known as the Metropolis algorithm) [32], the variational energy can be estimated as

$$
E \simeq \frac { 1 } { M } \sum _ { i = 1 } ^ { M } E _ { \mathrm { l o c } } ( R _ { i } ) ,\tag{4}
$$

where

$$
E _ { \mathrm { l o c } } ( R ) \equiv { \frac { H \psi _ { \mathrm { T } } ( R , \alpha ) } { \psi _ { \mathrm { T } } ( R , \alpha ) } }\tag{5}
$$

is the so-called local energy.

Even though appending a simple Jastrow correlation function to the trial WF enables to recover most of the dynamic correlation effects [33], we are considering the shadow wave function (SWF) of Kalos and coworkers [34, 35], as our trial WF. Its main advantage is that it allows to accurately describe localized and delocalized phases within the same functional form [36]. Hence, it is possible to use the same wave function for describing both insulating and metallic electronic structures. In addition, it even admits to compute inhomogeneous systems [37–39]. Finally, the SWF has additional advantageous properties, such as for instance that many-body correlations are taken into account and that it obeys a strong similitude with the exact ground state WF [40].

## A. Shadow Wave Function

The SWF formalism allows to systematically improve an arbitrary trial WF ψT by applying the imaginarytime propagator $e ^ { - \tau H }$ that projects ψT 6⊥ ψGS onto the

ground state WF $\psi _ { \mathrm { G S } }$ . In order to demonstrate this, let us decompose the trial WF into

$$
\psi _ { \mathrm { T } } = \sum _ { n = 0 } ^ { + \infty } c _ { n } \phi _ { n } ,\tag{6}
$$

where $\phi _ { n }$ are the eigenfunctions of the Schr¨odinger equation and $c _ { n }$ the corresponding expansion coefficients. Employing the imaginary-time propagator onto $\psi _ { \mathrm { T } }$ , we obtain

$$
e ^ { - \tau H } \psi _ { \mathrm { T } } = \sum _ { n = 0 } ^ { + \infty } c _ { n } e ^ { - \tau E _ { n } } \phi _ { n } .\tag{7}
$$

The projector $e ^ { - \tau H }$ implies that all excited components are exponentially decaying [41], so that eventually the ground state energy $E _ { 0 }$ is projected out, i.e.

$$
\operatorname * { l i m } _ { \tau  \infty } e ^ { - \tau H } \psi _ { \mathrm { T } } = \operatorname * { l i m } _ { \tau  \infty } \sum _ { n = 0 } ^ { + \infty } c _ { n } e ^ { - \tau E _ { n } } \phi _ { n } \propto \phi _ { 0 } .\tag{8}
$$

From this it follows that $\psi _ { \mathrm { T } } ( R )$ can be systematically improved by

$$
\begin{array} { l } { { \displaystyle e ^ { - \tau H } \psi _ { \mathrm { T } } ( R ) = \langle R | e ^ { - \tau H } | \psi _ { \mathrm { T } } \rangle } } \\ { { \displaystyle ~ = \int d S \langle R | e ^ { - \tau H } | S \rangle \langle S | \psi _ { \mathrm { T } } \rangle , } } \end{array}\tag{9a}
$$

(9b)

where we have introduced an integral over a complete set of Dirac deltas |Si and omitted the inessential normalization factor. Assuming that $\tau \ll 1$ , we now use the Trotter formula to approximate

$$
e ^ { - \tau ( K + V ) } \sim e ^ { - \frac { \tau } { 2 } V } e ^ { - \tau K } e ^ { - \frac { \tau } { 2 } V } ,\tag{10}
$$

where K and V are the operators corresponding to the kinetic and potential energies, respectively [42]. Using the identity

$$
\langle x | e ^ { - \tau K } | y \rangle = \frac { e ^ { - \frac { ( x - y ) ^ { 2 } } { 4 \tau } } } { a } ,\tag{11}
$$

where a is a normalization factor, the eventual expression for the improved trial WF reads as

$$
e ^ { - \tau H } \psi _ { \mathrm { T } } ( R ) = e ^ { - \frac { \tau } { 2 } V ( R ) } \int d S e ^ { - \frac { \tau } { 2 } V ( S ) } e ^ { - \frac { ( R - S ) ^ { 2 } } { 4 \tau } } \langle S | \psi _ { \mathrm { T } } \rangle .\tag{12}
$$

Yet, throughout our derivation we have assumed that $\tau \ll 1$ , which causes that the imaginary-time propagation is rather short and the trial WF only slightly improved. In order to elongate the propagation in imaginary-time and to solve the Schr¨odinger equation exactly, the described procedure needs to be applied repeatedly, which eventually results in a formalism rather similar to the path-integral approach [43, 44]. However, there is no explicit importance sampling in path-integral MC methods [45]. Thus, following our original objective to find an improved and computational efficient trial WF, we rather truncate the projection after one step and refine the obtained functional form variationally. In other words, instead of approaching the limit $\tau  0$ , we substitute τ by a variational parameter C in the gaussian term. Furthermorer, we interpret the exponential $e ^ { - V ( R ) }$ as the Jastrow correlation factor $J _ { \mathrm { p } } ( { \bar { R } } )$ for the protons and likewise $e ^ { - V ( S ) }$ as the corresponding two-body correlation term $J _ { \mathrm { s } } ( S )$ for the shadows. The definition $\langle S | \psi _ { \mathrm { T } } \rangle = \psi _ { \mathrm { T } } ( S )$ entails that the original trial WF has to be evaluated on the shadow coordinates ${ \cal S } \equiv ( \mathbf { s } _ { 1 } , \mathbf { s } _ { 2 } , \ldots , \mathbf { s } _ { N } )$ . The latter is particularly important for the term that determines the symmetry of the SWF, which corresponds to a product of orbitals for a bosonic and a Slater-Determinant (SD) for a fermionic system, respectively [46]. As a consequence, any trial WF ψT can be systematically improved by shadow formalism. The resulting SWF for a bosonic system then reads as

$$
\psi _ { \mathrm { S W F } } ( R ) = J _ { \mathrm { p } } ( R ) \int d S e ^ { - C \sum _ { i = 1 } ^ { N } ( { \bf r } _ { i } - { \bf s } _ { i } ) ^ { 2 } } J _ { \mathrm { s } } ( S ) \psi _ { \mathrm { T } } ( S ) ,\tag{13}
$$

where exp $\Big ( - C \sum _ { i = 1 } ^ { N } ( { \bf r } _ { i } - { \bf s } _ { i } ) ^ { 2 } \Big ) = \boldsymbol \Xi _ { e s }$ is the kernel that connects the electronic coordinates with the associated shadows by means of a gaussian term. From the discussion above, it is apparent that the SWF can also be thought of as an one-step variational path-integral [47].

## B. Shadow Wave Function for Fermionic Systems

Since electrons are spin-1/2 fermions, Fermi-Dirac statistics dictates that the WF must obey the antisymmetry requirement to comply with the Pauli exclusion principle. Thus, a fermionic version of the SWF requires dealing with antisymmetric functions that are changing its sign upon interchanging any two like-spin particles, but whose nodes are inherently unknown.

The most natural way to devise an antisymmetrized SWF is to introduce a SD for each of the spins as a function of $S ,$ i.e. det $\big ( \phi _ { \alpha } ( \mathbf { s } _ { \beta } ^ { \uparrow } ) \big )$ and det $\big ( \phi _ { \alpha } ( \mathbf { s } _ { \beta } ^ { \downarrow } ) \big )$ , where $\phi _ { \alpha }$ are single-particle orbitals that are typically determined by mean-field theories, such as HF or DFT. This results in the so-called Fermionic Shadow Wave Function (FSWF)

$$
\psi _ { \mathrm { F S W F } } ( R ) = J _ { \mathrm { e e } } ( R ) J _ { \mathrm { e p } } ( R , Q ) \int d S e ^ { - C ( R - S ) ^ { 2 } } J _ { \mathrm { s e } } ( S , R )
$$

$$
\times ~ J _ { \mathrm { s p } } ( S , Q ) \operatorname * { d e t } ( \phi _ { \alpha } ( \mathbf { s } _ { \beta } ^ { \uparrow } ) ) \operatorname * { d e t } ( \phi _ { \alpha } ( \mathbf { s } _ { \beta } ^ { \downarrow } ) ) ,\tag{14}
$$

where α and $\beta$ are the row and column indexes of the SDs for the spin-up and spin-down electrons, $J _ { \mathrm { s e } } ( S , R )$ the electron-shadow and $J _ { \mathrm { s p } } ( S , Q )$ the shadow-proton Jastrow correlation factor [40, 48–50], while $\begin{array} { r l } { Q } & { { } \equiv } \end{array}$ $( \pmb { q } _ { 1 } , \pmb { q } _ { 2 } , \dots , \pmb { q } _ { M } )$ are the coordinates of all M protons. However, the FSWF is plagued by a sign problem [49– 51], which differs from the infamous fermion sign problem of projection QMC methods such as Green’s function or diffusion MC [52, 53], but limits its applicability to relatively small systems.

A simple ansatz to circumvent the sign problem is the Antisymmetric Shadow Wave Function (ASWF)

$$
\begin{array} { l } { { \displaystyle \psi _ { \mathrm { A S W F } } ( R ) = J _ { \mathrm { e e } } ( R ) J _ { \mathrm { e p } } ( R , Q ) \operatorname* { d e t } ( \phi _ { \alpha } ( { \bf r } _ { \beta } ^ { \uparrow } ) ) \operatorname* { d e t } ( \phi _ { \alpha } ( { \bf r } _ { \beta } ^ { \downarrow } ) ) } } \\ { { \displaystyle \qquad \times \int d S e ^ { - C ( R - S ) ^ { 2 } } J _ { \mathrm { s e } } ( S , R ) J _ { \mathrm { s p } } ( S , Q ) , \quad ( 1 5 ) } } \end{array}
$$

where $\operatorname* { d e t } ( \phi _ { \alpha } ( \mathbf { r } _ { \beta } ^ { \uparrow } ) )$ and det $( \phi _ { \alpha } ( \mathbf { r } _ { \beta } ^ { \downarrow } ) )$ are SDs as a function of the electronic coordinates [54]. Even though the ASWF already includes many-body correlation effects of any order, the FSWF is superior since it accounts not only for symmetric, but, in addition, also for backflow correlation effects [55, 56].

## C. Trial Wave Functions

We now introduce the trial wave functions that we have employed in our calculations. In particular, the so-called Jastrow-Slater (JS) WF consists of a single SD that is multiplied by a simple Jastrow correlation factor to recover most of the dynamic correlation effects [33, 57, 58]:

$$
\psi _ { \mathrm { J S } } ( R ) \equiv \mathrm { d e t } ( \phi _ { \alpha } ( \mathbf { r } _ { \beta } ^ { \uparrow } ) ) \mathrm { d e t } ( \phi _ { \alpha } ( \mathbf { r } _ { \beta } ^ { \downarrow } ) ) J _ { \mathrm { e e } } ( R ) J _ { \mathrm { e p } } ( R , Q ) ,\tag{16}
$$

where $J _ { \mathrm { e e } }$ and $J _ { \mathrm { e p } }$ are the Jastrow correlation factors $\begin{array} { r } { J = e ^ { - \sum _ { i , j } u ( r _ { i j } ) } } \end{array}$ for the electron-electron and electronproton interactions, whereas $u ( r _ { i j } )$ is a two-body pseudopotential.

For the latter, here we have chosen the Yukawa-Jastrow pseudopotential for $J _ { \mathrm { e e } }$ and $J _ { \mathrm { e p } }$ , respectively, which is defined as

$$
u _ { \mathrm { Y U K } } ( r ) \equiv A \frac { 1 - e ^ { - F r } } { r } ,\tag{17}
$$

where A and F are both variational parameters. The Yukawa-Jastrow pseudopotential is able to satisfy Kato’s cusp condition from the outset, since

$$
u _ { \mathrm { Y U K } } ( r ) \xrightarrow { r  0 } A F - \frac { A F ^ { 2 } } { 2 } r + \mathcal { O } ( r ^ { 2 } ) .\tag{18}
$$

Nevertheless, we have not utilized the cusp condition to fix one of the two parameters, but instead have determined both of them by means of the modified stochastic reconfiguration (SR) algorithm [59], as detailed in section III.

Moving our attention to the orbitals employed in the SD, we have considered four type of orbitals:

1. simple plane wave (pw):

$$
e ^ { i \mathbf { k } _ { i } \mathbf { r } _ { i } } ,
$$

where ki are k-vectors in the Fermi sphere. More details about its actual implementation to include finite size effects are provided in subsection IV C.

2. DFT, computed by the PWscf code of the Quantum Espresso suite of programs [60]. In particular, the Perdew-Burke-Ernzerhof (PBE) generalized gradient approximation to the exact exchangecorrelation functional was used together with the bare Coulomb potential and an associated PW cutoff of just 8 Ry [61]. In order to accurately sample the first Brillouin zone, a dense k-point mesh with at least $5 ^ { 3 }$ special points was utilized [62]. Again, more details are duly appropriated in subsection IV C.

3. 1s, corresponding to the lowest energy solution of the Schr¨odinger equation for an isolated hydrogen atom and is parametrized by the corresponding proton position:

$$
\phi _ { \mathrm { 1 s } } ( \mathbf { r } , \mathbf { q } ) = e ^ { - \gamma | \mathbf { r } - \mathbf { q } | } ,
$$

where γ is a variational parameter.

4. bi-atomic, defined as

$$
\psi _ { \mathrm { b i - a t o m i c } } ( \mathbf { r } , \mathbf { q } _ { 1 } , \mathbf { q } _ { 2 } ) = \phi _ { \mathrm { 1 s } } ( \mathbf { r } , \mathbf { q } _ { 1 } ) + \phi _ { \mathrm { 1 s } } ( \mathbf { r } , \mathbf { q } _ { 2 } ) ,
$$

where $\mathbf { q } _ { 1 }$ and $\mathbf { q } _ { 2 }$ are the positions of the protons of the same $H _ { 2 }$ molecule.

## III. COMPUTATIONAL DETAILS

In the following we are investigating a system comprising of $N = 1 2 8$ hydrogen atoms as specified by the Hamiltonian

$$
\begin{array} { l } { { \displaystyle { \cal H } = - \sum _ { i = 1 } ^ { N } \hbar ^ { 2 } \frac { \nabla _ { i } ^ { 2 } } { 2 m _ { e } } - \sum _ { I = 1 } ^ { N } \hbar ^ { 2 } \frac { \nabla _ { I } ^ { 2 } } { 2 N _ { I } } - \sum _ { i , I = 1 } ^ { N } \frac { K _ { C } } { \left| \mathbf { r } _ { i } - \mathbf { q } _ { I } \right| } } } \\ { { \displaystyle ~ + \sum _ { i < j } \frac { K _ { C } } { \left| \mathbf { r } _ { i } - \mathbf { r } _ { j } \right| } + \sum _ { I < J } \frac { K _ { C } } { \left| \mathbf { q } _ { I } - \mathbf { q } _ { J } \right| } , } } \end{array}\tag{19}
$$

where $K _ { C } = 1 / ( 4 \pi \epsilon _ { 0 } )$ is the Coulomb constant and $\epsilon _ { \mathrm { 0 } }$ the electric free space permittivity.

For the sake of simplicity, we have confined ourselves to the hcp and bcc phases as representatives for the insulating molecular and metallic atomic phases of solid hydrogen, respectively. To simulate an extended solid, 3-dimensional periodic boundary conditions (pbc) were deployed throughout, whereas the volume of the corresponding unit cell was determined by the Wigner-Seitz radius $r _ { s } = \sqrt [ 3 ] { 3 / ( 4 \pi \rho ) }$ , with $\rho$ being the particle density.

The electronic Schr¨odinger equation is approximately solved by VMC in conjunction with the various trial WFs described above using the HswfQMC code [63]. Since it is well known that conducting a QMC calculation by displacing all particles concurrently from a flat distribution entails a rather strong autocorrelation, here we have elected to use single-particle moves instead. This is to say that $\pmb { r } _ { l } ^ { \mathrm { n e w } } = \pmb { r } _ { l } ^ { \mathrm { o l d } } + \bar { \Delta } ( \eta _ { 1 } , \eta _ { 2 } , \eta _ { 3 } )$ , where l is the index of the moved electron and $\Delta$ the corresponding magnitude of the displacement, while $\eta _ { i }$ are random numbers from the interval $[ - 1 / 2 , + 1 / 2 ]$ . Whereas efficiently recomputing the Jastrow correlation factor after a single particle move efficiently is relatively straightforward, this is not the case for the update of the SD. Following Ceperley et al. [64],

$$
\mathrm { S D } ^ { \mathrm { n e w } } = \mathrm { S D } ^ { \mathrm { o l d } } \sum _ { j } \left( A ^ { - 1 } \right) _ { j l } ^ { \mathrm { o l d } } A _ { l j } ^ { \mathrm { n e w } } ,\tag{20}
$$

where A is the matrix that generates the SD, i.e $\operatorname* { d e t } ( A ) =$ SD. Similarly, also the inverse matrix $\left( A ^ { - 1 } \right)$ can be conveniently updated by means of

$$
\left\{ \begin{array} { l l l } { \left( A ^ { - 1 } \right) _ { i l } ^ { \mathrm { n e w } } = } & { \left( A ^ { - 1 } \right) _ { i l } ^ { \mathrm { o l d } } \frac { \mathrm { S D } ^ { \mathrm { o l d } } } { \mathrm { S D } ^ { \mathrm { n e w } } } } \\ { \left( A ^ { - 1 } \right) _ { i j } ^ { \mathrm { n e w } } = } & { \left( A ^ { - 1 } \right) _ { i j } ^ { \mathrm { o l d } } - \left( A ^ { - 1 } \right) _ { i l } ^ { \mathrm { o l d } } \frac { \mathrm { S D } ^ { \mathrm { o l d } } } { \mathrm { S D } ^ { \mathrm { n e w } } } } \\ { \times } & { \sum _ { s } \left( A ^ { - 1 } \right) _ { s j } ^ { \mathrm { o l d } } A _ { l s } ^ { \mathrm { n e w } } , } \end{array} \right.\tag{21}
$$

with $j \neq l$ . At the beginning of each VMC simulation, we set $\Delta$ so as to realize an acceptance rate of $\sim 5 0 \%$ Moreover, in order to reduce the autocorrelations, $3 N / 2$ single-particle moves were attempted between every successive evaluation of the estimators.

Even though the high-dimensional integral of Eq. (1) can be efficiently computed using the $\mathrm { { M } ( R \bar { T } ) ^ { 2 } }$ algorithm, it is nevertheless essential to determine the optimal variational parameters α that minimizes the variational energy. For that purpose we utilize the recently proposed modified SR algorithm [59], originally proposed by Sorella [65]. Specifically, the SR method prescribes that the variational parameters are varied according to

$$
\delta \alpha _ { l } = \lambda \sum _ { k = 1 } ^ { n } f _ { k } \left( s ^ { - 1 } \right) _ { k l } ,\tag{22}
$$

where

$$
\left\{ \begin{array} { l } { s _ { l k } = \langle O _ { k } O _ { l } \rangle - \langle O _ { l } \rangle \langle O _ { k } \rangle } \\ { f _ { k } = \langle H \rangle \langle O _ { k } \rangle - \langle O _ { k } H \rangle } \\ { O _ { k } = \frac { \partial } { \partial \alpha _ { k } } \ln ( \psi _ { \mathrm { T } } ) } \end{array} \right\}\tag{23}
$$

and $\langle \cdot \rangle \equiv \langle \psi _ { \mathrm { T } } | \cdot | \psi _ { \mathrm { T } } \rangle$ Once the gradient in the variational parameters space, which minimizes the variational energy has been computed, the step length λ along this direction needs to be identified. Since determining the new direction $\delta \alpha$ is computational approximately equally expensive than calculating the the variational energy, it is convenient to start with a rather small value for λ and continuously adjusting it on the fly during the optimization.

## IV. VARIATIONAL MONTE CARLO FOR EXTENDED SYSTEMS

When dealing with extended systems, special care is required to accurately consider pbc and single-electron finite size effects.

## A. Periodic Coordinates

If computed in its straightforward fashion, the Yukawa-Jastrow, as any other slowly decaying Jastrow correlation factor, leads to a spurious bias in the kinetic energy. Therefore, all contributions that are originating from the periodic images of the unit cell must be taken explicitly into account in order to avoid discontinuities in the derivatives of the WF when the particle distances switch from one closest image to the other. Needless to say that this approach is computationally relatively time consuming and a more economic strategy very desirable.

However, before presenting our solution to this effect, let us start by introducing a particular useful test to verify if all correlations are correctly taken into account. To that extent, the expression for the kinetic energy (for simplicity we consider the kinetic contribution of only one particle j) is integrated by parts

$$
\begin{array} { l } { { \displaystyle \int _ { \Omega } d R \psi ^ { * } ( R ) \nabla _ { j } ^ { 2 } \psi ( R ) = \sum _ { x _ { \alpha } } \int _ { \Omega } d R \psi ^ { * } ( R ) \left( \frac { \partial ^ { 2 } } { \partial x _ { \alpha } ^ { 2 } } \right) \psi ( R ) } }  \\ { { \displaystyle = \sum _ { x _ { \alpha } } \int _ { \Omega } d R ^ { \tilde { x } _ { \alpha } } \int _ { - L _ { x _ { \alpha } } / 2 } ^ { + L _ { x _ { \alpha } } / 2 } d x _ { \alpha } \psi ^ { * } ( R ) \frac { \partial ^ { 2 } } { \partial x _ { \alpha } ^ { 2 } } \psi ( R ) } }  \\ { { \displaystyle = \sum _ { x _ { \alpha } } \left\{ \int _ { \Omega } d R ^ { \tilde { x } _ { \alpha } } \left[ \psi ^ { * } ( R ) \frac { \partial } { \partial x _ { \alpha } } \psi ( R ) \right] _ { - L _ { x _ { \alpha } } / 2 } ^ { + L _ { x _ { \alpha } } / 2 } \right. } } \\ { { \displaystyle \left. - \int _ { \Omega } d R \left( \frac { \partial } { \partial x _ { \alpha } } \psi ^ { * } ( R ) \right) \left( \frac { \partial } { \partial x _ { \alpha } } \psi ( R ) \right) \right\} } , \qquad { \displaystyle ( 2 4 } }  \end{array}
$$

where $\boldsymbol { x } _ { \alpha } = ( x _ { j } , y _ { j } , z _ { j } )$ , dRx¯α is the same as dR but excluding the infinitesimal element $d x _ { \alpha } , \Omega$ represents the domain of integration, i.e. the simulation cell, while $L _ { x _ { \alpha } }$ is the length of the edge of Ω along the $x _ { \alpha }$ axis. However, at the presence of periodic boundary conditions, the WF $\psi ( R )$ and also its derivatives are required to be periodic, meaning that they are invariant with respect to particle translations $\mathbf { v } = ( n _ { x } L _ { x } , n _ { y } L _ { y } , n _ { z } L _ { z } )$ , where $n _ { x } , \ n _ { y }$ and $n _ { z }$ are all integers. From this follows that the term

$$
\left[ \psi ^ { * } ( R ) \frac { \partial } { \partial x _ { \alpha } } \psi ( R ) \right] _ { - L _ { x _ { \alpha } } / 2 } ^ { + L _ { x _ { \alpha } } / 2 }\tag{25}
$$

vanishes, which leads to a modified Jackson-Feenberg (JF) kinetic energy expression [66]

$$
E _ { \mathrm { J F } } = \hbar ^ { 2 } \sum _ { j = 1 } ^ { N } \frac { 1 } { 2 m _ { j } } \int _ { \Omega } d R \nabla _ { j } \psi ^ { * } ( R ) \cdot \nabla _ { j } \psi ( R ) .\tag{26}
$$

As a consequence, the equivalence of the Pandharipande-Bethe (PB)

$$
E _ { \mathrm { P B } } = - \hbar ^ { 2 } \sum _ { i = 1 } ^ { N } \frac { 1 } { 2 m _ { j } } \int _ { \Omega } d R \psi ^ { * } ( R ) \nabla _ { j } ^ { 2 } \psi ( R )\tag{27}
$$

and Jackson-Feenberg (JF) expressions for the kinetic energy is a necessary, but not sufficient condition for the required periodic properties of the WF. Thus, in all of our calculations we have computed both expressions and explicitly verified that both are indeed identical, within the corresponding statistical uncertainties.

<!-- image-->  
FIG. 1. Integration of $\mathbf { r } _ { j }$ within a box with periodic boundary conditions, from the point of view of particle i. The continuous black line represents the simulation box, whereas the dotted black line denotes the effective volume of integration for the distance between the particles i and j.

However, at the presence of additional correlation terms, such as the Jastrow, $\operatorname { E q } .$ 25 must be correctly interpreted since the inter-particle distances are computed using the closest periodic image. In fact, even though the particle coordinate $\mathbf { r } _ { j }$ is confined to the unit cell of volume $\mathcal { V } = L _ { x } L _ { y } L _ { z }$ , the distance $\mathbf { r } _ { i j }$ between the particles i (assumed as fixed) and $j$ does not range within $\sqrt { \left( \mathbf { r } _ { j } - \mathbf { r } _ { i } \right) ^ { 2 } }$ , but always within the box of volume V centered on particle i. This concept is illustrated in Fig. 1. Therefore, it is possible and convenient to fix the origin at the position of the i-th particle that is considered. As a consequence, in the following we will set $\mathbf { r } _ { i } = 0$ , so that $\mathbf { r } = ( x , y , z ) \equiv \mathbf { r } _ { i j } = \mathbf { r } _ { j }$ and $r = | \mathbf { r } | = { \sqrt { x ^ { 2 } + y ^ { 2 } + z ^ { 2 } } }$

Let us demonstrate the JF test by showing that the Yukawa-Jastrow violates it. For that purpose we consider the simple case of only two interacting particles, i.e.

$$
\begin{array} { r } { J ( r ) = e ^ { - \frac { A \left( 1 - e ^ { - F r } \right) } { r } } . } \end{array}\tag{28}
$$

Its first derivative along the x axis reads as

$$
\begin{array} { r l } & { \frac { \partial J ( r ) } { \partial x } = \frac { \partial e ^ { - \frac { A \left( 1 - e ^ { - F r } \right) } { r } } } { \partial x } } \\ & { \qquad = \frac { \partial e ^ { - \frac { A \left( 1 - e ^ { - F r } \right) } { r } } } { \partial r } \frac { \partial r } { \partial x } } \\ & { \qquad = e ^ { - \frac { A \left( 1 - e ^ { - F r } \right) } { r } } \left( A \frac { 1 - e ^ { - F r } } { r ^ { 2 } } - A F \frac { e ^ { - F r } } { r } \right) \frac { x } { r } } \\ & { \qquad = \frac { \partial J ( r ) } { \partial r } \frac { x } { r } . } \end{array}\tag{9}
$$

It is then apparent that the first derivative is not continuous at $x = \pm L / 2$ , which is the border between its two closest periodic images. As a consequence,

$$
\operatorname* { l i m } _ { \varepsilon \to 0 } \left[ J ( r ) { \frac { \partial } { \partial x } } J ( r ) \right] _ { - ( L _ { x } / 2 ) + \varepsilon } ^ { + ( L _ { x } / 2 ) - \varepsilon } = J ( r ) { \frac { \partial J ( r ) } { \partial r } } { \frac { L _ { x } } { r } } \neq 0 .\tag{30}
$$

Therefore, the JF and PB kinetic energies differ since the term in Eq. 25 does not vanish. Yet, if $\frac { 1 } { r } \frac { \partial J } { \partial r }$ is small enough at $\begin{array} { r } { x = \pm \frac { L } { 2 } } \end{array}$ , the difference is negligible.

An even deeper understanding can be obtained by means of the distribution theory. In fact,

$$
\begin{array} { c } { \displaystyle \frac { \partial ^ { 2 } } { \partial x ^ { 2 } } J ( r ) = \displaystyle \frac { x } { r } \left( \frac { \partial ^ { 2 } J ( r ) } { \partial r ^ { 2 } } \frac { x } { r } + \frac { \partial J ( r ) } { \partial r } \frac { 1 } { r } - \frac { \partial J ( r ) } { \partial r } \frac { 1 } { r ^ { 2 } } \right) } \\ { \displaystyle - \left( \frac { \partial J ( r ) } { \partial x } \frac { x } { r } \right) 2 \delta \left( x - \frac { L _ { x } } { 2 } \right) , } \end{array}\tag{31}
$$

since

$$
\int _ { L _ { x } / 2 - \varepsilon } ^ { L _ { x } / 2 + \varepsilon } \frac { \partial } { \partial x } \left( \frac { \partial J ( r ) } { \partial x } \right) = - \frac { \partial J ( r ) } { \partial x } \frac { L _ { x } } { r } .\tag{32}
$$

In other words, the discontinuity in the first derivative entails a Dirac delta in the second derivative. Obviously, this artifact must be circumvented in order to avoid a bias in the computation of the kinetic energy.

Since a discontinuity in the first derivative not only affects the validity of the JF expression but also the PB one, the kinetic energy contribution provided by the Yukawa-Jastrow is biased. Mathematically, the problem can be eliminated by enforcing a smooth change between the closest periodic images. Physically, all of this originates from the fact that the simulation box is not large enough to “contain” all correlations between the particles.

A straightforward solution to remedy the latter is inspired by the Ewald summation technique [67]. More specifically, the Jastrow is decomposed into a quickly and a slowly decaying part, which are computed separately in real and reciprocal k-space, respectively. However, this method requires a summation over the whole momentum space, which is computationally rather demanding.

An alternative approach, which is not only more elegant and simpler, but at the same time also more efficient, is due to Attaccalite and Sorella and results from exploiting Periodic Coordinates (PC) [68]. As the name suggests, the only modification required is to substitute the original coordinates by

<!-- image-->  
FIG. 2. Comparison between the Yukawa Jastrow correlation factor for two interacting particles exp $\left( - { \frac { A ( 1 - \exp ( - F r ) ) } { r } } \right)$ (black line) and its modified version as obtained by employing Periodic Coordinates (dashed line). We have set $y = z = 0 ,$ $L = 1 0$ and $A = F = 1$ , respectively.

$$
x ^ { \prime } = { \frac { L } { \pi } } \sin \left( { \frac { \pi x } { L } } \right) ,\tag{33a}
$$

$$
y ^ { \prime } = { \frac { L } { \pi } } \sin \left( { \frac { \pi y } { L } } \right) ,\tag{33b}
$$

$$
z ^ { \prime } = { \frac { L } { \pi } } \sin \left( { \frac { \pi z } { L } } \right) ,\tag{33c}
$$

and hence evaluate the distances via

$$
r ^ { \prime } = \frac { L } { \pi } \sqrt { \sin ^ { 2 } \left( \frac { \pi x } { L } \right) + \sin ^ { 2 } \left( \frac { \pi y } { L } \right) + \sin ^ { 2 } \left( \frac { \pi z } { L } \right) } .\tag{34}
$$

The employment of Periodic Coordinates enforces the correct periodicity of the WF. For example, the first derivative

$$
\begin{array} { c } { \displaystyle \frac { \partial J ( \boldsymbol { r } ^ { \prime } ) } { \partial x } = \frac { \partial J ( \boldsymbol { r } ^ { \prime } ) } { \partial \boldsymbol { r } ^ { \prime } } \frac { \partial \boldsymbol { r } ^ { \prime } } { \partial x ^ { \prime } } \frac { \partial x ^ { \prime } } { \partial x } } \\ { \displaystyle = \frac { \partial J ( \boldsymbol { r } ^ { \prime } ) } { \partial \boldsymbol { r } ^ { \prime } } \frac { x ^ { \prime } } { \boldsymbol { r } ^ { \prime } } \cos \left( \frac { \pi x } { L } \right) , } \end{array}\tag{35}
$$

is continuous in $x = \pm L / 2$ , i.e. on the borders of the simulation box. The same also holds for all higher order derivatives. The consequential modifications of the Yukawa Jastrow are illustrated in Fig. 2.

To demonstrate the effectiveness of PC, we have calculated the kinetic energy using the JS-pw trial WF for two different systems, each consisting of 16 hydrogen atoms. The results of the atomic bcc (atm-bcc) and the molecular hcp (mol-hcp) phases of solid hydrogen including the corresponding Wigner-Seitz radii are shown in Table I. As can be extracted by comparing $E _ { \mathrm { k i n } }$ with $E _ { \mathrm { J F } }$ , the aforementioned spurious bias can be completely eliminated by the use of PC with a only negligible additional computational cost. Nevertheless, we find it important to remark that employing PC leads to a somewhat modified Yukawa Jastrow, which may slightly violate the electron-electron and electron-proton cusp conditions [69]. However, the accuracy of employed Jastrow in conjunction with PC can be easily checked by means of the variational principle. We have explicitly verified that in practice the latter bias is generally tiny.

TABLE I. Kinetic energies (in Ry) for the atm-bcc and molhcp phases of solid hydrogen as obtained with and without PC.
<table><tr><td rowspan=1 colspan=1></td><td rowspan=1 colspan=1>without PC</td><td rowspan=1 colspan=1>with PC</td></tr><tr><td rowspan=1 colspan=1>atm-bcc $r _ { s } = 1 . 3 1$ </td><td rowspan=1 colspan=1> $\overline { { E _ { \mathrm { k i n } } = 5 . 5 7 6 6 ( 5 ) } }$  $E _ { \mathrm { J F } } = 2 . 9 8 4 1 ( 3 2 )$ </td><td rowspan=1 colspan=1> $\overline { { E _ { \mathrm { k i n } } = 1 . 5 4 8 0 ( 6 ) } }$  $E _ { \mathrm { J F } } = 1 . 5 4 5 4 ( 1 6 )$ </td></tr><tr><td rowspan=1 colspan=1>mol-hcp $r _ { s } = 2 . 6 1$ </td><td rowspan=1 colspan=1> $\overline { { E _ { \mathrm { k i n } } = 2 . 2 4 2 8 ( 6 ) } }$  $E _ { \mathrm { J F } } = 2 . 1 2 5 2 ( 1 0 )$ </td><td rowspan=1 colspan=1> $\overline { { E _ { \mathrm { k i n } } = 1 . 0 3 0 7 ( 1 1 ) } }$  $E _ { \mathrm { J F } } = 1 . 0 2 9 0 ( 5 )$ </td></tr></table>

## B. SWF Kernel Truncation

If the variational parameter C of the SWF kernel is small, the simulation box is typically not large enough to constrain each particle to its associated shadows s1 and s2 within the limit $L / 2$ . As before, this entails a bias in the kinetic energy, as can be seen by the difference between $E _ { \mathrm { k i n } } ~ = ~ 1 . 6 2 4 ( 4 )$ Ry and $E _ { \mathrm { J F } } ~ = ~ 1 . 4 7 8 ( 5 )$ Ry, respectively [70].

In order to eliminate this shortcoming, the kernel must be modified so that it vanishes for $| \mathbf { r } - \mathbf { s } |  L / 2$ . An appropriate choice for the modified kernel reads as

$$
\Xi _ { e s } ( R , S ) = \prod _ { i = 1 } ^ { N } e ^ { - \varsigma ( | \mathbf { r } _ { i } - \mathbf { s } _ { i } | ) } ,\tag{36}
$$

where

$$
\varsigma ( x ) = \left\{ \begin{array} { l l } { C x ^ { 2 } } & { \mathrm { ~ i f ~ } x \leq \frac { L } { 6 } } \\ { \alpha _ { 0 } + \frac { \alpha _ { 1 } } { x - L / 2 } } & { \mathrm { ~ i f ~ } \frac { L } { 6 } < x \leq \frac { L } { 2 } - \varepsilon } \\ { \beta _ { 0 } + \beta _ { 1 } x ^ { n } } & { \mathrm { ~ i f ~ } x > \frac { L } { 2 } - \varepsilon } \end{array} \right. ,\tag{37}
$$

with

$$
\begin{array} { r } { \left\{ \begin{array} { l l } { \varepsilon = \frac { L } { n + 1 } } \\ { \beta 1 = - \frac { 2 \alpha _ { 1 } } { \varepsilon ^ { 3 } n \left( n - 1 \right) \left( \frac { L } { 2 } - \varepsilon \right) ^ { n - 2 } } } \\ { \beta _ { 0 } = \alpha _ { 0 } - \frac { \alpha _ { 1 } } { \varepsilon } - \beta _ { 1 } \left( l - \varepsilon \right) ^ { n } } \end{array} \right. } \end{array}\tag{38}
$$

and $n \ \geq \ 2$ Our simulations have suggested that a suitable choice is $\ n \ = \ 1 2$ The modification introduced by Eq. (37) are illustrated in Fig. 3. The corresponding kinetic energies are $E _ { \mathrm { k i n } } = 1 . 5 3 8 ( 6 )$ Ry and $E _ { \mathrm { J F } } = 1 . 5 3 3 ( 6 )$ Ry, which demonstrates that the proposed SWF kernel truncation method completely alleviates the aforementioned limitation.

<!-- image-->  
FIG. 3. Illustration of the SWF kernel truncation method prescribed by Eq. (37) with $C = 0 . 5 4 2$ , for 16 hydrogen atoms in the metallic atm-bcc phase at $r _ { s } = 1 . 3 1$ $( L / 2 \simeq 2 . 6 6 ~ \mathrm { a _ { 0 } } )$ .

## C. Twist Averaged Boundary Conditions

As already alluded to previously, the application of pbc do not automatically result in an accurate description of an infinite system. In fact, identical simulations but using distinct values for N may entail rather different results. As a consequence, these effects are generally referred to as finite-size effects, which can be minimized by the usage of so-called Twist Averaged Boundary Conditions (TABC) [71]. The origin of these finite-size effects are that the embedded k-vectors do not well represent an infinite system, since in general a discrete grid of points cannot reproduce the whole Fermi sphere (see Fig. 4). The TABC method, which allows to bypass this limitation by means of an integration over the Fermi sphere, prescribes a recurrent random shift

$$
\mathbf { v } _ { \mathrm { t w i s t } } = \frac { 2 \pi } { L } \left( \eta _ { 1 } , \eta _ { 2 } , \eta _ { 3 } \right)\tag{39}
$$

of the k-grid, where ηi are random numbers sampled in the range $[ - \frac { 1 } { 2 } , \frac { 1 } { 2 } ]$ Within the context of TABC, the translation vector $\mathbf { v } _ { \mathrm { t w i s t } }$ is referred to as twist. The corresponding integral over k-space is computed by means of MC. As can be seen in Fig. 5, the application of TABC results in an accelerated convergence to the thermodynamic limit. Beyond solely reducing finite-size effects, employing TABC also permits calculations, where the number of particles M are distinct from magic numbers[72], without spurious drift and anisotropy effects.

The eventual algorithm for a VMC simulation of a 3D unpolarized system employing TABC with $N _ { \mathrm { t w i s t } }$ twists reads as follows:

1. Determine the smallest magic number n that is larger than $N / 2 ;$

<!-- image-->  
FIG. 4. Two-dimensional cross section of the momenta k of $\mathrm { ~ a ~ } N \times N \mathrm { ~ S D _ { p w } ~ }$ matrix, for several values of N . The dotted circles delineate the Fermi sphere.

<!-- image-->  
FIG. 5. Electronic kinetic energy of solid hydrogen at $r _ { s } = $ 1.31 as computed using $\mathrm { S D } _ { \mathrm { p w } }$ , where the number of particles M are magic numbers VI.

2. Find the first n Fermi k-vectors, yielding $\Gamma _ { 0 } ~ =$ $\left\{ \mathbf { k } _ { 1 } , \mathbf { k } _ { 2 } , \ldots , \mathbf { k } _ { n } \right\}$ ;

3. Generate $\mathbf { v } _ { \mathrm { t w i s t } }$ as described in Eq. (39);

4. $K = \Gamma _ { 0 } + \left\{ \mathbf { v } _ { \mathrm { t w i s t } } \right\} _ { n } ;$

5. Sort the k-vectors in K by increasing magnitude, and then use the first $N / 2$ k-vectors to build up $\mathrm { S D _ { p w } }$ ;

6. Perform $M _ { \mathrm { r e l a x } }$ relaxation steps;

<!-- image-->  
FIG. 6. Kinetic energy of solid hydrogen with 54 and 66 particles for the twists $\begin{array} { r } { \mathbf { v } _ { \mathrm { t w i s t } } = \frac { 2 \pi } { L } \eta ( 1 , 0 , 0 ) } \end{array}$ 7 $\begin{array} { r } { \mathbf { v } _ { \mathrm { t w i s t } } = \frac { 2 \pi } { L } \eta ( 1 , 1 , 0 ) } \end{array}$ and $\begin{array} { r } { \mathbf { v } _ { \mathrm { t w i s t } } = \frac { 2 \pi } { L } \eta ( 1 , 1 , 1 ) } \end{array}$ , respectively.

7. Sample $M / N _ { \mathrm { t w i s t } }$ points and accumulate the estimators of the observables of interest (normally the kinetic and potential energies);

8. Repeat the points $3 – 7 \ N _ { \mathrm { t w i s t } }$ times.

The $M _ { \mathrm { r e l a x } }$ relaxation steps of point 6 are essential to prevent the emergence of a bias in the calculation. Even though it is possible to circumvent this step by submitting the twist to the acceptance/refuse process of the $\mathrm { { M } ( R T ) ^ { 2 } }$ algorithm, we have not exploited this possibility, since the number of relaxation steps is small and its computational cost negligible.

However, as a consequence of the twist, a momentum in the external shell, which initially was not included in $\mathrm { S D } _ { \mathrm { p w } } ,$ may indeed have a lower magnitude than the employed ones. This is to say that such a momentum actually replaces the one with the actual highest magnitude. Therefore, in step 1 of the just outlined algorithm, more k-vectors than strictly necessary to generate $\mathrm { S D _ { p w } }$ are considered and eventually selected as described in point 5. The corresponding kinetic energies generated by this method are reported in Fig. 6.

In the following we present our extension of the TABC approach to $\mathrm { S D } _ { \mathrm { D F T } }$ . In fact, the DFT method itself also suffers from finite-size effects, which requires to sum over contributions from different K-points in the first Brillouin zone. The simplest grid consists of just one point, denoted as $\Gamma _ { 0 } .$ which corresponds to the Fermi gas momenta. In order to reduce finite-size effects within DFT, it is essential to consider multiple K-points to yield a more accurate averaged estimate of the aforementioned integral, similarly to TABC technique. Typically, the Kpoint grids are generated using the Monkhorst and Pack construction scheme [62]. Due to the fact that each $K \mathfrak { - }$ point has an associated weight, instead of summing over all weighted configurations, we propose here to adopt the

TABC approach with a probability proportional to their weight. In other words, we average over all K-points, while making the most of importance sampling.

The implementation of the modified TABC method for $\mathrm { S D } _ { \mathrm { D F T } }$ can be summarized by the following instructions:

1. Conduct a DFT plane-wave calculation with an energy cutoff $E _ { c t f }$ in order to obtain nK solutions, one for each K-points $K _ { i }$ and its associated weight wi;

2. Sample each K-point $K _ { j }$ with probability

$$
P _ { j } = \frac { w _ { j } } { \sum _ { l = 1 } ^ { n _ { \mathrm { K } } } w _ { l } }
$$

and employ its associated solutions in the $\operatorname { S D } _ { \mathrm { { D F T } } } ;$

3. Perform $M _ { \mathrm { r e l a x } }$ relaxation steps;

4. Sample $M / N _ { \mathrm { t w i s t } }$ points and accumulate the estimators;

5. Repeat the points $2 \mathrm { - 4 } ~ N _ { \mathrm { t w i s t } }$ times.

The results, as obtained employing the modified TABC method in conjunction with a JS-DFT trial ${ \mathrm { W F } } ,$ are reported in Fig. 7. As can be seen, the convergence with respect to $E _ { c t f }$ is much slower for the metallic atm-bcc than for the insulating mol-hcp phase of solid hydrogen, where as few as 10 Ry is adequate. Moreover, in all cases $n _ { K } = 5$ is sufficient to consider all finite size effects for N larger than 16. Nevertheless, since the accumulated statistics for each K-point contribute to the overall average, the total computational cost is essentially independent from $n _ { K }$

The effectiveness of the modified TABC approach as a function of N is demonstrated in Fig. 8. As can be seen, the TABC provide a quicker convergence to the thermodynamic limit especially in the case of the metallic atm-bcc phase that obeys rather large finite size effects.

## V. RESULTS AND DISCUSSION

To demonstrate the predictive power of the SWF in general and the ASWF-DFT trial WF in particular, we investigate the metal-insulator-transition (MIT) from the metallic atm-bcc to the insulating mol-hcp phase of solid hydrogen. The corresponding results using the conventional JS trial WF are shown in Fig. 9. Not surprisingly, using the JS-DFT trial WF, the variational energies are throughout more favorable than the ones obtained by the JS-pw trial WF. However, while the latter are in reasonable good agreement with the former for the metallic atm-bcc phase, the $\mathrm { J S - p w }$ trial WF fails to describe the insulating mol-hcp phase. In general, the results of the JS-pw and JS-DFT trial WFs are deviating from each other with increasing distance between the monomers that implies with larger multireference character. Interestingly, we find that especially for large monomer separation the rather simple JS-bi-atomic and JS-1s trial

<!-- image-->

<!-- image-->  
FIG. 7. Variational energy of solid hydrogen in the atm-bcc phase at $r _ { s } ~ = ~ 1 . 3 1$ and mol-hcp phase at $r _ { s } ~ = ~ 2 . 6 1$ as a function of $E _ { \mathrm { c t f } }$ and $n _ { K }$ . The energies were calculated for $N = 1 6$ using the JS-DFT trial WF.

WFs are in fact even more accurate than the JS-DFT results. Considering its simplicity, the JS-1s trial WF performs relatively well for both of the considered phases. However, as can be seen in Fig. 10, the increased accuracy of the ASWF with respect to the JS-type WFs is rather limited. Although, the improvement is noticeable in the case of the JS-pw trial WF, for the more accurate JS-DFT approach it renders inessential. This is to say that the observed improvement in the employed WF is nearly entirely due to the application of DFT to construct the SD, which subsequently is not further enhanced by the present shadow formalism. The latter suggests that the eventual DFT-based trial WFs are already very accurate.

In order to the determine the transition pressure of the MIT for the various trial WF investigated here, in Fig. 11 the energies the metallic atm-bcc and the insulating molhcp phases are shown as a function of $r _ { s } .$ . Using the common tangent construction, we find an MIT pressure of 12 GPa for the JS-pw trial WF, which is even lower than predicted by Wigner and Huntington back in 1935 [1, 73]. Applying the more accurate ASWF formalism instead of the plain JS trial WF, the MIT pressure slightly increases to and 45 GPa. However, as before, substituting the pw orbitals within the SD by those of a meanfield DFT calculation, results not only in a substantially reduced variational energy, but also in a dramatically increased MIT pressure. Specifically, employing the JS-DFT trial WF results in a transition pressure of 395 GPa, while the usage of the present ASWF transformation increases the MIT pressure to even 520 GPa, which is still beyond the largest pressures experimentally realized so far at low temperature. Therefore, although the variational energy is only slightly improved by the ASWF when using DFT orbitals in the SD, the impact on the transition pressure is rather large. Moreover, the present results immediately suggest the general trend that more accurate the employed trial WF, the higher the resulting MIT pressure. In fact, despite the simplicity of the underlying JS-type trial WF, the present ASWF-DFT results compares relatively favorable with recent state-ofthe-art finite-temperature QMC calculations using much more sophisticated trial WF [74–78]. Nevertheless, it is important to note that the here considered solid phases of insulating molecular and metallic atomic hydrogen are not the energetically most favorable structures known to date and as such only qualitative representatives of the MIT [7]. Furthermore, the possible existence of a quantum fluid phase at zero temperature, which is consistent with a maximum in the melting curve [4, 8, 9, 11, 79], is

<!-- image-->

<!-- image-->  
FIG. 8. Variational energy for solid hydrogen in the atmbcc phase at $r _ { s } = 1 . 3 1$ and mol-hcp phase at $r _ { s } = 2 . 6 1$ as a function of N. The energies were calculated for $N = 1 6$ using the JS-DFT trial WF.

<!-- image-->

<!-- image-->  
FIG. 9. Variational energies of the metallic atm-bcc and the insulating mol-hcp phases of solid hydrogen using various JS-type trial WFs.

neglected.

## VI. CONCLUSION

In conclusion, we have extended the ASWF to periodic large-scale systems made up fermions. For that purpose, we have exploited an improved SR scheme to efficiently optimize the employed ASWF [59], and combined it with enhanced PC and TABC techniques. To demonstrate the predictive power of this approach, we investigated the MIT of solid hydrogen at very high pressure. In particular we found that the ameliorated accuracy of the ASWF results in a significantly increased transition pressure of 520 GPa.

## ACKNOWLEDGMENTS

The authors would like to thank the Graduate School of Excellence MAINZ for financial support and Markus Holzmann for useful comments. The Gauss Center for Supercomputing (GCS) is kindly acknowledged for providing computing time through the John von Neumann Institute for Computing (NIC) on the GCS share of the supercomputer JUQUEEN at the J¨ulich Supercomputing Centre (JSC).

<!-- image-->

<!-- image-->

FIG. 10. Variational energies of the metallic atm-bcc and the insulating mol-hcp phases of solid hydrogen using the ASWF and JS-type trial WFs.  
<!-- image-->

<!-- image-->

<!-- image-->

<!-- image-->  
FIG. 11. The MIT between the metallic atm-bcc and the insulating mol-hcp phases of solid hydrogen using the ASWF and JS-type trial WFs.

[1] E. Wigner and H. B. Huntington, J. Chem. Phys. 3, 764 (1935).

[2] A. Alavi, M. Parrinello, and D. Frenkel, Science 269, 1252 (1995).

[3] N. W. Ashcroft, Phys. Rev. Lett. 21, 1748 (1968).

[4] S. A. Bonev, E. Schwegler, T. Ogitsu, and G. Galli, Nature (London) 431, 669 (2004).

[5] I. F. Silvera, Rev. Mod. Phys. 52, 393 (1980).

[6] H.-k. Mao and R. J. Hemley, Rev. Mod. Phys. 66, 671 (1994).

[7] J. M. McMahon, M. A. Morales, C. Pierleoni, and D. M. Ceperley, Rev. Mod. Phys. 84, 1607 (2012).

[8] S. Scandolo, Proc. Nat. Acad. Sci. USA 100, 3051 (2003).

[9] S. Deemyad and I. F. Silvera, Phys. Rev. Lett. 100, 155701 (2008).

[10] I. F. Silvera and S. Deemyad, Low Temperature Physics 35, 318 (2009), ISSN 10906517.

[11] M. I. Eremets and I. A. Trojan, JETP Lett. 89, 174 (2009).

[12] M. I. Eremets and I. A. Troyan, Nature Mater. 10, 927 (2011).

[13] D. E. Ramaker, L. Kumar, and F. E. Harris, Phys. Rev. Lett. 34, 812 (1975).

[14] T. W. Barbee, M. L. Cohen, and J. L. Martins, Phys. Rev. Lett. 62, 1150 (1989).

[15] T. W. Barbee, III, A. Carcia, and M. L. Cohen, Nature (London) 340, 369 (1989).

[16] C. F. Richardson and N. W. Ashcroft, Phys. Rev. Lett. 78, 118 (1997).

[17] C. J. Pickard and R. J. Needs, Nature Phys. 3, 473 (2007).

[18] S. Azadi and T. D. K¨uhne, JETP Lett. 95, 449 (2012).

[19] S. Azadi, W. M. C. Foulkes, and T. D. K¨uhne, New Journal of Physics 15, 113005 (2013).

[20] R. Singh, S. Azadi, and T. D. K¨uhne, Phys. Rev. B 90, 014110 (2014).

[21] R. O. Jones and O. Gunnarsson, Rev. Mod. Phys. 61, 689 (1989).

[22] W. Kohn, Rev. Mod. Phys. 71, 1253 (1999).

[23] W. M. C. Foulkes, L. Mitas, R. J. Needs, and G. Rajagopal, Rev. Mod. Phys. 73, 33 (2001).

[24] A. L¨uchow, WIREs Comput. Mol. Sci. 1, 388 (2011).

[25] J. Kolorenc and L. Mitas, Rep. Prog. Phys. 74, 026502 (2011).

[26] B. M. Austin, D. Y. Zubarev, and W. A. Lester, Chem. Rev. 112, 263 (2012).

[27] W. L. McMillan, Phys. Rev. 138, A442 (1965).

[28] M. H. Kalos and P. A. Whitlock, Monte Carlo Methods (Wiley-VCH, Weinheim, 2008).

[29] D. P. Landau and K. Binder, A Guide to Monte Carlo Simulations in Statistical Physics (Cambridge University Press, Cambridge, 2013).

[30] J. A. Pople, Rev. Mod. Phys. 71, 1267 (1999).

[31] T. Helgaker, P. Jorgensen, and J. Olsen, Molecular Electronic-Structure Theory (Wiley, Chichester, 2013).

[32] N. Metropolis, A. W. Rosenbluth, M. N. Rosenbluth, A. H. Teller, and E. Teller, J. Chem. Phys. 21, 1087 (1953).

[33] R. Jastrow, Phys. Rev. 98, 1479 (1955).

[34] S. Vitiello, K. Runge, and M. H. Kalos, Phys. Rev. Lett. 60, 1970 (1988).

[35] L. Reatto and G. L. Masserini, Phys. Rev. B 38, 4516 (1988).

[36] F. Pederiva, A. Ferrante, S. Fantoni, and L. Reatto, Phys. Rev. Lett. 72, 2589 (1994).

[37] F. Pederiva, G. V. Chester, S. Fantoni, and L. Reatto, Phys. Rev. B 56, 5909 (1997).

[38] F. Operetto and F. Pederiva, Phys. Rev. B 69, 024203 (2004).

[39] L. Dandrea, F. Pederiva, S. Gandolfi, and M. H. Kalos, Phys. Rev. Lett. 102, 255302 (2009).

[40] M. H. Kalos and L. Reatto, in Progress in Computational Physics of Matter, edited by L. Reatto and F. Manghi (World Scientific, Singapore, 1995).

[41] Note1, if some energy eigenvalues $E _ { n }$ are negative, the corresponding term is exponentially increasing instead of decaying. Nevertheless, it is always possible to add an appropriately chosen constant energy-shift to the Hamiltonian H, so that all excited components are again exponentially decaying.

[42] H. F. Trotter, Proc. Am. Math. Soc. 10, 545 (1959).

[43] R. P. Feynman and A. R. Hibbs, Quantum Mechanics and Path Integrals (McGraw-Hill, New York, 1965).

[44] H. Kleinert, Path Integrals in Quantum Mechanics, Statistics, Polymer Physics, and Financial Markets (World Scientific, Singapore, 2009).

[45] D. Ceperley and B. Alder, Science 231, 555 (1986).

[46] J. C. Slater, Phys. Rev. 34, 1293 (1929).

[47] D. M. Ceperley, Rev. Mod. Phys. 67, 279 (1995).

[48] F. Pederiva and G. V. Chester, J. Low Temp. Phys. 113, 741 (1998).

[49] F. Calcavecchia, F. Pederiva, and T. D. K¨uhne, Journal of Unsolved Questions 1, 13 (2011).

[50] F. Calcavecchia, F. Pederiva, M. H. Kalos, and T. D. K¨uhne, Phys. Rev. E 90, 053304 (2014).

[51] F. Calcavecchia and M. Holzmann, arXiv:1601.01558 (2016).

[52] M. H. Kalos, D. Levesque, and L. Verlet, Phys. Rev. A 9, 2178 (1974).

[53] D. M. Ceperley and B. J. Alder, Phys. Rev. Lett. 45, 566 (1980).

[54] F. Pederiva, S. A. Vitiello, K. Gernoth, S. Fantoni, and L. Reatto, Phys. Rev. B 53, 15129 (1996).

[55] R. P. Feynman, Phys. Rev. 94, 262 (1954).

[56] R. P. Feynman and M. Cohen, Phys. Rev. 102, 1189 (1956).

[57] A. Bijl, Physica 7, 869 (1940).

[58] R. B. Dingle, Philos. Mag. 40, 573 (1949).

[59] F. Calcavecchia and T. D. K¨uhne, Europhys. Lett. 110, 20011 (2015).

[60] P. Giannozzi et al., J. Phys.: Condens. Matter 21, 5502 (2009).

[61] J. P. Perdew, K. Burke, and M. Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996).

[62] H. J. Monkhorst and J. D. Pack, Phys. Rev. B 13, 5188 (1976).

[63] Note2, https://github.com/francesco086/HswfQMC.

[64] D. Ceperley, G. V. Chester, and M. H. Kalos, Phys. Rev. B 16, 3081 (1977).

[65] S. Sorella, Phys. Rev. B 71, 241103 (2005).

[66] Note3, the original Jackson-Feenberg expression reads as:

$$
{ \frac { \hbar ^ { 2 } } { 2 } } \sum _ { j = 1 } ^ { N } { \frac { 1 } { 2 m _ { j } } } \int _ { \Omega } d R \left( \nabla _ { j } \psi ^ { * } ( R ) \nabla _ { j } \psi ( R ) - \psi ^ { * } ( R ) \nabla _ { j } ^ { 2 } \psi ( R ) \right)
$$

[67] V. Natoli and D. M. Ceperley, Journal of Computational Physics 117, 171 (1995), ISSN 0021-9991.

[68] C. Attaccalite, Ph.D. thesis, SISSA Trieste, Italy (2005).

[69] T. Kato, Comm. Pure Appl. Math. 10, 151 (1957).

[70] Note4, the following estimated kinetic energies have been computed for 16 hydrogen atoms in the metallic atmbcc phase at $r _ { s } ~ = 1 . 3 1$ using the ASWF-pw trial WF. The employed variational parameters are: $A _ { e e } ^ { \uparrow \uparrow } = 0 . 4 2 3 $ $A _ { e e } ^ { \uparrow \downarrow } = \bar { 0 . 8 2 9 } , F _ { e e } ^ { \uparrow \uparrow } = 2 . 5 6 8 , \bar { F } _ { e e } ^ { \uparrow \downarrow } = 1 . 8 3 4 , A _ { e p } ^ { \uparrow \uparrow } = - 7 4 . 9 3 0 ,$ $A _ { e p } ^ { \uparrow \downarrow } = - 6 8 . 1 9 1 , F _ { e p } ^ { \uparrow \uparrow } = 0 . 2 3 1 , F _ { e p } ^ { \uparrow \downarrow } = 0 . 2 4 2 , C = 0 . 5 4 2$ $A _ { s s } ^ { \uparrow \uparrow } = 2 . 4 0 0 , A _ { s s } ^ { \uparrow \downarrow } = 2 . 1 1 2 , F _ { s s } ^ { \uparrow \uparrow } = 5 . 5 0 8 , F _ { s s } ^ { \uparrow \downarrow } = 1 9 . 0 3 9$ $A _ { s p } ^ { \uparrow \uparrow } ~ = ~ 2 . 4 0 0 , ~ A _ { s p } ^ { \uparrow \downarrow } ~ = ~ 2 . 1 1 2 , ~ F _ { s p } ^ { \uparrow \uparrow } ~ = ~ 5 . 5 0 8 ~ \mathrm { a n d } ~ F _ { s p } ^ { \uparrow \downarrow } ~ =$ 19.039, respectively.

[71] C. Lin, F. H. Zong, and D. M. Ceperley, Phys. Rev. E 64, 016702 (2001).

[72] Note5, magic numbers are those that close the Fermi momenta shell in a simple cubic box. For a three dimensional system these are 1, 7, 19, 27, 33, 57, 81, 93, 123, 147, 171, 179, 203, 251 . . . . .

[73] W. J. Nellis, High Pressure Research 33, 369 (2013).

[74] M. A. Morales, C. Pierleoni, E. Schwegler, and D. M. Ceperley, Proc. Nat. Acad. Sci. USA 107, 12799 (2010).

[75] E. Liberatore, M. A. Morales, D. M. Ceperley, and C. Pierleoni, Mol. Phys. 109, 3029 (2010).

[76] G. Mazzola, S. Yunoki, and S. Sorella, Nature Comm. 5, 3487 (2014).

[77] G. Mazzola and S. Sorella, Phys. Rev. Lett. 114, 105701 (2015).

[78] C. Pierleoni, M. A. Morales, C. Rillo, M. Holzmann, and D. M. Ceperley, Proc. Nat. Acad. Sci. USA 113, 4953 (2016).

[79] F. Datchi, P. Loubeyre, and R. LeToullec, Phys. Rev. B 61, 6535 (2000).