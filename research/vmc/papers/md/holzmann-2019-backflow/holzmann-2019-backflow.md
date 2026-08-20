# Orbital–dependent backflow wave functions for real–space quantum Monte Carlo

Markus Holzmann

Univ. Grenoble Alpes, CNRS, LPMMC, 3800 Grenoble, France and

Institut Laue Langevin, BP 156, F-38042 Grenoble Cedex 9, France

Saverio Moroni

CNR-IOM DEMOCRITOS, Istituto Officina dei Materiali,

and SISSA Scuola Internazionale Superiore di Studi Avanzati, Via Bonomea 265, I-34136 Trieste, Italy

We present and motivate an efficient way to include orbital dependent many–body correlations in trial wave function of real–space Quantum Monte Carlo methods for use in electronic structure calculations. We apply our new orbital–dependent backflow wave function to calculate ground state energies of the first row atoms using variational and diffusion Monte Carlo methods. The systematic overall gain of correlation energy with respect to single determinant Jastrow-Slater wave functions is competitive with the best single determinant trial wave functions currently available. The computational cost per Monte Carlo step is comparable to that of simple backflow calculations.

## I. INTRODUCTION

The fermion sign problem in general prevents electronic quantum Monte Carlo (QMC) calculations from determining unbiased ground–state properties within a controlled precision and only polynomial increasing computational cost in the number of particles. Real–space QMC methods1 avoid the sign problem through the fixed–node approximation, solving the Schr¨odinger equation with Dirichlet boundary conditions on the nodes of a trial function Ψ. While fixed–node results are often accurate, the quest for reducing the systematic error incurred has prompted generalizations of the standard Jastrow–Slater (JS) wave function. Better wave functions are obtained replacing the Slater determinant by a multideterminant expansion2, antisymmetrized geminal product $( \mathrm { A G P } ) ^ { 3 }$ or Pfaffian (PF)4. As an alternative or in addition, backflow (BF) transformations5 can be applied to the particles’ coordinates. All these variations include correlation effects in the nodal structure of Ψ, which in turn determines the accuracy of the fixed–node approximation.

In this paper we introduce a way of including electron correlations in the antisymmetric factor of Ψ improving the nodal structure of strongly inhomogeneous systems. Whereas in previous BF wave functions the ith particle’s coordinate $\mathbf { r } _ { i }$ in the argument of the nth single–particle orbital is substituted by the BF–transformed coordinate qi [e.g. given by Eq. (6) below],

$$
\phi _ { n } ( { \bf r } _ { i } ) \longrightarrow \phi _ { n } [ { \bf q } _ { i } ( X ) ] ,\tag{1}
$$

where X specifies the configuration of the system (e.g. electronic and nuclear coordinates), we instead replace each orbital by two or more orbitals coupled via BF correlations in the amplitudes,

$$
\phi _ { n } ( \mathbf { r } _ { i } ) \longrightarrow \phi _ { n } ^ { ( 1 ) } ( \mathbf { r } _ { i } ) + \left[ \mathbf { q } _ { i } ( X ) - \mathbf { r } _ { i } \right] \cdot \nabla \phi _ { n } ^ { ( 2 ) } ( \mathbf { r } _ { i } ) .\tag{2}
$$

Here $\phi _ { n } ^ { ( a ) } , a = 1 , 2 , . . .$ , denote reoptimized orbitals of the same spatial symmetry as $\phi _ { n }$ , specific to the nth orbital.

Thus, the same BF transformation $\mathbf { q } _ { i } ( X )$ affects differently the various orbitals describing the antisymmetric part of Ψ. We call ”orbital backflow” (OBF) this way of using the transformed coordinates.

The OBF functional form is motivated in Sec. II using the local energy method6,7 for a single–determinant wave function, and applied to the first row atoms in Sec. III, where it proves competitive with inhomogeneous backflow (IBF)8,9, AGP3,10 and PF4,11 wave functions.

## II. ORBITAL BACKFLOW TRIAL WAVE FUNCTION

We briefly outline how normal and orbital backflow may emerge naturally from approximating a generalized Feynman-Kac path integral formula. We are merely interested in possible functional forms, suitable for numerical evaluation, so that most of the approximations in this section are driven more by the need of simplifcation than by mathematical rigour. Thus, anticipating the eventual optimization of the functional parameters of any resulting trial wave function, we use variational freedom already in intermediate simplification steps, to modify some of the detailed expressions into plausible functional forms suggested by physical intuition. The notation $\widetilde f ( \cdot )$ will be used to indicate changes of an explicit function f (·) due to parameter optimization.

The ratio between the exact ground–state wave function Φ(R) and a trial wave function $\Psi _ { 0 } ( R )$ not orthgonal to Φ is7,12,13

$$
\frac { \Phi ( R ) } { \Psi _ { 0 } ( R ) } \propto \langle e ^ { - \int _ { 0 } ^ { \infty } E _ { L } ( R ( t ) ) d t } \rangle ,\tag{3}
$$

where $\begin{array} { l } { R ~ = ~ \left( \mathbf { r } _ { 1 } , \ldots , \mathbf { r } _ { N } \right) } \end{array}$ are the coordinates of the N particles, and the brackets denote the average over all random walks R(t) starting at R generated by the importance–sampled Green’s function. The local energy method6,7 uses an analytic approximation of Eq. (3) to give an explicit expression for an improved wave function Ψ in terms of $\Psi _ { 0 }$ and its local energy $E _ { L } ( R ) \ =$ $\langle R | H | \Psi _ { 0 } \rangle / \langle R | \Psi _ { 0 } \rangle$ ,

$$
\frac { \Phi ( R ) } { \Psi _ { 0 } ( R ) } \approx e ^ { - \langle \int _ { 0 } ^ { \tau } E _ { L } ( R ( t ) ) d t \rangle } \approx e ^ { - \tau \widetilde { E } _ { L } ( R ) } \equiv \frac { \Psi ( R ) } { \Psi _ { 0 } ( R ) } .\tag{4}
$$

The approximations underlying Eq. (4) are the truncation of the cumulant expansion at first order over a finite projection time τ, and the assumption that the random walk average of time integrals of $E _ { L } [ R ( t ) ]$ merely reproduces the same functional form of the local energy, but with a smoother R dependence in the relevant phasespace region. The resulting expression $\widetilde { E } _ { L } ( R )$ in the exponent of the improved wave function is therefore given by a functional expression similar to $E _ { L } ( R )$ containing modified/optimized pseudopotentials and orbitals.

We take $\Psi _ { 0 }$ as a simple wave function with a Jastrow factor $e ^ { - U ( R ) }$ and a Hartree product of single–particle orbitals $\phi _ { n } ( \mathbf { r } _ { i } )$ (the antisymmetrization being applied afterwards, on the improved wave function Ψ). The modified local energy $\widetilde { E } _ { L } ( R )$ then contains terms proportional to $\nabla _ { i } \widetilde { U } ( R ) \cdot \nabla _ { i }$ i ln $\phi _ { n } ( \mathbf { r } _ { i } )$ . Specializing further to a two–body Jastrow factor $\begin{array} { r } { \dot { U ( R ) } = \sum _ { i < j } u ( \dot { r } _ { i j } ) } \end{array}$ , Eq. (4) suggests that the one–particle orbitals in the Slater determinant of the improved wave function Ψ are given by

$$
\widetilde { \phi } _ { n } ( \mathbf { r } _ { i } ) \exp \left[ \sum _ { j \neq i } \frac { \widetilde { u } ^ { \prime } } { r _ { i j } } ( \mathbf { r } _ { i } - \mathbf { r } _ { j } ) \cdot \nabla _ { i } \ln \widetilde { \phi } _ { n } ( \mathbf { r } _ { i } ) \right] .\tag{5}
$$

When ln $\phi _ { n }$ is linear in $\mathbf { r } _ { i } , \mathrm { e . g . } - i \mathbf { k } _ { n } \cdot \mathbf { r } _ { i }$ for plane waves of wave vector $\mathbf { k } _ { n }$ describing homogeneous systems, we recover the familiar case of Eq. (1) with the usual backflow transformation

$$
\mathbf { q } _ { i } = \mathbf { r } _ { i } + \sum _ { j \neq i } \eta ( r _ { i j } ) ( \mathbf { r } _ { i } - \mathbf { r } _ { j } ) ,\tag{6}
$$

where $\eta = \widetilde { u } ^ { \prime } / r _ { i j }$

Whereas the cumulant expansion in the local energy method guarantees the extensivity of the logarithm of Ψ for extended systems, this approximation may poorly describe modifications of strongly inhomogeneous, localized orbitals. Local modifications of orbitals may better be captured by keeping only the linear term of the exponential of Eq. (5). By further choosing different modified orbitals $\phi _ { n } ^ { ( a ) }$ for each $n ,$ to improve the variational flexibility of our trial wave function, the OBF form of Eq. (2) is obtained. In our case, $\mathbf { q } _ { i }$ remains a simple backflow coordinate with homogeneous two–body correlations of the form given by Eq. (6).

Let us stress that the sequence of approximations made to simplify Eq. (3) are rather crude and remain on a heuristic level. However, the procedure is not aimed to directly obtain accurate expressions, but to suggest flexible functional forms for the trial function suitable for approximating the exact ground state within polynomial computational cost. The quality of the resulting functional form is determined a posteriori for specific systems after optimization of the radial function η and the modified orbitals $\phi _ { n } ^ { ( a ) }$ involved.

## III. CASE STUDY OF THE FIRST ROW ATOMS

To benchmark the accuracy of the OBF wave function we have calculated the energies of all–electron first row atoms with variational Monte Carlo (VMC) and fixed– node diffusion Monte Carlo (DMC) for a trial wave function represented by the product of a Jastrow factor and a single–determinant per spin component composed from backflow improved orbitals according to the transformation (2).

In particular, s orbitals now obtain the following form

$$
\phi _ { n } ^ { s } ( \mathbf { r } _ { i } , \mathbf { q } _ { i } ) = \chi _ { n } ^ { ( 1 ) } ( r _ { i } ) + ( \mathbf { q } _ { i } - \mathbf { r } _ { i } ) \cdot \mathbf { r } _ { i } \chi _ { n } ^ { ( 2 ) } ( r _ { i } ) ,\tag{7}
$$

where $\mathbf { q } _ { i }$ is given by Eq. (6) using different η functions for like– and unlike– spin electrons expressed as locally piecewise–quintic Hermite interpolants $( \mathrm { L P Q H I } ) ^ { 1 4 }$ , and $\chi _ { j } ^ { ( \alpha ) }$ are radial functions expanded in a basis of Slater type orbitals15. The p orbitals read

$$
\begin{array} { c } { { \phi _ { n } ^ { p _ { \alpha } } ( \mathbf { r } _ { i } , \mathbf { q } _ { i } ) = r _ { i } ^ { \alpha } \chi _ { n } ^ { ( 1 ) } ( r _ { i } ) + ( q _ { i } ^ { \alpha } - r _ { i } ^ { \alpha } ) \chi _ { n } ^ { ( 2 ) } ( r _ { i } ) } } \\ { { + r _ { i } ^ { \alpha } ( \mathbf { q } _ { i } - \mathbf { r } _ { i } ) \cdot \mathbf { r } _ { i } \chi _ { n } ^ { ( 3 ) } ( r _ { i } ) } } \end{array}\tag{8}
$$

where α is the cartesian component required in the nth orbital. Instead of using $[ \partial \chi _ { n } ^ { ( 2 ) } ( r ) / \partial r ] / r$ as suggested by Eq. (2), we introduced a third independent radial function $\chi _ { n } ^ { ( 3 ) } ( r )$ for increased variational freedom.

Implementation of OBF is rather straightforward by considering both the particles’ coordinates $\mathbf { r } _ { i }$ and the renormalized BF coordinates $\mathbf { q } _ { i }$ as independent variables of the modified orbitals on the right–hand side of (2). Gradient and Laplacian of the trial wave functions are then obtained by applying the chain rule in the same way as for standard $\mathrm { \bar { B F } } ^ { \bar { 6 } }$ . Compared to a direct inclusion of orbital–dependent BF correlations through different coordinate transformations for different orbitals, the computational cost of our OBF wave function thus maintains the overall $N ^ { 3 }$ scaling of standard BF, with a small additional cost of less than a factor 2. The increased number of independent terms in each orbital can be dealt with by modern optimization techniques16,17.

The symmetric Jastrow factor of our case study on first row atoms contains an electron–electron term $\begin{array} { r } { \prod _ { i < j } \exp [ - u ( r _ { i j } ) ] } \end{array}$ with different pseudopotentials u for like and unlike spins, an electron–nucleus term $\begin{array} { r } { \prod _ { i } \exp [ - w ( r _ { i } ) ] } \end{array}$ , and electron–electron–nucleus correlations

$$
\prod _ { i \neq j } \exp \{ - [ \xi _ { 0 } ( r _ { i } ) \xi _ { 0 } ( r _ { j } ) - \xi _ { 1 } ( r _ { i } ) \xi _ { 1 } ( r _ { j } ) { \bf r } _ { i } \cdot { \bf r } _ { j } ] \} .\tag{9}
$$

All radial functions u, w, $\xi _ { 0 }$ , and $\xi _ { 1 }$ are expressed as LPQHI. The variational parameters (58 for Li and Be, 67 for the other atoms) are optimized by minimization of the variational energy16. The resulting VMC and DMC energies obtained are listed in Table I.

<!-- image-->  
TABLE I. Energies in Hartree a.u. of the first row atoms obtained with VMC and fixed–node DMC using the OBF wave function.

Energies of the first row atoms have been calculated by several authors using a variety of different trial wave functions beyond the simple JS form providing useful comparisons. In Figs. 1 and 2 our OBF data from Table I, indicated by full red circles, are compared with selected results from the literature, as indicated by the labels in the body of the figures with the reference in brackets (unpublished calculations10 using AGP and PF, and an earlier AGP result from Ref. 3; the IBF energies from Ref. 9 for VMC in Fig. 1 and from Ref. 8 for DMC in Fig. 2; a PF calculation4 and its version (PFBF) with IBF included, and a general PF form dubbed STU11 which encompasses both singlet and triplet pairing, as well as unpaired orbitals). For AGP and IBF, subsequent calculations with the same kind of wave function found lower energies on account of more aggressive optimization and/or use of extended basis sets. We also show by empty symbols in Fig. 1 the VMC JS result from the respective sources. Orbital dependent Jastrow correlations applied to the oxygen atom18 have not led to significant improvement compared to the common JS trial wave function within VMC and DMC.

The most pertinent and systematic comparison is possible between IBF and OBF. Both are ways to include backflow effects in inhomogeneous systems, where standard backflow of the form of Eq. (1) with the simple coordinate transformation of Eq. (6) does not significantly lower the energies. Within IBF, nuclear coordinates are included inside the standard BF transformation through atom–specific electron–nucleus7,19 and electron– electron–nucleus8 terms. Instead, OBF only uses the basic BF transformation with homogeneous electronelectron term, but introduces an orbital–specific dependence through the modified orbitals of Eq. (2). We mention that a yet different BF wave function, featuring iterative coordinate renormalization20, gives excellent results for both homogeneous and inhomogeneous strongly correlated systems21. However, it becomes less beneficial as correlations weaken, providing only marginal improvements for the first row atoms.

Whereas OBF obtains slightly less correlation energy than IBF9 within VMC (see Fig. 1), a small gain relative to IBF is obtained by OBF at the level of fixed-node DMC, except for lithium and neon where they are very close. We attribute the qualitatively different behavior of VMC and DMC to a better parametrization of the symmetric Jastrow factor of the IBF in Ref.9 compared to the Jastrow factor of the present work, clearly visible in the difference in the respective bare Slater-Jastrow (JS) data (see Fig. 1). We further note that the DMC results for single–determinant IBF are only provided by Ref. 8. In Ref. 9, the IBF wave function was optimized much better, lowering the VMC energies particularly in the case of beryllium and boron, but DMC values have not been provided. It is natural to expect that the better optimized IBF wave function will also lower the corresponding DMC values, thus reducing the rather large difference between IBF and OBF of those two atoms shown in Fig. 2. For the other atoms, however, the variational quality of single– determinant IBF wave functions given in Refs. 8 and 9 is very similar, and the DMC results of Ref. 8 shown in Fig. 2 should be representative of a well–optimized IBF. Overall, it seems fair to conclude that the OBF nodes tend to provide a slightly better description than those of IBF.

<!-- image-->  
FIG. 1. Fraction of correlation energy recovered in VMC for the first row atoms of orbital backflow wave functions (OBF) compared to previous results using different kinds of single–determinant wave function, namely antisymmetrized geminal product (AGP)3,10, Pfaffian (PF)4, inhomogeneous backflow (IBF)8,9, Pfaffian including inhomogeneous backflow (PFBF)4, and a general Pfaffian form dubbed STU11. Empty symbols denote the bare Slater-Jastrow result from the respective sources. Small horizontal shifts of some data are added for clarity. The Hartree–Fock and estimated exact energies are taken from Table I of Ref. 8.

<!-- image-->  
FIG. 2. Fraction of correlation energy recovered in fixed–node DMC for the first row atoms of orbital backflow wave functions (OBF) compared to previous results of different kinds of single–determinant wave functions. Notations are the same as those of Fig. 1.

Pairing wave functions (AGP, PF or STU) overcome some topological inadequacies22 that Hartree-Fock nodes share with single–determinant wave functions, including OBF and IBF. The excellent result –particularly in DMC– obtained for beryllium with the AGP wave function is due to its multi-determinant character3 with respect to the nearly–degenerate p orbitals. However, for the heavier atoms, single–determinant wave functions with OBF or IBF provide essentially as good energies as pairing wave functions. The further, non–negligible gain obtained by inclusion of IBF in a PF wave function4 (see the PFBF points in Figs. 1 and 2) suggests that pairing and backflow correlations improve complementary aspects of the wave function, at least to some extent.

So far, the discussion and our comparison in Figs. 1 and 2 has been restricted exclusively to single– determinant trial wave functions. For small atoms, nearly exact energies can be retrieved in QMC multideterminant expansions with a modest number of terms9,23,24. However, the same accuracy cannot be maintained for heavier atoms or molecules with an affordable number of determinants, whereas backflow wave functions should improve the accuracy fairly independent of the number of electrons with only polynomial increasing computational cost.

We briefly mention that a simple linear extrapolation of the JS and OBF energies against the corresponding variancies of the energy, as succesfully used in strongly correlated homogeneous quantum fluids20,21, does not provide any systematic improvement for the first row atoms’ energies. The discrete nature of the density of states for the electrons in the nuclear potential seems to considerably shrink down the region of validity for such extrapolations.

## IV. DISCUSSION

Orbital–dependent backflow wave functions provide a simple and efficient way of introducing and tuning a physically appealing orbital dependence in many–body correlations. We have shown that the resulting gain in energy for first row atoms is competitive with inhomogeneneous backflow wave functions9, the best currently available single–determinant trial wave function for electronic structure of atoms and molecules. For more complex systems, the orbital dependence of OBF presents an appealing alternative to the atom–specific IBF transformation, and may be better suited to study orbital–selective phenomena in strongly correlated systems25.

Variational flexibility of OBF is added by two main ingredients, a coordinate renormalization, Eq. (6), and an orbital modification, Eq. (2). The latter could be used without the former, for instance replacing the backflow coordinate by the fluctuation of a local dipole or by the wave vector of a density fluctuation in an extended system in the scalar product with the gradient. In such a version, OBF would correspond to an earlier representation26 of backflow correlations used in lattice models27.

Generically, OBF provides a modification of orbitals enlarging the functional flexibility of trial wave functions suited for standard real space QMC methods. It can directly be combined with IBF including further dependency on the electron–nucleus distances or electron– electron–nucleus in the backflow coordinates, and extended to iterated backflow wave functions, as well as in the use of Pfaffian and multi–determinant trial wave function. Efficient optimization among all possible combinations may request improved optimization strategies28.

Despite obvious limitations, the accuracy reached by real–space QMC methods should be sufficient to tackle and provide new insights to the role of correlation in electronic structure. Flexible trial wave functions capturing different aspects of correlation put up a frame to estimate and reduce the bias of the underlying trial wave function. Still, reliabe estimates and/or control of the the systematic error of fixed–node QMC involving large number of electrons remains challenging.

1 J. Kolorenˇc. and L. Mitas, Rep. Prog. Phys. 74, 026502 (2011).

2 C. Filippi and C. J. Umrigar, J. Chem. Phys. 105, 213 (1996).

3 M. Casula, C. Attaccalite, and S. Sorella, J. Phys. Chem. 121, 7110 (2004).

4 M. Bajdich, L. Mitas, L. K. Wagner and K. E. Schmidt, Phys. Rev. B 77, 115112 (2008).

5 K. E. Schmidt, M. A. Lee, M. H. Kalos, and G. V. Chester, Phys. Rev. Lett. 47, 807 (1981).

6 Y. Kwon, D. M. Ceperley and R. M. Martin, Phys. Rev. B 48, 12037 (1993).

7 M. Holzmann, D. M. Ceperley, C. Pierleoni and K. Esler, Phys. Rev. E 68, 046707 (2003).

8 M. D. Brown, J. R. Trail, P. L´opez R´ıos, and R. J. Needs, J. Chem. Phys. 126, 224110 (2007).

9 P. Seth, P. L´opez R´ıos, and R. J. Needs J. Chem. Phys. 134, 084105 (2011).

10 S. Sorella, private communication.

11 M. Bajdich, L. Mitas, G. Drobn´y, L. K. Wagner, and K. E. Schmidt, Phys. Rev. Lett. 96, 130201 (2006).

12 K. S. Liu, M. H. Kalos, and G. V. Chester, Phys. Rev. A 10, 303 (1974).

13 M. Caffarel and P. Claverie, J. Chem. Phys. 88, 1088 (1988).

14 V. Natoli and D. M. Ceperley, J. Comput. Phys. 117, 171 (1995).

## ACKNOWLEDGMENT

We thank the Fondation NanoSciences (Grenoble) and the CNRS, INP for support.

15 E. Clementi and C. Roetti, At. Data Nucl. Data Tables 14, 177 (1974).

16 S. Sorella, M. Casula, and D. Rocca, J. Chem. Phys. 127, 014105 (2007).

17 C. J. Umrigar, J. Toulouse, C. Filippi, S. Sorella, and R. G. Hennig, Phys. Rev. Lett. 98, 110201 (2007).

18 T. Bouab¸ca, B. Bra¨ıda, and M. Caffarel, J. Chem. Phys. 133, 044111 (2010).

19 P. L´opez R´ıos, A. Ma, N. D. Drummond, M. D. Towler, and R. J. Needs, Phys. Rev. E 74, 066701 (2006).

20 M. Taddei, M. Ruggeri, S. Moroni, and M. Holzmann, Phys. Rev. B 91, 115106 (2015).

21 M. Ruggeri, S. Moroni, and M. Holzmann, Phys. Rev. Lett. 120, 205302 (2018).

22 K. M. Rasch and L. Mitas, Chem. Phys. Lett. 528, 59 (2012).

23 M. A. Morales, J. McMinis, B. K. Clark, J. Kim, and G. E. Scuseria, J. Chem. Theory Comput. 8 (7), 2181 (2012).

24 J. Toulouse and C. J. Umrigar, J. Chem. Phys. 128, 174101 (2008).

25 See, e.g., V. I. Anisimov, I. A. Nekrasov, D. E. Kondakov, T. M. Rice and M. Sigrist, Eur. Phys. J. B 25, 191-201 (2002).

26 L. F. Tocchio, F. Becca, A. Parola, and S. Sorella, Phys. Rev. B 78, 041101(R) (2008).

27 Di Luo and B K. Clark, arXiv:1807.10770.

28 D. Kochkov and B. K. Clark, arXiv:1811.12423.