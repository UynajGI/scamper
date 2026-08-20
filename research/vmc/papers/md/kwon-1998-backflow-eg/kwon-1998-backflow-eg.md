# Effects of Backflow Correlation in the Three-Dimensional Electron Gas: Quantum Monte Carlo Study

Yongkyung Kwon Department of Physics and Center for Advanced Materials and Devices, Kon-Kuk University, Seoul 143-701, Korea

D. M. Ceperley Department of Physics and National Computational Science Alliance, University of Illinois, Urbana, Illinois 61801

Richard M. Martin Department of Physics and Materials Research Laboratory, University of Illinois, Urbana, Illinois 61801 (July 2, 2018)

## Abstract

The correlation energy of the homogeneous three-dimensional interacting electron gas is calculated using the variational and fixed-node diffusion Monte Carlo methods, with trial functions that include backflow and three-body correlations. In the high density regime (rs ≤ 5) the effects of backflow dominate over those due to three-body correlations, but the relative importance of the latter increases as the density decreases. Since the backflow correlations vary the nodes of the trial function, this leads to improved energies in the fixednode diffusion Monte Carlo calculations. The effects are comparable to those found for the two-dimensional electron gas, leading to much improved variational energies and fixed-node diffusion energies equal to the release-node energies of Ceperley and Alder within statistical and systematic errors.

71.10.Ca, 71.15.Nc, 71.15.Pd

## I. INTRODUCTION

The homogeneous electron gas in three dimensions is the simplest model to study the effects of correlation between electrons in metals [1]. Its correlation energy, defined as the total ground-state energy minus the Hartree-Fock energy, has been used to give the exchange-correlation potential in density functional calculations with the Local Density Approximation [2,3]. The possible phases which this simple system can display are prototypes for understanding interacting electrons in extended matter [4–8].

The theoretical study of the interacting electron gas began with Bloch [4] who discovered by using the Hartree-Fock approximation that the system would favor a ferromagnetic liquid state over the normal paramagnetic state at low electron densities. Wigner [5] first calculated the correlation energy of the homogeneous electrons at high density limit, using the secondorder perturbation theory. He also pointed out that for sufficiently low densities the electrons would become localized and form an ordered array. After calculating the correlation energy of this electron solid with the Wigner-Seitz approximation [1,8], he proposed an interpolation formula for the correlation energy in a wide range of densities having his high- and lowdensity limits. The development of the field-theoretic approaches in 1950s led to various approximate methods [8] to calculate the ground-state properties of the electron gas. Among them, Gell-Mann and Brueckner [9] summed the ring diagrams to compute the correlation energy in the high-density limit. The dielectric function formalism, especially with the self-consistent treatment of the screening process introduced by Singwi, Tosi, Land, and Sj¨olander [10], gave more accurate ground-state properties at a wider range of densities.

On the other hand, various stochastic numerical methods known collectively as quantum Monte Carlo (QMC) have been developed to compute the properties of a quantum many-body system such as the electron gas. Ceperley [11] first applied variational Monte Carlo (VMC) to calculate a much more accurate upper bound to the ground-state energy of the electron gas than is given by Hartree-Fock. More accurate correlation energies were computed by Ceperley and Alder [12] with the diffusion Monte Carlo (DMC) method which projects the true ground state of a many-body system from a trial state. Even though the DMC method gives the exact ground-state energy for a system of many bosons, it has a serious difficulty in treating fermion systems, because fermion wave functions must be antisymmetric under particle exchanges [13]. In order to address this problem, Ceperley and Alder [12] developed the released-node method. Its only limitation is that the statistical fluctuations can grow rapidly at large projection time. So the statistical noise can dominate the signal before converging to the ground state.

In this work we use the fixed-node method [14,12,15], where the nodal surface of the exact ground-state wave function is approximated by that of the trial wave function. We adopt the approach of systematically improving the fixed-node DMC results by using a trial function with better nodes, analogous to our work on the two dimensional electron gas [16]. Unlike the released-node method this method is stable and does not have the convergence problem. It gives the best upper bound to the exact energy consistent with the assumed nodes.

Ceperley and Alder used the Slater-Jastrow trial wave function in both released-node [12] and fixed-node calculations [17], which consists of the Slater determinant of single-body orbitals and products of two-body correlation functions. According to their released-node calculation, the electron gas could exhibit three different phases at zero temperature, the paramagnetic and ferromagnetic liquids and the Wigner crystal, depending on its density. More recently, Ortiz and Ballone [18] reported a new fixed-node DMC calculations with the Slater-Jastrow wave function. Their correlation energies were found to be, as expected, smaller in magnitude than Ceperley and Alder’s released-node calculations, especially at high metallic densities. In this paper, we perform a fixed-node DMC calculation using a trial function with backflow and three-body correlations in addition to two-body correlation. Our calculations in the two-dimensional electron gas [16] showed that the inclusion of the backflow correlation in a trial state greatly improved the Slater-Jastrow fixed-node results. We found that, although the Slater-Jastrow wave function accounts for most of the correlation energy of the electron gas, the remaining errors that are in the Slater-Jastrow function are mostly due to backflow and three-body correlations. This is similar to what has been observed on calculations of the other strongly correlated system of fermions, liquid 3He [19]. Our fixed-node DMC results in three dimensions will be compared with Ceperley and Alder’s released-node results.

All ground-state properties of the electron gas at zero magnetic fields are determined only by the dimensionless density parameter $r _ { s } = a / a _ { 0 }$ , where $a _ { 0 }$ is the Bohr radius, $\textstyle a = { \bigl ( } { \frac { 3 } { 4 \pi \rho } } { \bigr ) } ^ { 1 / 3 }$ is the radius of a sphere which encloses one electron on the average and $\rho$ is the number density. With energy units of Rydbergs (Ry) and the length units of $a ,$ the Hamiltonian of the electron gas is

$$
H = - \frac { 1 } { r _ { s } ^ { 2 } } \sum _ { i = 1 } ^ { N } \nabla _ { i } ^ { 2 } + \frac { 2 } { r _ { s } } \sum _ { i < j } \frac { 1 } { | { \bf r } _ { i } - { \bf r } _ { j } | } ~ + ~ c o n s t a n t ,\tag{1}
$$

where the constant is the term due to the uniform background of opposite charge. We consider the density range of $1 \leq r _ { s } \leq 2 0$ , where Ceperley and Alder found the system in the normal liquid phase. We do not consider spin polarized or superconducting states.

## II. METHODOLOGY

In a VMC calculation, one estimates the properties of a quantum state, by assuming a trial wave function $\Psi _ { T } ( R )$ with the correct symmetry, where $R = \left( \mathbf { r } _ { 1 } , \mathbf { r } _ { 2 } , \cdots , \mathbf { r } _ { N } \right)$ is a 3Ndimensional vector representing the positions of N particles. With a set of configurations $\{ R _ { i } \}$ sampled with a probability density proportional to $\Psi _ { T } ^ { 2 } ( R )$ , the variational energy is just the average of local energies, $E _ { L } ( R _ { i } ) = H \Psi _ { T } ( R _ { i } ) / \Psi _ { T } ( R _ { i } )$ . This method can give a good upper bound to the exact energy if the trial state is accurate as it is for the homogeneous electron gas [11].

Even more accurate ground-state properties of a many-body system can be obtained with the DMC method, where the Schr¨odinger equation is solved by treating it as a diffusion equation [12]. The solution of the Schr¨odinger equation in imaginary time $t , - \partial | \Phi \rangle / \partial t =$ $( { \hat { H } } - E _ { T } ) { \dot { | \Phi } } \rangle$ , can be expressed in terms of the exact energy eigenvalues $E _ { i }$ and eigenstates $\phi _ { i } \colon$

$$
\Phi ( R , t ) = \sum _ { i } c _ { i } \exp [ - t ( E _ { i } - E _ { T } ) ] \phi _ { i } ( R ) .\tag{2}
$$

At sufficiently long times only the ground state $\phi _ { 0 }$ survives in Eq. (2), if $\Phi ( R , 0 )$ is not orthogonal to it. In order to implement this idea with a stochastic procedure, we consider the real-space representation of the Schr¨odinger equation:

$$
\frac { \partial f ( R , t ) } { \partial t } = \frac { 1 } { r _ { s } ^ { 2 } } \sum _ { i = 1 } ^ { N } \nabla _ { i } \cdot ( \nabla _ { i } f - f \nabla _ { i } \ln \Psi _ { T } ^ { 2 } ) - \left( E _ { L } ( R ) - E _ { T } \right) f \ ,\tag{3}
$$

where $f ( R , t ) = \Phi ( R , t ) \Psi _ { T } ( R )$ . Note that the Schr¨odinger equation is multiplied by the trial wave function $\Psi _ { T } ( R )$ Eq. (3) can be viewed as a diffusion equation in a 3N-dimensional space with the density of diffusing particles $f ( R , t )$ . Its second term imposes a drift and the final term gives rise to a branching process by which the sampled configurations converge to the lowest-energy state. The initial ensemble of configurations $\{ R \}$ with probability density $f ( R , 0 ) = \Psi _ { T } ^ { 2 } ( R )$ is evolved forward in time by the above diffusion equation and reaches the equilibrium distribution $f ( R , \infty ) = \phi _ { 0 } ( R ) \Psi _ { T } ( R )$ at large enough t. From this distribution of random walks, the exact ground-state energy $E _ { 0 } = \langle \phi _ { 0 } | \hat { H } | \Psi _ { T } \rangle / \langle \phi _ { 0 } | \Psi _ { T } \rangle$ can be estimated as the average of the local energy: $E _ { L } ( R ) = H \Psi _ { T } ( R ) / \Psi _ { T } ( R )$

The diffusion equation formulation described above requires for implementation that the population density $f ( R , t )$ be non-negative. For Bose systems, this is not a problem since their ground-state wave functions can be chosen to be non-negative. However, fermion wave functions are antisymmetric, change sign, and have nodes. This leads to the famous sign problem [13] in the QMC calculations of Fermi systems. The apparent limitation of the diffusion analogy in this case can be dealt with by treating positive and negative regions separately. One easy way to accomplish this is not to allow diffusion between these two regions, which corresponds to the f ixed-node approximation [14,12]. If $\Psi _ { T }$ were to have the exact nodes of the ground state, one could treat the fermion system immediately and exactly, since f would never change sign. Unfortunately, the exact location of the nodes in many-fermion systems is not known [20]. The fixed-node approximation is based on the requirement that $\phi _ { 0 } ( R ) \Psi _ { T } ( R )$ be non-negative. The fixed-node DMC energy is an upper bound to the exact energy; the best upper bound with the given nodes, and usually lies well below the variational energy [20].

Another way to deal with the sign problem is to use the released-node method [12] which puts no constraints on the nodal structure of the true ground-state wave function. In this method, there is a population of positive random walks which give positive contributions to any average, and a population of negative walks with negative contributions. Whenever a random walk diffuses across the nodes of the trial function, the sign of its contribution changes. Even though it can be shown that the difference population converges to the antisymmetric fermion ground state, it does not proceed without problems. Since both the positive and negative populations grow geometrically with a large number of projections, the statistical fluctuations in the average increase exponentially [12]. So for this method to be successful the diffusion process needs to converge to the ground state before the fluctuations become large. As the system size gets larger, the fluctuations grow. Hence, the fixed-node method is more useful for systems with many fermions.

## III. MONTE CARLO CALCULATIONS

## A. Trial Wave Function

In all QMC methods mentioned above, a good trial function is very important for accurate results. The convergence time in a released-node calculation can be reduced with a better trial function while the nodes of a trial function determine the ultimate accuracy of a fixed-node calculation. The usual choice of a trial function is of the Slater-Jastrow type

$$
\Psi _ { T } ( R ) = \mathrm { d e t } ( \varphi _ { m n } ) \exp [ - \sum _ { i < j } ^ { N } u ( r _ { i j } ) ] ,\tag{4}
$$

where $\varphi _ { m n } = e ^ { i \mathbf { k } _ { m } \cdot \mathbf { r } _ { n } }$ for a homogeneous liquid phase. The nodes are determined by only the Slater determinant. We use the two-body correlation function $u ( r )$ that minimizes the variational energy in the Random Phase Approximation [21,11]. With this trial wave function, Ceperley calculated the ground-state properties of the electron gas, using the VMC [11], the fixed-node [17] and the released-node DMC method [12].

In order to improve the nodes, we consider a more complicated trial function which includes backflow and three-body correlations [16]. Our wave function has the form of

$$
\Psi _ { T } ( R ) = \operatorname * { d e t } ( e ^ { i { \bf k } _ { i } \cdot { \bf x } _ { j } } ) \exp [ - \sum _ { i < j } ^ { N } \widetilde { u } ( r _ { i j } ) - \frac { \lambda _ { T } } { 2 } \sum _ { l = 1 } ^ { N } { \bf G } ( l ) \cdot { \bf G } ( l ) ] ,\tag{5}
$$

where $\mathbf { x } _ { i } ^ { \prime } \mathrm { s }$ are quasiparticle coordinates defined as:

$$
{ \bf x } _ { i } = { \bf r } _ { i } + \sum _ { j \neq i } ^ { N } \eta ( r _ { i j } ) \left( { \bf r } _ { i } - { \bf r } _ { j } \right) ,\tag{6}
$$

$$
{ \bf G } ( l ) = \sum _ { i \neq l } ^ { N } \xi ( r _ { l i } ) \left( { \bf r } _ { l } - { \bf r } _ { i } \right) ,\tag{7}
$$

and

$$
\tilde { u } ( r ) = u ( r ) - \lambda _ { T } \xi ^ { 2 } ( r ) r ^ { 2 } .\tag{8}
$$

In addition to the two-body correlation, this trial function includes the three-body correlation, $\mathbf { G } ( l ) \cdot \mathbf { G } ( l )$ , and the state-dependent correlation, $\mathbf { k } \cdot ( \mathbf { r } _ { i } - \mathbf { r } _ { j } ) \eta ( r _ { i j } )$ , which incorporates the hydrodynamic backflow [22]. We call $\xi ( r )$ the “three-body correlation function” and $\eta ( r )$ the “backflow correlation function”. Note that it is the backflow correlation which makes the nodes of the wave function different from those of the Slater-Jastrow trial function.

Our calculations are done for N electrons in a cube with periodic boundary conditions. The Ewald method [23] is used for the Coulomb potential and the two-body correlation $u ( r )$ to minimize size effects. The higher-order correlation functions, $\eta ( r )$ and $\xi ( r )$ , are required to go to zero smoothly at a cutoff distance $r _ { c }$ set to half the side of the simulation cell we use:

$$
f ( r ) \longrightarrow f ( r ) + f ( 2 r _ { c } - r ) - 2 f ( r _ { c } ) .\tag{9}
$$

The backflow and the three-body correlation function are parametrized as

$$
\eta ( r ) = \lambda _ { B } { \frac { 1 + s _ { B } r } { r _ { B } + w _ { B } r + r ^ { 4 } } } ,\tag{10}
$$

and

$$
\xi ( r ) = \exp [ - ( r - r _ { T } ) ^ { 2 } / w _ { T } ^ { 2 } ] \ .\tag{11}
$$

This functional form for $\eta ( r )$ satisfies the long-range behavior $( \sim 1 / r ^ { 3 } )$ in three dimensions predicted by the local-energy method of Ref. [16]. It should be noted that the optimized $\eta ( r )$ goes to zero rapidly at the edge of the simulation box. Our three-body correlation has the same form as used for liquid 3He in Ref. [24].

In order to optimize our higher-order correlation functions, we minimize the variance of the local energy [25], defined by

$$
V _ { \Psi _ { T } } = \frac { \int d R \Psi _ { T } ^ { 2 } ( R ) \ : ( \ : E _ { L } ( R ) - E _ { v } \ : ) ^ { 2 } } { \int d R \Psi _ { T } ^ { 2 } ( R ) } .\tag{12}
$$

If our trial function $\Psi _ { T }$ were an exact eigenfunction of the Hamiltonian, the variance would be zero. Because the variance is a non-linear function of the parameters we cannot be certain that we have achieved converged results for this class of trial functions.

The optimum variational parameters that we have obtained as a function of density, are given in Table I. Fig. (1) shows the effect (in the logarithm of the wave function) of two electrons a distance r apart coming from the three-body term. There is a very strong density dependence. The effect is almost negligible at $r _ { s } = 1$ but as large as 10% at $r _ { s } \ge 1 0$ . Negative values imply that electron configurations in which the “forces” (coming from $\xi ( r ) )$ are not “balanced” are slightly enhanced. Fig. (2) shows the magnitude of the displacement of the quasiparticle coordinate caused by an electron a distance r away. This is an estimate of the distance that the free-fermion nodal surfaces are displaced by backflow. The strongest effects are observed when two electrons are very close, for distances less than the average nearest neighbor distance which is 2 in the units we have used. We also note that the backflow potential is attractive for $r _ { s } \ge 1 0$ . We expect that the displacement of the quasiparticle coordinates is on the order of 0.01 a. Assuming that this is the case, the released-node calculation with a relatively short projection time should be able to correct the nodal surfaces from those of a free fermion trial function.

The main difficulties in the use of the backflow wave function is that firstly, the implementation is considerably more complex than for the Slater-Jastrow form, and secondly, because update formulas cannot be used to speed up single particle moves, Monte Carlo moves are of updating all particles simultaneously. Details and fuller discussion of the algorithm are given in Ref. [16].

## B. Ground State Energy

We first calculated the ground-state energy of the system with N = 54 electrons at a density range of $1 \leq r _ { s } \leq 2 0$ by both VMC and fixed-node DMC methods. Table II shows the results obtained from the improved trial wave functions in Eq. (5) as well as the Slater-Jastrow wave functions. It can be seen that both VMC and fixed-node calculations with the trial functions including backflow correlation improve significantly the Slater-Jastrow results at all densities considered. However, the three-body correlation is found to have minimal effect for $r _ { s } \le 5$ , which corresponds to typical metallic densities.

Fig. (3) shows the effects of backflow and three-body correlations on the correlation energy that is missed by the Slater-Jastrow wave function both from the variational and the fixed-node calculation. In the following discussion, our best results (the backflow fixednode energies) are assumed to be exact. We will examine this assumption at the end of this section. At high densities of $r _ { s } ~ \le ~ 5$ , the effect due to the three-body correlation is negligible and the backflow effect is dominant. However, as the density decreases, the three-body effect increases while the backflow effect decreases. We can conclude from the trends of Fig. (3) that at the density where Wigner crystallization occurs, estimated to be $r _ { s } \sim 1 0 0$ by Ceperley and Alder [12], the effect in the energy of the three-body term will be much larger than the backflow term. This is consistent with the expectation that backflow correlation is energetically less important as electrons are localized by strong correlation at low densities. Note, however, that the actual effect of the backflow correlation on the wave function decreases with density, as shown in Fig. (2).

The combined effects of both higher-order correlations in the variational wave function account for 60% to 80% of the correlation energy missing in the Slater-Jastrow function. At high densities $( r _ { s } \le 5 )$ , this variational energy is shown to be roughly as good as the Slater-Jastrow fixed-node DMC energy, which captures 70 to 80% of the missing correlation energy throughout our density range.

The backflow and three-body effects in the electron gas discussed above are very similar to the situation in two dimensions. See Fig. (4) of Ref. [16]. The only notable difference is that at the lowest density considered $( r _ { s } = 2 0 )$ , the backflow effect is more important than the static three-body correlation in three dimensions while two correlations have virtually equal importance in two dimensions. This can be understood in terms of the increased importance of correlations in lower dimensions for the same value of $r _ { s } ;$ for example, this is reflected in the fact that Wigner crystalization occurs at smaller $r _ { s }$ in two dimensions than in three dimensions [26].

Fig. (4) shows the correlation energies missing from the Slater-Jastrow wave function and from the three-body and backflow wave function divided by the kinetic energy. Since the kinetic energy operator does not commute with the Hamiltonian, the kinetic energy cannot be computed directly with the distribution $\phi _ { 0 } ( R ) \Psi _ { T } ( R )$ sampled through the diffusion process in Eq. (3). It has been estimated by making an extrapolation between the VMC and the DMC results [13]. It is clear from the figure that the missing correlation energy from both types of trial wave functions becomes a smaller fraction of the kinetic energy at higher densities.

Since our calculation has been done on the system with a finite number of electrons, we extrapolate the energies to the thermodynamic limit to compare with other calculations. We follow the extrapolation scheme based upon the Fermi liquid theory [27,26], which assumes that the energy per particle for a finite system with the periodic boundary condition is related to the bulk energy by

$$
E _ { N } = E _ { \infty } + b _ { 1 } ( r _ { s } ) \Delta T _ { N } + b _ { 2 } ( r _ { s } ) \frac { 1 } { N } .\tag{13}
$$

Here, $E _ { N } \ ( E _ { \infty } )$ is the total energy per electron of the finite (infinite) system and $\Delta T _ { N }$ is the free particle kinetic energy differences between two systems. We determine the parameters $E _ { \infty } , b _ { 1 } .$ and $b _ { 2 }$ by a least-squares fit to VMC calculations with Slater-Jastrow trial functions at different values of $N = 5 4 , 6 6 , 1 1 4 , 1 6 2 , 2 4 6$ . In Table III are shown the energies, fitted parameters, and the $\chi ^ { 2 }$ value of the fit. The reasonable values of $\chi ^ { 2 }$ show that the Fermi liquid theory completely explains the size dependence of the energy to statistical accuracy of the VMC energies over this range of particle numbers. To extract the extrapolated threebody and backflow DMC energy for the infinite system, $E _ { \infty } ^ { 3 B F - D M C }$ , we did the DMC runs only at $N = 5 4$ whose results are shown in Table II and then use the parameters determined from VMC to get $E _ { \infty } ^ { 3 B F - D M C }$ . It is assumed that the size dependences for the VMC (SJ) and the DMC (3BF) results are the same. This assumption needs to be tested in future calculations. The same procedure was successfully applied to assess the finite-size effects in our previous QMC calculation for the two-dimensional electrons [16].

One can see in Table III that our extrapolated backflow fixed-node energies are lower, even if the differences are small, than Ceperley and Alder’s released-node results as well as Ortiz and Ballone’s Slater-Jastrow fixed-node results. Our present results show that the calculations of Ceperley and Alder only got approximately half of the Slater-Jastrow fixednode error with their released-node procedure due to computer limitations at that time. Considering that a fixed-node energy is an upper bound to the true ground-state energy, this validates our assertion that our backflow fixed-node results are accurate.

Since the fixed-node results depend only on the nodal structures of the trial functions used, one can speculate that the nodal locations of the backflow wave function are fairly close to those of the exact ground state. Without more investigation, we cannot quantify this statement, because there is not a simple relationship between nodal locations and fixed-node energy. The accuracy of the backflow nodes was also shown in our previous released-node (transient-estimate) calculation for the two-dimensional electron gas [28].

Although comparison with well-converged exact results is the best method of assessing the accuracy of a fixed-node result for the energy, in the remainder of this section we develop two other methods that require only the VMC and fixed-node DMC energies. Both methods rely on the fact that the errors in the variational energy $E _ { V M C }$ , the variance of the local energy V , and the fixed-node energy $E _ { F N }$ should all be quadratic in the difference between a trial function and the true ground state. Thus, as a trial function is significantly improved in going from a two-body level (Slater-Jastrow) to a three-body level (backflow and threebody), one can estimate the exact energy by the relative improvements of the variational energy relative to the variance and the fixed-node energy.

The variances of the local energy (Eq. (12)) for the various trial wave functions, are given in Table II and plotted in Fig. (5) at $r _ { s } = 1 0$ . As can be seen, the variance decreases roughly proportional to the drop in energy for the four trial functions considered. The dotted line in Fig. (5) represents a linear fit and the triangle our best (backflow) fixed-node energy. There is no fundamental reason why the energy and variance for general trial functions would have a linear relationship. However, in practice this relation is often observed [16]. The observed linear relationship both validates our optimization procedure and provides an independent estimate of the exact energy.

Shown in Table II is our estimate of the error of the computed backflow fixed-node energy obtained from the energies and variances of the best (backflow + three-body) and worst (Slater-Jastrow) trial functions, which is based on the following assumption:

$$
\frac { V ^ { ( k ) } } { E _ { V M C } ^ { ( k ) } - E _ { 0 } } = c o n s t a n t .\tag{14}
$$

$\epsilon _ { V }$ in Table II is the difference between this extrapolation $E _ { 0 }$ and our best fixed-node energy. We extrapolated using only the results from the best and worst trial functions to minimize the extrapolation error. There is the Temple lower bound [29] to the ground-state energy which involves the energy and the variance. However, it is not useful for many-body systems. Because our procedure is not rigorous, there is no guarantee that the estimate will lie below our computed best fixed-node result. In fact at $r _ { s } = 1 0$ the estimate lies above it. Our next extrapolation procedure does not have this problem.

In going from the two-body to three-body level, one can also assume that the nodal positions improve at the same rate as the variational energy so that we can assume:

$$
\frac { E _ { F N } ^ { ( k ) } - E _ { 0 } } { E _ { V M C } ^ { ( k ) } - E _ { 0 } } = c o n s t a n t .\tag{15}
$$

Using this equation with our best and worst energies in both VMC and DMC calculations, we determine $E _ { 0 }$ and hence the error in the backflow fixed-node energy is shown as $\epsilon _ { F N }$ in Table II. Again this procedure has no fundamental validity since it is possible to improve the variational energy without affecting the nodes by improving the bosonic correlations. This estimate shows that considerably larger corrections might be expected from exact calculations, from 0.6mRy at $r _ { s } = 1$ to 0.1mRy at $r _ { s } = 2 0$

It can be seen in Table II that the estimated fixed-node errors $( \epsilon _ { V }$ and $\epsilon _ { F N } )$ are smaller at all densities considered than the energy improvements due to the nodal change from the Slater-Jastrow function to the backflow wave function.

## IV. CONCLUSION

We have studied the correlation energy of the interacting three-dimensional electron gas, using VMC and fixed-node DMC calculations including the three-body and the backflow correlation. The additional correlation energy due to backflow is dominant over the threebody effect in the high density regime but the relative importance of the former decreases as the density is reduced. This is the same trend as was found for the two-dimensional electron gas [16] except that the importance of backflow is more significant in higher dimensions, especially at low densities. This is due to the fact that in two dimensions the effects of interactions are larger than in three dimensions at a given $r _ { s }$ and other effects tend to dominate more over the effects of backflow.

The variational wave function with backflow and three-body correlations is a large improvement over the Slater-Jastrow function. We find that these higher-order correlations account for 60 to 80% of the remaining correlation energy beyond the Slater-Jastrow variational results. Since backflow changes the nodes, the fixed-node DMC results are also significantly improved. The fixed-node method based upon the Slater-Jastrow nodes is found to capture no more than 80% of the remaining correlation energy.

After making a careful finite-size analysis, we have compared our backflow fixed-node energies with Ceperley and Alder’s released-node results. These two independent calculations using different methods are found to give nearly identical results within statistical and systematic errors. From a linear extrapolation to zero variance of the local energy, we find further evidence that our backflow fixed-node results are very close to the true ground-state energy.

For future work, we conclude that one should be able to use the much improved wave functions, better released-node methods [30], with more size-dependence studies and full utilization of current computer hardware to achieve an order of magnitude more accurate results for the energy of the electron gas than was done nearly two decades ago.

## V. ACKNOWLEDGEMENT

This work has been supported by the Korea Science and Engineering Foundation under grant 96-0207-045-2 and through its SRC program, and by the National Science Foundation under grant DMR 94-224-96.

## REFERENCES

[1] D. Pines, Elementary Excitations in Solids (Addison-Wesley, New York, 1963).

[2] W. Kohn and L. J. Sham, Phys. Rev. 140, A1133 (1965).

[3] R. G. Parr and W. Yang, Density-Functional Theory of Atoms and Molecules (Oxford University Press, New York, 1989).

[4] F. Bloch, Z. Phys. 57, 545 (1929).

[5] E. Wigner, Phys. Rev. 46, 1002 (1934); Trans. Faraday Soc. 34, 678 (1938).

[6] A. W. Overhauser, Phys. Rev. Lett. 3, 414 (1959); ibid. 4, 462 (1960); Phys. Rev. 128, 1437 (1962).

[7] K. Moulopoulos and N. W. Ashcroft, Phys. Rev. Lett. 69, 2555(1992).

[8] G. D. Mahan, Many-Particle Physics (Plenum, New York, 1991).

[9] M. Gell-Mann and K. A. Brueckner, Phys. Rev. 106, 364 (1957).

[10] K. S. Singwi, M. P. Tosi, R. H. Land, and A. Sj¨olander, Phys. Rev. 176, 589 (1968).

[11] D. M. Ceperley, Phys. Rev. B 18, 3126 (1978).

[12] D. M. Ceperley and B. J. Alder, Phys. Rev. Lett. 45, 566 (1980).

[13] K. Schmidt and M. H. Kalos, in Applications of the Monte Carlo method in Statistical Physics, edited by K. Binder (Springer-Verlag, Berlin, 1984).

[14] J. B. Anderson, J. Chem. Phys. 63, 1499 (1975); 65, 4121 (1976)

[15] P. J. Reynolds, D. M. Ceperley, B. J. Alder, and J. W. A. Lester, J. Chem. Phys. 77, 5593 (1982).

[16] Y. Kwon, D. M. Ceperley, and R. M. Martin, Phys. Rev. B 48, 12037 (1993).

[17] D. M. Ceperley, in Recent Progress in Many Body Theories ed. J. Zabolitzky, Lecture Notes in Physics, 142, Springer-Verlag (1981).

[18] G. Ortiz and P. Ballone, Europhys. Lett. 23, 7 (1993); Phys. Rev. B 50, 1391 (1994).

[19] S. Moroni, S. Fantoni and G. Senatore, Phys. Rev. B 52 13547 (1995).

[20] D. M. Ceperley, J. Stat. Phys. 63, 1237 (1991).

[21] T. Gaskell, Proc. Phys. Soc. London 77, 1182(1961).

[22] R. P. Feynman and M. Cohen, Phys. Rev. 102, 1189 (1956).

[23] P. Ewald, Ann. Phys. 64, 253 (1921).

[24] R. M. Panoff and J. Carlson, Phys. Rev. Lett. 62, 1130 (1989).

[25] C. J. Umrigar, K. G. Wilson, and J. W. Wilkins, Phys. Rev. Lett. 60, 1719 (1988).

[26] B. Tanatar and D. M. Ceperley, Phys. Rev. B 39, 5005 (1989).

[27] D. M. Ceperley and B. J. Alder, Phys. Rev. B 36, 2092 (1987).

[28] Y. Kwon, D. M. Ceperley, and R. M. Martin, Phys. Rev. B 53, 7376 (1996).

[29] G. Temple, Proc. R. Soc. London, Ser. A 119, 22 (1928).

[30] M. Caffarel and D. M. Ceperley, J. Chem. Phys. 97, 8415 (1992).

## FIGURES

FIG. 1. The three-body contribution to the logarithm of the wave function for a pair of electrons separated by a distance r. Solid line, $r _ { s } = 1$ ; dotted, $r _ { s } = 5 ;$ dot-dashed, $r _ { s } = 1 0 \mathrm { : }$ ; long dashed, $r _ { s } = 2 0$

FIG. 2. The change in quasi-electron coordinate due to a pair of electrons separated by a distance r. Solid line, $r _ { s } = 1$ ; dotted, $r _ { s } = 5 ;$ dot-dashed, $r _ { s } = 1 0 \mathrm { : }$ ; long dashed, $r _ { s } = 2 0$

FIG. 3. Effects of three-body and backflow correlations as a function of the density of the system. The vertical axis shows $\Delta E / \Delta E _ { S J } = ( E - E _ { D M C } ^ { 3 B F } ) / ( E _ { V } ^ { S J } - E _ { D M C } ^ { 3 B F } )$ ,that is, top axis corresponds to the Slater-Jastrow variational energy $E _ { V } ^ { S J }$ and bottom axis to the fixed-node DMC energy $E _ { D M C } ^ { 3 B F }$ , calculated with our best trial function including three-body and backflow correlations. The diamonds show the effect of only three-body correlation, the circles the effect of only backflow and the squares represent the combined effect of both correlations. Finally, the filled triangles show the result using the fixed-node DMC method with free-fermion nodes of the Slater-Jastrow function.

FIG. 4. The energy missing from the Slater-Jastrow wave function (⋄) and from the three-body and backflow wave function (•) divided by the kinetic energy as a function of the density parameter $r _ { s } .$ The vertical axis shows $( E _ { V } - E _ { D M C } ^ { 3 B F } ) / < T >$

FIG. 5. Variational energy versus the variance of local energy at $r _ { s } = 1 0$ . Each point • represents one variational calculation: from higher to lower energies, the Slater-Jastrow, three-body, backflow, and $( \mathrm { b a c k f l o w ~ + ~ t h r e e { - } b o d y } )$ results. The filled triangle represents our backflow fixed-node result and the dotted line shows a linear fit through • points. The statistical errors of the data are smaller than the sizes of the symbols.

TABLES  
TABLE I. Optimized variational parameters of three-body and backflow correlation functions for $N = 5 4$
<table><tr><td> $r _ { s }$ </td><td> $\lambda _ { B }$ </td><td> $s _ { B }$ </td><td> $r _ { B }$ </td><td> $w _ { B }$ </td><td> $\lambda _ { T }$ </td><td> $r _ { T }$ </td><td> $w _ { T }$ </td></tr><tr><td>1.0</td><td>0.025</td><td>0.395</td><td>0.210</td><td>0.689</td><td>0.006</td><td>0.293</td><td>0.949</td></tr><tr><td>5.0</td><td>0.105</td><td>0.158</td><td>0.180</td><td>0.670</td><td>-0.060</td><td>0.286</td><td>1.176</td></tr><tr><td>10.0</td><td>0.959</td><td>-0.672</td><td>0.247</td><td>3.788</td><td>-0.258</td><td>0.257</td><td>0.918</td></tr><tr><td>20.0</td><td>1.249</td><td>-0.938</td><td>0.275</td><td>3.787</td><td>-0.255</td><td>0.252</td><td>0.911</td></tr></table>

TABLE II. VMC and fixed-node (FN) DMC energies E and the variances of the local energy V obtained with various trial wave functions for $N = 5 4$ (SJ: the Slater-Jastrow function, 3BD: three-body correlation, BF: backflow correlation). The energies are in units of $R y$ per electron and the variances in units of $r _ { s } ^ { 4 } ( R y / \mathrm { e l e c t r o n } ) ^ { 2 }$ . Also shown are our estimations of the fixed-node error ǫ in the backflow fixed-node DMC calculation.
<table><tr><td></td><td> $r _ { s } = 1 . 0$ </td><td> $r _ { s } = 5 . 0 \AA$ </td><td> $r _ { s } = 1 0 . 0$ </td><td> $r _ { s } = 2 0 . 0$ </td></tr><tr><td> $\overline { { E _ { V M C } ^ { S J } } }$ </td><td>1.0669(6)</td><td>-0.15558(7)</td><td>-0.10745(2)</td><td>-0.06333(1)</td></tr><tr><td> $_ { \mathrm { \it { F } } } { \it { S } } \dot { \cal { J } } + \dot { 3 } \breve { B } { \cal { D } }$  EV M CI</td><td>1.0663(5)</td><td>-0.15569(5)</td><td>-0.10773(2)</td><td>-0.06348(1)</td></tr><tr><td> $_ { F } { \dot { S } } { \dot { J } } { \dot { + } } B F$  EV M C</td><td>1.0617(4)</td><td>-0.15729(5)</td><td>-0.10829(2)</td><td>-0.06365(1)</td></tr><tr><td> $E _ { V M C } ^ { S J + 3 B D + B F }$ </td><td>1.0613(4)</td><td>-0.15735(5)</td><td>-0.10835(2)</td><td>-0.06378(2)</td></tr><tr><td> $E _ { F N } ^ { S J }$ </td><td>1.0619(4)</td><td>-0.15734(3)</td><td>-0.10849(2)</td><td>-0.06388(1)</td></tr><tr><td> $E _ { F N } ^ { S J + 3 B D + B F }$ </td><td>1.0601(2)</td><td>-0.15798(4)</td><td>-0.10882(2)</td><td>-0.06403(1)</td></tr><tr><td> $V ^ { S J }$ </td><td></td><td></td><td></td><td></td></tr><tr><td> $V ^ { S J + 3 B D }$ </td><td>0.0213(4) 0.0205(4)</td><td>0.0266(3) 0.0229(4)</td><td>0.074(1)</td><td>0.189(3)</td></tr><tr><td> $V ^ { S J + B F }$ </td><td></td><td></td><td>0.054(2)</td><td>0.144(3)</td></tr><tr><td> $V ^ { S J + 3 B D + B F }$ </td><td>0.0054(3) 0.0053(2)</td><td>0.0069(2) 0.0066(2)</td><td>0.027(1) 0.026(1)</td><td>0.111(2)</td></tr><tr><td></td><td></td><td></td><td></td><td>0.079(2)</td></tr><tr><td> $\epsilon _ { V }$ </td><td>0.0007(6)</td><td>-0.00005(8)</td><td>0.00002(5)</td><td>0.00007(4)</td></tr><tr><td> $\epsilon _ { F N }$ </td><td>0.0006(4)</td><td>0.00036(8)</td><td>0.00027(5)</td><td>0.00013(3)</td></tr></table>

TABLE III. Size dependence in the Slater-Jastrow VMC method of normal electron liquid at $1 \leq r _ { s } \leq 2 0$ and $\chi ^ { 2 } { \cdot } \mathrm { f i t }$ parameters. Also shown are the extrapolated DMC energies at an infinite system $( E _ { \infty } ^ { S J - D M C }$ and $E _ { \infty } ^ { 3 B F - D M C } )$ , Ceperley and Alder’s released-node result $\left( \mathrm { C A } ^ { * } \right)$ , and Ortiz and Ballone’s Slater-Jastrow fixed-node result $\mathrm { ( O B ^ { * * } ) }$
<table><tr><td></td><td> $r _ { s } = 1 . 0$ </td><td> $r _ { s } = 5 . 0 \AA$ </td><td> $r _ { s } = 1 0 . 0$ </td><td> $r _ { s } = 2 0 . 0$ </td></tr><tr><td></td><td> $\overline { { N = 5 4 } }$  1.0669(6)</td><td>-0.15558(7)</td><td>-0.10745(2)</td><td>-0.06333(1)</td></tr><tr><td> $N = 6 6$ </td><td>1.1496(5)</td><td>-0.15166(4)</td><td>-0.10637(2)</td><td>-0.06303(1)</td></tr><tr><td> $N = 1 1 4$ </td><td>1.2079(5)</td><td>-0.14867(3)</td><td>-0.10552(2)</td><td>-0.06278(1)</td></tr><tr><td> $N = 1 6 2$ </td><td>1.1162(4)</td><td>-0.15238(3)</td><td>-0.10642(1)</td><td>-0.06270(1)</td></tr><tr><td> $N = 2 4 6$ </td><td>1.1938(3)</td><td>-0.14886(3)</td><td>-0.10548(1)</td><td>-0.06275(1)</td></tr><tr><td> $E _ { \infty } ^ { S J - V M C }$ </td><td>1.1795(4)</td><td>-0.14914(3)</td><td>-0.10549(2)</td><td>-0.06273(1)</td></tr><tr><td> $b _ { 1 } ( r _ { s } )$ </td><td>1.096(6)</td><td>1.18(1)</td><td>1.21(2)</td><td>1.22(3)</td></tr><tr><td> $b _ { 2 } ( r _ { s } )$ </td><td>-1.16(5)</td><td>-0.134(4)</td><td>-0.051(2)</td><td>-0.0181(7)</td></tr><tr><td> $\chi ^ { 2 }$ </td><td>1.20</td><td>1.29</td><td>2.22</td><td>2.26</td></tr><tr><td></td><td></td><td></td><td></td><td></td></tr><tr><td> $E _ { \infty } ^ { S J - D M C }$ </td><td>1.1744(4)</td><td>-0.15094(4)</td><td>-0.10654(2)</td><td>-0.06329(1)</td></tr><tr><td> $E _ { \infty } ^ { 3 \bar { B } F - D M C }$ </td><td>1.1726(2)</td><td>-0.15158(5)</td><td>-0.10687(2)</td><td>-0.06344(1)</td></tr><tr><td> $\mathrm { C A } ^ { * }$ </td><td>1.174(1)</td><td>-0.1512(1)</td><td>-0.10675(5)</td><td>-0.06329(3)</td></tr><tr><td> $\mathrm { O B ^ { * * } }$ </td><td>1.181(1)</td><td>-0.1514(3)</td><td>-</td><td></td></tr></table>

<!-- image-->

<!-- image-->

<!-- image-->

<!-- image-->

<!-- image-->