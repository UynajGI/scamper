NATIONAL   
RESOURCE   
FOR COMPUTATION   
IN CHEMISTRY

λ i THE GROUND STATE OF THE ELECTRON GAS BY A STOCHASTIC METHOD

MASTER

fi - D.M. CeperTey and B. J. Alder

May 1980

LAWRENCE BERKELEY LABORATORY 15 UNIVERSITY OF CALIFORNIA

OISTRIBUTION OF THIS DOCUMENT IS UNLIMITER

THE GROUND STATE. OF THE ELECTRON GAS \* BY A STOCHASTIC METHOD {4

D. M. Ceperley

National Resource for Computation In Chemistry Lawrence Berkeley Laboratory University of California Berkeley, CA 94720 .

B. J. Alder

Lawrence Livermore Laboratory University of California Livermore; CA 94550

May 1980

ABSTRAct: We have used an exact stochastic simulation of the Scbroedinger equation for charged Bosons and Fermions to calculate the correlation energies, to locate the transitions to their respective crystal phases at zero temperature within l0%, and to establish the stability at intermediate densities of a ferromagnetic fluid of electrons.

## ,Disclamer

This bo  U Necher the Uni Sutas Governmond ror y agancy thero ne y  their mploy mkss an werranty, eapres or implied, or gumes any lngal liability or recorabjlity tor the saatacy. comodlesaneas. or usetuines of any Intormenion, spargtus, product, or procie dactomd, or esents t nti pidReleherin t c crc podu pc cey  eadeerk lctu bom not negamarty censtituse or imoly hs endoramant, recommendation, ar favariny by the Unitad Su Gotr yha Th npinon utohe not neury sate or re thom ol the Unitd Sam Govenmnt or wny agincy phero.

\*This research was supported in part by the Office of Basic Energy Sciences of the U. S. Department of Energy under Concract No. W-7405-ENG-48 and by the National Science Foundation umder Grant No. CHE-7721305.

DTRIBUTION FT CUMENT I TE

The properties of the ground state of the electron gas, also

referred to as the Fermion one component plasma and jellium, have

rigorously only been established in the limit of high densitiesk where

the system approaches a perfect gas and at low density' where the

electrons crystallize. Furthermore, Hartree-Fock calcula ions', and

variational calculations" suggest that at intermediate dersities, the

spin aligned state of the electrons will be more stable than the normal,

unpolarized state. Precise calculations of this many-Fermion system are

required to establish the regions of stability of the various phases

because of the small energy differences among them. This note outlines

a Monte Carlo method, that if run long enough on a computer, can give as

precise a solution for the grourd state of a given Fermion system as

desired.

In practice, the precision of such a calculation is limited to about

two orders of magnitude smaller than that of an approximate trial wave

function that is introduced as an importance function in the Monte Carlo

process. That the introduction of such an importance function is

essential, was previously demonstrated for the many-Boson problem.2 The

extension of this Boson calculation to Fermions requires dealing with

the probability density of a random walk cannot be chosen everywhere

positive, and unless prevented the random walk will always converge to

the all positive, Boson ground state. It is demonstrated here, for the

electron system, that before the effect of this inherent instability

becomes serious, it is possible to extract the properties of the lowest

antisymmetric state. A more general procedure which removes the effects

of the instability has yet to be perfected.

The solution of the Fermion problem was carried out in two steps.

In the first step the nodes, the places where the trial function

vanishes, act as fixed absorbing barriers to the diffusion process.

Inside a connected nodal region the wavefunction is everywhere positive

and vanishes at the boundaries. With these boundary conditions, the

Fermion problem is equivalent to a Boson problem. The energy calculated

is an upper bound to the exact Fermion ground state energy and generally

very close to it. In principle one could next vary the nodal locations

to obtain the best upper bound, by for example, varying the functions

used as elements in the Slater determinant of the trial wave function.

In practice the highly dimensional nocal surfaces are difficult to

parameterize in a systematic fashion.

The second step, called 'nodal relaxation', begins with the

population of walks from the 'fixed-node' approximation. In this second

procedure, if a random walk strays across the node of the trial function

it is not terminated, but the sign of its contribution to any average is

reversed. At any stage of the random walk there is a population of

positive walks (those that remained in the same nodal region or crossed

an even number of nodes) and a population of negative walks (those that

crossed an odd number of nodes). The importance function used in this process is the absolute value of the trial function. It can be easily shown that the difference population converges to the antisymmetric eigenfunction. However both the positive and negative populations grow geometrically with a rate equal to the difference between the Fermi and Bose energies. If the relaxation time from the fixed-node distribution times this energy difference is less than unity, the Fermion energy can be reliably extracted. We have found that for the electron gas this condition is satisfied if the nodes of the Hartree-Fock wavefunction are used.

Our simulation method is a simpler, though approximate, version of the Greens function Monte Carlo method of Kalos et. ${ \mathsf { a l . } } ^ { 5 }$ A trial wavefunction $\Psi _ { \mathbf { T } } ^ { \mathbf { \alpha } } ( \mathbb { R } )$ of the Bijl-Jastrow-Slater $\pm \tt y p e ^ { 4 }$ and an ensemble of about l0o systems are selected from a variational Monte Carlo calculation, where R represents the 3N spatial coordinates of the systems of N electrons. Let the probability density of finding a random walk in $\tt R d R ^ { 3 N }$ at time t be given by $\mathtt { f } ( \mathtt { R } , \mathtt { t } ) \mathtt { d } \mathtt { R } ^ { 3 \mathtt { N } }$ . Then the value of f at t=0 is given by $\big | \Psi _ { \mathbf { T } } ( { \boldsymbol { \mathsf { R } } } ) \big | ^ { 2 }$ properly normalized. The diffusion equation for f(R,t) is:

$$
\frac { \partial \mathbf { f } } { \partial \mathbf { t } } = \frac { \hbar ^ { 2 } } { 2 m } [ \sum _ { \mathbf { i } = 1 } ^ { \mathbf { N } } \nabla _ { \mathbf { i } } ^ { 2 } \mathbf { f } -  \vec { \nabla } _ { \mathbf { i } } ( \mathbf { f } \vec { \nabla } _ { \mathbf { i } } \mathbf { 1 } \cdot \vec { \mathbf { l } } \cdot \vec { \mathbf { l } } ) \mathbf { \nabla } _ { \mathbf { T } } | ^ { 2 } ] - [ \frac { \hbar \Psi _ { \mathbf { i } } } { \Psi _ { \mathbf { T } } } -  \mathbf { \nabla } _ { \mathbf { E } } \mathbf { _ { r e f } } ] \textbf { f }\tag{1}
$$

where H is the Hamiltonian

$$
\textrm { \textbf { H } } = \frac { \pi ^ { 2 } } { 2 \pi } \sum _ { 1 = 1 } ^ { \textrm { N } } \nabla _ { 1 } ^ { 2 } - \sum _ { 1 < 1 } e ^ { 2 } / \mathbf { r _ { \textrm { i j } } }\tag{2}
$$

It is easily verified that for large times, $\begin{array} { r l r } { \mathbf { f } ( \mathbf { r } , \mathbf { t } ) } & { = } & { \Psi _ { \mathbf { T } } \boldsymbol { \dot { \Phi } } _ { 0 } \exp ( - \mathbf { t } ( \mathbf { E } _ { \mathbf { r e f } } - \mathbf { E } _ { o } ) ) } \end{array}$

where $\Xi _ { \circ }$ and $\$ 0$ are the exact eigenvalue and eigenfunction. The

above equation for f(R,t) has a simple interpretation as a stochastic

process. Each member of the ensemble of systems undergoes i) random

diffusion caused by the zero point motion, ii) biasing or drift by the

trial quantum force, $\nabla 1 _ { \overline { { \mathbf { n } } } } | \Psi _ { \overline { { \mathbf { T } } } } | ^ { 2 }$ , and iii) branching with probability

given by the difference between the local trlal energy, $\Xi _ { \overline { { \mathbf { \updownarrow } } } } \approx \mathbb { \updownarrow } \ \Psi _ { \mathbf { \ T } } / \ \Psi _ { \mathbf { \updownarrow } }$

and the arbitrarily chosen refernce energy, $\mathbf { E } _ { \mathbf { r e f } } .$ , By "branching", it

is meant that a particular system is either eliminated from the ensemble

(if the local enargy is greater than the reference energy) or duplicated

in the ensemble (otherwise). A steady state population of the ensemble

requires that the reference energy equal the lowest eigenvalue. This is

one way of determining the eigenvalue.

The trial wavefunction employed in the present calculations are

identical with those used in an earlier Monte Carlo variational

calculation.4 This trial function is a product of two-body

correlation factors times a Slater determinant of single particle

orbitals. The two body correlation factors are chosen such that they

remove exactly the singularities in the local energy when two electrons

approach each other, thus reducing tremendously the variance of the

estimate of the ground state energy. For the fluid phase the single

particle orbitals are plane waves, with the wave vector lying within the

Fermi sea. For the polarized state, where there is only one spin for

each spatial state, as opposed to two for the normal unpolarized state, the Fermi wavevector has been increased to allow. for twice as many spatial orbitals. In the crystal phase, the orbitals are Gaussians centered around body centered cubic lattice sites with a width chosen variationally.

Fig. 1 shcws that the relaxation from the unpolarized nodes to the

ground state is rapid with a small lowering of the energy. A less

accurate trial wave function with different nodes obtained from a linear

combination of polarized and unpolarized Slater determinants is

nevertheless shown to lead to similar energies with a somewhat larger

relaxation time. This shows the insensitivity of results to the

original location of the nodes. Since at all densities the relaxation

from the Hartree-Fock nodes was rapid, the ground state energy of the

electron gas by the method employed could be obtained with very little

uncertainty.

The largest uncertainty in the results is in fact due to the number

dependence. Due to the high accuracy of the results derived from

employing a good trial wave function and the consequent small

established for systems ranging from 38 to 246 particles is an order of

magnitude larger than the statistical error. Extrapolation to infinite

particle results was carried out at each density on the basis of ${ \mathfrak { E } } ( { \mathfrak { N } } ) \ =$

$$
\begin{array} { r l r l r l r l r l } { \bf { E } _ { 0 } } & { { } + } & { } & { \bf { E } _ { 1 } / \tt N } & { { } + } & { \bf { E } _ { 2 } } & { } & { \Delta _ { \tt N } } \end{array}
$$

$$
\begin{array} { r l } { { \tt E } _ { \tt o } , } & { { } { \tt E } _ { 1 } } \end{array}
$$

$$
\mathtt { E } _ { 2 }
$$

were empirically determined from the simulations. The $\mathtt { E } _ { 1 }$ term arises

from the potential energy and is due to the correlation between a particle and its images in the periodically extended space that is used in the Ewald summation $\mathtt { p r o c e d u r e } ^ { 4 }$ to eliminate the major surface effects. The $\mathtt { E } _ { 2 }$ tarm comes from the discrate nature of the Fermi sea for finite systems, and $\spadesuit$ is the size dependence of an ideal Fermi system at the same density. That term is absent for Bosons. In addition the energies have been extrapoiated to zero time step by empirically establishing the validity of linear extrapolation. This correction is quite small, on the order of the statistical error for the time steps used. However this correction can be completely avoided by using an integral formulation of eq.(1).5

The results for the energy of the plasma in four different phases is given in Table I. These energies multiplied by $\tau _ { \mathsf { s } } ^ { 2 }$ are plotted in Fig. 2 relative to the lowest Boson state. Multiplying by $\tau _ { \mathrm { s } } ^ { 2 }$ corresponds to holding the density fixed and increasing the charge. Plotted in this manner the minute differences in energy at low density can be more clearly seen. The Boson system undergoes $\mathtt { w i g n e r } ^ { 6 }$ crystallization at $\mathrm { ~ \bf ~ r ~ } _ { \mathrm { ~ s ~ } } = 1 6 0 \mathrm { ~ \bf ~ \pm ~ } 1 0$ , The Fermion system has two phase transitions, crystallization at $\mathrm { ~ \bf ~ r ~ } _ { \mathrm { ~ \bf ~ s ~ } } = 1 0 0 \mathrm { ~ \bf ~ \pm ~ } 2 0$ and depolarization at ${ \bf { r } } _ { \textbf { s } } = { \bf { \sigma } } _ { 7 5 } \pm { \bf { \sigma } } 5$ . The difference in energy between a Boson crystal and a Fermion crystal is less than $1 . 0 ~ \mathrm { ~ x ~ } ~ 1 0 ^ { - 6 } \mathrm { ~ R ~ }$ at $\mathrm { ~ \bf ~ r ~ } _ { \mathrm { ~ \bf ~ s ~ } } = \mathrm { ~ \bf ~ l 0 ^ { \wedge ~ } ~ }$ . The energies of the three Fermion states are sufficiently close in the low density regime that still more accurate calculations on larger systems would be desirable to confirm these results.

In Table II the correlation energy for the unpolarized Fermi fluid, that is the ground state energy relative to the Hartree-Fock energy, is compared to that of several other theories in the metallic density range. The correlation energies are very similar for all methods. The coupled-cluster7,8 formalism give the most accurate results. It is seen that a variational integral equation theory, the Fermi hypernetted chain, gives energies below the present results, indicating that the approximations employed have compromised the variational principle.

Finally, Table III displays the differences between the pair product variational results, the fixed-node results and the final energies. Although the Bijl-Jastrow-Slater results are quite accurate, the error is different for the different phases, changing their relative stability. This demonstrates how essential it is to perform exact simulations to reliably calculate phase transitions densitites.

The authors would like to thank M. H. Kalos for numerous useful discussions and for inspiring the present work. We thank Mary Ann Mansigh for computational assistance.

## References

1M. Gell-Mann and K.A. Brueckner, Phys. Rev. 106, 364 (1957).

2W. J. Carr Phys. Rev. 122, 1437 (1961).

3F. Bloch, Z. Phys. 57, 549 (1929).

4D, Ceperley, Phys. Rev. B18, 3126 (1978).

5M.H. Kalos, D. Levesque, L. Verlet, Phys. Rev. A9, 2178 (1974). D. M. Ceperley and M. H. KAlos "Monte Carlo Methods in Statistical Physics" pg. 145, Springer-Verlag (1979).

E. P. Wigner, Phys. Rev. 46, 1002 (1934). Trans. FaraDay Soc. 34, 678 (1938).

7D. L. Freeman, Phys. Rev. Bl5, 5513 (1977).

8p. F. Bishop, K. H. Luhrman, Phys. Rev. B (in press).

9L. J. Lantto, P. J. Simens, NuCl. Phys. A317, 55 (1979); L. J. Lantto, Nucl. Phys. A (in press).

10p. Vashishta and K. S. Singwi, Phys. Rev. B6, 875 (1972).

<table><tr><td rowspan=1 colspan=1>rs</td><td rowspan=1 colspan=1>EpMF</td><td rowspan=1 colspan=1>EFMF</td><td rowspan=1 colspan=1>EBF</td><td rowspan=1 colspan=1>EBCC</td></tr><tr><td rowspan=1 colspan=1>1.0</td><td rowspan=1 colspan=1>1.174(1)</td><td rowspan=1 colspan=1></td><td rowspan=1 colspan=1></td><td rowspan=1 colspan=1></td></tr><tr><td rowspan=1 colspan=1>2.0</td><td rowspan=1 colspan=1>0.0041(4)</td><td rowspan=1 colspan=1>0.2517(6)</td><td rowspan=1 colspan=1>-0.4531(1)</td><td rowspan=1 colspan=1></td></tr><tr><td rowspan=1 colspan=1>5.0</td><td rowspan=1 colspan=1>-0.1512(1)</td><td rowspan=1 colspan=1>-0.1214(2)</td><td rowspan=1 colspan=1>-0.21663(6)</td><td rowspan=1 colspan=1></td></tr><tr><td rowspan=1 colspan=1>10.0</td><td rowspan=1 colspan=1>-0.10675(5)</td><td rowspan=1 colspan=1>-0.1013(1)</td><td rowspan=1 colspan=1>-0.12150(3)</td><td rowspan=1 colspan=1></td></tr><tr><td rowspan=1 colspan=1>20.0</td><td rowspan=1 colspan=1>-0.06329(3)</td><td rowspan=1 colspan=1>-0.06251(3)</td><td rowspan=1 colspan=1>-0.06666(2)</td><td rowspan=1 colspan=1></td></tr><tr><td rowspan=1 colspan=1>50.0</td><td rowspan=1 colspan=1>-0.02884(1)</td><td rowspan=1 colspan=1>-0.02878(2)</td><td rowspan=1 colspan=1>-0.02927(1)</td><td rowspan=1 colspan=1>-0.02876(1)</td></tr><tr><td rowspan=1 colspan=1>100.0</td><td rowspan=1 colspan=1>-0.015321(5)</td><td rowspan=1 colspan=1>-0.015340(5)</td><td rowspan=1 colspan=1>-0.015427(4)</td><td rowspan=1 colspan=1>-0.015339(3)</td></tr><tr><td rowspan=1 colspan=1>130.0</td><td rowspan=1 colspan=1></td><td rowspan=1 colspan=1></td><td rowspan=1 colspan=1>-0.012072(4)</td><td rowspan=1 colspan=1>-0.012037(2)</td></tr><tr><td rowspan=1 colspan=1>200.0</td><td rowspan=1 colspan=1></td><td rowspan=1 colspan=1></td><td rowspan=1 colspan=1>-0.008007(3)</td><td rowspan=1 colspan=1>-0.008035(1)</td></tr></table>

The ground state energy of the charged Fermi and Bose systems. The density parameter, rs, is the Wigner sphere radius in units of Bohr radii. The energies are Rydbergs and the digits in parenthesis represent the error bar in the last decimal place. The four phases are: paramagnetic or unpolarzed Fermi fluid (PMF); the ferromagnetic or polarized Fermi fluid (FMF); the Bose fluid (BF); and the Bose crystal with a BCC lattice.

<table><tr><td rowspan=1 colspan=1>rs</td><td rowspan=1 colspan=1>EMC</td><td rowspan=1 colspan=1>εcCl</td><td rowspan=1 colspan=1>εCC2</td><td rowspan=1 colspan=1>εDE</td><td rowspan=1 colspan=1>EFHNC</td></tr><tr><td rowspan=1 colspan=1>1.0</td><td rowspan=1 colspan=1>0.121(1)</td><td rowspan=1 colspan=1>0.118</td><td rowspan=1 colspan=1>0.123</td><td rowspan=1 colspan=1>0.112</td><td rowspan=1 colspan=1>0.138</td></tr><tr><td rowspan=1 colspan=1>2.0</td><td rowspan=1 colspan=1>0.0902(4)</td><td rowspan=1 colspan=1>0.0884</td><td rowspan=1 colspan=1>0.0917</td><td rowspan=1 colspan=1>0.089</td><td rowspan=1 colspan=1>0.098</td></tr><tr><td rowspan=1 colspan=1>5.0</td><td rowspan=1 colspan=1>0.0563(1)</td><td rowspan=1 colspan=1>0.0567</td><td rowspan=1 colspan=1>0.0568</td><td rowspan=1 colspan=1>0.058</td><td rowspan=1 colspan=1>0.058</td></tr><tr><td rowspan=1 colspan=1>10 0</td><td rowspan=1 colspan=1>0.03722(5)</td><td rowspan=1 colspan=1>0.03888</td><td rowspan=1 colspan=1></td><td rowspan=1 colspan=1></td><td rowspan=1 colspan=1>0.037</td></tr></table>

Caption Comparison of the correlation energy with other theories. EMC is the correlation energy from this calculation with the parenthesis representing the error bar in the last decimal place. ECCl and εCC2 8 are the first and second order of the coupled cluster or (eg) theory. EE is the correlation energy in the dielectric formulation and EFHNC is the Fermi-hypernetted chain correlation energy.

<table><tr><td rowspan=1 colspan=1> $\pmb { r _ { s } }$ </td><td rowspan=1 colspan=2> $\delta _ { \overline { { \tt P } } \overline { { \tt E } } }$         $\Upsilon _ { \mathbf { p u p } }$ </td><td rowspan=1 colspan=1> $\delta _ { \overline { { \mathbf { F M F } } } }$ </td><td rowspan=1 colspan=1> $\Upsilon _ { \underline { { \mathbf { F M E } } } }$ </td><td rowspan=1 colspan=1> $\delta _ { \tt B F }$ </td><td rowspan=1 colspan=1> $\delta _ { \mathtt { B C C } }$ </td></tr><tr><td rowspan=1 colspan=1>2</td><td rowspan=1 colspan=1>40</td><td rowspan=1 colspan=1>6</td><td rowspan=1 colspan=1>11.0</td><td rowspan=1 colspan=1>--</td><td rowspan=1 colspan=1>12.0</td><td rowspan=1 colspan=1></td></tr><tr><td rowspan=1 colspan=1>5</td><td rowspan=1 colspan=1>17</td><td rowspan=1 colspan=1>2</td><td rowspan=1 colspan=1>7.2</td><td rowspan=1 colspan=1>--</td><td rowspan=1 colspan=1>6.8</td><td rowspan=1 colspan=1></td></tr><tr><td rowspan=1 colspan=1>10</td><td rowspan=1 colspan=1>11</td><td rowspan=1 colspan=1>→</td><td rowspan=1 colspan=1>6.5</td><td rowspan=1 colspan=1>1.8</td><td rowspan=1 colspan=1>5.1</td><td rowspan=1 colspan=1></td></tr><tr><td rowspan=1 colspan=1>20</td><td rowspan=1 colspan=1>6.7</td><td rowspan=1 colspan=1>0.7</td><td rowspan=1 colspan=1>3.0</td><td rowspan=1 colspan=1>1.0</td><td rowspan=1 colspan=1>3.3</td><td rowspan=1 colspan=1></td></tr><tr><td rowspan=1 colspan=1>50</td><td rowspan=1 colspan=1>2.9</td><td rowspan=1 colspan=1>0.31</td><td rowspan=1 colspan=1>1.6</td><td rowspan=1 colspan=1>0.25</td><td rowspan=1 colspan=1>1.7</td><td rowspan=1 colspan=1>2.0</td></tr><tr><td rowspan=1 colspan=1>100</td><td rowspan=1 colspan=1>1.7</td><td rowspan=1 colspan=1></td><td rowspan=1 colspan=1>1.2</td><td rowspan=1 colspan=1></td><td rowspan=1 colspan=1>1.2</td><td rowspan=1 colspan=1>0.41</td></tr><tr><td rowspan=1 colspan=1>130</td><td rowspan=1 colspan=1>-*•-</td><td rowspan=1 colspan=1></td><td rowspan=1 colspan=1>--</td><td rowspan=1 colspan=1></td><td rowspan=1 colspan=1>1.1</td><td rowspan=1 colspan=1>0.30</td></tr></table>

The error in the variational approximation in $1 0 ^ { - 4 }$ Rydbergs for four different phases. $\delta = E _ { v } - \mathrm { ~ \bf ~ E } _ { 0 }$ (the difference between the Jastrow trial function and the exact ground state energy). $\begin{array} { r l } { \gamma = { } } & { { } \tt E _ { F N } - \tt E _ { 0 } } \end{array}$ (the difference between the 'fixed-node' energy with plane wave nodes and the exact ground state energy).

<!-- image-->  
FIG. l. The energy in Rydbergs per particle of a 38 electron system at the density ${ \pmb { \tau } } _ { \textbf { s } } = \ { \bf \epsilon } _ { 1 0 }$ versus diffusion time (in inverse Rydbergs) from removal of the fixed-nodes. The lower curve is the relaxation of an ensemble of $1 . 6 ~ \mathrm { ~ \bf ~ x ~ } ~ 1 0 ^ { 4 }$ systems from the nodes of the unpolarized determinant of plane waves. The upper curve is the relaxation of l.0 x ${ \mathsf { 1 0 } } ^ { 5 }$ systems from the nodes of a linear combination of polarized and unpolarized determinants.

<!-- image-->  
FIG. 2. The energy of the four phases studied relative to that of the lowest Boson state times $\mathbf { r } _ { : } ^ { 2 }$ in Rydbergs versus $\mathbf { r _ { s } }$ in Bohr radii. Below $\mathrm { ~ \pmb { r } ~ } _ { \mathsf { s } } = \mathrm { ~ \pmb { \imath } 6 0 ~ }$ the Bose fluid is the most stable phase, while above, the Wigner crystal is most stable. The energies of the polarized and unpolarized Fermi fluid are seen to intersect of $\yen 123$ The polarized (ferromagnetic) Fermi fluid will be stable between $\yen 123$ and $\mathbf { r _ { s } } = 1 0 0$