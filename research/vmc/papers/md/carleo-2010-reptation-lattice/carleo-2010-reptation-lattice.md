# Reptation quantum Monte Carlo for lattice Hamiltonians with a directed-update scheme

Giuseppe Carleo, Federico Becca, Saverio Moroni, and Stefano Baroni

SISSA – Scuola Internazionale Superiore di Studi Avanzati and

DEMOCRITOS National Simulation Center,

Istituto Officina dei Materiali del CNR

Via Bonomea 265, I-34136, Trieste, Italy

(Dated: November 12, 2018)

We provide an extension to lattice systems of the reptation quantum Monte Carlo (RQMC) algorithm, originally devised for continuous Hamiltonians. For systems affected by the sign problem, a method to systematically improve upon the so-called fixed-node approximation is also proposed. The generality of the method, which also takes advantage of a canonical worm algorithm scheme to measure off-diagonal observables, makes it applicable to a vast variety of quantum systems and eases the study of their ground-state and excited-states properties. As a case study, we investigate the quantum dynamics of the one-dimensional Heisenberg model and we provide accurate estimates of the ground-state energy of the two-dimensional fermionic Hubbard model.

## I. INTRODUCTION

The path-integral formulation of quantum mechanics is the foundation of many numerical methods that allow one to study with great accuracy the rich physics of interacting quantum systems. At finite temperature, a pathintegral Monte Carlo (PIMC) technique for continuous systems has been developed and applied by Ceperley and Pollock.1,2 Recently, this approach has been renovated in a new class of methods known as worm algorithms.3,4 The zero-temperature counterparts of the PIMC algorithm are the reptation quantum Monte Carlo (RQMC)5 and the path-integral ground state (PIGS) methods,6 which have been demonstrated useful to simulate coupled electron-ion systems,7 as well as to infer spectral properties from imaginary time dynamics.8 A number of important physical problems –particularly in the fields of strongly correlated fermions and cold atoms– can be fruitfully modeled by lattice Hamiltonians. A first application of path-integral techniques to (boson) lattice models was proposed by Krauth et al. in 1991.9 Few other attempts to apply PIMC to lattice models have been made ever since, with a recent application of the RQMC idea to the quantum dimer model Hamiltonian.10 In this paper, we propose a new method that generalizes and improves the approach of Ref. 10 in several ways. Our method is based on continuous-time random walks and is therefore unaffected by time-step errors. Inspired by the work of Syljuasen and Sandvik11 and Rousseau,12 we adopt a generalization of the bounce algorithm of Pierleoni and Ceperley,7 called directed updates, which helps reducing the correlation time in path sampling. We also introduce a worm-algorithm based method to calculate pure expectation values of arbitrary off-diagonal observables, which are generally out of the scope of existing lattice ground-state methods.

The resulting algorithm naturally applies to fermions, using the fixed-node approximation. A technique to improve systematically upon this approximation is proposed, based on the calculation of a few moments of the

Hamiltonian. Our methodology is demonstrated by a few case studies on the one-dimensional Heisenberg and the two-dimensional fermion Hubbard models.

This paper is organized as follow: in Sec. II we present the general formalism of ground-state PIMC for lattice models; in Sec. III our implementation of the RQMC algorithm on a lattice is presented. In particular, we give a detailed account of the above mentioned directed update technique (Sec. III A) and of the continuous-time propagator (Sec. III B); in Sec. III C, we introduce an extension of the algorithm to cope with off-diagonal observables, while in Sec. III D a further extension to systems affected by sign problems is presented, including a strategy to improve systematically upon the fixed-node approximation. Sec. IV contains a few case applications, including the simulation of the spectral properties and spin correlations of the one-dimensional Heisenberg model and the calculation of the ground-state energies of the fermionic Hubbard model with a significantly better accuracy than that achieved by the fixed-node approximation. In Sec. V we finally draw our conclusions.

## II. GENERAL FORMALISM

Let us consider a generic lattice Hamiltonian $\hat { H }$ and a complete and orthogonal basis set, whose states are denoted by x . Given the generic wave function Ψ , its amplitude on the configuration |xi will be denoted by Ψ(x), namely Ψ $\mathsf { \Pi } ^ { \prime } ( x ) = \langle \bar { x } | \Psi \rangle$ . The exact ground-state wave function $\left| { \Psi _ { 0 } } \right.$ can be obtained by the imaginary time evolution of a given variational state $\vert \Psi _ { V } \rangle$ :

$$
\vert \Psi _ { 0 } \rangle \propto \operatorname* { l i m } _ { \beta \to \infty } \vert \Psi _ { \beta } \rangle ,\tag{1}
$$

where $| \Psi _ { \beta } \rangle \equiv e ^ { - \beta \hat { H } } | \Psi _ { V } \rangle$ , provided that the variational state is non-orthogonal to $\vert \Psi _ { 0 } \rangle , \mathrm { i . e . , } \langle \Psi _ { V } \vert \Psi _ { 0 } \rangle \neq 0$ 0. Then, the ground-state expectation value of a quantum opera-

tor $\hat { O }$ can be obtained by

$$
\langle \hat { O } \rangle = \operatorname* { l i m } _ { \beta \to \infty } \frac { \langle \Psi _ { \beta } | \hat { O } | \Psi _ { \beta } \rangle } { \langle \Psi _ { \beta } | \Psi _ { \beta } \rangle } .\tag{2}
$$

A practical computational scheme can be conveniently introduced by considering a path-integral representation of the imaginary time evolution. To such a purpose, we split the total imaginary time $\beta$ into M slices of “duration” $\tau = \beta / M$ , in such a way that the value of the evolved wave function on a generic many-body state of the system reads

$$
\Psi _ { \beta } ( x _ { 0 } ) = \sum _ { x _ { 1 } \ldots x _ { M } } \prod _ { i = 1 } ^ { M } G _ { x _ { i - 1 } x _ { i } } ^ { \tau } \Psi _ { V } ( x _ { M } ) ,\tag{3}
$$

where we have introduced the imaginary time propagators

$$
\begin{array} { r } { G _ { x _ { i - 1 } x _ { i } } ^ { \tau } = \langle x _ { i - 1 } | e ^ { - \tau \hat { H } } | x _ { i } \rangle . } \end{array}\tag{4}
$$

Within this approach, it is easy to write expectation values of operators $\hat { O }$ that are diagonal in the chosen basis |xi, i.e., $\langle x | \hat { O } | y \rangle = O ( x ) \delta _ { x , y }$ . In fact, in this case we have that:

$$
\langle \hat { O } \rangle = \operatorname* { l i m } _ { \beta \to \infty } \frac { \sum _ { \mathbf { X } } \Pi ^ { \beta } ( \mathbf { X } ) O ( x _ { M } ) } { \sum _ { \mathbf { X } } \Pi ^ { \beta } ( \mathbf { X } ) } ,\tag{5}
$$

where the summation is extended to all possible imaginary time paths $\mathbf { X } \equiv \{ x _ { 0 } , x _ { 1 } , \ldots , x _ { 2 M } \}$ , and $\Pi ^ { \beta } ( \mathbf { X } )$ is given by:

$$
\Pi ^ { \beta } ( \mathbf { X } ) = \Psi _ { V } ( x _ { 0 } ) \left[ \prod _ { i = 1 } ^ { 2 M } G _ { x _ { i - 1 } x _ { i } } ^ { \tau } \right] \Psi _ { V } ( x _ { 2 M } ) .\tag{6}
$$

The ground-state energy can be conveniently obtained by means of the mixed average:

$$
E _ { 0 } = \operatorname* { l i m } _ { \beta  \infty } \frac { \sum _ { \mathbf { X } } \Pi ^ { \beta } ( \mathbf { X } ) E _ { L } ( x _ { 0 } ) } { \sum _ { \mathbf { X } } \Pi ^ { \beta } ( \mathbf { X } ) } ,\tag{7}
$$

where $E _ { L } ( x ) = \langle x | \hat { H } | \Psi _ { V } \rangle / \langle x | \Psi _ { V } \rangle$ is the so-called local energy.

Besides the static (i.e., equal-time) correlation functions, this formalism allows one to calculate also dynamical correlations in imaginary time $\begin{array} { r l } { C _ { A B } ( { \boldsymbol { \mathcal { T } } } ) } & { { } = } \end{array}$ $\langle \hat { A } ( \mathcal { T } ) \hat { B } ( 0 ) \rangle$ that can be computed as

$$
C _ { A B } ( \mathcal { T } ) = \operatorname* { l i m } _ { \beta  \infty } \frac { \sum _ { \mathbf { X } } \Pi ^ { \beta } ( \mathbf { X } ) A ( x _ { n } ) B ( x _ { m } ) } { \sum _ { \mathbf { X } } \Pi ^ { \beta } ( \mathbf { X } ) } ,\tag{8}
$$

where $x _ { n }$ and $x _ { m }$ are two coordinates of the path such that $( m - n ) \tau = \tau$

## III. REPTATION QUANTUM MONTE CARLO

A probabilistic interpretation of the previous expectation values (5), (7), and (8) can be immediately recovered whenever $\Pi ^ { \beta } ( \mathbf { X } ) \geq 0$ for all configurations X. Indeed, in this case, $\Pi ^ { \beta } ( \mathbf { X } )$ can be interpreted as a probability distribution that may be readily sampled by using Monte Carlo algorithms. This fact allows ground-state expectation values to be calculated exactly, within statistical errors.

The basic idea of the RQMC algorithm is to sample the distribution probability $\Pi ^ { \beta } ( \mathbf { X } )$ by using a Markov process with simple moves. Given the configuration $\mathbf { X } \equiv$ $\left\{ x _ { 0 } , x _ { 1 } , \dots , x _ { 2 M } \right\}$ , a new configuration is proposed in two possible ways: either $\mathbf { X } _ { L } \equiv \{ x _ { T } , x _ { 0 } , \dots , x _ { 2 M - 1 } \}$ (which we call “left move”) or $\mathbf { X } _ { R } \equiv \{ x _ { 1 } , \ldots , x _ { 2 M } , x _ { T } \}$ , (which we call “right move”). In both cases, xT is a new configuration proposed according to a suitable transition probability $R ^ { \tau } ( x  x _ { T } )$ , where x stays for x0 $\left( { { x } _ { 2 M } } \right)$ when the left (right) move is considered. Such “sliding moves” are depicted in Fig. 1. Ideally, the transition probability should guarantee the minimum possible statistical error on the desired observables and, to such a purpose, it has been proved useful to consider the propagator with importance sampling, i.e., $\tilde { G } _ { x y } ^ { \tau } \ = \ G _ { x y } ^ { \bar { \tau } } \Psi _ { V } ( y ) / \Psi _ { V } ( x )$ and take the following transition probability

$$
R ^ { \tau } ( x  y ) = \frac { \tilde { G } _ { x y } ^ { \tau } } { w ( x ) } ,\tag{9}
$$

where

$$
w ( x ) = \sum _ { x ^ { \prime } } \tilde { G } _ { x x ^ { \prime } } ^ { \tau }\tag{10}
$$

represents the normalization factor. The explicit form of $R ^ { \tau } ( x  x _ { T } )$ will be discussed in more detail in Sec. III B. The proposed configuration $\mathbf { X } _ { d }$ (where $d = L$ or R) is accepted or rejected according to the usual Metropolis algorithm, where the acceptance rate is given by:

$$
A = \operatorname* { m i n } \{ 1 , { \frac { \Pi ^ { \beta } ( { \bf X } _ { d } ) R ^ { \tau } ( x _ { T }  x ) } { \Pi ^ { \beta } ( { \bf X } ) R ^ { \tau } ( x  x _ { T } ) } } \} .\tag{11}
$$

In this way, a sequence of configurations $\mathbf { X } ^ { k }$ is generated, k being the (discrete) time index of the Markov chain.

In order to reduce the auto-correlation time of the observables it is convenient to make several consecutive sliding moves along the same imaginary time direction.5 To such a purpose, a recent development called “bounce” algorithm has been proposed.7 Although the bounce algorithm sampling procedure does not fulfill the microscopic detailed balance, the equilibrium probability $\Pi ^ { \beta } ( { \bf x } )$ is correctly sampled.7 The RQMC algorithm with bounce moves can be then summarized in the following steps:

1. For the current direction of the move and for the present configuration $\mathbf { X } ^ { k }$ , propose $x _ { T }$ according to the transition probability $R ^ { \tau } ( x  x _ { T } )$ , where $x =$ $x _ { 0 } ^ { k }$ if d = L and $x = x _ { 2 M } ^ { k }$ if $d = R$

2. Given the form of the acceptance ratio A of Eq. (11), accept the proposed configuration accord-

<!-- image-->  
Figure 1: Pictorial representation of the “sliding moves” along the right imaginary time direction. In the new configuration (bottom), a new head for the reptile is generated from the old configuration (top) and the tail is discarded.

ing to the probability

$$
A _ { L } = \operatorname* { m i n } \left\{ 1 , \frac { w ( x _ { 0 } ^ { k } ) } { w ( x _ { 2 M - 1 } ^ { k } ) } \right\} ,\tag{12}
$$

if d = L, or with probability

$$
A _ { R } = \operatorname* { m i n } \left\{ 1 , \frac { w ( x _ { 2 M } ^ { k } ) } { w ( x _ { 1 } ^ { k } ) } \right\} ,\tag{13}
$$

if d = R.

3. If the move is accepted, update the path configurations according to $\bar { \mathbf { X } } ^ { k + \bar { 1 } } \mathbf { \Lambda } = \mathbf { X } _ { d }$ and continue along the same direction, otherwise $\mathbf { X } ^ { k + 1 } \mathbf { \Phi } = \mathbf { \Phi } \mathbf { X } ^ { k }$ and change direction.

4. Go back to 1.

## A. Directed updates

At this point we introduce a novel alternative sampling approach, which generalizes the bounce idea while strictly fulfilling the detailed balance condition. Such a scheme, which is largely inspired by the loop algorithm methods devised for the stochastic series expansion11,13 and for the worm algorithm,12,14 allows one to choose the time direction in a purely Markovian way, i.e., independently of the previous history.

In our algorithm, a Markov step consists of many simple consecutive “sliding moves”, whose number is not fixed a-priori but is determined by a certain probability (see below). The actual Monte Carlo step takes place at the end of few consecutive updates along the currently chosen direction. In the examples below, we denote the number of these sliding moves between two Monte Carlo steps by s.

At the beginning of each Markov step we choose a direction d according to the probability $\dot { \mathbf { \zeta } } P ( \mathbf { X } ^ { k } , d )$ , whose form will be specified later. Assuming that the right direction has been chosen, we propose a new configuration xT , according to the transition probability $R ^ { \tau } { \bar { ( x _ { 2 M } ^ { k } \to } }$ $x _ { T } )$ and the configuration labels are shifted according to $\mathbf { X } ^ { k + 1 } \ = \ \{ x _ { 1 } ^ { k } , \ldots , x _ { 2 M } ^ { k } , x _ { T } \}$ with $x _ { 2 M } ^ { k + 1 } ~ = ~ x _ { T }$ At this point, we continue updates along this direction with probability $K ( \mathbf { X } ^ { k + 1 } , \dot {  } )$ , or stop with probability $[ 1 - \bar { K } ( \mathbf { X } ^ { k + 1 } , \bar {  } ) ]$ ]. If it has been decided to continue the updates, then a new configuration is generated according to $R ^ { \tau } ( x _ { 2 M } ^ { k + 1 }  x _ { T } )$ and the labels of the configuration are again shifted according to $\mathbf { X } ^ { k + 2 } = \{ x _ { 1 } ^ { k + 1 } , \ldots , x _ { 2 M } ^ { k + 1 } , x _ { T } \}$ The Markov step finishes after s consecutive updates along the right direction only when $K ( \mathbf { X } ^ { k + s } ,  ) < \xi _ { s } ,$ where $\xi _ { s }$ is a random number uniformly distributed in [0, 1). At this point a Metropolis test should be done, in order to accept or reject the sequence of intermediate s sliding moves:

$$
A = \operatorname* { m i n } \left\{ 1 , { \frac { q ( \mathbf { X } ^ { k + s } ) } { q ( \mathbf { X } ^ { k } ) } } \right\} ,\tag{14}
$$

where (see Appendix A)

$$
\begin{array} { l } { { q ( { \bf X } ) ~ = ~ \frac { P ( { \bf X } , \left. ) } { 1 - K ( { \bf X } , \right. ) } w ( x _ { 2 M - 1 } ) } } \\ { { ~ = ~ \frac { P ( { \bf X } , \right. ) } { 1 - K ( { \bf X } , \left. ) } w ( x _ { 1 } ) . } } \end{array}\tag{15}
$$

However, in order to avoid time-consuming restorations of the original configuration, it is preferable to accept all the moves, while keeping track of the residual weight $q ( \mathbf { X } )$ . This is possible since A only depends upon initial and final configurations, so that, given that all the intermediate moves are accepted, the sampled distribution probability is $\Pi ^ { \beta } ( { \mathbf { X } } ) \times \mathbf { \bar { \mu } } _ { q } ( { \mathbf { X } } )$ . The contribution of the current configuration to statistical averages must be then weighted by the factor $1 / q ( \mathbf { X } )$

To proceed to the next Markov step, a new direction d is chosen according to $P ( \mathbf { X } ^ { k + s } ,  )$ and $P ( \mathbf { X } ^ { k + s } ,  )$ and the updates are carried along the extracted new direction.

Let us now show the actual expressions for the aforementioned probabilities. In Appendix $\mathrm { A } ,$ it is demonstrated that the detailed balance is satisfied if one chooses the probabilities for the directions as

$$
P ( { \bf X } ,  ) = \frac { 1 } { 1 + a ( { \bf X } ) } ,\tag{16}
$$

$$
P ( { \bf X } ,  ) = \frac { a ( { \bf X } ) } { 1 + a ( { \bf X } ) } ,\tag{17}
$$

where

$$
a ( \mathbf { X } ) = \frac { w ( x _ { 2 M - 1 } ) } { w ( x _ { 1 } ) } \frac { 1 - K ( \mathbf { X } , \left. ) } { 1 - K ( \mathbf { X } , \right. ) } ,\tag{18}
$$

which is positive and, therefore, guarantees that the above defined quantities are well defined probabilities, i.e., $0 \leq P ( \mathbf { X } ,  ) \leq 1$ and $0 \leq P ( \mathbf { X } ,  ) \leq 1$ , with the additional property that $P ( \mathbf { X } , \longleftrightarrow ) + P ( \mathbf { X } , \to ) = 1$

Regarding the probabilities to continue the updates along the current direction, we have a substantial freedom of choice, provided that the condition $\begin{array} { r } { \frac { K ( \mathbf { X } , \left. ) } { K ( \mathbf { X } , \right. ) } = } \end{array}$ $\frac { w ( x _ { 1 } ) } { w ( x _ { 2 M - 1 } ) }$ , is satisfied (see Appendix A). In this paper we have adopted

$$
K ( { \bf X } ,  ) = \alpha \operatorname* { m i n } \{ 1 , b ( { \bf X } ) \} ,\tag{19}
$$

$$
K ( { \bf X } ,  ) = \alpha \operatorname* { m i n } \{ 1 , { \frac { 1 } { b ( { \bf X } ) } } \} ,\tag{20}
$$

where we have defined:

$$
b ( \mathbf { X } ) = \frac { w ( x _ { 1 } ) } { w ( x _ { 2 M - 1 } ) } .\tag{21}
$$

and $0 \textless \alpha \textless 1$ is an arbitrary parameter of the algorithm, which controls the average number of consecutive updates along the same direction.

Summarizing, the RQMC algorithm with directed updates consists of a sequence of Markov steps determined by the following rules:

1. Choose a time direction d according to the probabilities of Eqs. (16) and (17).

2. Propose a new configuration xT according to the transition probability $R ^ { \tau } ( x  x _ { T } )$ , where $x = x _ { 0 } ^ { k }$ if $d = L$ and $x = x _ { 2 M } ^ { k }$ if d = R.

3. Shift the configuration indexes according to $\mathbf { X } ^ { k + 1 } = \{ x _ { T } , x _ { 0 } ^ { k } , \ldots , x _ { 2 M - 1 } ^ { k } \}$ if d = L or $\mathbf { X } ^ { k + 1 } =$ $\{ x _ { 1 } ^ { k } , \ldots , x _ { 2 M } ^ { k } , x _ { T } \} { \mathrm { ~ i f ~ } } d = R$

4. According to the probability $K ( \mathbf { X } ^ { k } ,  )$ or $K ( \mathbf { X } ^ { k } ,  )$ , decide whether keep moving in the same direction or change direction. In the former case, go to 2, otherwise go to 5.

5. The Markov step ends here and the current configuration carries the weight $1 / q ( \mathbf { X } ^ { k + s } )$ , where s is the number of intermediate moves along the direction chosen.

The relationship between the directed update scheme and the bounce algorithm is further elucidated in $\mathrm { A p \mathrm { - } }$ pendix B, where general considerations about the efficiency of the algorithms are presented.

## B. Continuous-time propagator

One of the most striking differences between the original formulation of the RQMC on the continuum and the present formulation on the lattice is the lack of the discretization error appearing in the Trotter decomposition of the propagator. Indeed it is easier to carry the propagation in continuous imaginary time on a lattice,15 than on the continuum.16 To such a purpose, let us consider the limit of an infinitesimal imaginary time $\epsilon ,$ for which the transition probability of Eq. (9) can be written as

$$
R ^ { \epsilon } ( x  y ) \simeq \frac { \delta _ { x y } - \epsilon \Psi _ { V } ( y ) H _ { x y } / \Psi _ { V } ( x ) } { 1 - \epsilon E _ { L } ( x ) }\tag{22}
$$

$$
\simeq \delta _ { x y } \left[ 1 + \epsilon E _ { L } ( x ) \right] - \epsilon \left[ H _ { x y } \frac { \Psi _ { V } ( y ) } { \Psi _ { V } ( x ) } \right] + o ( \epsilon ^ { 2 } ) ,\tag{23}
$$

where $E _ { L } ( x )$ is the previously defined local energy and $H _ { x , y } = \langle x | \dot { H } | y \rangle$ denotes the matrix element of the Hamiltonian. Whenever $\Psi _ { V } ( y ) H _ { x y } / \Psi _ { V } ( x )$ is non positive for all x and $y ,$ this equation takes the form of a continuoustime Markov process, whose analytical properties are well known. In particular, the probability distribution for the “waiting $\mathrm { t i m e } ^ { \mathrm { 3 } } \tau _ { w }$ in a given state x, i.e., the average time that the system spends in the state x before making an off-diagonal transition to another state $y \neq x$ , is exactly known, namely $P ( \tau _ { w } ; x ) = \exp \{ - \tau _ { w } \left[ H _ { x x } - E _ { L } ( x ) \right] \}$ . As a consequence, the finite-time propagator $R ^ { \tau } ( x  y )$ can be directly sampled, giving rise to a succession of a certain number n of consecutive transitions $x  z _ { 1 }  z _ { 2 } $ $\cdots  y .$ with corresponding waiting times $\tau _ { w } ( z _ { i } )$ (such that $\begin{array} { r } { \sum _ { i } \tau _ { w } ( z _ { i } ) = \tau ) } \end{array}$ . The normalization of the whole process is

$$
w ( x ) = \exp \left[ - \sum _ { i } \tau _ { w } ( z _ { i } ) E _ { L } ( z _ { i } ) \right] ,\tag{24}
$$

where the waiting times are extracted according to the exponential probability $P ( \tau _ { w } ; z _ { i } )$ . The transitions between the intermediate configurations are done according to the off-diagonal elements of Eq. (23), $\mathrm { i . e . , } \ z _ { i + 1 }$ is chosen with probability proportional to $[ - \Psi _ { V } ( \dot { z } _ { i + 1 } ) H _ { z _ { i } z _ { i + 1 } } / \Psi _ { V } ( z _ { i } ) ]$

## C. Off-diagonal observables

The formalism so-far developed allows one to successfully compute pure ground-state expectation values of operators that are diagonal in the local basis x, with the expectation values of off-diagonal operators restricted to the so-called mixed averages.5,6,15 Nonetheless, it is often of great interest to remove such a limitation (whose result is biased by the quality of the variational wave function) and a dedicated sampling strategy has to be devised in order to cope with such a need. In the following, we show that a relatively easy modification of the sampling scheme can accomplish this task, providing us with a general tool to compute ground-state averages of operators that are non local in the chosen basis x.

Let us consider an arbitrary off-diagonal observable $\hat { \mathcal { O } } _ { : }$ which can be in turn considered as the summation of many observables we are interested in, i.e., $\begin{array} { r } { \hat { \mathcal { O } } = \sum _ { d } \hat { \mathcal { O } } ^ { ( d ) } } \end{array}$ For example, we can imagine these operators to be the components of the one-body density matrix at a given distance, $\begin{array} { r } { \hat { \mathcal { O } } ^ { ( d ) } = \sum _ { \langle r , r ^ { \prime } \rangle _ { d } } b _ { r } ^ { \dagger } \dot { b } _ { r ^ { \prime } } } \end{array}$ with the summation extended to all lattice coordinates at a fixed distance d.

In the spirit of Refs 12,14 we introduce a wormoperator defined by

$$
\mathcal { W } _ { x , y } = \delta _ { x , y } + \lambda \mathcal { O } _ { x , y } ,\tag{25}
$$

where λ is a positive constant, and consider the extended configuration space spanned by the paths

$$
\begin{array} { l } { { \displaystyle \Pi _ { \mathcal { W } } ^ { \beta } ( { \bf X } ) ~ = ~ \Psi _ { V } \left( x _ { 0 } \right) \times \prod _ { i = 1 } ^ { L } { G _ { x _ { i - 1 } x _ { i } } ^ { \tau } \times \mathcal { W } _ { x _ { L } x _ { R } } \times } } } \\ { { \displaystyle ~ \times \prod _ { i = R + 1 } ^ { 2 M + 1 } { G _ { x _ { i - 1 } x _ { i } } ^ { \tau } \times \Psi _ { V } \left( x _ { 2 M + 1 } \right) } . } } \end{array}\tag{26}
$$

The extended paths are broken in two (space)- discontinuous pieces by the worm operator, which is placed at an imaginary-time $0 \le \tau _ { L R } \le \beta$ Therefore, paths contain $2 ( M + 1 )$ coordinates, including $x _ { L }$ and $x _ { R }$ that refer to the same imaginary time $\tau _ { L R }$

The configuration space spanned by Eq. (26) is clearly larger than the one spanned by Eq. (6), which is recovered whenever $x _ { L } = x _ { R }$ , i.e., when the worm operator is diagonal.

The pure ground-state expectation value of the operator Oˆ is conveniently written in terms of the extended paths as

$$
\langle \hat { \mathcal { O } } \rangle = \frac { 1 } { \lambda } \operatorname* { l i m } _ { \beta \to \infty } \frac { \sum _ { \mathbf { X } } \Pi _ { \mathcal { W } } ^ { \beta } ( \mathbf { X } ) \times \Theta ( x _ { L } \neq x _ { R } ) } { \sum _ { \mathbf { X } } \Pi _ { \mathcal { W } } ^ { \beta } ( \mathbf { X } ) \times \Theta ( x _ { L } = x _ { R } ) } ,\tag{27}
$$

where $\Theta ( C ) \neq 0$ whenever condition C is satisfied. The modulus of Eq. (26) can be in turn interpreted as a probability distribution and stochastically sampled by means of the elementary sliding moves considered before. Indeed, whenever the worm operator is far from the ends of the imaginary-time paths, the sampling scheme remains unchanged. In this case, a move along direction d will generate a new head (or tail) for the reptile according to $R ^ { \tau } ( x  x _ { T } )$ while shifting the worm position of τ. On the other hand, whenever the worm operator reaches the ends of the reptile, a new worm configuration is proposed on the other side; in analogy with the previous analysis, new configurations are generated according to a transition probability

$$
R ^ { W } ( x \to y ) = \frac { 1 } { \bar { w } ( x ) } \left| { \mathcal W } _ { x y } \frac { \psi _ { V } ( y ) } { \psi _ { V } ( x ) } \right| ,\tag{28}
$$

where $\bar { w } ( x )$ is the normalization factor. Due to the particular form of the matrix elements (25), the transition probability will lead either to diagonal configurations $\mathbf { \Phi } ( x \ \mathbf { \Phi } = \ \ - \ y )$ or to off-diagonal configurations $( x ~ \neq ~ y )$ , thus generating continuous and discontinuous paths. The relative probability for diagonal and off-diagonal configurations depends on the value of λ that can be tuned in order to reach a balanced sampling frequency for the different sectors of the extended paths. In order to exemplify the worm updates, let us consider the case in which $d \ = \ R$ and a configuration $\begin{array} { r } { \Psi _ { V } \left( x _ { 0 } \right) \mathcal { W } _ { x _ { 0 } x _ { 1 } } \left[ \prod _ { i = 2 } ^ { 2 M + 1 } G _ { x _ { i - 1 } x _ { i } } ^ { \tau } \right] \Psi _ { V } \left( x _ { 2 M + 1 } \right) } \end{array}$ , after a sliding update in the right direction, we will have $\begin{array} { r } { \Psi _ { V } \left( x _ { 1 } \right) \left[ \prod _ { i = 2 } ^ { 2 M + 1 } G _ { x _ { i - 1 } x _ { i } } ^ { \tau } \right] \mathcal { W } _ { x _ { 2 M + 1 } x _ { T } } \Psi _ { V } \left( x _ { T } \right) } \end{array}$ , where xT is proposed according to the transition probability $R ^ { \hat { \mathcal { W } } } ( x _ { 2 M + 1 } \to x _ { T } )$ (see Fig. 2). In analogy with the previous case, the acceptance factor for the bounce moves reads $\begin{array} { r } { \bar { A } _ { R } = \operatorname* { m i n } \left\{ 1 , \frac { \bar { w } ( x _ { 2 M + 1 } ) } { \bar { w } ( x _ { 1 } ^ { k } ) } \right\} } \end{array}$

<!-- image-->  
Figure 2: Pictorial representation of the “sliding moves” along the right imaginary-time direction when the worm operator sits at the tail of the reptile. In the new configuration (bottom), a new head for the reptile is generated from the old configuration (top), the old tail configuration is discarded and the worm discontinuity is moved to the “neck” of the reptile.

Summarizing, the RQMC with worm-updates consists of the following steps:

1. For the current direction of the move d and for the present configuration $\mathbf { X } ^ { k }$ consider the wormoperator position $\tau _ { L R }$

2. If the worm is not at the ends of the reptile (i.e., $\tau _ { L R } \neq 0$ when $d = L$ and $\tau _ { L R } \neq \beta$ when $d = R )$ go to step (a), otherwise go to step (b).

(a) Propose a new configuration $x _ { T }$ according to the transition probability $R ^ { \tau } ( x  x _ { T } )$ , where $x = x _ { 0 } ^ { k }$ if $d = L$ and $x = x _ { 2 M + 1 } ^ { k } { \mathrm { ~ i f ~ } } d = R$ . The new configuration is accepted with probability

$$
{ \cal A } _ { L } = \operatorname* { m i n } \left\{ 1 , \frac { w ( x _ { 0 } ^ { k } ) } { w ( x _ { 2 M } ^ { k } ) } \right\} ,\tag{29}
$$

if d = L, or with probability

$$
A _ { R } = \operatorname* { m i n } \left\{ 1 , \frac { w ( x _ { 2 M + 1 } ^ { k } ) } { w ( x _ { 1 } ^ { k } ) } \right\} ,\tag{30}
$$

if $d = R$ . In the proposed state $\mathbf { X } _ { d } .$ , all the configuration labels are shifted in the d direction, determining in turn a shift of the worm operator of a time interval ±τ, depending on d.

(b) Propose a new configuration xT according to the worm transition probability $R ^ { \mathcal { W } } ( x  x _ { T } )$ , where $x = x _ { L } ^ { k } = x _ { 0 } ^ { k }$ if $d = L$ and $x = x _ { R } ^ { k } =$ $x _ { 2 M + 1 } ^ { k } \mathrm { i f } d = \mathbf { \bar { \Gamma } }$ . Accept the new configuration with probability

$$
\bar { A } _ { L } = \operatorname* { m i n } \left\{ 1 , \frac { \bar { w } ( x _ { 0 } ^ { k } ) } { \bar { w } ( x _ { 2 M } ^ { k } ) } \right\} ,\tag{31}
$$

if $d = L$ , or with probability

$$
\begin{array} { r } { \bar { A } _ { R } = \operatorname* { m i n } \left\{ 1 , \frac { \bar { w } ( x _ { 2 M + 1 } ^ { k } ) } { \bar { w } ( x _ { 1 } ^ { k } ) } \right\} , } \end{array}\tag{32}
$$

if $d = R$ . In the proposed state $\mathbf { X } _ { d } ,$ all the configuration labels are shifted in the d direction, and the worm operator is moved from the head (tail) to the tail (head) of the reptile, depending on d.

3. If the move is accepted, update the path configurations according to $\bar { \mathbf { X } } ^ { k + \bar { 1 } } \mathbf { \Lambda } = \mathbf { X } _ { d }$ and continue along the same direction, otherwise $\mathbf { X } ^ { k + 1 } \mathbf { \Phi } = \mathbf { \Phi } \mathbf { X } ^ { k }$ and change direction.

4. Go back to 1.

This scheme samples the probability density associated to the modulus of Eq. (26), and the expectation values of the individual components $\hat { \mathcal { O } } _ { d }$ can be recast as statistical averages over such a probability distribution, while keeping track of the overall sign of the extended paths. In particular the best estimate of the ground-state expectation values is obtained when the worm is in the central part of the path, at $\tau _ { L R } = \beta / 2$ , leading to

$$
\begin{array} { l l l } { { \displaystyle \langle \hat { \mathcal { O } } ^ { ( d ) } \rangle = \frac { \sum _ { \mathbf { X } } \Pi _ { \mathcal { W } } ^ { \beta } ( { \bf X } ) \times \Theta \left( \mathcal { O } _ { x _ { L } x _ { R } } ^ { ( d ) } \neq 0 , \tau _ { L R } = \frac { \beta } { 2 } \right) } { \sum _ { \mathbf { X } } \Pi _ { \mathcal { W } } ^ { \beta } ( { \bf X } ) \times \Theta \left( x _ { L } = x _ { R } , \tau _ { L R } = \frac { \beta } { 2 } \right) } \ ~ } } \\ { { \displaystyle = \frac { 1 } { \lambda } \frac { \left. \Theta \left( \mathcal { O } _ { x _ { L } x _ { R } } ^ { ( d ) } \neq 0 \right) \times \mathrm { s i g n } \left[ \Pi _ { \mathcal { W } } ^ { \beta } ( { \bf X } ) \right] \right. _ { \mathrm { O D } } ^ { \mathrm { c e n t e r } } } { N _ { D } ^ { \mathrm { c e n t e r } } } , ~ \ ~ ( 3 ) } } \end{array}\tag{3}
$$

where $\langle \dots \rangle _ { \mathrm { O D } } ^ { \mathrm { c e n t e r } }$ denotes statistical averages over the offdiagonal distribution $\textstyle \left| \Pi _ { \mathcal { W } } ^ { \beta } ( \mathbf { X } ) \right| \Theta ( x _ { L } \ \neq \ x _ { R } , \tau _ { L R } \ = \ \frac { \beta } { 2 } )$ and $N _ { D } ^ { \mathrm { c e n t e r } }$ is the number of configurations sampled with a diagonal worm operator in the center of the paths.

## D. Tackling the sign problem

When the probability distribution of Eq. (6) is not positive defined, as is generally the case with fermions, the probabilistic interpretation of the imaginary time paths breaks down. This circumstance, which is known as “sign problem”, originates whenever $\Psi _ { V } ( y ) H _ { x y } / \Psi _ { V } ( x ) > 0$ for some element $x \ \neq \ y .$ In this case, it is not possible to have polynomial algorithms that are able to obtain an exact solution of the problem, which would imply to sample correctly the resulting signs. Therefore, approximated schemes are welcome and often adopted, the most widespread one being the so-called fixed-node (FN) approximation. For lattice systems, this approach relies on the definition of an effective Hamiltonian, which depends parametrically on the nodal structure of a given variational wave function $\Psi _ { V } ( x ) = \langle x | \Psi _ { V } \rangle . ^ { 1 7 }$ The matrix elements of the FN Hamiltonian are defined as

$$
H _ { x y } ^ { \mathrm { f n } } = { \left\{ \begin{array} { l l } { H _ { x x } + \nu _ { \mathrm { s f } } ( x ) } & { { \mathrm { i f ~ } } x = y } \\ { H _ { x y } } & { { \mathrm { i f ~ } } \Psi _ { V } ( y ) H _ { x y } \Psi _ { V } ( x ) \leq 0 } \\ { 0 } & { { \mathrm { i f ~ } } \Psi _ { V } ( y ) H _ { x y } \Psi _ { V } ( x ) > 0 } \end{array} \right. }\tag{34}
$$

where the sign-flip potential is νsf(x) = $\begin{array} { r l } { \sum _ { u : \mathrm { s f } } \Psi _ { V } ( y ) H _ { x y } / \bar { \Psi } _ { V } ( \bar { x } ) } & { { } \bar { } } \end{array}$ the sum being extended to all the sign-flip states defined by the condition $\Psi _ { V } ( y ) H _ { x y } \Psi _ { V } ( x )$ > 0. With such a choice, the transition matrix $R ^ { \epsilon } ( x  y )$ of Eq. (23) is always positive definite and the summation of Eq. (3) is now restricted –which results in the FN approximation– to a region of the Hilbert space in which imaginary time paths are positive definite. Therefore, within the FN approximation, the ground-state wave function $| \Psi ^ { \mathrm { f n } } \rangle$ of $\hat { H } ^ { \mathrm { f i n } }$ can be stochastically sampled without any sign problem. Moreover, it is easy to show that the FN approximation becomes exact whenever the signs of the exact ground state are known. Most importantly, it has been proven17 that the FN ground-state energy $\hat { E } ^ { \mathrm { f n } } = \langle \hat { H } ^ { \mathrm { f n } } \rangle$ gives a rigorous upper-bound to the exact ground-state one and improves the pure variational results.

At this point, we introduce a straightforward, although computationally expensive, way to improve further the FN energy. Our strategy amounts to compute the expectation values of arbitrary powers of the original Hamiltonian $\hat { H }$ on the FN ground state $| \Psi _ { \mathrm { f n } } \rangle$ , namely

$$
L _ { k } = \frac { \langle \Psi _ { \mathrm { f n } } | \hat { H } ^ { k } | \Psi _ { \mathrm { f n } } \rangle } { \langle \Psi _ { \mathrm { f n } } | \Psi _ { \mathrm { f n } } \rangle } .\tag{35}
$$

The FN ground state can be expanded in the basis set of the eigenstates of Hˆ as $| \Psi _ { \mathrm { f n } } \rangle = \hat { \gamma } _ { 0 } | \Psi _ { 0 } \rangle + \gamma _ { 1 } | \Psi _ { 1 } \rangle + \gamma _ { 2 } | \Psi _ { 2 } \rangle +$ . . . and $L _ { k } = \gamma _ { 0 } ^ { 2 } E _ { 0 } ^ { k } + \gamma _ { 1 } ^ { 2 } E _ { 1 } ^ { k } + \gamma _ { 2 } ^ { 2 } E _ { 2 } ^ { k } + . . . , \mathrm { w i t h } \sum _ { i } \gamma _ { i } ^ { 2 } = 1$

Since very often the FN wave function has a considerable overlap with only few low-energy states, the knowledge of the first few moments of the Hamiltonian are enough to approximately reconstruct both the coefficients $\gamma _ { i }$ and the energies $E _ { i }$ . To such a purpose, let us consider a typical situation in which only the first 2n moments of the Hamiltonian have been numerically calculated and are therefore known. We can then truncate the expansion for $L _ { k }$ to the order n  1 having a closed system of 2n equations

$$
L _ { k } = \sum _ { i = 0 } ^ { n - 1 } \gamma _ { i , n } ^ { 2 } E _ { i , n } ^ { k } ,\tag{36}
$$

for $k = 0 , \ldots 2 n - 1$ that can be solved for the unknowns $\gamma _ { i , n }$ and $E _ { i , n }$ . In the limit of large n, the approximated

$E _ { 0 , n }$ converges to the exact ground-state energy. Moreover, we verified that $E _ { 0 , n } \ \geq \ E _ { 0 }$ , as a result of a connection between the solutions of the $\operatorname { E q . }$ (36) and the Lanczos procedure written in terms of the moments of the Hamiltonian.18

The Hamiltonian moments are off-diagonal operators and can, in principle, measured according to the sampling procedure detailed in Sec. III C. In the present implementation we are able to achieve sufficient statistical accuracy only for the first moment of the Hamiltonian, i.e., $L _ { 1 } = \langle \hat { H } \rangle$ i, while higher moments are too noisy. Yet, to our knowledge our algorithm is the only one that allows the calculation of the expectation value of the original Hamiltonian $\hat { H }$ . This is known17 to be a better upper bound than the expectation value of the FN Hamiltonian accessible with other zero-temperature algorithms.

Although we are not currently in position to measure directly the Hamiltonian moments $L _ { k }$ we have a controlled access to the mixed averages

$$
L _ { k } ^ { \mathrm { m i x } } = \frac { \langle \Psi _ { \mathrm { f n } } | \hat { H } ^ { k } | \Psi _ { V } \rangle } { \langle \Psi _ { \mathrm { f n } } | \Psi _ { V } \rangle } ,\tag{37}
$$

which present optimal statistical uncertainty. Moreover, an improved estimate of the ground-state energy based on the knowledge of the first few moments $L _ { k } ^ { \mathrm { m i x } }$ can be obtained solving a system of equation similar to Eq. (36) that leads to the approximate ground-state energies $\dot { E } _ { i , n } ^ { \mathrm { m i x } }$ . Unfortunately, the proof that $E _ { i , n } ^ { \mathrm { m i x } } \geq E _ { 0 }$ for $n > 1 ,$ is far from being trivial, requiring a generalization of the already non-trivial upper bound for $n = 1$ described in Ref. 17. Nonetheless, we have numerically verified that, in all the cases treated in this paper (where $E _ { 0 }$ is a-priori known), the condition $E _ { i , n } ^ { \mathrm { m i x } } \geq E _ { 0 }$ is always verified. We are then led to conjecture that this may always be the case.

## IV. RESULTS

## A. Low-energy excitations and spin correlations of the Heisenberg model

Hereafter, we present a simple application of the previous ideas to sign-problem free spin Hamiltonians. Let us consider the one-dimensional quantum Heisenberg model

$$
{ \hat { H } } = J \sum _ { i } { \hat { \mathbf { S } } } _ { i } \cdot { \hat { \mathbf { S } } } _ { i + 1 } ,\tag{38}
$$

where $\hat { \mathbf { S } } _ { i } = \left( \hat { S } _ { i } ^ { x } , \hat { S } _ { i } ^ { y } , \hat { S } _ { i } ^ { z } \right)$ is the spin $1 / 2$ operator on the site i and $J \overset { \cdot } { > } 0$ is the nearest-neighbor super-exchange coupling.

The total number of sites is denoted by L and periodicboundary conditions are assumed. This model can be solved exactly by using the so-called Bethe ansatz technique.19 Information on the excitation spectrum can be obtained from the dynamical structure factor

<!-- image-->  
Figure 3: Lowest-energy excitations as a function of the wavevector q for an $L = 2 0$ Heisenberg chain. The energies are extracted from the dynamical structure factor $S ( q , \omega )$ and are compared to exact results by the Lanczos method.

<!-- image-->  
Figure 4: The same as in Fig. 3 for $L = 8 0$ . Exact results are given by Bethe ansatz.

$$
S ( q , \omega ) = \int d t \langle \hat { S } _ { q } ^ { z } ( t ) \hat { S } _ { - q } ^ { z } ( 0 ) \rangle e ^ { i \omega t } ,\tag{39}
$$

where $\hat { S } _ { q } ^ { z } ( t ) = 1 / \sqrt { L } \sum _ { j } \hat { S } _ { j } ^ { z } ( t ) e ^ { i q j }$ is the Fourier transform of time-evolved spin projection on the z-axis. By introducing a complete set of eigenstates of the Hamiltonian $\left| { \Psi _ { n } } \right.$ with eigenvalues $E _ { n } ,$ we have that

$$
S ( q , \omega ) = \sum _ { n \neq 0 } | \langle \Psi _ { 0 } | \hat { S } _ { q } ^ { z } | \Psi _ { n } \rangle | ^ { 2 } \delta ( \omega - \omega _ { n } ) ,\tag{40}
$$

where $\omega _ { n } = ( E _ { n } - E _ { 0 } )$ . In the thermodynamic limit, the spin-1 states form a branch, which is very similar to spin waves in standard ordered systems, although no long-range order is found in one dimension.

Imaginary time correlation functions of arbitrary (diagonal) operators can be efficiently evaluated via Eq. (8).

<!-- image-->  
Figure 5: Ground-state expectation value of the spin-spin correlation function $\mathcal { C } ( d )$ for the Heisenberg model on a 80-site chain.

This fact allows us to have a direct access to $S ( q , T ) =$ $\langle \hat { S } _ { q } ^ { z } ( T ) \hat { S } _ { - q } ^ { z } ( 0 ) \rangle$ . This imaginary time correlation function can be then analytically continued, by using the Maximum-Entropy method,20 in order to have a reasonable numerical estimate for the dynamical structure factor of Eq. (40).

Before presenting the results, let us mention that we consider the following Jastrow state as a variational wave function:21,22

$$
\left| \Psi _ { V } \right. = \exp \left[ \sum _ { i , j } v _ { i j } \hat { S } _ { i } ^ { z } \hat { S } _ { j } ^ { z } \right] \left| F M \right.\tag{41}
$$

where $| F M \rangle$ is the ferromagnetic state along the x direction, for which $\langle x | F M \rangle$ does not depend upon the spin configuration and the variational parameters $v _ { i j }$ are optimized by using the method of Ref. 23.

In Fig. 3, we show the results for a small $L = 2 0$ system, where exact diagonalizations are possible by using the Lanczos method. We report the energy excitations $\Delta E ( q ) = E _ { q } - E _ { 0 }$ for the lowest state with $S = 1$ and fixed momentum q. In this case a perfect agreement between our RQMC results and the exact ones is found. Moreover, also on larger systems a very good accuracy is possible (see Fig. 4), demonstrating the performances of our numerical algorithm.

In order to exemplify the potentialities of the scheme outlined in III C, we conclude this part of the results devoted to the Heisenberg model showing the ground-state expectation value of the spin-spin correlation at distance d

$$
\mathcal { C } ( d ) \ = \ \frac { 1 } { L } \sum _ { i } \left( \hat { \bf S } _ { i } \cdot \hat { \bf S } _ { i + d } \right) .\tag{42}
$$

The desired observable is used as a worm operator and the value of the correlation function at the various distances is computed by means of the estimator of Eq. (33).

In Fig. 5, we show the expectation value of $\mathcal C ( d )$ for a 80-site one-dimensional lattice. In this case we are able to achieve very good statistics for the off-diagonal observable, with a relatively negligible computational effort, when compared to the evaluation of the ground-state expectation value of other diagonal observables.

## B. Ground-state properties of the fermionic Hubbard model

As an example of the application of the RQMC to signproblem affected Hamiltonians, we present some results for the fermionic Hubbard model on a square lattice, defined by:

$$
\hat { H } = - t \sum _ { { \left. i , j \right. } , \sigma } \hat { c } _ { { i , \sigma } } ^ { \dagger } \hat { c } _ { { j , \sigma } } + { h . c . } + U \sum _ { i } \hat { n } _ { { i , \uparrow } } \hat { n } _ { { i , \downarrow } } ,\tag{43}
$$

where $\langle \dots \rangle$ indicate nearest-neighbor sites, $\hat { c } _ { i , \sigma } ^ { \dag } \left( \hat { c } _ { i , \sigma } \right)$ creates (destroys) an electron on the site i with spin σ, and $\hat { n } _ { i , \sigma } = \hat { c } _ { i , \sigma } ^ { \dagger } \hat { c } _ { i , \sigma }$ . As a variational state we consider

$$
\left| \Psi _ { V } \right. = \exp \left[ \sum _ { i , j } v _ { i j } \hat { n } _ { i } \hat { n } _ { j } \right] \left| F S \right.\tag{44}
$$

where $| F S \rangle$ is the non-interacting Fermi sea and the Jastrow factor involves density-density correlations. The variational parameters $v _ { i j }$ entering in the Jastrow factor may be optimized again by minimizing the variational energy with the method of Ref. 23. In order to avoid open shells in $| F S \rangle$ , we consider 45-degrees tilted lattices with $L = 2 \times l ^ { 2 }$ sites, such that both the half-filled case and selected holes-doped cases are closed shells.

<!-- image-->  
Figure 6: Ground-state energy for the fermionic Hubbard model at half filling on a 18-sites tilted-square lattice. The energy difference $\Delta E = E _ { \mathrm { e x a c t } } - E$ is computed with distinct approximations described in the text.

<table><tr><td> $U / t$ </td><td> $E ^ { \mathrm { f n } }$ </td><td>h Hi</td><td> $\underline { { { E } _ { 0 , 2 } ^ { \mathrm { m i x } } } }$ </td></tr><tr><td>4</td><td>−42.850(1)</td><td>43.16(1)</td><td>−43.282(1)</td></tr><tr><td>5</td><td>−36.364(1)</td><td>-36.51(1)</td><td>−37.052(1)</td></tr><tr><td>6</td><td>−31.885(1)</td><td>-32.17(1)</td><td>−32.640(1)</td></tr><tr><td>7</td><td>−28.318(1)</td><td>-28.66(1)</td><td>−29.022(1)</td></tr><tr><td>8</td><td>−25.382(1)</td><td>-25.62(1)</td><td>−26.056(1)</td></tr></table>

Table I: Ground-state energy as a function of the Hubbard U repulsion on the 50-site lattice at half filling.

<table><tr><td> $N$ </td><td> $E ^ { \mathrm { f n } }$ </td><td>h hH</td><td> $\underline { { E _ { 0 , 2 } ^ { \mathrm { m i x } } } }$ </td><td> $E _ { A F }$ </td></tr><tr><td>50</td><td>−42.850(1)</td><td>−43.16(1)</td><td>−43.282(1)</td><td>−43.983(1)</td></tr><tr><td>42</td><td>−53.402(1)</td><td>−53.57(1)</td><td>−53.769(1)</td><td>−54.001(1)</td></tr><tr><td>26</td><td>−55.4325(1)</td><td>−55.63(1)</td><td>−55.6112(1)</td><td>−55.782(1)</td></tr><tr><td>18</td><td>−50.4127(1)</td><td>−50.50(1)</td><td>−50.4383(1)</td><td>−50.474(1)</td></tr></table>

Table II: Ground-state energy as a function of the number of electrons N for Hubbard repulsion $U / t = 4$ on a 50-site lattice. The numerically exact results obtained by the Auxiliary-Field Monte Carlo method $E _ { A F }$ are also shown for comparison.25

Let us start by showing the results for 18 electrons on 18 sites, where Lanczos diagonalizations are possi-$\mathrm { b l e . ^ { 2 4 } }$ In Fig. 6, we report our results for the groundstate energy. The FN approach gives rather accurate results for small values of $U / t ,$ i.e., $U / t \lesssim 4$ , where $( E _ { \mathrm { e x a c t } } - E ^ { \mathrm { f n } } ) / E _ { \mathrm { e x a c t } } \lesssim 0 . 0 1$ . By increasing the on-site interaction, the FN approach becomes worse and worse. This fact is due to the choice of the variational wave function that does not contain antiferromagnetic order. Remarkably, a considerable improvement may be obtained by considering the pure expectation value of the Hamiltonian, which is systematically lower than the FN energy, as demonstrated in Ref. 17 and now accessible within our algorithm. Further improvements to the FN energy can be obtained upon considering few (up to three) higher moments of the Hamiltonian measured as mixed-averages, see Fig. 6. The scheme based upon the Hamiltonian moments (described in Sec. III D) allows us to reach a great accuracy for the ground state energy, with a residual error almost independent of $U / t$ . Indeed, in this way we have $( E _ { \mathrm { e x a c t } } - \bar { E } ) / E _ { \mathrm { e x a c t } } \lesssim 0 . 0 0 2$ up to $U / t = 8$

This approach remains very effective also for larger systems, even though the variational wave function loses accuracy by increasing the cluster size (because the ground state has antiferromagnetic order in the thermodynamic limit, while the variational state is paramagnetic). In Table I, we report the ground-state energy for 50 sites for the half-filled case, while in Table II we report the ground-state energies for selected cases at finite holedoping, where numerically exact results (for moderate values of U and moderate lattice sizes) can be obtained by the Auxiliary-Field Monte Carlo method.25

## V. CONCLUSIONS

In this paper we have provided an efficient and general formulation of the reptation quantum Monte Carlo technique on lattice models. In particular, we showed an alternative sampling approach which generalizes the bounce algorithm, previously introduced to reduce autocorrelation time of the observables. Our scheme allows one to choose the time direction in a purely Markovian way. In addition, the average number of consecutive moves along the time directions may be optimized by a fine tuning of a certain parameter that has been expressly introduced in the transition probabilities. We reported benchmarks for two different models with pure bosonic and fermionic degrees of freedom, by showing to what extent it is possible to have accurate results both on the ground state and low-energy excitations. The introduction of a general method to compute ground-state expectation values of arbitrary off-diagonal observables also constitutes an important achievement, which will ease the study of relevant properties such as Bose-Einstein condensation and superconductivity phenomena in strongly interacting models. In addition, the possibility to directly measure the pure ground-state expectation values may open the way to a better optimization of the correlated wave function associated to the ground-state of an effective Hamiltonian which is not the FN one.

## Acknowledgments

It is a pleasure to acknowledge here precious discussions with S. Sorella and A. Parola. We also acknowledge support from CINECA and COFIN 07.

## Appendix A: Derivation of the probabilities for the directed-update scheme

In this Appendix we give a detailed derivation of the probabilities for the directed updates. The detailed balance condition guarantees that the given probability distribution $\Pi ^ { \beta } ( \mathbf { X } )$ is sampled if transitions from an initial state $\mathbf { X } ^ { k }$ to a final state $\mathbf { \bar { X } } ^ { k + s }$ differing for s intermediate updates are accepted according to:

$$
A ^ { s } = \operatorname* { m i n } \{ 1 , { \frac { \Pi ^ { \beta } ( { \bf X } ^ { k + s } ) } { \Pi ^ { \beta } ( { \bf X } ^ { k } ) } } { \frac { { \cal T } ^ { s } ( { \bf X } ^ { k + s }  { \bf X } ^ { k } ) } { { \cal T } ^ { s } ( { \bf X } ^ { k }  { \bf X } ^ { k + s } ) } } \} ,\tag{A1}
$$

$\mathcal { T } ^ { s }$ being the overall transition probability between the two states. Let us first consider the case when $s = 1$ and fix the right direction $d = R$ (a similar derivation can be obtained for $d = L )$ . In this case, the transition probability from the initial state to the final state reads

$$
\begin{array} { r c l } { { \mathcal { T } ^ { 1 } ( { \bf X } ^ { k }  { \bf X } ^ { k + 1 } ) ~ = ~ { \cal P } ( { \bf X } ^ { k } ,  ) \times { \cal R } ^ { \tau } ( x _ { 2 M } ^ { k }  x _ { 2 M } ^ { k + 1 } ) \times } } \\ { { } } & { { } } & { { } } \\ { { \mathrm {  ~ \times ~ } [ 1 - K ( { \bf X } ^ { k + 1 } ,  ) ] , ~ \mathrm {  ~ \Lambda ~ } \mathrm {  ~ \Lambda ~ } ( \mathrm { A } 2 ) } } \end{array}
$$

namely, it is the product of the probability of having chosen the right direction, times the transition probability for the new tail of the reptile, times the probability of stopping the updates after one intermediate step. The inverse transition probability instead reads

$$
\begin{array} { r l r } { \mathcal { T } ^ { 1 } ( \mathbf { X } ^ { k + 1 } \to \mathbf { X } ^ { k } ) ~ = ~ } & { { } } & { P ( \mathbf { X } ^ { k + 1 } ,  ) \times R ^ { \tau } ( x _ { 0 } ^ { k + 1 } \to x _ { 0 } ^ { k } ) \times } \\ { ~ } & { { } } & { \times ~ [ 1 - K ( \mathbf { X } ^ { k } ,  ) ] , ~ \mathrm { ( A 3 ) } } \end{array}
$$

which can be obtained reversing the time directions and considering transitions from the head of the reptile instead that from the tail. Therefore, the acceptance factor reads as

$$
\begin{array} { r c l } { A ^ { 1 } = \displaystyle \operatorname* { m i n } \left\{ 1 , \frac { 1 - K ( { \mathbf { X } } ^ { k } , \left. ) } { P ( { \mathbf { X } } ^ { k } , \right. ) } \times \right. } \\ { \displaystyle \left. \times \frac { w ( x _ { 2 M - 1 } ^ { k + 1 } ) } { w ( x _ { 1 } ^ { k } ) } \times \frac { P ( { \mathbf { X } } ^ { k + 1 } , \left. ) } { 1 - K ( { \mathbf { X } } ^ { k + 1 } , \right. ) } \right\} . } \end{array}\tag{A4}
$$

For two intermediate transitions instead the transition probabilities are

$$
\begin{array} { r c l } { { \mathcal { T } ^ { 2 } ( { \bf X } ^ { k }  { \bf X } ^ { k + 2 } ) ~ = ~ { \cal P } ( { \bf X } ^ { k } ,  ) \times { \cal R } ^ { \tau } ( x _ { 2 M } ^ { k }  x _ { 2 M } ^ { k + 1 } ) \times } } \\ { { } } & { { \times ~ { \cal K } ( { \bf X } ^ { k + 1 } ,  ) \times { \cal R } ^ { \tau } ( x _ { 2 M } ^ { k + 1 }  x _ { 2 M } ^ { k + 2 } ) \times } } \\ { { } } & { { \times ~ [ 1 - K ( { \bf X } ^ { k + 2 } ,  ) ] , } } & { { ( \mathrm { A 5 } ) } } \end{array}
$$

and

$$
\begin{array} { r l r } { { \mathcal { T } } ^ { 2 } ( \mathbf { X } ^ { k + 2 } \to \mathbf { X } ^ { k } ) } & { = } & { P ( \mathbf { X } ^ { k + 2 } ,  ) \times R ^ { \tau } ( x _ { 0 } ^ { k + 2 } \to x _ { 0 } ^ { k + 1 } ) \times } \\ & { \times } & { K ( \mathbf { X } ^ { k + 1 } ,  ) \times R ^ { \tau } ( x _ { 0 } ^ { k + 1 } \to x _ { 0 } ^ { k } ) \times } \\ & { \times } & { [ 1 - K ( \mathbf { X } ^ { k } ,  ) ] , \quad \quad \quad \quad \quad \quad \mathrm { ( A 6 ) } } \end{array}
$$

leading to the acceptance factor

$$
\begin{array} { r l } & { A ^ { 2 } = \operatorname* { m i n } \bigg \{ 1 , \frac { 1 - K ( \mathbf { X } ^ { k } , \left. ) } { P ( \mathbf { X } ^ { k } , \right. ) } \times \frac { K ( \mathbf { X } ^ { k + 1 } , \left. ) } { K ( \mathbf { X } ^ { k + 1 } , \right. ) } \times } \\ & { \times \frac { w ( x _ { 2 M - 1 } ^ { k + 1 } ) } { w ( x _ { 1 } ^ { k + 1 } ) } \times \frac { P ( \mathbf { X } ^ { k + 2 } , \left. ) } { 1 - K ( \mathbf { X } ^ { k + 2 } , \right. ) } \times \frac { w ( x _ { 2 M - 1 } ^ { k + 2 } ) } { w ( x _ { 1 } ^ { k } ) } \bigg \} } \end{array}\tag{A7) .}
$$

The generalization to generic s intermediate steps is straightforward and can be written as

$$
\begin{array} { l } { { \displaystyle { \cal A } ^ { s } = \operatorname* { m i n } \left\{ 1 , \frac { 1 - K ( { \bf X } ^ { k } , \left. ) } { P ( { \bf X } ^ { k } , \right. ) } \times \frac { P ( { \bf X } ^ { k + s } , \left. ) } { 1 - K ( { \bf X } ^ { k + s } , \right. ) } \times \right. } } \\ { { \displaystyle \left. \times \frac { w ( x _ { 2 M - 1 } ^ { k + s } ) } { w ( x _ { 1 } ^ { k } ) } \times \left[ \prod _ { l = 1 } ^ { s - 1 } \frac { K ( { \bf X } ^ { k + l } , \left. ) } { K ( { \bf X } ^ { k + l } , \right. ) } \times \frac { w ( x _ { 2 M - 1 } ^ { k + l } ) } { w ( x _ { 1 } ^ { k + l } ) } \right] \right\} . } } \end{array}\tag{A8}
$$

To find a simple solution for the unknown probabilities, we first impose a cancellation for the intermediate acceptance factors, namely

$$
\frac { K ( \mathbf { X } , \left. ) } { K ( \mathbf { X } , \right. ) } = \frac { w ( x _ { 1 } ) } { w ( x _ { 2 M - 1 } ) } ,\tag{A9}
$$

this condition is satisfied by Eqs. (19) and (20). Then, we notice that the acceptance factor can be written only in terms of the final and the initial states as

$$
\begin{array} { r } { A ^ { s } = \operatorname* { m i n } \left\{ 1 , \frac { q ( \mathbf { X } ^ { k + s } , \left. ) } { q ( \mathbf { X } ^ { k } , \right. ) } \right\} . } \end{array}\tag{A10}
$$

Further, we can impose the two factors q to be independent on the direction, i.e., the condition $q ( \mathbf { X } ,  ) =$ $q ( \mathbf { X } ,  ) = q ( \mathbf { X } )$ , which is satisfied if

$$
\begin{array} { r l r } { \frac { P ( \mathbf { X } , \left. ) } { 1 - K ( \mathbf { X } , \right. ) } \ \times \ w ( x _ { 2 M - 1 } ) = } \\ { \ } & { = \ \frac { P ( \mathbf { x } , \right. ) } { 1 - K ( \mathbf { x } , \left. ) } \times w ( x _ { 1 } ) . } \end{array}\tag{A11}
$$

Since the two time directions are mutually exclusive, it is also true that $P ( \mathbf { X } , \longleftrightarrow ) + P ( \mathbf { X } , \to ) = 1$ , which allows us to solve Eq. (A11) and obtain Eqs. (16) and (17). The same reasoning can be repeated for the left direction and, due to imposed homogeneity for the probabilities, it can be checked that the detailed balance is satisfied for the left direction too.

## Appendix B: Bounce algorithm, directed updates, and efficiency

In this Appendix we comment on the relationship between the directed-update scheme and the bounce algorithm. If $\alpha = 1$ is taken in Eqs. (19) and (20), then after s updates along the direction $d ,$ at the end of the Markov step $P ( \mathbf { X } ^ { k + s } , \bar { d } ) = 0$ , i.e., the next Markov step will be taken in the opposite direction, just like the bounce algorithm. Although the two algorithms are similar in this particular limit, there is an important difference which eventually leads to a different computational efficiency. In order to elucidate this point and to show the α-dependence of the efficiency of the directed updates, we have done a systematic comparison of the two algorithms.

In particular, we have compared the efficiency of the directed updates with the bounce algorithm for a onedimensional Heisenberg model. The computational efficiency is generally defined as

$$
\mathcal { E } = \frac { 1 } { \sigma _ { O } ^ { 2 } T } ,\tag{B1}
$$

where $\sigma _ { O } ^ { 2 }$ is the square of the statistical error associated to a given observable after a given computational time T . In Fig. 7, we show the ratio between the directed-update scheme efficiency over the bounce algorithm efficiency, for the measurement of the ground-state energy of a onedimensional chain.

We notice that the two sampling schemes have comparable performances, being both based on a similar approach. As anticipated, it clearly emerges from Fig. 7 that the two algorithms do not have exactly the same behavior at $\alpha = 1$ , the maximum efficiency of the directed updates being reached for lower values of α. This feature is due to the fact that when α is very close to its saturation value, then a single Markov step can consist of a conspicuous number of individual “sliding moves”. Even if this situation leads to a fast decorrelation of configurations it also leads to a rarefaction of the possibility to measure the desired observables, which can eventually take place only at the end of the Markov step and not during the individual moves. This leads to a worse efficiency if compared to the bounce algorithm, where measurements can be in principle done after every sliding move.

<!-- image-->  
Figure 7: Relative efficiency of the directed update scheme and the bounce algorithm. The measured quantity is the ground-state energy of the one-dimensional Heisenberg model on a chain of size $L = 8 0$ sites.

1 D.M. Ceperley and E.L. Pollock, Phys. Rev. Lett. 56, 351 (1986).

2 D.M. Ceperley, Rev. Mod. Phys. 67, 279 (1995).

3 N.V. Prokof’ev, B.V. Svistunov, and I.S. Tupitsyn, Phys. Lett. A 238, 253 (1998).

4 M. Boninsegni, N. Prokof’ev, and B. Svistunov, Phys. Rev. Lett. 96, 070601 (2006).

5 S. Baroni and S. Moroni, Phys. Rev. Lett. 82, 4745 (1999).

6 A. Sarsa, K.E. Schmidt, and W.R. Magro, Journal of Chemical Physics 113, 1366 (2000).

7 C. Pierleoni and D.M. Ceperley, ChemPhysChem 6, 1872 (2005).

8 G. Carleo, S. Moroni, and S. Baroni, Phys. Rev. B 80, 094301 (2009).

9 W. Krauth, N. Trivedi, and D. Ceperley, Phys. Rev. Lett. 67, 2307 (1991).

10 O.F. Syljuasen, Phys. Rev. B 73, 245105 (2006).

11 O.F. Syljuasen and A.W. Sandvik, Phys. Rev. E 66, 046701 (2002).

12 V.G. Rousseau, Phys. Rev. E 78, 056707 (2008).

13 A.W. Sandvik, Phys. Rev. B 59, R14157 (1999).

14 S.M.A. Rombouts, K. Van Houcke, and L. Pollet, Phys.

In conclusion, the performances of the two algorithms are very close, although some advantages may arise from the use of the directed-updates. We further notice that the purely Markovian approach introduced in this paper could be slightly more efficient in cases where the number of rejected configurations by the bounce algorithm is substantial whereas all the generated configurations are accepted in the directed update scheme.

Rev. Lett 96, 180603 (2006)

15 S. Sorella and L. Capriotti, Phys. Rev. B 61, 2599 (2000).

16 K.H. Schmidt, P. Niyaz, A. Vaught, and M.A. Lee, Phys. Rev. E 71, 016707 (2005).

17 D.F.B. ten Haaf, H.J.M. van Bemmel, J.M.J. van Leeuwen, W. van Saarloos, and D.M. Ceperley, Phys. Rev. B 51, 13039 (1995).

18 R.R. Whitehead and A. Watt, J. Phys.G: Nucl. Phys. 4, 835 (1978).

19 See for example, T. Giamarchi, Quantum Physics in One Dimension (Oxford University Press, Oxford, 2004).

20 J.E. Gubernatis, M. Jarrell, R.N. Silver, and D.S. Silvia, Phys. Rev. B 44, 6011 (1991).

21 E. Manousakis, Rev. Mod. Phys. 63, 1 (1991).

22 F. Franjic and S. Sorella, Prog. Theor. Phys. 97, 399 (1997).

23 S. Sorella, Phys. Rev. B 71, 241103 (2005).

24 F. Becca, A. Parola, and S. Sorella, Phys. Rev. B 61, 16287(R) (2000).

25 S. Sorella, private communication.