# iterateKT
Solver for iterative solutions to general Omnes-Khuri-Treiman problems.
That is, solutions to the system of coupled integral equations involving any number of single-variable analytic functions of the form:
```math
    F_i(s) = P_{n-1}(s) + \frac{s^n}{\pi} \int ds^\prime \, \frac{\text{disc }F_i(s^\prime)}{s^{\prime n} \, (s^\prime - s)}
```
satisfying the unitarity condition
```math
    \text{disc }F_i(s) =  \sin\delta_i(s) \, e^{-i\delta_i(s)} \left[ F_i(s) + \sum_{j} \int dt \,  K_{ij}(s,t) \,  F_j(t) \right] ~.
```
For maximum flexibility, the code only requires specifying the elastic phase shift $\delta_i(s)$ and kernel functions $K_{ij}(s,t)$ of each isobar. Things such as isospin and/or helicity amplitudes can be built outside of the core iterative functionality by combining isobars into a full amplitude.

##  INSTALLATION

Compilation of the base library requires only [ROOT](https://root.cern.ch/) (tested with version 6.24) with [*MathMore*](https://root.cern.ch/mathmore-library) and [Boost C++](https://www.boost.org/) (version $\geq$ 1.68) libraries.

To install, clone normally and use:
```bash
cd iterateKT
mkdir build && cd build
cmake ..
cmake --build . --target install
```
This will create the core library `/lib/libITERATEKT.so` with the linkable library as well as ROOT dictionary (.pcm) files. 

Additionally a [scripting executable](./src/cling/iterateKT.cpp) will be installed into `/bin/iterateKT` which short-cuts loading the libraries into an interactive ROOT session and running a .cpp file as an interpreted script.   This executable requires the environment variable `ITERATEKT` to be set to the top-level directory in order to find auxilary files. This can be done as such:
```bash
export ITERATEKT=/path/to/iterateKT # for bash
setenv ITERATEKT /path/to/iterateKT # for csh
```

##  USAGE
The compiled executable pipes an analysis script, relevent header files, and the compiled library into ROOT's cling interpeter to run. 
This set up mimics a Python-like environment without requiring recompilation of the whole library when changes are made to amplitude files. To run a script located in the bin directory simply run 
```bash
iterateKT my_script.cpp
```
or add the bin directory to $PATH to call `iterateKT` from any directory. 

The main classes of interest are:
- [`kinematics`](./src/kinematics.hpp) This contains all relevant information regarding the masses of particles involved. So far, the three final state particles must have the same mass. 
- [`amplitude`](./src/amplitude.hpp) Acts as a container class which specifies how different isobars contribute to a specific process. This allows a user to implement different isospin/helicity amplitudes and calculate decay widths.
- [`isobar`](./src/isobar.hpp) This is the main physics object as it reconstructs two-particle subsystems in terms of basis functions after arbitrary iterations of the KT equations.

A typical script may look like this
```c++
// Specify decay masses
kinematics kin = new_kinematics(m_decay, m_final_state);

// Specify amplitude structure (quantum numbers)
amplitude  amp = new_amplitude<my_amplitude>(kin, "My Decay Process");

// Specify each isobar and number of subtractions
amp->add_isobar<my_first_isobar> (n_subtractions);
amp->add_isobar<my_second_isobar>(m_subtractions);

// Iterate the KT equations
amp->iterate(N_iterations);

// Access all isobars
std::vector<isobar> isobars = amp->get_isobars();
// or an individual one
isobar first_isobar = amp->get_isobar(id_of_first_isobar);

// Evaluate the ith iteration of the jth basis function
print(first_isobar->basis_function(i, j, s+IEPS));

// Evaluate the full amplitude
print(amp->evaluate(s, t, u));
```

### Virtual functions
As illustrated above, `isobar` is a pointer to an instance of an abstract template class ( `raw_isobar`). The following virtual functions which must be implemented by the user in a derived class in order to specify the physics case of interest:

##### `uint raw_isobar::id()`
Each isobar $F_i(s)$ needs to be assigned an integer $i$ to identify it, which should be unique among the set of isobars considered.

##### `double raw_isobar::phase_shift(double s)`
The elastic phase shift $\delta_i(s)$ fully determines the Omnes function $\Omega_i(s)$ and therefore the initial guess for each isobar.

##### `complex raw_isobar::ksf_kernel(uint j, complex s, complex t)` and `uint raw_isobar::singularity_power()`
The kernel function $K_{ij}(s,t)$ which enters in the inhomogeneity of the KT equations. In order to avoid kinematic singularities, we actually specify the KSF kernel defined by
```math
    \hat{K}_{ij}(s,t) = \kappa^{n_i+1} \, K_{ij}(s,t) ~,
```
in terms of the Kacser function $\kappa$. The function `ksf_kernel(j, s, t)` then specifies $\hat{K}_{ij}(s,t)$ and `singularity_power()` returns the exponent $n_i$ (note the total power is $n_i+1$ with one factor always coming from the Jacobian of the angular integral).

The above are sufficient if one is only interested in finding the basis functions which solve the KT equations. One may also combine isobars together using `amplitude` in the form:
```math
\mathcal{A}(s,t,u) = \sum_i \left[P^i_s(s,t,u) \, F_i(s) + P^i_t(s,t,u) \, F_i(t) + P^i_u(s,t,u)\, F_i(u) \right] ~,
```
for arbirary complex $s$, $t$, and $u$. The function $P_s^i$ is specified by overriding  `raw_amplitude::prefactor_s(uint i, complex s, complex t, complex u)` and analogous functions for $P_t^i$ and $P_u^i$ (i.e. `prefactor_t` and `prefactor_u`). These provide any barrier factors, isospin coefficients, and angular structure which are irrelevant for individual isobars. 

