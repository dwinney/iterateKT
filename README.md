# iterateKT
Solver for iterative solutions to general Omnes-Khuri-Treiman problems.
That is, solutions to a system of coupled integral equations for any number of single-variable analytic functions of the form:
```math
    F_i(s) = P_{n-1}(s) + \frac{s^n}{\pi} \int ds^\prime \, \frac{\text{disc }F_i(s^\prime)}{s^{\prime n} \, (s^\prime - s)}
```
satisfying the unitarity condition
```math
    \text{disc }F_i(s) =  \sin\delta_i(s) \, e^{-i\delta_i(s)} \left[ F_i(s) + \sum_{j} \int dt \,  K_{ij}(s,t) \,  F_j(t) \right] ~.
```
For maximum flexibility, the code is agnostic to any quantum numbers and only requires specifying the elastic phase shift $\delta_i(s)$ and kernel functions $K_{ij}(s,t)$ for each isobar considered.

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
- [`kinematics`](./src/kinematics.hpp) which contains all relevant information regarding masses of particles involved. So far the three final state particles must have the same mass. 
- [`amplitude`](./src/amplitude.hpp) which contains the information of how many isobars contribute to a specific process and how they are combined. This is where a user can implement different isospin/helicity amplitudes. 
- [`isobar`](./src/isobar.hpp) which contains all of the dynamical information of a single two-body subsystem. 

A typical script may look like this
```c++
// Specify decay masses
kinematics kin = new_kinematics(m_decay, m_pi);

// Specify amplitude structure (quantum numbers)
amplitude  amp = new_amplitude<my_amplitude>(kin, "My Decay Process");

// Specify each isobar and number of subtractions
amp->add_isobar<my_first_isobar> (n_subtractions);
amp->add_isobar<my_second_isobar>(m_subtractions);

// Iterate the KT equations N times
amp->iterate(n_iterations);

// Access all isobars
std::vector<isobar> isobars = amp->get_isobars();
// or an individual one
isobar first_isobar = amp->get_isobar(kId);

// Evaluate the ith iteration of the jth basis function
complex first_isobar->basis_function(i, j, s+IEP);
```

### Virtual functions
As illustrated above, `amplitude` and `isobar` are pointers to instances of an abstract template classes (`raw_amplitude` and `raw_isobar` respectively). These contain the following virtual functions which must be specified by the user in a derived class:

##### `raw_amplitude::prefactor_s(uint i, complex s, complex t, complex u)`

The full amplitude can be evaluated from all the isobars at arbitrary $s$, $t$, and $u$ by evaluating:
```math
\mathcal{A}(s,t,u) = \sum_i \left[P^i_s(s,t,u) \, F_i(s) + P^i_t(s,t,u) \, F_i(t) + P^i_u(s,t,u)\, F_i(u) \right] ~.
```
The function $P_s^i$ is specified by the `prefactor_s`. Similarly $P_t^i$ and $P_u^i$ given by `prefactor_t` and `prefactor_u`. These provide the kinematic prefactors, isospin coefficients, and angular structure to the full amplitude which are irrelevant for individual isobars. 

##### `raw_isobar::id()` and `raw_isobar::name()`
Each isobar needs to be assigned an integer identifier and a string name with which it differentiate different isobars which get looped over. The int id is used internally while name provides a more human-readable identifier for command-line messages.

##### `raw_isobar::phase_shift(double s)`
The elastic phase shift $\delta(s)$ of a given isobar. This determes the Omnes function $\Omega(s)$ and therefore the initial guess for each isobar.

##### `raw_isobar::singularity_power()` and `raw_isobar::ksf_kernel(uint j, complex s, complex t)`
The last thing that must be specified is the kernel function $K_{ij}(s,t,u)$ which enters in the inhomogeneity of the KT equations. In order to avoid kinematic singularities, we instead specify the KSF kernel defined by
```math
    \hat{K}_{ij}(s,t) = \kappa^{n+1} \, K_{ij}(s,t) ~,
```
in terms of the Kacser function $\kappa$. The function `ksf_kernel(j, s, t)` then specifies $\hat{K}_{ij}(s,t)$ and `singularity_power()` returns the exponent $n$ (note the total power is $n+1$ with one factor always coming from the Jacobian of the angular integral).
