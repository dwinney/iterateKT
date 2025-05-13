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

The driving term, $P_{n-1}(s)$, parameterizes the left-hand cuts associated with the production of the 3-body system. This function should therefore not contain right-hand cuts and is bounded by $s^{n-1}$. Traditionally, this is simply a polynomial of order $n-1$ but the code allows arbitrary functions to incorporate production effects. 

Note that convergence of the KT equations is not guaranteed! This may depend on the number of isobars, number of subtractions, masses and quantum numbers considered.

##  INSTALLATION

Compilation of the base library requires only [CMake](https://cmake.org/) (version $\geq$ 3.30), [ROOT](https://root.cern.ch/) (tested with version 6.24) with [*MathMore*](https://root.cern.ch/mathmore-library), and [Boost C++](https://www.boost.org/) (version $\geq$ 1.68).

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

The usual CMake installation may not work if you are running macOS. To circumvent this, the shell script `macos_build.sh` is included to get things working on a Mac. 

##  USAGE
The compiled executable pipes an analysis script, relevent header files, and the compiled library into ROOT's cling interpeter to run. 
This set up mimics a Python-like environment without requiring recompilation of the whole library when changes are made to amplitude files. To run a script located in the bin directory simply run 
```bash
iterateKT my_script.cpp
```
or add the bin directory to $PATH to call `iterateKT` from any directory. 

The classes of interest are:
- [`kinematics`](./src/kinematics.hpp) contains all relevant information regarding the masses of particles involved and kinematic quantities. So far, the three final state particles must have the same mass. 
- [`amplitude`](./src/amplitude.hpp) acts as a container class which specifies how different isobars contribute to a specific process and how to combine them to a full amplitude in terms of all Mandelstam variables.
- [`isobar`](./src/isobar.hpp) is the main physics object as it reconstructs two-particle subsystems in terms of basis functions after arbitrary iterations of the KT equations.

A typical script may look like this
```c++
// Specify decay masses
kinematics kin = new_kinematics(m_decay, m_final_state);

// Specify amplitude structure (quantum numbers)
amplitude  amp = new_amplitude<my_amplitude>(kin);

// Specify each isobar and number of subtractions
// Total number of basis functions (per isobar) will be (i+j+k)
amp->add_isobar<first_isobar> (i, id::first_isobar);
amp->add_isobar<second_isobar>(j, id::second_isobar);
amp->add_isobar<third_isobar> (k, id::third_isobar);

// Iterate the KT equations
amp->iterate(N);

// Access all isobars
std::vector<isobar> isobars = amp->get_isobars();
// or an individual one
isobar first_isobar = amp->get_isobar(id::first_isobar);

// Evaluate the lth basis function above and below cut
print("above", first_isobar->basis_function(l, s+IEPS));
print("below", first_isobar->basis_function(l, s-IEPS));

// Evaluate the full amplitude
print(amp->evaluate(s, t, u));
```

### Virtual functions
As illustrated above, `isobar` is a pointer to an instance of an abstract template class ( `raw_isobar`). The following virtual functions which must be implemented by the user in a derived class in order to specify the physics case of interest:

##### `double raw_isobar::phase_shift(double s)`
The elastic phase shift $\delta_i(s)$ fully determines the Omnes function $\Omega_i(s)$ and therefore the initial guess for each isobar.

##### `complex raw_isobar::ksf_kernel(uint j, complex s, complex t)` and `uint raw_isobar::angular_momentum()`
The kernel function $K_{ij}(s,t)$ which enters in the inhomogeneity of the KT equations. In order to avoid kinematic singularities, we actually specify the KSF kernel defined by
```math
    \hat{K}_{ij}(s,t) = \kappa^{2j_i+1} \, K_{ij}(s,t) ~,
```
in terms of the Kacser function $\kappa$. The function `ksf_kernel(j, s, t)` then specifies $\hat{K}_{ij}(s,t)$ and `angular_momentum()` returns the exponent $j_i$ which is specified by the spin-projection of the 2-body state (note the total power is $2j_i+1$ with one factor always coming from the Jacobian of the angular integral).

### Amplitudes
The above are sufficient if one is only interested in finding the basis functions which solve the KT equations. One may also combine isobars together using `amplitude` in the form:
```math
\mathcal{A}(s,t,u) = \sum_i \left[P^i_s(s,t,u) \, F_i(s) + P^i_t(s,t,u) \, F_i(t) + P^i_u(s,t,u)\, F_i(u) \right] ~,
```
for arbitrary complex $s$, $t$, and $u$. The function $P_s^i$ is specified by overriding  `raw_amplitude::prefactor_s(uint i, complex s, complex t, complex u)` and analogous functions for $P_t^i$ and $P_u^i$ (i.e. `prefactor_t` and `prefactor_u`). These can be used to provide any barrier factors, isospin coefficients, or angular structure which are irrelevant to solving the KT equations. 

From here one may calculate the double-differential decay width using `raw_amplitude::differential_width(double s, double t)`:
```math
\frac{dΓ}{ds\,dt} = \frac{1}{(2\pi)^3 \, 32 \, M^3} \frac{1}{\mathcal{N}} \, \left|\mathcal{A}(s,t,u)\right|^2 ~,
```
where $\mathcal{N}$ is a numerical factor specificied by `raw_amplitude::combinatorial_factor()` and can be used to add constants related to indentical particles and/or averaging over initial-state helicities. Single differential or fully integrated widths may also be accessed with `raw_amplitude::differential_width(double s)` and `raw_amplitude::width()`.  

### Plotting and Fitting
Many utilities are available to effectively fit amplitudes to data and plot the results. 
See documentation in [`fitter.hpp`](./src/fitter.hpp) and [`plotter.hpp`](./src/plotter.hpp) for details or the example scripts in [`/scripts`](./scripts/).