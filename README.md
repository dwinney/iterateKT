# iterateKT
Solver for iterative solutions to general Omnes-Khuri-Treiman problems.
That is, solutions to the system of coupled integral equations involving any number of single-variable analytic functions of the form:
```math
    F_i(s) = P_{n-1}(s) + \frac{s^n}{\pi} \int ds^\prime \, \frac{\text{disc }F_i(s^\prime)}{s^{\prime n} \, (s^\prime - s)}
```
satisfying the unitarity condition
```math
    \text{disc }F_i(s) =  \sin\delta_i(s) \, e^{-i\delta_i(s)} \left[ F_i(s) + \sum_{j} \int dt \, K_{ij}(s,t) \, F_j(t) \right] ~.
```
For maximum flexibility, the code only requires specifying the elastic phase shift $\delta_i(s)$ and kernel functions $K_{ij}(s,t)$ of each isobar. Things such as isospin and/or helicity amplitudes can be built outside of the core iterative functionality by combining isobars into a full amplitude.

The driving term, $P_{n-1}(s)$, parameterizes the left-hand cuts associated with the production of the 3-body system. This function should therefore not contain right-hand cuts and is bounded by $s^{n-1}$. Traditionally, this is simply a polynomial of order $n-1$ but the code allows arbitrary functions to incorporate production effects. 

Note that convergence of the KT equations is not guaranteed! This may depend on the number of isobars, number of subtractions, masses and quantum numbers considered.

## Quick Start

If you have the prerequisites installed, you can build and run the project with:
```bash
git clone <repository-url>
cd iterateKT
mkdir build && cd build
cmake ..
make
make install
cd ..
./bin/iterateKT scripts/your_script.cpp
```

## Prerequisites

Before building, ensure you have:
- **CMake** (version ≥ 3.16)
- **ROOT** (tested with version 6.24) with MathMore library
- **Boost C++** (version ≥ 1.68) with system, filesystem, and math components

## Installation

To install, clone normally and use:
```bash
cd iterateKT
mkdir build && cd build
cmake ..
make
make install
```
This will create the core library in `lib/` (relative to project root) with the linkable library (`.so` on Linux, `.dylib` on macOS) as well as ROOT dictionary (`.pcm`) files. 

Additionally, a [scripting executable](./src/cling/iterateKT.cpp) will be installed into `bin/iterateKT` which shortcuts loading the libraries into an interactive ROOT session and running a .cpp file as an interpreted script.

The CMake build system automatically detects your operating system and handles platform-specific requirements. 

## USAGE
The compiled executable pipes an analysis script, relevant header files, and the compiled library into ROOT's cling interpreter to run. 
This setup mimics a Python-like environment without requiring recompilation of the whole library when changes are made to amplitude files. To run a script, you can either:
```bash
# Run from the project root directory
./bin/iterateKT my_script.cpp

# Or add the bin directory to your PATH
export PATH=$PATH:$(pwd)/bin
iterateKT my_script.cpp
``` 

The classes of interest are:
- [`kinematics`](./src/kinematics.hpp) contains all relevant information regarding the masses of particles involved and kinematic quantities. So far, the three final state particles must have the same mass. 
- [`amplitude`](./src/amplitude.hpp) acts as a container class which specifies how different isobars contribute to a specific process and how to combine them to a full amplitude in terms of all Mandelstam variables.
- [`isobar`](./src/isobar.hpp) is the main physics object as it reconstructs two-particle subsystems in terms of basis functions after arbitrary iterations of the KT equations.

A typical script may look like this:
```c++
// Specify decay masses (example: eta -> 3pi)
kinematics kin = new_kinematics(0.548, 0.140); // eta mass, pion mass

// Specify amplitude structure (quantum numbers)
amplitude amp = new_amplitude<my_amplitude>(kin);

// Specify each isobar and number of subtractions
// Total number of basis functions (per isobar) will be (i+j+k)
amp->add_isobar<first_isobar>(2, id::first_isobar);   // 2 subtractions
amp->add_isobar<second_isobar>(1, id::second_isobar); // 1 subtraction
amp->add_isobar<third_isobar>(1, id::third_isobar);   // 1 subtraction

// Iterate the KT equations (example: 3 iterations)
amp->iterate(3);

// Access all isobars
std::vector<isobar> isobars = amp->get_isobars();
// or an individual one
isobar first_isobar = amp->get_isobar(id::first_isobar);

// Evaluate the 0th basis function above and below cut
double s = 0.5; // example Mandelstam variable
double IEPS = 1e-6; // small imaginary part for analytic continuation
print("above", first_isobar->basis_function(0, s + IEPS));
print("below", first_isobar->basis_function(0, s - IEPS));

// Evaluate the full amplitude
print(amp->evaluate(s, 0.3, 0.2)); // example s, t, u values
```

### Virtual functions
As illustrated above, `isobar` is a pointer to an instance of an abstract template class (`raw_isobar`). The following virtual functions which must be implemented by the user in a derived class in order to specify the physics case of interest:

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
for arbitrary complex $s$, $t$, and $u$. The function $P_s^i$ is specified by overriding `raw_amplitude::prefactor_s(uint i, complex s, complex t, complex u)` and analogous functions for $P_t^i$ and $P_u^i$ (i.e. `prefactor_t` and `prefactor_u`). These can be used to provide any barrier factors, isospin coefficients, or angular structure which are irrelevant to solving the KT equations. 

From here one may calculate the double-differential decay width using `raw_amplitude::differential_width(double s, double t)`:
```math
\frac{dΓ}{ds\,dt} = \frac{1}{(2\pi)^3 \, 32 \, M^3} \frac{1}{\mathcal{N}} \, \left|\mathcal{A}(s,t,u)\right|^2 ~,
```
where $\mathcal{N}$ is a numerical factor specified by `raw_amplitude::combinatorial_factor()` and can be used to add constants related to identical particles and/or averaging over initial-state helicities. Single differential or fully integrated widths may also be accessed with `raw_amplitude::differential_width(double s)` and `raw_amplitude::width()`.  

### Plotting and Fitting
Many utilities are available to effectively fit amplitudes to data and plot the results. 
See documentation in [`fitter.hpp`](./src/fitter.hpp) and [`plotter.hpp`](./src/plotter.hpp) for details or the example scripts in [`/scripts`](./scripts/).

## Troubleshooting

**Common build issues:**

- **ROOT not found**: Ensure ROOT is installed and `root-config` is in your PATH
- **Boost not found**: Install Boost with `brew install boost` (macOS) or `sudo apt-get install libboost-all-dev` (Ubuntu)
- **CMake version too old**: Update CMake to version 3.16 or higher
- **Library loading errors**: Ensure the build completed successfully and `make install` was run

**Runtime issues:**

- **Script not found**: Check that the script path is correct relative to where you run `iterateKT`
- **Header files not found**: Ensure the project was built and installed correctly