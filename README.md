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
- [`kinematics`](./src/kinematics.hpp) which contains all relevant information regarding masses of particles involved. So far only decays into three equal mass particles is implemented. 
- [`amplitude`](./src/amplitude.hpp) which contains the information of how many isobars contribute to a specific process and how they are combined. This is where a user can implement different isospin/helicity amplitudes. 
- [`isobar`](./src/isobar.hpp) which contains all of the dynamical information of a single two-body subsystem. 