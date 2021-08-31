<div align="center"><img src='logo.png' alt='Cait.jl Logo'></img></div>

# Project Cait

Code for stochastically simulating oscillatory systems in multicellular models of zebrafish. Current iterations include stable CPU and GPU versions written in Julia. There is also an experimental low-level version in CUDA C. Code by Tazmilur Saad and Ahmet Ay.

## Installation 

This package is not available on Julia's registries and must be installed manually. You can clone the repository using:

```
git clone https://github.com/syedtaz/project-cait.git
```

Currently, this package has been tested for Julia 1.5.0 but it should work on the latest releases. Please see this [link](https://github.com/syedtaz/project-cait.git) on how to install Julia. Once installed, you can launch the program by entering the Julia REPL and loading the package using

```
julia> include("src/simulator.jl")
```

This will load all of the functions of the package into the REPL scope. If there are dependencies that are not present on your system Julia's package manager will download and precompile them. Since Julia is just-in-time compiled, the first run of any function in the package will trigger compilation and thus will be relatively slower than subsequent runs.

Conside looking into [Revise.jl](https://github.com/timholy/Revise.jl) and documentation online on how to manage Julia REPL sessions and compilation times.

## Quickstart

The following code will create a 1-cell model with 100 seconds of simulation time using the default reactants and rate parameters.

```
julia> model = nrm(1,100)
```

You can plot this using Plots.jl. For example, the following code snippet will plot the cell contents vs time in the Julia REPL which is useful for plotting in remote environments.

```
using Plots
unicodeplots()
plot(model.cells[1].Time, model.cells[1].levels)
```

> Consult the documentation for more details on the algorithm and the underlying implementation.
