# Usage

There are two main ways of using this package. You can either write your own Julia driver code (or any language that can inferface with Julia) in either a separate script or use an interactive environment such as the Julia REPL or Jupyter/Pluto. In this case, we are assuming you are using the REPL, but there should be minimal differences in using a notebook.

Initiate the Julia REPL by calling julia from the terminal or launching the Julia application. If you are using the CLI, please see this [link](https://docs.julialang.org/en/v1/manual/environment-variables/) to ensure that your PATH variables are set up properly.

The program assumes that you are in the `/project-cait/` directory. You can either `cd` from the REPL itself or launch it from that directory. You will then have to load the files by running

```
> include("src/simulator.jl")
```

## Compilation

Julia is just-in-time compiled -- which means that the system will only compile code when needed. Usually, this means that only parts of a program that you use are compiled when you call them. So the very first time you run any function from this package (in a fresh Julia session) it will take some time to compile but every use afterwards (until you close the session) will be much faster. Conside looking into [Revise.jl](https://github.com/timholy/Revise.jl) and related documentation online on how to manage Julia REPL sessions and compilation times.

## Dependencies

When you run the `include` code, Julia will load all of the externed `.jl` into the current scope. If there are scripts that require other packages, Julia will download and install the required dependencies from the project.toml if you run `activate .` and `instantiate` from inside the Julia [package manager](https://docs.julialang.org/en/v1/stdlib/Pkg/). However, I recommend just installing the dependencies manually for now. 

You can enter the `pkg` mode in the repl by typing `]` followed by a space. From there you should run the following commands:

```
> add DataStructures
> add Parameters
> add JLD2
> add StaticArrays
```

This should install and precompile the dependencies. It is probably a good idea to install them into a separate environment (similiar to Python's `virtualenv` or `conda` environments). See [here](https://pkgdocs.julialang.org/v1.2/environments/) for details.

## Implementation 

The `automata.jl` file defines two structs `Model` and `Cell`. The `Model` can contain a vector of cells as well as an array of functions that map to every reaction and an array of reaction effects. The code comes with a convenient function to instantiate the model with identical cells, but there is nothing stopping you from creating an initial state of your choosing. It also comes with the reaction-function mapping (`proprules`) and effects (`stoichiometry`) but it is possible to change these reaction dynamics or add your own.

The `Cell` contains all the necessary elements for for tracking its internal state (reactants, reaction rates, etc) as well as those for the Next Reaction Method such as internal time vectors, propensity vectors, priority queues for delayed and non-delayed reactions as well. For plotting, the `levels` vector contains the concentration of whichever reactant you want to measure (set to mh1 by default) and the `Time` vector contains the corresponding timestamp.

The dependency graph for the associated reactions have been generated and stored in `.jld2` file which is loaded automatically. As of right now there is no way to dynamically generate the dependency graph. You can use the `JLD2` package to open and decrypt the graph if you want to.

The `nrm` file contains modular functions for the Next Reaction Method algorithm and `simulator` contains driver code that runs the algorithm itself. Beyond that, please see the associated papers and docstrings in the source files: they are pretty straightforward and map the pseudocode in the paper to Julia code. 

## Running

The main function in `simulator` is `nrm`. This function takes in the number of cells you want, the simulation time in seconds and a vector of reaction rates. The defaults are set to one cell for 0.05s. It returns a `Model` containing cells, and each cell contains their own internal vectors containing the levels of the reactans at each timestep. If you want to simulate two cells for 100s using the default rates, run

```
model = nrm(2,100)
```

## Plotting

Ensure you have the standard plotting library installed, and load it into your session by running

```
> using Plots
```

You can then access the cells in your model and then plot the fields you want. For example, the following code will print the graph for the first cell in the model.

```
> plot(model.cells[1].Time, model.cells[1].levels)
```

If you are running this on AWS servers, you need to set the backend to unicodeplots() in order to plot the graph on the terminal. On macOS, the plots are generated in GTK+ by default. 