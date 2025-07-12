# GasSteadyGeneralNetworkPartitioning.jl 

This package implements an algorithm that first sets up a general partition of the underlying gas network into subnetworks and then  solves the steady state gas flow problem by  Newton iterations where each Newton iteration involves the solution of the steady-state gas flow problem for each pipeline subnetwork to obtain the solution for the full network. The algorithm allows use of both ideal and non-ideal equations of state.

``GasSteadyGeneralNetworkPartitioning.jl`` is not a registered Julia package. Hence installation of the package should be done as follows:

```julia 
using Pkg
Pkg.add("https://github.com/me10240/GasSteadyGeneralNetworkPartitioning.jl.git")
```

For the API usage, users are referred to the ``test/`` and the ``examples/`` directories.

## Note
This is different from the hierarchical partitioning technique implemented in ``GasSteadyHierarchicallNetworkPartitioning.jl`` ([link](https://github.com/me10240/GasSteadyHierarchicalNetworkPartitioningSim.jl.git)).

