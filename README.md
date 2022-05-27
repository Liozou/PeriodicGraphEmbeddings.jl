# PeriodicGraphEmbeddings

[![Build Status](https://github.com/Liozou/PeriodicGraphEmbeddings.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Liozou/PeriodicGraphEmbeddings.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://liozou.github.io/PeriodicGraphEmbeddings.jl/dev)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A Julia package for manipulating periodic graph embeddings in euclidean spaces, such as
those representing atoms in a crystal.

Partly wraps the [spglib](https://github.com/spglib/spglib) library to provide symmetry
detection for 3-periodic graphs.

See also:

- [PeriodicGraphs.jl](https://github.com/Liozou/PeriodicGraphs.jl) for the all that
  relates to the periodic graph itself, irrespective of its euclidean embedding.
- [CrystalNets.jl](https://github.com/coudertlab/CrystalNets.jl) for a dependent package
  specialized on crystal nets.
