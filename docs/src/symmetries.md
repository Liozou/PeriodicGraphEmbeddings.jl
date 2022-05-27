# Symmetries

Symmetry detection is provided by the [spglib](https://github.com/spglib/spglib) library,
wrapped in helper functions detailed below.

## Manual

The main function is [`find_symmetries`](@ref) which returns a [`SymmetryGroup3D`](@ref):

```@docs
find_symmetries
SymmetryGroup3D
PeriodicGraphEmbeddings.PeriodicSymmetry3D
```

## Space group database API

```@docs
PeriodicGraphEmbeddings.SPACE_GROUP_HALL
PeriodicGraphEmbeddings.SPACE_GROUP_HM
PeriodicGraphEmbeddings.SPACE_GROUP_FULL
PeriodicGraphEmbeddings.SPACE_GROUP_IT
PeriodicGraphEmbeddings.HALL_SYMBOLS
```

## Internal API

```@docs
find_hall_number
PeriodicGraphEmbeddings.SpglibDataset
get_symmetry_equivalents
PeriodicGraphEmbeddings.get_spglib_dataset
check_valid_symmetry
```
