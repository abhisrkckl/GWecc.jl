GWecc.jl
========

[![CI](https://github.com/abhisrkckl/GWecc.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/abhisrkckl/GWecc.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/abhisrkckl/GWecc.jl/branch/main/graph/badge.svg?token=ZLIW2MYZ7L)](https://codecov.io/gh/abhisrkckl/GWecc.jl)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![arXiv](https://img.shields.io/badge/arXiv-2002.03285-b31b1b.svg)](https://arxiv.org/abs/2002.03285)
[![arXiv](https://img.shields.io/badge/arXiv-2210.11454-b31b1b.svg)](https://arxiv.org/abs/2210.11454)

`GWecc.jl` computes Pulsar Timing Array signals due to gravitational waves from eccentric supermassive binary sources.
This is a partial Julia rewrite of the C++ package [`GWecc`](https://github.com/abhisrkckl/GWecc).

This code is based on Susobhanan et al. 2020 and Susobhanan 2022. If you use this code in your work please cite the original articles

* Abhimanyu Susobhanan, Achamveedu Gopakumar, George Hobbs, and Stephen Taylor, 2020, *"Pulsar timing array signals induced by black hole binaries in relativistic eccentric orbits"*, Physical Review D, 101(4), 043022, DOI:[10.1103/PhysRevD.101.043022](https://doi.org/10.1103/PhysRevD.101.043022), arXiv:[2002.03285](https://arxiv.org/abs/2002.03285)
* Abhimanyu Susobhanan, 2022, *"Post-Newtonian-accurate pulsar timing array signals induced by inspiralling eccentric binaries: accuracy, computational cost, and single-pulsar search"*, Classical and Quantum Gravity, 40, 155014, DOI: [10.1088/1361-6382/ace234](https://doi.org/10.1088/1361-6382/ace234),  arXiv:[2210.11454](https://arxiv.org/abs/2210.11454)

See `citation.bib` for bibtex entries.

Installation
------------
I suggest installing this package in a new `conda` environment to avoid conflicts with existing package installations.
I also suggest installing `enterprise-pulsar` via conda before installing `GWecc.jl` to make sure that `tempo2`, which is
a dependency of `enterprise-pulsar`, is installed properly.

```
> conda create -n gwecc python=3.10
> conda activate gwecc
> conda install -c conda-forge enterprise-pulsar
> conda install -c conda-forge julia
```

`GWecc.jl` can be installed by typing the following in a Julia REPL:

```
> julia
julia> import Pkg
julia> Pkg.add(url="https://github.com/abhisrkckl/GWecc.jl.git")
julia> exit()
```

The Python wrapper can be installed by typing

```
> pip install git+https://github.com/abhisrkckl/GWecc.jl.git
```

Note that the Julia package should be installed before installing the Python interface.

Usage
-----
`GWecc.jl` is intended to be used used together with [`ENTERPRISE`](https://github.com/nanograv/enterprise/) package to search for eccentric supermassive binary sources and with [`libstempo`](https://github.com/vallis/libstempo/) to simulate such sources. Examples of such usage is given in the `examples/` directory (work in progress...). 
