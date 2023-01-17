GWecc.jl
========

[![CI](https://github.com/abhisrkckl/GWecc.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/abhisrkckl/GWecc.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/abhisrkckl/GWecc.jl/branch/main/graph/badge.svg?token=ZLIW2MYZ7L)](https://codecov.io/gh/abhisrkckl/GWecc.jl)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


`GWecc.jl` computes pulsar TOA delays due to gravitational waves from eccentric supermassive binary sources.
This is a partial Julia rewrite of the C++ package [`GWecc`](https://github.com/abhisrkckl/GWecc).

This code is based on Susobhanan et al. 2020 and Susobhanan 2022. If you use this code in your work please cite the original articles

* Abhimanyu Susobhanan, Achamveedu Gopakumar, George Hobbs, and Stephen Taylor, 2020, *"Pulsar timing array signals induced by black hole binaries in relativistic eccentric orbits"*, Physical Review D, 101, 4, 043022, DOI:[10.1103/PhysRevD.101.043022](https://doi.org/10.1103/PhysRevD.101.043022), arXiv:[2002.03285](https://arxiv.org/abs/2002.03285)
* Abhimanyu Susobhanan, 2022, *"Post-Newtonian-accurate pulsar timing array signals induced by inspiralling eccentric binaries: accuracy and computational cost"*, arXiv e-prints, arXiv:[2210.11454](https://arxiv.org/abs/2210.11454)
 
Installation
------------
`GWecc.jl` can be installed by typing the following in a Julia REPL:

```
> import Pkg
> Pkg.add(url="https://github.com/abhisrkckl/GWecc.jl.git")
```

The Python wrapper can be installed by typing

```
$ pip install git+https://github.com/abhisrkckl/GWecc.jl.git
```

Note that the Julia package should be installed first.

Usage
-----
`GWecc.jl` is intended to be used used together with [`ENTERPRISE`](https://github.com/nanograv/enterprise/) package to search for eccentric supermassive binary sources and with [`libstempo`](https://github.com/vallis/libstempo/) to simulate such sources. Examples of such usage is given in the `examples/` directory (work in progress...). 
