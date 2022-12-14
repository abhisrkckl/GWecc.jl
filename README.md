GWecc.jl
========

`GWecc.jl` computes pulsar TOA delays due to gravitational waves from eccentric supermassive binary sources.
This is a partial Julia rewrite of the C++ package [`GWecc`](https://github.com/abhisrkckl/GWecc).

This code is based on Susobhanan et al. 2020 and Susobhanan 2022. If you use this code in your work please cite the original articles

* Abhimanyu Susobhanan, Achamveedu Gopakumar, George Hobbs, and Stephen Taylor, 2020, *"Pulsar timing array signals induced by black hole binaries in relativistic eccentric orbits"*, Physical Review D, 101, 4,  043022, DOI: 10.1103/PhysRevD.101.043022, arXiv:2002.03285
* Abhimanyu Susobhanan, 2022, *"Post-Newtonian-accurate pulsar timing array signals induced by inspiralling eccentric binaries: accuracy and computational cost"*, arXiv e-prints, arXiv:2210.11454
 
Installation
------------
`GWecc.jl` can be installed by typing the following in a Julia REPL:

```
> import Pkg
> Pkg.add(url="https://github.com/abhisrkckl/GWecc.jl.git")
```

The `juliacall` package is also required to call `GWecc.jl` from Python. It can be installed by typing

```
$ pip install juliacall
```

Usage
-----
`GWecc.jl` is intended to be used used together with [`ENTERPRISE`](https://github.com/nanograv/enterprise/) package to search for eccentric supermassive binary sources and with [`libstempo`](https://github.com/vallis/libstempo/) to simulate such sources. Examples of such usage is given in the `examples/` directory (Work in progress). 
