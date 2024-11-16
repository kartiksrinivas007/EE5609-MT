## Linear System over Galois Field

---

This contains analysis of the solution space of a linear system over a four element Galois Field and implements an efficient gaussian elimination algorithm to solve it.
It also deduces whether the snumber of solutions is degenerate and reports the rank of the linear transform. Details are in `assignment-1.pdf`

The code is written in [Julia](https://julialang.org/).

### Prerequisites
---

Ijulia and Jupyter-Notebook needs to be installed.
Please install the `NbInclude` package using Julia to import Julia Notebook files

```bash
julia
import Pkg
Pkg.add("NBInclude")
Pkg.add("JLD2")
```
For running test scripts

```bash
julia evaluation.jl
```

