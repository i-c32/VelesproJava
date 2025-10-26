# VelesproJava
Quatum chemistry developed in Java.
## Input
```hocon
Velespro {
  water {
    coord = [
      {element:"H", x: 0.8668118,    y: 0.6014357,    z:-0.0000000}
      {element:"H", x:-0.8668118,    y: 0.6014357,    z:-0.0000000}
      {element:"O", x: 0.0000000,    y:-0.0757918,    z: 0.0000000}
    ]
    method = "HF"
    basis set = "STO-3G"
    charge = 0
    multiplicity = 1
    num_atom = 3
    option = {
      bohr = false
    }
  }
}
```

## Methodology
### Integrals
#### Overlap
M. Hô and J. M. Hernández-Pérez, “Evaluation of Gaussian Molecular Integrals,” The Mathematica Journal, 2012. dx.doi.org/doi:10.3888/tmj.14-3.