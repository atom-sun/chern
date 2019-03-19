# chern
Chern number calculator for 2D materials.


In 

* Mathematica
* python (to appear)

## Usage

### Mathematica version

The Chern calculator is generated in chern.nb notebook. After then, a chern.m 
is generated for the next time use. One can either

* Evaluate chern.nb each time in current (working) kernal, or
* import chern from chern.m in the user's working notebook

to get the chern calculator in the working memory.

To be noticed:

* There are two args to be given in the chern[hk, discretized], one for the
Hamiltonian and one for the discretization in calculation (e.g., discretized=16
stands for 16-by-16 discretization in calculation).

* The Hamiltonian should be two-dimensional, in the sense that it should be either
    - Function that can accept two arguments (..#1..#2..&), or
    - Expression in terms of {kx, ky} (not a function, but when {kx, ky} replaced
    by numbers it should reduces to matrix with all entry numeric)

    These two form of Hamiltonian can both be accepted by our calculator.

* The Hamiltonian is not limited to two-band models (example: graphene, haldane etc.),
    but it should be Hermittian and periodic in (2\pi-\prod-\pi) brillioun zone.

* The Hamiltonian should not contain other unevaluated Symbols except {kx, ky},
    which means after we assign value to {kx, ky}, the Hamiltonian should be
    numeric in all its entries.

See Example [haldane_honeycomb](Mathematica/test_haldane.nb)
        
## Reference

The algorithm is given by [chern 2005](
            https://journals.jps.jp/doi/10.1143/JPSJ.74.1674) this paper. 
Credits to them.


## Authors

* Ning Sun
    - email: ningsun.atom@gmail.com
