# Chern Calculator
Chern number calculator for 2D materials.

In 
* Mathematica (v11.3)
* python (python3.7.2)

## Usage
### Mathematica version
The Chern calculator is generated in chern.nb notebook. After then, a chern.m 
is generated for the next time use. One can either
* Evaluate chern.nb each time in current (working) kernal, or
* import chern from chern.m in the user's working notebook

to get the chern calculator in the working memory. Then use in two ways:

* chern[hk]
* chern[hk, *discretized*]

To be noticed:

* It accepts one required arg + one optional arg in the calculator, i.e. 
chern[hk, *discretized*]
    - required: hk, the Hamiltonian
    - optional: *discretized*, the discretization to use in calculation, default to 16.
      (e.g., discretized=16 stands for 16-by-16 discretization of the Brillioun zone).

* The Hamiltonian should be two-dimensional, in the sense that it should be either
    - Function that can accept two arguments (..#1..#2..&), or
    - Expression in terms of {kx, ky} (not a function, but when {kx, ky} replaced
    by numbers it should reduces to matrix with all entries numeric).

    These two form of Hamiltonian can both be accepted by our calculator.

* The Hamiltonian is not limited to two-band models (example: graphene, haldane etc.),
    but it should be Hermittian and periodic in ($2\pi\times2\pi$) brillioun zone.

* The Hamiltonian should not contain other unevaluated Symbols except {kx, ky},
    which means that after we assign values to {kx, ky}, the Hamiltonian should be
    numeric in all its entries.

* It returns the list of Chern number for all bands, which should be **EXACT
INTEGERs**. If it's not, the model possibly be gapless --- which is, by the 
theory, not allowed to calculate the Chern for gapless bands in the first place.

#### Examples
* See example [easy_chern_insulator](examples/easy.m)
* See example [haldane_honeycomb](examples/haldane.m)

### python version
#### Installing and getting started
```bash
git clone https://github.com/atom-sun/chern.git
cd python/chern
pip install .
python examples/haldane.py
```
See also [haldane](examples/haldane.ipynb).

#### Uninstalling
```bash
pip uninstall chern
```
 
#### Examples
* Haldane honeycomb 2band model
```bash
python example/haldane.py
```
or
```bash
python examples/haldane.py --m 0.5 --phi 0
```
to specify a set of $$\{m, \phi\}$$ values of haldane model. 
Also in [example:haldane](examples/haldane.ipynb).
        
## Reference
The algorithm is given by [chern 2005](
            https://journals.jps.jp/doi/10.1143/JPSJ.74.1674) this paper. 
Credits to them.

## Authors
* Ning Sun
    - email: ningsun.atom@gmail.com
