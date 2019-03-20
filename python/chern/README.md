# Chern

Chern number calculator for 2D materials.


## Installing and getting started

```bash
git clone https://github.com/atom-sun/chern.git
cd python/chern
pip install .
python examples/haldane.py
```

See also [example](examples/example.ipynb).


## Uninstalling

```bash
pip uninstall chern
```
    
 
## Examples

**Haldane honeycomb 2band model**

```bash
python example/haldane.py
```

or

```bash
python examples/haldane.py --m 0.5 --phi 0
```

to specify a set of $$\{m, \phi\}$$ values of haldane model. 
Also in [example:haldane](examples/example.ipynb).


## Reference

The algorithm is given by [chern 2005](
            https://journals.jps.jp/doi/10.1143/JPSJ.74.1674) this paper. 
Credits to them.


## Authors

* Ning Sun
    - email: ningsun.atom@gmail.com
