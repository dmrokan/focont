# Focont

**Static output feedback and fixed order controller design package for Python**

## Static output feedback (SOF)

SOF is the simplest feedback controller structure. More precisely, it redirects the system output to the system input after multiplying by a constant gain matrix. You can find brief information about SOF in [this page](https://otomatik.art/content/static-feedback-calculator).

The algorithm implemented in this package can calculate a stabilizing SOF gain which also minimizes the $\mathcal{H}_2$ norm of the closed loop system. The resulting controller is comparable to the result obtained by linear quadratic regulators (LQR) with respect to the impulse response energy of the closed loop system.

However, this algorithm works when some sufficient conditions are satisfied. If the problem parameters (listed below) is not appropriate, the algorithm fails and prints an error message. Please see the [article](https://journals.sagepub.com/doi/abs/10.1177/0142331220943071), and the [PhD thesis](http://hdl.handle.net/11693/54900), for detailed information and analysis.

The algorithm is mainly developed for discrete time systems, but it may also compute similar SOF gains for continuous time systems when this algorithm is applied to the zero-order hold (ZOH) discretized versions with a sufficiently large sampling frequency.

Furthermore, the algorithm can be used to calculate fixed-order controllers. Please, check [tests](./tests/test_01.py) for examples and [docs](./docs/focont.md) for detailed information.


## Installation

```
python3 -m venv venv
source venv/bin/activate
pip install -e '.[dev]'
pytest
```

### Also,

It can be installed via pip from `pypi`.
```
pip install focont
```

## Example

```python
from focont import foc, system

pdata = system.load(json_or_mat_filename)
foc.solve(pdata)
foc.print_results(pdata)
```

You can find json and mat file examples in the `/tests` directory.

