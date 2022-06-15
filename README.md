# Focont

**Static output feedback and fixed order controller design package for Python**

## Static output feedback (SOF)

The SOF is the simplest feedback controller structure. It basically feedbacks
the system output to the system input after multiplying a constant gain matrix.

This package can calculate a stabilizing SOF gain which also optimizes the ![H2](/doc/h2.gif)
norm of the closed loop system.

However, this algorithm works under sufficient conditions. If the problem
parameters (listed below) is not appropriate, the algorithm fails and
prints an error message.

(See the article, https://journals.sagepub.com/doi/abs/10.1177/0142331220943071 ,
and the PhD thesis, http://hdl.handle.net/11693/54900 , for detailed
information and analysis)

The algorithm is purposedly developed for discrete time systems, but it also works
for continuous time systems when the SOF is calculated for the zero-order hold
discretized version with a sufficiently large sampling frequency.

Furthermore, the algorithm can be used to calculate fixed-order controllers.
This need additional entries in the described JSON or mat file.

Please, vizit [API docs](/doc/focont.md) for the detailed information.


## Installation

Create and activate a virtual environment (this step is not required but recommended).
Then,
```
pip install -r requirements.txt
```
or
```
pip install numpy scipy
```
Finally,
```
pip install .
```

### Also,

It can be installed with pip from `pypi`.
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

