# Focont

**Static output feedback and fixed order controller design package for Python**

## Static output feedback (SOF)

The SOF is the simplest feedback controller structure. It basically feedbacks
the system output to the system input after multiplying a constant gain matrix.

This package can calculate a stabilizing SOF gain which also optimizes the ![H2](/doc/h2.gif)
norm of the closed loop system.

(See the article, https://journals.sagepub.com/doi/abs/10.1177/0142331220943071 ,
and the PhD thesis, http://hdl.handle.net/11693/54900 , for detailed
information and analysis)

The algorithm is purposedly developed for discrete time systems, but it also works
for continuous time systems when the SOF is calculated for the zero-order hold
discretized version with a sufficiently large sampling frequency.

System definition can be provided by a JSON or MATLAB mat file. The required
entries in the file are:

* **A**: (2d array) System matrix
* **B**: (2d array) Input matrix
* **C**: (2d array) Output matrix.

If you want to adjust const function weights, these should also be provided:

* **Q**: (2d array) Cost function weight for the state variables (similar to LQR problem)
* **R**: (2d array) Cost function weight for the system input.

If (A, B, C) defines a continuous time system, these should also be provided:

* **type**: (string) 'C' (which means continuous, 'D' means discrete)
* **Ts**: (number) Sampling period

## Fixed order controller

Furthermore, the algorithm can be used to calculate fixed-order controllers.
This need additional entries in the described JSON or mat file;

* **structure**: (string) 'FO' (meaning fixed-order or 'SOF' for SOF controller)
* **controller\_order**: (integer) The order of controller (by default it is chosen as the number of inputs)

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

## Example

```python
from focont import foc, system

pdata = system.load(json_or_mat_filename)
foc.solve(pdata)
foc.print_results(pdata)
```

You can find json and mat file examples in the `/tests` directory.

