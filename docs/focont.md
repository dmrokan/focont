# focont package

## Submodules

## focont.accessories module

### *exception* focont.accessories.FocontError(message: str = 'An error occured.')

Bases: `Exception`

General exception object

### focont.accessories.convert_to_lti(A: ndarray[Any, dtype[Any]], B: ndarray[Any, dtype[Any]], C: ndarray[Any, dtype[Any]], D: ndarray[Any, dtype[Any]] = array([], shape=(0, 0), dtype=float64), t: str = 'D')

Create MIMO `scipy.signal.lti` or `scipy.signal.dlti` from state space matrices.

* **Parameters:**
  * **A** – State matrix $A in mathbb{R}^{n times n}$
  * **B** – Input matrix $B in mathbb{R}^{n times m}$
  * **C** – State matrix $C in mathbb{R}^{r times n}$
  * **D** – Input to output matrix $D in mathbb{R}^{r times m}$
  * **t** – Discrete (`"D"`) or continuous (`"C"`) time
* **Returns:**
  2D (r by n) list of `scipy.signal.lti` or `scipy.signal.dlti` instances.

### focont.accessories.freq_response(A: ndarray[Any, dtype[Any]], B: ndarray[Any, dtype[Any]], C: ndarray[Any, dtype[Any]], D: ndarray[Any, dtype[Any]], N: int, xscale: str = 'log')

Calculate the frequency response of discrete time MIMO LTI system.

* **Parameters:**
  * **A** – State matrix
  * **B** – Input matrix
  * **C** – Output matrix
  * **D** – Input to output matrix
  * **N** – Calculate the frequency response at `N` distinct points.
  * **xscale** – Frequency is in logarithmic (`"log"`) or linear (`"lin"`) scale
* **Returns:**
  A tuple of 2D complex `numpy.array`

The second entry in the returned tuple has a 2D list of frequency responses
and the first entry has a 2D list of corresponding frequencies.

### focont.accessories.h2_norm(lti_mimo: List[List[Any]])

Calculate the energy of impulse response of the LTI system.

* **Parameters:**
  **lti_mimo** – 2D list of LTI system in scipy.signal.lti or scipy.signal.dlti structure.
* **Returns:**
  2D list of impulse responses in numpy.array structre

### focont.accessories.is_stable(type: str, evals: ndarray[Any, dtype[Any]])

Check if the eigenvalues satisfies the system stability condition.

* **Parameters:**
  * **type** – Discrete (`"D"`) or continuous (`"C"`) time system
  * **evals** – Vector of eigenvalues of state matrix $A$
* **Returns:**
  `True` if stable

System is stable if

> - $|evals| < 1$ for discrete time
> - $Re{evals} < 0$ for continuous time

### focont.accessories.is_symmetric(a: ndarray[Any, dtype[Any]], atol: float = 1e-05, rtol: float = 1e-08)

Check if square matrix `a` is symmetric.

* **Parameters:**
  * **a** – Square numpy matrix
  * **atol** – Absolute tolerance
  * **rtol** – Relative tolernace
* **Returns:**
  `True` if symmetric.

Check numpy.allclose <https://numpy.org/devdocs/reference/generated/numpy.allclose.html#numpy.allclose>
for more information.

### focont.accessories.message(msg: str, indent: int = 0)

### focont.accessories.warning(msg: str, indent: int = 0)

## focont.foc module

### focont.foc.bode(pdata: Dict[str, Any], loop: str, N: int = 256, xscale: str = 'log', i: int = -1, j: int = -1)

Calculate closed or open loop frequency response of the system defined in the problem data.

* **Parameters:**
  * **pdata** – Problem data
  * **loop** – Open loop (`"open"`) or closed loop (`"closed"`) frequency response
  * **N** – calculate at `N` distinct frequency points
  * **xscale** – Frequency is in logarithmic (`"log"`) or linear (`"lin"`) scale
  * **i** – Calculate for `i` th output
  * **j** – Calculate for `j` th input
* **Returns:**
  A tuple of 2D complex `numpy.array`

The second entry in the returned tuple has a 2D list of frequency responses
and the first entry has a 2D list of corresponding frequencies.

### focont.foc.get_closed_loop_system(pdata: Dict[str, Any], i: int = -1, j: int = -1)

Returns the closed loop system in SciPy discrete LTI system form.

* **Parameters:**
  * **pdata** – Problem data structure.
  * **i** – Controller output index for the MIMO controller.
  * **j** – Controller input index for the MIMO controller.
* **Returns:**
  SciPy (discrete) LTI system representation.

Returns an m by r Python array when i or j is not provided.
The ith row and jth column of the return value gives the discrete LTI
system from jth input to the ith output.

### focont.foc.get_controller(pdata: Dict[str, Any], i: int = -1, j: int = -1)

Returns the controller in SciPy discrete LTI system form.

* **Parameters:**
  * **pdata** – Problem data structure.
  * **i** – Controller output index for the MIMO controller.
  * **j** – Controller input index for the MIMO controller.

Returns an m by r Python array when i or j is not provided.
The ith row and jth column of the return value gives the discrete LTI
system from jth input to the ith output.

### focont.foc.h2_improvement(pdata: Dict[str, Any])

Compares the $mathcal{H}_2$ norms of the closed loop
system obtained by the algortihm and the open loop system
if the open loop system is also stable.

* **Parameters:**
  **pdata** – Problem data structure.
* **Returns:**
  Ratio of the closed and open loop $mathcal{H}_2$ norms.

### focont.foc.norm(pdata: Dict[str, Any], cl: bool = True)

Calculates $mathcal{H}_2$ norm of the closed or open loop
MIMO system.

* **Parameters:**
  * **pdata** – Problem data structure.
  * **cl** – Calculate closed loop systems 2-norm if it is True.
* **Returns:**
  $mathcal{H}_2$ norm.

### focont.foc.print_results(pdata: Dict[str, Any])

Prints the resulting controller’s system matrices and closed loop system’s eigenvalues.

* **Parameters:**
  **pdata** – Problem data

### focont.foc.solve(pdata: Dict[str, Any])

Solves the SOF (static output feedback) or FOC (fixed order controller)
problem for the given LTI (discrete or continous) system by applying the
proposed solution method [1-2].

* **Parameters:**
  **pdata** – Python dictionary of problem parameters.

Controller is calculated by performing the following steps;

1. Find an appropriate realization of the LTI system.
2. Apply the approximate dyanmic programming (ADP) iterations to
   calculate the stabilizing controller which minimize a quadratic cost
   function similart to the well-known linear quadratic regulator (LQR)
   problem.

Please check [`focont.system.ProblemDataStructure`](#focont.system.ProblemDataStructure) for detailed
information on the problem data.

*NOTE*: Solution is appended to the input argument `pdata`.

[1]: Demir, O. and Özbay, H., 2020. Static output feedback stabilization
of discrete time linear time invariant systems based on approximate dynamic
programming. Transactions of the Institute of Measurement and Control,
42(16), pp.3168-3182.

[2]: Demir, O., 2020. Optimality based structured control of distributed
parameter systems (Doctoral dissertation, Bilkent University).

## focont.system module

### *class* focont.system.ProblemDataStructure

Bases: `object`

Problem data structure

Users can provide a file path or a Python dictionary to the
[`focont.system.load()`](#focont.system.load) function to define the problem.

If `problem_data` it is a `dict[str, Any]`, it must contain the definitions below.

#### A *: list[list[float]]*

System’s state matrix (required)

#### B *: list[list[float]]*

System’s input matrix (required)

#### C *: list[list[float]]*

System’s output matrix (required)

#### Ccont *: list[list[str]]*

The output matrix of controller.

Default value: $Q_{cont} = [ I ~~ 0 ]$

#### Dcont *: list[list[str]]*

The input to output matrix of controller. (required when `structre` is `"FO"`)

#### Q *: list[list[float]] | str*

Cost function weight on system’s state vector $x_t$

Variants:

> - Q must be symmetric semi-positive definite matrix.
> - `Q = "cI"` is equivalent to $Q = cI$

Default value: $Q = C^TC$

#### Q0 *: list[list[float]] | str*

Variants:

> - Q0 must be symmetric semi-positive definite matrix.
> - `Q0 = "cI"` is equivalent to $Q0 = cI$

Default value: $Q0 = I$

#### Q0cont *: list[list[str]]*

When the control structure is `"FO"`, the system state and input vectors are expanded by adding
new state variables. Therefore, the cost function weights `Q`, `Q0` and `R` should also
be expanded. Expanded versions are:

$$Q0_{extended} = diag{ Q_0, Q0_{cont} }$$

Default value: $Q0_{cont} = I$

#### Qcont *: list[list[str]]*

When the control structure is `"FO"`, the system state and input vectors are expanded by adding
new state variables. Therefore, the cost function weights `Q`, `Q0` and `R` should also
be expanded. Expanded versions are:

$$Q_{extended} = diag{ Q, Q_{cont} }$$

Default value: $Q_{cont} = I$

#### R *: list[list[float]] | str*

Cost function weight on system’s input vector $u_t$

Variants:

> - R must be symmetric positive definite matrix.
> - `R = "cI"` is equivalent to $R = cI$

Default value: $R = I$

#### R0 *: list[list[float]] | str*

Variants:

> - R0 must be symmetric positive definite matrix.
> - `R0 = "cI"` is equivalent to $R_0 = cI$

Default value: $R_0 = I$

#### Rcont *: list[list[str]]*

When the control structure is `"FO"`, the system state and input vectors are expanded by adding
new state variables. Therefore, the cost function weights `Q`, `Q0` and `R` should also
be expanded. Expanded versions are:

$$R_{extended} = diag{ R, R_{cont} }$$

Default value: $R_{cont} = I$

#### Ts *: float*

The discretization period in seconds.

Default value: 1e-2

#### eps_conv *: float*

Termination tolerance on the solver. Terminates if

$$||P-P_{pre}||_2^2/||P||_2^2 < eps_conv$$

Default value: 1e-12

#### max_iter *: int*

The maximum number of iterations allowed.

Default value: 1000000

#### structure *: str*

The resulting controller structure.

Variants:

> - “SOF”: Static output feedback controller.
> - “FO”: Fixed order controller.

#### type *: str*

System type, discrete (“D”) or continuous (“C”) time. If it is a continuous time
system it will be discretized first by applying ZOH discretization method.

Default value: “D”

#### zoh_calc_step *: int*

Number of time domain intervals used in ZOH discretization. ZOH is
calculated by a numerical integral when continuous time system’s
state matrix is singular.

Default value: 256

### focont.system.load(input_data: Dict[str, Any] | str)

Load Fixed Order Controller problem paramters from
a Python data structure or from a json, or mat file.

* **Parameters:**
  **input_data** – The source from which the problem parameters will be loaded.
* **Returns:**
  Dictionary of all required problem data to be provided foc.solve.

input_data can be json or mat filepath. In this case, file will
be read and problem parameters data structre will be created from
the json or mat file. Or, it can be dictionary of proplem parameters.

*NOTE*: Matrices must be Python array of array of floats with appropriate row
and column sizes (They are not numpy arrays!). Some matrices can be defined
as a string for ease of use. E.g:

C = ‘I’ or Q = ‘1e-2I’

They will be translated to numpy identity matrices, np.eye(n) and
1e-2 \* np.eye(n), where n is the dimension of LTI systems state vector.

Please check [`focont.system.ProblemDataStructure`](#focont.system.ProblemDataStructure) for detailed
information on the problem data.

### focont.system.load_from_json_file(json_filename: str)

### focont.system.load_from_mat_file(filename: str)

## Module contents
