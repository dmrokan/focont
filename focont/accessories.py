import numpy as np
from scipy import signal

from typing import (
    List,
    TypeVar,
    Any,
    Optional,
    Tuple,
    Self,
)

from numpy.typing import (
    ArrayLike,
    NDArray,
)

# TODO: Make this type equivalent to ndarray dtype
DType = Any
FOCArray = NDArray[DType]


def is_stable(type: str, evals: FOCArray) -> bool:
    """
    Check if the eigenvalues satisfies the system stability condition.

    :param type: Discrete (``"D"``) or continuous (``"C"``) time system
    :param evals: Vector of eigenvalues of state matrix $A$

    :return: ``True`` if stable

    System is stable if

      - $|evals| < 1$ for discrete time
      - $\\Re\\{evals\\} < 0$ for continuous time
    """

    if type == "C":
        if np.any(np.real(evals) >= 0):
            return False
        else:
            return True
    elif type == "D":
        if np.any(np.abs(evals) >= 1):
            return False
        else:
            return True
    else:
        raise FocontError("Undefined system type '{}'.".format(type))


def is_symmetric(a: NDArray[DType], atol: float = 1e-05, rtol: float = 1e-08) -> bool:
    """
    Check if square matrix ``a`` is symmetric.

    :param a: Square numpy matrix
    :param atol: Absolute tolerance
    :param rtol: Relative tolernace

    :return: ``True`` if symmetric.

    Check `numpy.allclose <https://numpy.org/devdocs/reference/generated/numpy.allclose.html#numpy.allclose>`
    for more information.
    """

    return np.allclose(a, a.T, rtol=rtol, atol=atol)


def h2_norm(lti_mimo: List[List[Any]]) -> float:
    """
    Calculate the energy of impulse response of the LTI system.

    :param lti_mimo: 2D list of LTI system in `scipy.signal.lti` or `scipy.signal.dlti` structure.

    :return: 2D list of impulse responses in `numpy.array` structre
    """

    impulse_responses: List[List[NDArray[DType]]] = [[]]
    max_t: int = 0

    lti_siso = lti_mimo[0][0]
    is_D: bool = lti_siso.dt is not None

    r: int = len(lti_mimo)
    m: int = len(lti_mimo[0])

    h2: float = 0

    for i in range(r):
        for j in range(m):
            lti = lti_mimo[i][j]

            _, y = lti.impulse()

            # NOTE: `impulse` function implementations in SciPy
            # are different for discret and continuous LTI.
            if is_D:
                y = y[0]

            impulse_responses[i] += [y]

            if y.shape[0] > max_t:
                max_t = y.shape[0]

        impulse_responses += [[]]

    for t in range(max_t):
        y_mimo = np.zeros((r, m))

        for i in range(r):
            for j in range(m):
                resp = impulse_responses[i][j]
                if len(resp) > t:
                    yt = resp[t]
                    y_mimo[i, j] = yt.item(0)

        mag = np.linalg.norm(y_mimo)

        h2 += mag.item(0)

    return h2


def convert_to_lti(
    A: FOCArray,
    B: FOCArray,
    C: FOCArray,
    D: FOCArray = np.zeros((0, 0)),
    t: str = "D",
) -> Any:
    """
    Create MIMO ``scipy.signal.lti`` or ``scipy.signal.dlti`` from state space matrices.

    :param A: State matrix $A \\in \\mathbb{R}^{n \\times n}$
    :param B: Input matrix $B \\in \\mathbb{R}^{n \\times m}$
    :param C: State matrix $C \\in \\mathbb{R}^{r \\times n}$
    :param D: Input to output matrix $D \\in \\mathbb{R}^{r \\times m}$
    :param t: Discrete (``"D"``) or continuous (``"C"``) time

    :return: 2D (r by n) list of ``scipy.signal.lti`` or ``scipy.signal.dlti`` instances.
    """

    n: int = A.shape[0]
    m: int = B.shape[1]
    r: int = C.shape[0]

    if D.shape[0] == 0:
        D = np.zeros((r, m))

    if t == "D":
        lti_func = signal.dlti
    elif t == "C":
        lti_func = signal.lti

    lti: List[List[Any]] = [[]]
    for i in range(m):
        for j in range(r):
            lti_ij = lti_func(
                A, B[:, j : j + 1], C[i : i + 1, :], D[i : i + 1, j : j + 1]
            )

            lti[i] += [lti_ij]

        lti += [[]]

    if not lti[-1]:
        del lti[-1]

    return lti


def freq_response(
    A: FOCArray, B: FOCArray, C: FOCArray, D: FOCArray, N: int, xscale: str = "log"
) -> Tuple[FOCArray, FOCArray]:
    """
    Calculate the frequency response of discrete time MIMO LTI system.

    :param A: State matrix
    :param B: Input matrix
    :param C: Output matrix
    :param D: Input to output matrix
    :param N: Calculate the frequency response at ``N`` distinct points.
    :param xscale: Frequency is in logarithmic (``"log"``) or linear (``"lin"``) scale

    :return: A tuple of 2D complex ``numpy.array``

    The second entry in the returned tuple has a 2D list of frequency responses
    and the first entry has a 2D list of corresponding frequencies.
    """

    f: FOCArray
    if xscale == "log":
        f = np.logspace(-3, 0, N)
    elif xscale == "lin":
        f = np.linspace(0, 1, N)

    n: int = A.shape[0]
    m: int = B.shape[1]
    r: int = C.shape[0]

    om: FOCArray = np.pi * f
    z: FOCArray = np.exp(1j * om)
    I: FOCArray = np.eye(n)

    resp: FOCArray = np.zeros((r, m, N)) + 1j * np.zeros((r, m, N))

    i: int = 0
    for zi in z:
        IA = zi * I - A

        IAinv = np.linalg.inv(IA)

        R = np.matmul(C, np.matmul(IAinv, B)) + D

        resp[:, :, i : i + 1] = R

        i += 1

    return f, resp


def message(msg: str, indent: int = 0) -> None:
    print("-" * indent + " " + msg)


def warning(msg: str, indent: int = 0) -> None:
    print("-" * indent + " WARNING: " + msg)


class FocontError(Exception):
    """General exception object"""

    def __init__(self: Self, message: str = "An error occured.") -> None:
        self.message = message
        super(FocontError, self).__init__(self.message)
