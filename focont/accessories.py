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
DType = Any  # TypeVar("DType", bound=np.float64)
FOCArray = NDArray[DType]


def is_stable(type: str, evals: ArrayLike) -> bool:
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


def is_symmetric(a: NDArray[DType], rtol: float = 1e-05, atol: float = 1e-08) -> bool:
    return np.allclose(a, a.T, rtol=rtol, atol=atol)


def h2_norm(lti_mimo: Any) -> float:
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
    """General exception class for focont."""

    def __init__(self: Self, message: str = "An error occured.") -> None:
        self.message = message
        super(FocontError, self).__init__(self.message)
