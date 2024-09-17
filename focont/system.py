import sys
import os
import json
from copy import deepcopy

from typing import (
    Dict,
    Any,
)

import numpy as np
import scipy as sp
from scipy.linalg import svdvals, expm

from . import accessories as acc
from .accessories import (
    FOCArray,
    FocontError,
)


PData = Dict[str, Any]


class ProblemDataStructure:
    """
    Problem data structure

    Users can provide a file path or a Python dictionary to the
    :py:meth:`focont.system.load` function to define the problem.

    If ``problem_data`` it is a ``dict[str, Any]``, it must contain the definitions below.
    """

    A: list[list[float]]
    """
    System's state matrix (required)
    """

    B: list[list[float]]
    """
    System's input matrix (required)
    """

    C: list[list[float]]
    """
    System's output matrix (required)
    """

    Q: list[list[float]] | str
    """
    Cost function weight on system's state vector :math:`x_t`

    Variants:

      - Q must be symmetric semi-positive definite matrix.
      - ``Q = "cI"`` is equivalent to :math:`Q = cI`

    Default value: :math:`Q = C^TC`
    """

    R: list[list[float]] | str
    """
    Cost function weight on system's input vector :math:`u_t`

    Variants:

      - R must be symmetric positive definite matrix.
      - ``R = "cI"`` is equivalent to :math:`R = cI`

    Default value: :math:`R = I`
    """

    type: str
    """
    System type, discrete ("D") or continuous ("C") time. If it is a continuous time
    system it will be discretized first by applying ZOH discretization method.

    Default value: "D"
    """

    Ts: float
    """
    The discretization period in seconds.

    Default value: 1e-2
    """

    max_iter: int
    """
    The maximum number of iterations allowed.

    Default value: 1000000
    """

    eps_conv: float
    """
    Termination tolerance on the solver. Terminates if

    .. math::

       ||P-P_{pre}||_2^2/||P||_2^2 < eps_conv

    Default value: 1e-12
    """

    structure: str
    """
    The resulting controller structure.

    Variants:

      - "SOF": Static output feedback controller.
      - "FO": Fixed order controller.
    """

    zoh_calc_step: int
    """
    Number of time domain intervals used in ZOH discretization. ZOH is
    calculated by a numerical integral when continuous time system's
    state matrix is singular.

    Default value: 256
    """

    Q0: list[list[float]] | str
    """
    Variants:

      - Q0 must be symmetric semi-positive definite matrix.
      - ``Q0 = "cI"`` is equivalent to :math:`Q0 = cI`

    Default value: :math:`Q0 = I`
    """

    R0: list[list[float]] | str
    """
    Variants:

      - R0 must be symmetric positive definite matrix.
      - ``R0 = "cI"`` is equivalent to :math:`R_0 = cI`

    Default value: :math:`R_0 = I`
    """

    Ccont: list[list[str]]
    """
    The output matrix of controller.

    Default value: :math:`Q_{cont} = \\[ I ~~ 0 \\]`
    """

    Dcont: list[list[str]]
    """
    The input to output matrix of controller. (required when :py:attr:`structre` is ``"FO"``)
    """

    Qcont: list[list[str]]
    """
    When the control structure is ``"FO"``, the system state and input vectors are expanded by adding
    new state variables. Therefore, the cost function weights ``Q``, ``Q0`` and ``R`` should also
    be expanded. Expanded versions are:

    .. math::

       Q_{extended} = \\textrm{diag}\\set{ Q, Q_{cont} }

    Default value: :math:`Q_{cont} = I`
    """

    Rcont: list[list[str]]
    """
    When the control structure is ``"FO"``, the system state and input vectors are expanded by adding
    new state variables. Therefore, the cost function weights ``Q``, ``Q0`` and ``R`` should also
    be expanded. Expanded versions are:

    .. math::

       R_{extended} = \\textrm{diag}\\set{ R, R_{cont} }

    Default value: :math:`R_{cont} = I`
    """

    Q0cont: list[list[str]]
    """
    When the control structure is ``"FO"``, the system state and input vectors are expanded by adding
    new state variables. Therefore, the cost function weights ``Q``, ``Q0`` and ``R`` should also
    be expanded. Expanded versions are:

    .. math::

       Q0_{extended} = \\textrm{diag}\\set{ Q_0, Q0_{cont} }

    Default value: :math:`Q0_{cont} = I`
    """


def _parse_matrix_definition(d: str, n: int) -> FOCArray | None:
    if type(d) == str:
        d = d.strip()
        ind = d.find("I")
        if ind > -1:
            d = d.strip("I").strip()
            if d:
                g = float(d)
            else:
                g = 1.0

            return g * np.eye(n)

    return None


# TODO: Check if B and C are full rank
def _validate_input(pdata: PData) -> None:
    # ** A
    if "A" not in pdata:
        raise RuntimeError("System matrix 'A' is not defined.")
    else:
        n = len(pdata["A"])
        if n < 1:
            raise RuntimeError("System matrix 'A' has dimension zero.")
        else:
            if n != len(pdata["A"][0]):
                raise RuntimeError("System matrix 'A' is not square.")

    pdata["A"] = np.array(pdata["A"])

    # ** B
    if "B" not in pdata:
        raise RuntimeError("System matrix 'B' is not defined.")
    else:
        n1 = len(pdata["B"])
        if n1 != n:
            raise RuntimeError(
                "System matrices 'A' and 'B' must " "have the same number of rows"
            )
        m = len(pdata["B"][0])
        if m < 1:
            raise RuntimeError("System matrix 'B' has dimension zero.")

    pdata["B"] = np.array(pdata["B"])

    # ** C
    if "C" not in pdata:
        raise RuntimeError("System matrix 'C' is not defined.")
    elif type(pdata["C"]) == list:
        r = len(pdata["C"])
        if r < 1:
            raise RuntimeError("System matrix 'C' has dimension zero.")
        n2 = len(pdata["C"][0])
        if n2 != n:
            raise RuntimeError(
                "System matrices 'A' and 'C' must " "have the same number of columns."
            )

        pdata["C"] = np.array(pdata["C"])
    elif type(pdata["C"]) == str:
        try:
            pdata["C"] = _parse_matrix_definition(pdata["C"], n)
        except:
            raise RuntimeError(
                "Matrix 'C={}' is an invalid definition.".format(pdata["C"])
            )

    # ** Q
    if "Q" not in pdata:
        raise RuntimeError("Cost function weight 'Q' is not defined.")
    elif type(pdata["Q"]) == list:
        n3 = len(pdata["Q"])
        if n3 < 1:
            raise RuntimeError("Cost function weight 'Q' has dimension zero.")
        elif n3 != n:
            raise RuntimeError(
                "Matrices 'A' and 'Q' must " "have the same number of columns and rows."
            )
        n4 = len(pdata["Q"][0])
        if n3 != n4:
            raise RuntimeError("Cost function weight 'Q' must be a square matrix.")

        pdata["Q"] = np.array(pdata["Q"])
    elif type(pdata["Q"]) == str:
        if pdata["Q"] == "default":
            pdata["Q"] = np.matmul(pdata["C"].T, pdata["C"])
        else:
            try:
                pdata["Q"] = _parse_matrix_definition(pdata["Q"], n)
            except:
                raise RuntimeError(
                    "Matrix 'Q={}' has an invalid definition.".format(pdata["Q"])
                )

    if not acc.is_symmetric(pdata["Q"]):
        raise RuntimeError("Cost function weight 'Q' is not symmetric.")

    # ** R
    if "R" not in pdata:
        raise RuntimeError("Cost function weight 'R' is not defined.")
    elif type(pdata["R"]) == list:
        n5 = len(pdata["R"])
        if n5 < 1:
            raise RuntimeError("Cost function weight 'R' has dimension zero.")
        elif n5 != m:
            raise RuntimeError(
                "Matrices 'B' and 'R' must " "have the same number of columns."
            )
        n6 = len(pdata["R"][0])
        if n5 != n6:
            raise RuntimeError("Cost function weight 'R' must be a square matrix.")

        pdata["R"] = np.array(pdata["R"])
    elif type(pdata["R"]) == str:
        try:
            pdata["R"] = _parse_matrix_definition(pdata["R"], m)
        except:
            raise RuntimeError(
                "Matrix 'R={}' has an invalid definition.".format(pdata["R"])
            )

    if not acc.is_symmetric(pdata["R"]):
        raise RuntimeError("Cost function weight 'R' is not symmetric.")

    # ** Q
    if "Q0" not in pdata:
        raise RuntimeError("Cost function weight 'Q0' is not defined.")
    elif type(pdata["Q0"]) == list:
        n3 = len(pdata["Q0"])
        if n3 < 1:
            raise RuntimeError("Cost function weight 'Q0' has dimension zero.")
        elif n3 != n:
            raise RuntimeError(
                "Matrices 'A' and 'Q0' must "
                "have the same number of columns and rows."
            )
        n4 = len(pdata["Q0"][0])
        if n3 != n4:
            raise RuntimeError("Cost function weight 'Q0' must be a square matrix.")

        pdata["Q0"] = np.array(pdata["Q0"])
    elif type(pdata["Q0"]) == str:
        if pdata["Q0"] == "default":
            C = pdata["C"]
            pdata["Q0"] = np.eye(C.shape[1])
        else:
            try:
                pdata["Q0"] = _parse_matrix_definition(pdata["Q0"], n)
            except:
                raise RuntimeError(
                    "Matrix 'Q0={}' has an invalid definition.".format(pdata["Q0"])
                )

    if not acc.is_symmetric(pdata["Q0"]):
        raise RuntimeError("Cost function weight 'Q0' is not symmetric.")

    # ** type
    if "type" not in pdata:
        raise RuntimeError(
            "System type must be set as discrete (D) " "or continuous (C)."
        )
    elif not (pdata["type"] == "C" or pdata["type"] == "D"):
        raise RuntimeError(
            "System type must be set as discrete (D) " "or continuous (C)."
        )


def _discretize(pdata: PData) -> None:
    A: FOCArray = pdata["A"]
    B: FOCArray = pdata["B"]
    Ts: FOCArray = pdata["Ts"]
    TsA: FOCArray = Ts * A
    Ad: FOCArray = expm(TsA)

    n: int = A.shape[0]
    rankA: int = np.linalg.matrix_rank(A)

    if rankA == n:
        I = np.eye(n)
        Ainv = np.linalg.inv(A)
        Mb = np.matmul(Ainv, Ad - I)

        Bd = np.matmul(Mb, B)
    else:
        zoh_calc_step = pdata["zoh_calc_step"]
        zoh_Ts = Ts / zoh_calc_step

        Mb = np.eye(n)
        tau = zoh_Ts
        while tau < Ts:
            Atau = expm(tau * A)
            Mb += zoh_Ts * Atau

            tau += zoh_Ts

        Bd = np.matmul(Mb, B)

    pdata["Ac"] = A
    pdata["Bc"] = B

    pdata["A"] = Ad
    pdata["B"] = Bd

    pdata["Aplant_discretized"] = Ad
    pdata["Bplant_discretized"] = Bd


def _fo_controller_structure(pdata: PData) -> None:
    from scipy.linalg import block_diag

    A: FOCArray = pdata["A"]
    B: FOCArray = pdata["B"]
    C: FOCArray = pdata["C"]
    Q: FOCArray = pdata["Q"]
    R: FOCArray = pdata["R"]
    Q0: FOCArray = pdata["Q0"]

    nc: int = pdata["controller_order"]
    Ccont: FOCArray = pdata["Ccont"]
    Dcont: FOCArray = pdata["Dcont"]
    Qcont: FOCArray = pdata["Qcont"]
    Rcont: FOCArray = pdata["Rcont"]
    Q0cont: FOCArray = pdata["Q0cont"]

    n: int = A.shape[0]
    m: int = B.shape[1]
    r: int = C.shape[0]

    A11: FOCArray = A + B @ Dcont @ C
    A12: FOCArray = B @ Ccont

    Aext: FOCArray = np.block([[A11, A12], [np.zeros((nc, n + nc))]])
    Bext: FOCArray = np.block([[np.zeros((n, nc))], [np.eye(nc)]])
    Cext: FOCArray = np.block([[np.zeros((nc, n)), np.eye(nc)], [C, np.zeros((r, nc))]])
    Qext: FOCArray = block_diag(Q, Qcont)
    Q0ext: FOCArray = block_diag(Q0, Q0cont)
    Rext: FOCArray = block_diag(Rcont)

    pdata["A"] = Aext
    pdata["B"] = Bext
    pdata["C"] = Cext
    pdata["Q"] = Qext
    pdata["R"] = Rext
    pdata["Q0"] = Q0ext


def _controller_structure(pdata: PData) -> None:
    s: str = pdata["structure"]

    if s == "FO":
        m = pdata["B"].shape[1]
        if "controller_order" not in pdata:
            pdata["controller_order"] = m

        nc = pdata["controller_order"]

        if "Ccont" not in pdata:
            Ccont = np.eye(m)
            if nc > m:
                Ccont = np.block([Ccont, np.zeros((m, nc - m))])

            pdata["Ccont"] = Ccont

        rC, cC = pdata["Ccont"].shape
        if cC != nc or rC != m:
            raise FocontError("The dimension of 'Ccont' must be {}x{}.".format(m, nc))

        r = pdata["C"].shape[0]
        if "Dcont" not in pdata:
            pdata["Dcont"] = np.zeros((m, r))

        rD, cD = pdata["Dcont"].shape
        if rD != m and cD != r:
            raise FocontError("The dimension of 'Dcont' must be {}x{}.".format(m, r))

        for mtx in ["Qcont", "Q0cont", "Rcont"]:
            if mtx not in pdata:
                pdata[mtx] = np.eye(nc)
            elif type(pdata[mtx]) == str:
                pdata[mtx] = _parse_matrix_definition(pdata[mtx], nc)

        _fo_controller_structure(pdata)
    elif s == "SOF":
        pass
    else:
        raise FocontError("Undefined controller structure '{}'.".format(s))


def load_from_json_file(json_filename: str) -> PData:
    with open(json_filename, "r") as fp:
        jobj: PData = json.load(fp)

    return jobj


def load_from_mat_file(filename: str) -> PData:
    from scipy.io import loadmat

    mat: PData = loadmat(filename)

    for var_name in mat:
        if hasattr(mat[var_name], "dtype") and str(mat[var_name].dtype).startswith(
            "<U"
        ):
            mat[var_name] = str(mat[var_name][0])
        elif var_name == "controller_order":
            mat[var_name] = int(mat[var_name].item(0))

    return mat


def load(input_data: PData | str) -> PData:
    """
    Load Fixed Order Controller problem paramters from
    a Python data structure or from a json, or mat file.

    :param input_data: The source from which the problem parameters will be loaded.
    :return: Dictionary of all required problem data to be provided `foc.solve`.

    `input_data` can be json or mat filepath. In this case, file will
    be read and problem parameters data structre will be created from
    the json or mat file. Or, it can be dictionary of proplem parameters.

    *NOTE*: Matrices must be Python array of array of floats with appropriate row
    and column sizes (They are not `numpy` arrays!). Some matrices can be defined
    as a string for ease of use. E.g:

    `C = 'I'` or `Q = '1e-2I'`

    They will be translated to `numpy` identity matrices, `np.eye(n)` and
    `1e-2 * np.eye(n)`, where `n` is the dimension of LTI systems state vector.

    Please check :py:class:`focont.system.ProblemDataStructure` for detailed
    information on the problem data.
    """

    pdata: PData = {}
    if type(input_data) == str:
        filename = input_data
        _, ext = os.path.splitext(filename)

        if ext == ".json":
            pdata = load_from_json_file(filename)
        elif ext == ".mat":
            pdata = load_from_mat_file(filename)
    else:
        assert isinstance(input_data, dict)
        pdata = deepcopy(input_data)

    if "type" not in pdata:
        pdata["type"] = "D"

    if "Q" not in pdata:
        pdata["Q"] = "default"

    if "R" not in pdata:
        pdata["R"] = "I"

    if "Q0" not in pdata:
        pdata["Q0"] = "default"

    if "Ts" not in pdata:
        pdata["Ts"] = 0.01

    if "max_iter" not in pdata:
        pdata["max_iter"] = int(1e6)
    else:
        pdata["max_iter"] = int(pdata["max_iter"])

    if "eps_conv" not in pdata:
        pdata["eps_conv"] = 1e-12

    if "zoh_calc_step" not in pdata:
        pdata["zoh_calc_step"] = 256

    if "structure" not in pdata:
        pdata["structure"] = "SOF"

    _validate_input(pdata)

    stable = acc.is_stable(pdata["type"], np.linalg.eigvals(pdata["A"]))
    pdata["open_loop_stable"] = stable

    pdata["Aplant"] = pdata["A"]
    pdata["Bplant"] = pdata["B"]
    pdata["Cplant"] = pdata["C"]
    pdata["Qplant"] = pdata["Q"]
    pdata["Rplant"] = pdata["R"]
    pdata["Q0plant"] = pdata["Q0"]

    if pdata["type"] == "C":
        _discretize(pdata)

    pdata["open_loop_lti"] = acc.convert_to_lti(pdata["A"], pdata["B"], pdata["C"])
    _controller_structure(pdata)

    pdata["pdata_initialized"] = True

    return pdata
