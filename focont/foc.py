import sys
import os
import json

from typing import (
    Any,
    List,
    Tuple,
)

import numpy as np
import scipy as sp
from scipy.linalg import svdvals
from scipy import signal

from .accessories import (
    FocontError,
    message,
    is_stable,
    warning,
    h2_norm,
    convert_to_lti,
    freq_response,
    FOCArray,
)

from .system import PData


def _adequate_real(pdata: PData) -> None:
    A: FOCArray = pdata["A"]
    B: FOCArray = pdata["B"]
    C: FOCArray = pdata["C"]
    Q: FOCArray = pdata["Q"]
    Q0: FOCArray = pdata["Q0"]
    R0: FOCArray = np.eye(np.shape(C)[0])

    S0: FOCArray = sp.linalg.solve_discrete_are(A.T, C.T, Q0, R0)

    rankS0: int = np.linalg.matrix_rank(S0)
    if rankS0 != S0.shape[0]:
        raise RuntimeError("An adequate realization could not be found.")

    # S0 = v*e/v
    e: FOCArray
    v: FOCArray
    e, v = np.linalg.eigh(S0)
    e = np.sqrt(e)
    einv: FOCArray = np.divide(1.0, e)
    e = np.diag(e)
    einv = np.diag(einv)
    vinv: FOCArray = np.linalg.inv(v)

    # Similarity transformation
    Ta: FOCArray = v @ e @ vinv
    Tainv: FOCArray = v @ einv @ vinv

    Aa: FOCArray = Tainv @ A @ Ta
    Ba: FOCArray = Tainv @ B
    Ca: FOCArray = C @ Ta
    Qa: FOCArray = Ta.T @ Q @ Ta

    pdata["original_A"] = A
    pdata["original_B"] = B
    pdata["original_C"] = C
    pdata["original_Q"] = Q

    pdata["A"] = Aa
    pdata["B"] = Ba
    pdata["C"] = Ca
    pdata["Q"] = Qa
    pdata["Ta"] = Ta

    n: int = Aa.shape[0]
    I: FOCArray = np.eye(n)

    Cinv: FOCArray = np.linalg.pinv(Ca)
    Pi_c: FOCArray = Cinv @ Ca
    Pi_cbar: FOCArray = I - Pi_c
    Aa_cbar: FOCArray = Aa @ Pi_cbar

    svals: FOCArray = svdvals(Aa_cbar, overwrite_a=True)

    if np.any(svals >= 1):
        raise FocontError(
            "Projected system matrix has singular values "
            "greater than one\nfocont could not find a solution."
        )

    pdata["Pi_c"] = Pi_c
    pdata["Pi_cbar"] = Pi_cbar
    pdata["A_cbar"] = Aa_cbar
    pdata["Cinv"] = Cinv


def _closed_loop_FO_structure(pdata: PData) -> None:
    if pdata["type"] == "C":
        A = pdata["Aplant_discretized"]
        B = pdata["Bplant_discretized"]
    elif pdata["type"] == "D":
        A = pdata["Aplant"]
        B = pdata["Bplant"]

    C: FOCArray = pdata["Cplant"]

    m: int = B.shape[1]
    r: int = C.shape[0]

    # Controler dimension
    nc: int = pdata["controller_order"]

    K: FOCArray = pdata["K"]
    Acont: FOCArray = K[:, :nc]
    Bcont: FOCArray = K[:, nc:]
    Ccont: FOCArray = pdata["Ccont"]
    Dcont: FOCArray = pdata["Dcont"]

    A11: FOCArray = A + B @ Dcont @ C
    A12: FOCArray = B @ Ccont
    A21: FOCArray = Bcont @ C
    A22: FOCArray = Acont

    Acl: FOCArray = np.block([[A11, A12], [A21, A22]])

    ev_cl: FOCArray = np.linalg.eigvals(Acl)

    if not is_stable("D", ev_cl):
        raise FocontError(
            "Resulting fixed order controller does not stabilize "
            "the discrete time (or discretized) system."
        )

    Bcl: FOCArray = np.block([[B], [np.zeros((nc, m))]])

    Ccl: FOCArray = np.block([C, np.zeros((r, nc))])

    closed_loop_lti: Any = convert_to_lti(Acl, Bcl, Ccl)
    pdata["closed_loop_lti"] = closed_loop_lti

    controller_lti: Any = convert_to_lti(Acont, Bcont, Ccont, Dcont)
    pdata["controller_lti"] = controller_lti


def _closed_loop_system(pdata: PData) -> None:
    K: FOCArray = pdata["K"]

    A: FOCArray = pdata["original_A"]
    B: FOCArray = pdata["original_B"]
    C: FOCArray = pdata["original_C"]

    Acld: FOCArray = A + B @ K @ C
    ev_cl_d: FOCArray = np.linalg.eigvals(Acld)

    if not is_stable("D", ev_cl_d):
        raise FocontError(
            "Closed loop discrete time system matrix " "'A+BKC' is not stable."
        )

    pdata["Acld"] = Acld
    pdata["ev_cl_d"] = ev_cl_d

    if pdata["structure"] == "SOF":
        if pdata["type"] == "C":
            Ac = pdata["Ac"]
            Bc = pdata["Bc"]

            Aclc = Ac + np.matmul(Bc, np.matmul(K, C))
            ev_cl_c = np.linalg.eigvals(Aclc)

            if not is_stable(pdata["type"], ev_cl_c):
                raise FocontError(
                    "Closed loop continuous time system matrix "
                    "'Ac+BcKCc' is not stable. "
                    "Try to decrease sampling period."
                )

            pdata["Aclc"] = Aclc
            pdata["ev_cl_c"] = ev_cl_c

            m = Bc.shape[1]
            r = C.shape[0]

            closed_loop_lti: List[List[Any]] = [[]]
            for i in range(r):
                for j in range(m):
                    closed_loop_lti_ij = signal.lti(
                        Aclc, Bc[:, j : j + 1], C[i : i + 1, :], np.zeros((1, 1))
                    )

                    closed_loop_lti[i] += [closed_loop_lti_ij]

                closed_loop_lti += [[]]

            if not closed_loop_lti[-1]:
                del closed_loop_lti[-1]

            pdata["closed_loop_lti"] = closed_loop_lti
        elif pdata["type"] == "D":
            A = pdata["original_A"]
            B = pdata["original_B"]
            C = pdata["original_C"]

            Acl = A + np.matmul(B, np.matmul(K, C))

            m = B.shape[1]
            r = C.shape[0]

            closed_loop_lti = [[]]
            for i in range(r):
                for j in range(m):
                    closed_loop_lti_ij = signal.dlti(
                        Acl, B[:, j : j + 1], C[i : i + 1, :], np.zeros((1, 1))
                    )

                    closed_loop_lti[i] += [closed_loop_lti_ij]

                closed_loop_lti += [[]]

            if not closed_loop_lti[-1]:
                del closed_loop_lti[-1]

            pdata["closed_loop_lti"] = closed_loop_lti
    elif pdata["structure"] == "FO":
        _closed_loop_FO_structure(pdata)
    else:
        raise FocontError(
            "Unknown controller structure '{}'".format(pdata["structure"])
        )


def _calculate_sof(pdata: PData) -> None:
    A: FOCArray = pdata["A"]
    A_cbar: FOCArray = pdata["A_cbar"]
    B: FOCArray = pdata["B"]
    C: FOCArray = pdata["C"]
    Q: FOCArray = pdata["Q"]
    R: FOCArray = pdata["R"]
    Cinv: FOCArray = pdata["Cinv"]
    Pi_cbar: FOCArray = pdata["Pi_cbar"]

    max_iter: int = pdata["max_iter"]
    eps_conv: float = pdata["eps_conv"]

    progress: int = 10
    progress_step: int = max_iter // 10

    inaccurate_result: bool = False
    converged: bool = False
    P: FOCArray = np.copy(Q)

    for i in range(max_iter):
        P_pre = np.copy(P)

        Rbar = B.T @ P @ B + R
        Rinv = np.linalg.inv(Rbar)

        M1 = P @ B @ Rinv @ B.T @ P
        M2 = A.T @ (P - M1) @ A
        M3 = A_cbar.T @ M1 @ A_cbar

        P = Q + M2 + M3
        normP = np.linalg.norm(P)
        dP = np.linalg.norm(P - P_pre) / normP

        if dP < eps_conv:
            if np.isnan(dP) or np.isinf(dP):
                raise FocontError("Iterations did not converge.")
            else:
                message("Iterations converged, a solution is found")
                converged = True
                break

        if not inaccurate_result and normP * eps_conv > 1e2:
            warning("Cost-to-go is so large. Results can be inaccurate.")
            inaccurate_result = True

        if i % progress_step == 0:
            message("Progress:\t{}%, dP={}".format(progress, dP))
            progress += 10

    if not converged:
        raise FocontError(
            "Max iteration is reached but did not converge.\n"
            "Increase 'max_iter' or 'eps_conv' and try again."
        )

    F: FOCArray = -Rinv @ B.T @ P @ A
    K: FOCArray = F @ Cinv

    pdata["P"] = P
    pdata["F"] = F
    pdata["K"] = K

    _closed_loop_system(pdata)


def _create_open_loop_lti(pdata: PData) -> None:
    lti_func: Any
    if pdata["type"] == "C":
        lti_func = signal.lti
    elif pdata["type"] == "D":
        lti_func = signal.dlti

    A: FOCArray = pdata["Aplant"]
    B: FOCArray = pdata["Bplant"]
    C: FOCArray = pdata["Cplant"]

    r: int = C.shape[0]
    m: int = B.shape[1]

    open_loop_lti: List[List[Any]] = []
    for i in range(r):
        open_loop_lti += [[]]

        for j in range(m):
            open_loop_lti_ij = lti_func(
                A, B[:, j : j + 1], C[i : i + 1, :], np.zeros((1, 1))
            )

            open_loop_lti[i] += [open_loop_lti_ij]

    pdata["open_loop_lti"] = open_loop_lti


def _validate(pdata: PData) -> None:
    if "pdata_initialized" in pdata and pdata["pdata_initialized"]:
        pass
    else:
        raise FocontError(
            "Please, call 'system.load' function with "
            "your initial data set to generate problem "
            "data structure."
        )


def _invalidate(pdata: PData) -> None:
    pdata["pdata_initialized"] = False


def solve(pdata: PData) -> None:
    """
    Solves the SOF (static output feedback) or FOC (fixed order controller)
    problem for the given LTI (discrete or continous) system by applying the
    proposed solution method [1-2].

    [1]: Demir, O. and Özbay, H., 2020. Static output feedback stabilization
    of discrete time linear time invariant systems based on approximate dynamic
    programming. Transactions of the Institute of Measurement and Control,
    42(16), pp.3168-3182.

    [2]: Demir, O., 2020. Optimality based structured control of distributed
    parameter systems (Doctoral dissertation, Bilkent University).

    :arg pdata dict: Python dictionary of problem parameters obtained from
    `system.load` function of `focont` library.

    Controller is calculated by performing the following steps;

    1. Find an appropriate realization of the LTI system.
    2. Apply the approximate dyanmic programming (ADP) iterations to
    calculate the stabilizing controller which minimize a quadratic cost
    function similart to the well-known linear quadratic regulator (LQR)
    problem.

    *NOTE*: Solution is appended to the input argument `pdata`.
    """

    _validate(pdata)

    _create_open_loop_lti(pdata)
    _adequate_real(pdata)
    _calculate_sof(pdata)

    _invalidate(pdata)


def print_results(pdata: PData) -> None:
    with np.printoptions(precision=4):
        if pdata["structure"] == "SOF":
            message("Stabilizing SOF gain:", indent=1)
            print(pdata["K"])
            message("Eigenvalues of the closed loop system:", indent=1)
            if pdata["type"] == "C":
                print(pdata["ev_cl_c"])
            elif pdata["type"] == "D":
                print(pdata["ev_cl_d"])
                message("|e|:")
                print(np.abs(pdata["ev_cl_d"]))
        elif pdata["structure"] == "FO":
            nc = pdata["controller_order"]
            K = pdata["K"]
            Acont = K[:, :nc]
            Bcont = K[:, nc:]
            message("Acont:", indent=1)
            print(Acont)
            message("Bcont:", indent=1)
            print(Bcont)


def get_controller(pdata: PData, i: int = -1, j: int = -1) -> Any:
    """
    Returns the controller in SciPy discrete LTI system form.

    :arg pdata dict: Problem data structure.
    :arg i int: Controller output index for the MIMO controller.
    :arg j int: Controller input index for the MIMO controller.

    Returns an `m` by `r` Python array when `i` or `j` is not provided.
    The ith row and jth column of the return value gives the discrete LTI
    system from jth input to the ith output.
    """
    if i == -1 or j == -1:
        return pdata["controller_lti"]
    elif isinstance(pdata["controller_lti"], list):
        try:
            return pdata["controller_lti"][i][j]
        except IndexError:
            m = len(pdata["controller_lti"])
            if m > 0 and isinstance(pdata["controller_lti"][0], list):
                n = len(pdata["controller_lti"][0])
            else:
                return None
            raise FocontError(
                f"Controller dimension is {m}x{n}." f"Can not access ({i}, {j})"
            )
    else:
        return None


def get_closed_loop_system(pdata: PData, i: int = -1, j: int = -1) -> Any:
    """
    Returns the closed loop system in SciPy discrete LTI system form.

    :arg pdata dict: Problem data structure.
    :arg i int: Controller output index for the MIMO controller.
    :arg j int: Controller input index for the MIMO controller.

    :return scipy.signal.lti: SciPy (discrete) LTI system representation.

    Returns an `m` by `r` Python array when `i` or `j` is not provided.
    The ith row and jth column of the return value gives the discrete LTI
    system from jth input to the ith output.
    """
    if i == -1 or j == -1:
        return pdata["closed_loop_lti"]
    elif isinstance(pdata["closed_loop_lti"], list):
        try:
            return pdata["closed_loop_lti"][i][j]
        except IndexError:
            m = len(pdata["closed_loop_lti"])
            if m > 0 and isinstance(pdata["closed_loop_lti"][0], list):
                n = len(pdata["closed_loop_lti"][0])
            else:
                return None
            raise FocontError(
                f"Closed loop system dimension is {m}x{n}." f"Can not access ({i}, {j})"
            )
    else:
        return None


def norm(pdata: PData, cl: bool = True) -> float:
    r"""
    Calculates $\mathcal{H}_2$ norm of the closed or open loop
    MIMO system.

    :arg pdata dict: Problem data structure.
    :arg cl object: Calculate closed loop norm if it is `True`.

    :return float: $\mathcal{H}_2$ norm.
    """
    result: float = np.inf

    if cl:
        result = h2_norm(pdata["closed_loop_lti"])
    else:
        result = h2_norm(pdata["open_loop_lti"])

    return result


def h2_improvement(pdata: PData) -> float:
    r"""
    Compares the $\mathcal{H}_2$ norms of the closed loop
    system obtained by the algortihm and the open loop system
    if the open loop system is also stable.

    :arg pdata dict: Problem data structure.

    :return float: Ratio of the closed and open loop $\mathcal{H}_2$ norms.
    """
    ol_stable: bool = pdata["open_loop_stable"]

    if not ol_stable:
        warning(
            "Open loop system is not stable. " "Can not calculate H2 norm improvement."
        )

        return 0.0

    ol_h2: float = norm(pdata, False)
    print("Open loop H2 norm: {}".format(ol_h2))

    cl_h2: float = norm(pdata, True)
    print("Closed loop H2 norm: {}".format(cl_h2))

    result: float = cl_h2 / ol_h2
    print("Improvement: {}".format(result))

    return result


def bode(
    pdata: PData, loop: str, N: int = 256, xscale: str = "log", i: int = -1, j: int = -1
) -> Tuple[FOCArray, FOCArray]:
    lti: Any
    if loop == "open":
        lti = pdata["open_loop_lti"]
    elif loop == "closed":
        lti = pdata["closed_loop_lti"]

    r: int = len(lti)
    m: int = len(lti[0])
    A: FOCArray = lti[0][0].A
    B: FOCArray = np.zeros((A.shape[0], m))
    C: FOCArray = np.zeros((r, A.shape[0]))
    D: FOCArray = np.zeros((r, m))

    for k in range(r):
        C[k : k + 1, :] = lti[k][0].C

        for l in range(m):
            D[k : k + 1, l : l + 1] = lti[k][l].D

    for l in range(m):
        B[:, l : l + 1] = lti[0][l].B

    if i >= 0:
        C = C[i : i + 1, :]
        D = D[i : i + 1, :]

    if j >= 0:
        B = B[:, j : j + 1]
        D = D[:, j : j + 1]

    return freq_response(A, B, C, D, N)
