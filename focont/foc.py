import sys
import os
import json

import numpy as np
import scipy as sp
from scipy.linalg import svdvals
from scipy import signal

from focont.accessories import (FocontError, message, is_stable, warning,
                                h2_norm)


def _adequate_real(pdata):
    A = pdata['A']
    B = pdata['B']
    C = pdata['C']
    Q = pdata['Q']
    Q0 = pdata['Q0']
    R0 = np.eye(np.shape(C)[0])

    S0 = sp.linalg.solve_discrete_are(A.T, C.T, Q0, R0)

    rankS0 = np.linalg.matrix_rank(S0)
    if rankS0 != S0.shape[0]:
        raise RuntimeError("An adequate realization could not be found.")

    # S0 = v*e/v
    e, v = np.linalg.eigh(S0)
    e = np.sqrt(e)
    einv = np.divide(1.0, e)
    e = np.diag(e)
    einv = np.diag(einv)
    vinv = np.linalg.inv(v)

    # Similarity transformation
    Ta = np.matmul(v, np.matmul(e, vinv))
    Tainv = np.matmul(v, np.matmul(einv, vinv))

    Aa = np.matmul(Tainv, np.matmul(A, Ta))
    Ba = np.matmul(Tainv, B)
    Ca = np.matmul(C, Ta)
    Qa = np.matmul(Ta.T, np.matmul(Q, Ta))

    pdata['original_A'] = A
    pdata['original_B'] = B
    pdata['original_C'] = C
    pdata['original_Q'] = Q

    pdata['A'] = Aa
    pdata['B'] = Ba
    pdata['C'] = Ca
    pdata['Q'] = Qa
    pdata['Ta'] = Ta

    n = Aa.shape[0]
    I = np.eye(n)

    Cinv = np.linalg.pinv(Ca)
    Pi_c = np.matmul(Cinv, Ca)
    Pi_cbar = I - Pi_c
    Aa_cbar = np.matmul(Aa, Pi_cbar)

    svals = svdvals(Aa_cbar, overwrite_a=True)

    if sum(svals >= 1) > 0:
        raise FocontError("Projected system matrix has singular values "
                                            "greater than one\nfocont could not find a solution.")

    pdata['Pi_c'] = Pi_c
    pdata['Pi_cbar'] = Pi_cbar
    pdata['A_cbar'] = Aa_cbar
    pdata['Cinv'] = Cinv


def _closed_loop_FO_structure(pdata):
    if pdata['type'] == 'C':
        A = pdata['Aplant_discretized']
        B = pdata['Bplant_discretized']
    elif pdata['type'] == 'D':
        A = pdata['Aplant']
        B = pdata['Bplant']

    C = pdata['Cplant']

    m = B.shape[1]
    r = C.shape[0]

    # Controler dimension
    nc = pdata['controller_order']

    K = pdata['K']
    Acont = K[:, :nc]
    Bcont = K[:, nc:]
    Ccont = pdata['Ccont']
    Dcont = pdata['Dcont']

    A11 = A + np.matmul(B, np.matmul(Dcont, C))
    A12 = np.matmul(B, Ccont)
    A21 = np.matmul(Bcont, C)
    A22 = Acont

    Acl = np.block([
        [ A11, A12 ],
        [ A21, A22 ]
    ])

    ev_cl = np.linalg.eigvals(Acl)

    if not is_stable('D', ev_cl):
        raise FocontError('Resulting fixed order controller does not stabilize '
                          'the discrete time (or discretized) system.')

    Bcl = np.block([
        [ B ],
        [ np.zeros((nc, m)) ]
    ])

    Ccl = np.block([ C, np.zeros((r, nc)) ])

    closed_loop_lti = [ [] ]
    for i in range(r):
        for j in range(m):
            closed_loop_lti_ij = signal.dlti(
                Acl, Bcl[:, j:j+1], Ccl[i:i+1, :], np.zeros((1, 1))
            )

            closed_loop_lti[i] += [ closed_loop_lti_ij ]

        closed_loop_lti += [ [] ]

    if not closed_loop_lti[-1]:
        del closed_loop_lti[-1]

    pdata['closed_loop_lti'] = closed_loop_lti

    controller_lti = [ [] ]
    for i in range(m):
        for j in range(r):
            controller_lti_ij = signal.dlti(
                Acont, Bcont[:, j:j+1], Ccont[i:i+1, :], Dcont[i:i+1, j:j+1]
            )

            controller_lti[i] += [ controller_lti_ij ]

        controller_lti += [ [] ]

    if not controller_lti[-1]:
        del controller_lti[-1]

    pdata['controller_lti'] = controller_lti


def _closed_loop_system(pdata):
    K = pdata['K']

    A = pdata['original_A']
    B = pdata['original_B']
    C = pdata['original_C']

    Acld = A + np.matmul(B, np.matmul(K, C))
    ev_cl_d = np.linalg.eigvals(Acld)

    if not is_stable('D', ev_cl_d):
        raise FocontError("Closed loop discrete time system matrix "
                          "'A+BKC' is not stable.")

    pdata['Acld'] = Acld
    pdata['ev_cl_d'] = ev_cl_d

    if pdata['structure'] == 'SOF':
        if pdata['type'] == 'C':
            Ac = pdata['Ac']
            Bc = pdata['Bc']

            Aclc = Ac + np.matmul(Bc, np.matmul(K, C))
            ev_cl_c = np.linalg.eigvals(Aclc)

            if not is_stable(pdata['type'], ev_cl_c):
                raise FocontError("Closed loop continuous time system matrix "
                                  "'Ac+BcKCc' is not stable. "
                                  "Try to decrease sampling period.")

            pdata['Aclc'] = Aclc
            pdata['ev_cl_c'] = ev_cl_c

            m = Bc.shape[1]
            r = C.shape[0]

            closed_loop_lti = [ [] ]
            for i in range(r):
                for j in range(m):
                    closed_loop_lti_ij = signal.lti(
                        Aclc, Bc[:, j:j+1], C[i:i+1, :], np.zeros((1, 1))
                    )

                    closed_loop_lti[i] += [ closed_loop_lti_ij ]

                closed_loop_lti += [ [] ]

            if not closed_loop_lti[-1]:
                del closed_loop_lti[-1]

            pdata['closed_loop_lti'] = closed_loop_lti
        elif pdata['type'] == 'D':
            A = pdata['original_A']
            B = pdata['original_B']
            C = pdata['original_C']

            Acl = A + np.matmul(B, np.matmul(K, C))

            m = B.shape[1]
            r = C.shape[0]

            closed_loop_lti = [ [] ]
            for i in range(r):
                for j in range(m):
                    closed_loop_lti_ij = signal.dlti(
                        Acl, B[:, j:j+1], C[i:i+1, :], np.zeros((1, 1))
                    )

                    closed_loop_lti[i] += [ closed_loop_lti_ij ]

                closed_loop_lti += [ [] ]

            if not closed_loop_lti[-1]:
                del closed_loop_lti[-1]

            pdata['closed_loop_lti'] = closed_loop_lti
    else:
        _closed_loop_FO_structure(pdata)


def _calculate_sof(pdata):
    A = pdata['A']
    A_cbar = pdata['A_cbar']
    B = pdata['B']
    C = pdata['C']
    Q = pdata['Q']
    R = pdata['R']
    Cinv = pdata['Cinv']
    Pi_cbar = pdata['Pi_cbar']

    max_iter = pdata['max_iter']
    eps_conv = pdata['eps_conv']

    progress = 10
    progress_step = max_iter // 10

    inaccurate_result = False
    converged = False
    P = np.copy(Q)
    for i in range(max_iter):
        P_pre = np.copy(P)

        Rbar = np.matmul(B.T, np.matmul(P, B)) + R
        Rinv = np.linalg.inv(Rbar)

        M1 = np.matmul(P, np.matmul(B, np.matmul(Rinv, np.matmul(B.T, P))))
        M2 = np.matmul(A.T, np.matmul(P - M1, A))
        M3 = np.matmul(A_cbar.T, np.matmul(M1, A_cbar))

        P = Q + M2 + M3
        normP = np.linalg.norm(P)
        dP = np.linalg.norm(P - P_pre) / normP

        if dP < eps_conv:
            if np.isnan(dP) or np.isinf(dP):
                raise FocontError('Iterations did not converge.')
            else:
                message('Iterations converged, a solution is found')
                converged = True
                break

        if not inaccurate_result and normP * eps_conv > 1e2:
            warning('Cost-to-go is so large. Results can be inaccurate.')
            inaccurate_result = True

        if i % progress_step == 0:
            message('Progress:\t{}%, dP={}'.format(progress, dP))
            progress += 10

    if not converged:
        raise FocontError('Max iteration is reached but did not converge.\n'
                                            "Increase 'max_iter' or 'eps_conv' and try again.")

    F = -np.matmul(Rinv, np.matmul(B.T, np.matmul(P, A)))
    K = np.matmul(F, Cinv)

    pdata['P'] = P
    pdata['F'] = F
    pdata['K'] = K

    _closed_loop_system(pdata)


def _create_open_loop_lti(pdata):
    if pdata['type'] == 'C':
        lti_func = signal.lti
    elif pdata['type'] == 'D':
        lti_func = signal.dlti

    A = pdata['Aplant']
    B = pdata['Bplant']
    C = pdata['Cplant']

    r = C.shape[0]
    m = B.shape[1]

    open_loop_lti = [ [] ]
    for i in range(r):
        for j in range(m):
            open_loop_lti_ij = lti_func(
                A, B[:, j:j+1], C[i:i+1, :], np.zeros((1, 1))
            )

            open_loop_lti[i] += [ open_loop_lti_ij ]

        open_loop_lti += [ [] ]

    if not open_loop_lti[-1]:
        del open_loop_lti[-1]

    pdata['open_loop_lti'] = open_loop_lti


def solve(pdata):
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
    _create_open_loop_lti(pdata)
    _adequate_real(pdata)
    _calculate_sof(pdata)


def print_results(pdata):
    with np.printoptions(precision=4):
        if pdata['structure'] == 'SOF':
            message('Stabilizing SOF gain:', indent=1)
            print(pdata['K'])
            message('Eigenvalues of the closed loop system:', indent=1)
            if pdata['type'] == 'C':
                print(pdata['ev_cl_c'])
            elif pdata['type'] == 'D':
                print(pdata['ev_cl_d'])
                message('|e|:')
                print(np.abs(pdata['ev_cl_d']))
        elif pdata['structure'] == 'FO':
            nc = pdata['controller_order']
            K = pdata['K']
            Acont = K[:, :nc]
            Bcont = K[:, nc:]
            message('Acont:', indent=1)
            print(Acont)
            message('Bcont:', indent=1)
            print(Bcont)


def get_controller(pdata, i=-1, j=-1):
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
        return pdata['controller_lti']
    else:
        return pdata['controller_lti'][i][j]


def get_closed_loop_system(pdata, i=-1, j=-1):
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
        return pdata['closed_loop_lti']
    else:
        return pdata['closed_loop_lti'][i][j]


def norm(pdata, cl=True):
    r"""
    Calculates $\mathcal{H}_2$ norm of the closed or open loop
    MIMO system.

    :arg pdata dict: Problem data structure.
    :arg cl [TODO:type]: Calculate closed loop norm if it is `True`.

    :return float: $\mathcal{H}_2$ norm.
    """
    result = 0.0

    if cl:
        result = h2_norm(pdata['closed_loop_lti'])
    else:
        result = h2_norm(pdata['open_loop_lti'])

    return result


def h2_improvement(pdata):
    r"""
    Compares the $\mathcal{H}_2$ norms of the closed loop
    system obtained by the algortihm and the open loop system
    if the open loop system is also stable.

    :arg pdata dict: Problem data structure.

    :return float: Ratio of the closed and open loop $\mathcal{H}_2$ norms.
    """
    ol_stable = pdata['open_loop_stable']

    if not ol_stable:
        warning('Open loop system is not stable'
                'Can not calculate H2 norm improvement.')

        return 0.0

    ol_h2 = norm(pdata, False)
    print('Open loop H2 norm: {}'.format(ol_h2))

    cl_h2 = norm(pdata, True)
    print('Closed loop H2 norm: {}'.format(cl_h2))

    result = cl_h2 / ol_h2
    print('Improvement: {}'.format(result))

    return result

