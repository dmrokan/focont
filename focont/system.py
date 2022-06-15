import sys
import os
import json
from copy import deepcopy

import numpy as np
import scipy as sp
from scipy.linalg import svdvals, expm

from focont import accessories as acc


def _parse_matrix_definition(d, n):
    if type(d) == str:
        d = d.strip()
        ind = d.find('I')
        if ind > -1:
            d = d.strip('I').strip()
            if d:
                g = float(d)
            else:
                g = 1.0

            return g * np.eye(n)

# TODO: Check if B and C are full rank
def _validate_input(pdata):
    # ** A
    if 'A' not in pdata:
        raise RuntimeError("System matrix 'A' is not defined.")
    else:
        n = len(pdata['A'])
        if n < 1:
            raise RuntimeError("System matrix 'A' has dimension zero.")
        else:
            if n != len(pdata['A'][0]):
                raise RuntimeError("System matrix 'A' is not square.")

    pdata['A'] = np.array(pdata['A'])

    # ** B
    if 'B' not in pdata:
        raise RuntimeError("System matrix 'B' is not defined.")
    else:
        n1 = len(pdata['B'])
        if n1 != n:
            raise RuntimeError("System matrices 'A' and 'B' must "
                                                  "have the same number of rows")
        m = len(pdata['B'][0])
        if m < 1:
            raise RuntimeError("System matrix 'B' has dimension zero.")

    pdata['B'] = np.array(pdata['B'])

    # ** C
    if 'C' not in pdata:
        raise RuntimeError("System matrix 'C' is not defined.")
    elif type(pdata['C']) == list:
        r = len(pdata['C'])
        if r < 1:
            raise RuntimeError("System matrix 'C' has dimension zero.")
        n2 = len(pdata['C'][0])
        if n2 != n:
            raise RuntimeError("System matrices 'A' and 'C' must "
                                                  "have the same number of columns.")

        pdata['C'] = np.array(pdata['C'])

    elif type(pdata['C']) == str:
        try:
            pdata['C'] = _parse_matrix_definition(pdata['C'], n)
        except:
            raise RuntimeError(
                "Matrix 'C={}' is an invalid definition.".format(pdata['C'])
            )

    # ** Q
    if 'Q' not in pdata:
        raise RuntimeError("Cost function weight 'Q' is not defined.")
    elif type(pdata['Q']) == list:
        n3 = len(pdata['Q'])
        if n3 < 1:
            raise RuntimeError("Cost function weight 'Q' has dimension zero.")
        elif n3 != n:
            raise RuntimeError("Matrices 'A' and 'Q' must "
                                                  "have the same number of columns and rows.")
        n4 = len(pdata['Q'][0])
        if n3 != n4:
            raise RuntimeError("Cost function weight 'Q' must be a square matrix.")

        pdata['Q'] = np.array(pdata['Q'])
    elif type(pdata['Q']) == str:
        if pdata['Q'] == 'default':
            pdata['Q'] = np.matmul(pdata['C'].T, pdata['C'])
        else:
            try:
                pdata['Q'] = _parse_matrix_definition(pdata['Q'], n)
            except:
                raise RuntimeError(
                    "Matrix 'Q={}' has an invalid definition.".format(pdata['Q'])
                )

    if not acc.is_symmetric(pdata['Q']):
        raise RuntimeError("Cost function weight 'Q' is not symmetric.")

    # ** R
    if 'R' not in pdata:
        raise RuntimeError("Cost function weight 'R' is not defined.")
    elif type(pdata['R']) == list:
        n5 = len(pdata['R'])
        if n5 < 1:
            raise RuntimeError("Cost function weight 'R' has dimension zero.")
        elif n5 != m:
            raise RuntimeError("Matrices 'B' and 'R' must "
                                                  "have the same number of columns.")
        n6 = len(pdata['R'][0])
        if n5 != n6:
            raise RuntimeError("Cost function weight 'R' must be a square matrix.")

        pdata['R'] = np.array(pdata['R'])
    elif type(pdata['R']) == str:
        try:
            pdata['R'] = _parse_matrix_definition(pdata['R'], m)
        except:
            raise RuntimeError(
                "Matrix 'R={}' has an invalid definition.".format(pdata['R'])
            )

    if not acc.is_symmetric(pdata['R']):
        raise RuntimeError("Cost function weight 'R' is not symmetric.")

    # ** Q
    if 'Q0' not in pdata:
        raise RuntimeError("Cost function weight 'Q0' is not defined.")
    elif type(pdata['Q0']) == list:
        n3 = len(pdata['Q0'])
        if n3 < 1:
            raise RuntimeError("Cost function weight 'Q0' has dimension zero.")
        elif n3 != n:
            raise RuntimeError("Matrices 'A' and 'Q0' must "
                                                  "have the same number of columns and rows.")
        n4 = len(pdata['Q0'][0])
        if n3 != n4:
            raise RuntimeError("Cost function weight 'Q0' must be a square matrix.")

        pdata['Q0'] = np.array(pdata['Q0'])
    elif type(pdata['Q0']) == str:
        if pdata['Q0'] == 'default':
            C = pdata['C']
            pdata['Q0'] = np.eye(C.shape[1])
        else:
            try:
                pdata['Q0'] = _parse_matrix_definition(pdata['Q0'], n)
            except:
                raise RuntimeError(
                    "Matrix 'Q0={}' has an invalid definition.".format(pdata['Q0'])
                )

    if not acc.is_symmetric(pdata['Q0']):
        raise RuntimeError("Cost function weight 'Q0' is not symmetric.")

    # ** type
    if 'type' not in pdata:
        raise RuntimeError("System type must be set as discrete (D) "
                                              "or continuous (C).")
    elif not (pdata['type'] == 'C' or pdata['type'] == 'D'):
        raise RuntimeError("System type must be set as discrete (D) "
                                              "or continuous (C).")


def _calculate_bw(pdata):
    om_list = np.logspace(-2, 5)

    A = pdata['A']
    B = pdata['B']
    C = pdata['C']

    I = np.eye(np.shape(A)[0])
    DCgain = None
    for om in om_list:
        s = om * 1j
        M = np.linalg.inv(s*I-A)
        P = np.matmul(C, np.matmul(M, B))

        max_sval = np.amax(svdvals(P, overwrite_a=True))

        if DCgain == None:
            DCgain = max_sval
        elif max_sval < 0.7079 * DCgain:
            return om / 2.0 / np.pi

    raise RuntimeError("BW could not be calculated.")


def _discretize(pdata):
    A = pdata['A']
    B = pdata['B']
    Ts = pdata['Ts']

    TsA = Ts * A
    Ad = expm(TsA)

    n = A.shape[0]
    rankA = np.linalg.matrix_rank(A)

    if rankA == n:
        I = np.eye(n)
        Ainv = np.linalg.inv(A)
        Mb = np.matmul(Ad - I, Ainv)

        Bd = np.matmul(Mb, B)
    else:
        zoh_calc_step = pdata['zoh_calc_step']
        zoh_Ts = Ts / zoh_calc_step

        Mb = np.eye(n)
        tau = zoh_Ts
        while tau < Ts:
            Atau = expm(tau * A)
            Mb += zoh_Ts * Atau

            tau += zoh_Ts

        Bd = np.matmul(Mb, B)

    pdata['Ac'] = A
    pdata['Bc'] = B

    pdata['A'] = Ad
    pdata['B'] = Bd

    pdata['Aplant_discretized'] = Ad
    pdata['Bplant_discretized'] = Bd


def _fo_controller_structure(pdata):
    from scipy.linalg import block_diag
    A = pdata['A']
    B = pdata['B']
    C = pdata['C']
    Q = pdata['Q']
    R = pdata['R']
    Q0 = pdata['Q0']

    nc = pdata['controller_order']
    Ccont = pdata['Ccont']
    Dcont = pdata['Dcont']
    Qcont = pdata['Qcont']
    Rcont = pdata['Rcont']
    Q0cont = pdata['Q0cont']

    n = A.shape[0]
    m = B.shape[1]
    r = C.shape[0]

    A11 = A + np.matmul(B, np.matmul(Dcont, C))
    A12 = np.matmul(B, Ccont)

    Aext = np.block([
        [ A11,               A12 ],
        [ np.zeros((nc, n + nc)) ]
    ])
    Bext = np.block([
        [ np.zeros((n, nc)) ],
        [ np.eye(nc)        ]
    ])
    Cext = np.block([
        [ np.zeros((nc, n)), np.eye(nc)        ],
        [ C,                 np.zeros((r, nc)) ]
    ])
    Qext = block_diag(Q, Qcont)
    Q0ext = block_diag(Q0, Q0cont)
    Rext = block_diag(Rcont)

    pdata['A'] = Aext
    pdata['B'] = Bext
    pdata['C'] = Cext
    pdata['Q'] = Qext
    pdata['R'] = Rext
    pdata['Q0'] = Q0ext


def _controller_structure(pdata):
    s = pdata['structure']
    if s == 'FO':
        m = pdata['B'].shape[1]
        if 'controller_order' not in pdata:
            pdata['controller_order'] = m

        nc = pdata['controller_order']

        if 'Ccont' not in pdata:
            Ccont = np.eye(m)
            if nc > m:
                Ccont = np.block([ Ccont, np.zeros((m, nc - m)) ])

            pdata['Ccont'] = Ccont

        rC, cC = pdata['Ccont'].shape
        if cC != nc or rC != m:
            raise FocontError(
                "The dimension of 'Ccont' must be {}x{}.".format(m, nc)
            )

        r = pdata['C'].shape[0]
        if 'Dcont' not in pdata:
            pdata['Dcont'] = np.zeros((m, r))

        rD, cD = pdata['Dcont'].shape
        if rD != m and cD != r:
            raise FocontError(
                "The dimension of 'Dcont' must be {}x{}.".format(m, r)
            )

        for mtx in [ 'Qcont', 'Q0cont', 'Rcont' ]:
            if mtx not in pdata:
                pdata[mtx] = np.eye(nc)
            elif type(pdata[mtx]) == str:
                pdata[mtx] = _parse_matrix_definition(pdata[mtx], nc)

        _fo_controller_structure(pdata)
    elif s == 'SOF':
        pass
    else:
        raise FocontError("Undefined controller structure '{}'.".format(s))


def load_from_json_file(json_filename):
    with open(json_filename, 'r') as fp:
        jobj = json.load(fp)

    return jobj


def load_from_mat_file(filename):
    from scipy.io import loadmat

    mat = loadmat(filename)

    for var_name in mat:
        if hasattr(mat[var_name], 'dtype') \
            and str(mat[var_name].dtype).startswith('<U'):
            mat[var_name] = str(mat[var_name][0])
        elif var_name == 'controller_order':
            mat[var_name] = int(mat[var_name])

    return mat


def load(input_data):
    r"""
    Load Fixed Order Controller problem paramters from
    a Python data structure or from a json, or mat file.

    :arg input_data dict_or_str: The source from which the problem
            parameters will be loaded.

    `input_data` can be json or mat filepath. In this case, file will
    be read and problem parameters data structre will be created from
    the json or mat file.

    *NOTE*: Matrices must be Python array of array of floats with appropriate row
    and column sizes (They are not `numpy` arrays!). Some matrices can be defined
    as a string for ease of use. E.g:
        `C = 'I'` or `Q = '1e-2I'`

        They will be translated to `numpy` identity matrices, `np.eye(n)` and
        `1e-2 * np.eye(n)`, where `n` is the dimension of LTI systems state vector.

    `focont` expects the following paramters:
        * `A`: System matrix of the LTI sytem ($A \in \mathbb{R}_{n \times n}$).
        * `B`: Input matrix ($B \in \mathbb{R}_{n \times m}$).
        * `C`: Output matrix (can be defined as a string, see the note above.)
                ($C \in \mathbb{R}_{r \times n}$).

        * `Q` (optional): Cost function weight for LTI system states (can be defined as a string.)
        `Q` must have the same dimension as `A` and must be symmetric and
        semi-positive definite. If it is not provided, its default value is 'I'.
        * `R` (optional): Cost function weight for LTI system's input. (can be a string.)
        `R` must be square and have the same number of columns as `B`. It must be 
        symmetric and positive definite. Its default value is 'I'.
        * `Q0' (optional): It has the same properties as `Q`, but it is used for
        calculating an appropriate realization of the LTI system as an intermediary
        step of the algorithm. Its default value is 'I'.
        * `type` (optional): It can be 'D' if the LTI system is discrete and 
        'C' if it is continuous. Its default value is 'D'.
        * `Ts` (optional): It is the sampling period used for ZOH discretization of
        the LTI system. Its default value is '0.01'.
        * `max_iter` (optional): Dynamic programming iterations limit. Its default
        value is '1e6'.
        * `eps_conv` (optional): Condition for convergence. If change in the cost-to-go
        function is smaller than this value, iterations will be terminated. Its
        default value is '1e-12'.
        * `zoh_calc_step` (optional): Max number of iterations used in ZOH discretization. Its
        default value is '256'.
        * `structure` (optional): It is 'SOF', if a static output feedback is wanted
        to be calculated. It is 'FO' if controller is dynamic.

        If controller structre is dynamic, then the paramters below can be provided.
            * `Ccont` (optional): Output matrix of the proposed dynamic controller.
            Its default value is '$I_{m \times m}$'.
            * `Dcont` (optional): Input to output gain of the controller. Its
            default value is '$0_{m \times r}$'.
            * `Qcont` (optional): Cost function weight on controller's state vector.
            Its default value is 'I'.
            * `Q0cont` (optional)
            * `Rcont` (optional): Cons function weight on controller's input vector.
            Its default value is 'I'.
    """
    if type(input_data) == str:
        filename = input_data
        _, ext = os.path.splitext(filename)

        if ext == '.json':
            pdata = load_from_json_file(filename)
        elif ext == '.mat':
            pdata = load_from_mat_file(filename)
    else:
        pdata = deepcopy(input_data)

    if 'type' not in pdata:
        pdata['type'] = 'D'

    if 'Q' not in pdata:
        pdata['Q'] = 'default'

    if 'R' not in pdata:
        pdata['R'] = 'I'

    if 'Q0' not in pdata:
        pdata['Q0'] = 'default'

    if 'Ts' not in pdata:
        pdata['Ts'] = 0.01

    if 'max_iter' not in pdata:
        pdata['max_iter'] = int(1e6)
    else:
        pdata['max_iter'] = int(pdata['max_iter'])

    if 'eps_conv' not in pdata:
        pdata['eps_conv'] = 1e-12

    if 'zoh_calc_step' not in pdata:
        pdata['zoh_calc_step'] = 256

    if 'structure' not in pdata:
        pdata['structure'] = 'SOF'

    _validate_input(pdata)

    stable = acc.is_stable(pdata['type'], np.linalg.eigvals(pdata['A']))
    pdata['open_loop_stable'] = stable

    pdata['Aplant'] = pdata['A']
    pdata['Bplant'] = pdata['B']
    pdata['Cplant'] = pdata['C']
    pdata['Qplant'] = pdata['Q']
    pdata['Rplant'] = pdata['R']
    pdata['Q0plant'] = pdata['Q0']

    if pdata['type'] == 'C':
        _discretize(pdata)

    _controller_structure(pdata)

    return pdata

