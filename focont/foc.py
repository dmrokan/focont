import sys
import os
import json

import numpy as np
import scipy as sp
from scipy.linalg import svdvals

from focont.accessories import FocontError, message, is_stable, warning


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
  m = pdata['Bplant'].shape[1]

  # Controler dimension
  nc = pdata['controller_order']

  K = pdata['K']
  Acont = K[:, :nc]
  Bcont = K[:, nc:]
  Ccont = pdata['Ccont']
  Dcont = pdata['Dcont']

  A = pdata['Aplant']
  B = pdata['Bplant']
  C = pdata['Cplant']

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


def _closed_loop_structure(pdata):
  s = pdata['structure']
  if s == 'FO':
    _closed_loop_FO_structure(pdata)
  elif s == 'SOF':
    pass
  else:
    raise FocontError("Undefined controller structure '{}'.".format(s))


def _closed_loop_system(pdata):
  K = pdata['K']

  A = pdata['original_A']
  B = pdata['original_B']
  C = pdata['original_C']

  Acld = A+np.matmul(B, np.matmul(K, C))
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

      Aclc = Ac+np.matmul(Bc, np.matmul(K, C))
      ev_cl_c = np.linalg.eigvals(Aclc)

      if not is_stable(pdata['type'], ev_cl_c):
        raise FocontError("Closed loop continuous time system matrix "
                          "'Ac+BcKCc' is not stable. "
                          "Try to decrease sampling period.")

      pdata['Aclc'] = Aclc
      pdata['ev_cl_c'] = ev_cl_c
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


def solve(pdata):
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
      m = pdata['Bplant'].shape[1]
      nc = pdata['controller_order']
      K = pdata['K']
      Acont = K[:, :nc]
      Bcont = K[:, nc:]
      message('Acont:', indent=1)
      print(Acont)
      message('Bcont:', indent=1)
      print(Bcont)

