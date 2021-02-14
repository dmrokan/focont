import numpy as np


def is_stable(type, evals):
  if type == 'C':
    if sum(np.real(evals) >= 0) > 0:
      return False
    else:
      return True
  elif type == 'D':
    if sum(np.abs(evals) >= 1) > 0:
      return False
    else:
      return True
  else:
    raise FocontError('Undefined system type \'{}\'.'.format(type))


def is_symmetric(a, rtol=1e-05, atol=1e-08):
  return np.allclose(a, a.T, rtol=rtol, atol=atol)


def message(msg, indent=0):
  print('-' * indent + ' ' + msg)


class FocontError(Exception):
  """General exception class for focont.
  """

  def __init__(self, message="An error occured."):
      self.message = message
      super(FocontError, self).__init__(self.message)

