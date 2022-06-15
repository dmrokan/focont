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


def h2_norm(lti_mimo):
    impulse_responses = [ [] ]
    max_t = 0

    lti_siso = lti_mimo[0][0]
    is_D = True if lti_siso.dt else False

    r = len(lti_mimo)
    m = len(lti_mimo[0])

    h2 = 0

    for i in range(r):
        for j in range(m):
            lti = lti_mimo[i][j]

            _, y = lti.impulse()

            # NOTE: `impulse` function implementations in SciPy
            # are different for discret and continouse LTI.
            if is_D:
                y = y[0]

            impulse_responses[i] += [ y ]

            if y.shape[0] > max_t:
                max_t = y.shape[0]

        impulse_responses += [ [] ]

    for t in range(max_t):
        y_mimo = np.zeros((r, m))

        for i in range(r):
            for j in range(m):
                yt = 0
                resp = impulse_responses[i][j]
                if len(resp) > t:
                    yt = resp[t]

                y_mimo[i, j] = yt


        mag = np.linalg.norm(y_mimo)

        h2 += mag

    return h2


def message(msg, indent=0):
    print('-' * indent + ' ' + msg)


def warning(msg, indent=0):
    print('-' * indent + ' WARNING: ' + msg)


class FocontError(Exception):
    """General exception class for focont.
    """

    def __init__(self, message="An error occured."):
            self.message = message
            super(FocontError, self).__init__(self.message)

