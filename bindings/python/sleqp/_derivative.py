import numpy as np


EPS = np.finfo(float).eps

eps_twopoint = EPS**(1/2)
eps_threepoint = EPS**(1/3)

def perturbation(x, eps):
  return eps * np.maximum(np.abs(x), 1.)


def findiff_twopoint(f):

  def evaluate(x0, args=(), f0=None):
    x0 = np.atleast_1d(x0)
    h = perturbation(x0, eps_twopoint)

    if f0 is None:
      f0 = f(x0, *args)

    n = x0.size
    m = f0.size

    xd = np.copy(x0)

    df = np.empty((m, n))

    for j in range(n):
      xd[j] += h[j]

      fd = f(xd, *args)

      df[:, j] = (fd - f0) / h[j]

      xd[j] = x0[j]

    return df

  return evaluate


def findiff_threepoint(f):

  def evaluate(x0, args=(), f0=None):
    x0 = np.atleast_1d(x0)
    h = perturbation(x0, eps_threepoint)

    n = x0.size
    m = f0.size

    xd = np.copy(x0)

    df = np.empty((m, n))

    for j in range(n):
      xd[j] += h[j]

      fpos = f(xd, *args)

      xd[j] = x0[j]

      xd[j] -= h[j]

      fneg = f(xd, *args)

      xd[j] = x0[j]

      df[:, j] = (fpos - fneg) / (2 * h[j])

    return df

  return evaluate


def findiff_cs(f):
  def evaluate(x0, args=(), f0=None):
    x0 = np.atleast_1d(x0)
    h = perturbation(x0, eps_twopoint)

    n = x0.size
    m = f0.size

    df = np.empty((m, n))

    xd = x0.astype('complex')

    im = 1.j

    for j in range(n):
      xd[j] += h[j]*im

      df[:, j] = np.imag(f(xd)) / h[j]

      xd[j] = x0[j]

    return df

  return evaluate


def derivative(jac):
  def evaluate(x0, args=(), f0=None):
    return jac(x0)

  return evaluate
