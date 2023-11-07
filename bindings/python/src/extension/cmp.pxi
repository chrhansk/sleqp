# cython: language_level=3

cpdef double inf():
    """
    Numerical value of infinity used
    by SLEQP
    """
    return csleqp.sleqp_infinity()
