#!/usr/bin/python
#cython: language_level=3

cpdef double inf():
  return csleqp.sleqp_infinity()
