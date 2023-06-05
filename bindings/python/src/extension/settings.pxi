#cython: language_level=3

from enum import Enum


class _PropType(Enum):
  Bool = 0
  Integer = 1
  Enumerated = 2
  Real = 3


class _NullConverter:
  def from_py(self, val):
    return val

  def to_py(self, val):
    return val


class _EnumConverter:
  def __init__(self, enum_type):
    self.enum_type = enum_type

  def from_py(self, val):
    try:
      return val.value
    except AttributeError:
      pass

    try:
      return self.enum_type[val]
    except KeyError:
      raise AttributeError('Invalid attribute value \'{0}\''.format(val))

  def to_py(self, val):
    return self.enum_type(val)


class _Prop:
  def __init__(self, value, prop_type, converter):
    self.value = value
    self.prop_type = prop_type
    self.converter = converter

  def to_py(self, val):
    return self.converter.to_py(val)

  def from_py(self, val):
    return self.converter.from_py(val)

  @staticmethod
  def boolean(value):
    return _Prop(value, _PropType.Bool, _NullConverter())

  @staticmethod
  def integer(value):
    return _Prop(value, _PropType.Integer, _NullConverter())

  @staticmethod
  def enumerated(value, enum_type):
    return _Prop(value, _PropType.Enumerated, _EnumConverter(enum_type))

  @staticmethod
  def real(value):
    return _Prop(value, _PropType.Real, _NullConverter())


cdef dict opt_prop_map = {
  # Boolean properties
  'perform_newton_step':   _Prop.boolean(csleqp.SLEQP_SETTINGS_BOOL_PERFORM_NEWTON_STEP),
  'global_penalty_resets': _Prop.boolean(csleqp.SLEQP_SETTINGS_BOOL_GLOBAL_PENALTY_RESETS),
  'perform_soc':           _Prop.boolean(csleqp.SLEQP_SETTINGS_BOOL_PERFORM_SOC),
  'use_quadratic_model':   _Prop.boolean(csleqp.SLEQP_SETTINGS_BOOL_USE_QUADRATIC_MODEL),
  'enable_restoration':    _Prop.boolean(csleqp.SLEQP_SETTINGS_BOOL_ENABLE_RESTORATION_PHASE),
  'enable_preprocessor':   _Prop.boolean(csleqp.SLEQP_SETTINGS_BOOL_ENABLE_PREPROCESSOR),
  'lp_resolves':           _Prop.boolean(csleqp.SLEQP_SETTINGS_BOOL_LP_RESOLVES),

  # Integer properties
  'num_quasi_newton_iterates': _Prop.integer(csleqp.SLEQP_SETTINGS_INT_NUM_QUASI_NEWTON_ITERATES),
  'max_newton_iterations':     _Prop.integer(csleqp.SLEQP_SETTINGS_INT_MAX_NEWTON_ITERATIONS),
  'num_threads':               _Prop.integer(csleqp.SLEQP_SETTINGS_INT_NUM_THREADS),

  # SLEQP_SETTINGS_INT_FLOAT_WARNING_FLAGS,
  # SLEQP_SETTINGS_INT_FLOAT_ERROR_FLAGS,

  # Enumerated properties
  'deriv_check':          _Prop.enumerated(csleqp.SLEQP_SETTINGS_ENUM_DERIV_CHECK, DerivCheck),
  'hessian_eval':         _Prop.enumerated(csleqp.SLEQP_SETTINGS_ENUM_HESS_EVAL, HessianEval),
  'dual_estimation_type': _Prop.enumerated(csleqp.SLEQP_SETTINGS_ENUM_DUAL_ESTIMATION_TYPE, DualEstimationType),
  'initial_tr_choice':    _Prop.enumerated(csleqp.SLEQP_SETTINGS_ENUM_INITIAL_TR_CHOICE, InitialTRChoice),
  'bfgs_sizing':          _Prop.enumerated(csleqp.SLEQP_SETTINGS_ENUM_BFGS_SIZING, Sizing),
  'tr_solver':            _Prop.enumerated(csleqp.SLEQP_SETTINGS_ENUM_TR_SOLVER, TRSolver),
  'polishing_type':       _Prop.enumerated(csleqp.SLEQP_SETTINGS_ENUM_POLISHING_TYPE, PolishingType),
  'step_rule':            _Prop.enumerated(csleqp.SLEQP_SETTINGS_ENUM_STEP_RULE, StepRule),
  'linesearch':           _Prop.enumerated(csleqp.SLEQP_SETTINGS_ENUM_LINESEARCH, LineSearch),
  'parametric_cauchy':    _Prop.enumerated(csleqp.SLEQP_SETTINGS_ENUM_PARAMETRIC_CAUCHY, ParametricCauchy),
  'aug_jac_method':       _Prop.enumerated(csleqp.SLEQP_SETTINGS_ENUM_AUG_JAC_METHOD, AugJacMethod),

  'zero_eps':           _Prop.real(csleqp.SLEQP_SETTINGS_REAL_ZERO_EPS),
  'eps':                _Prop.real(csleqp.SLEQP_SETTINGS_REAL_EPS),
  'deriv_perturbation': _Prop.real(csleqp.SLEQP_SETTINGS_REAL_DERIV_PERTURBATION),
  'deriv_tol':          _Prop.real(csleqp.SLEQP_SETTINGS_REAL_DERIV_TOL),
  'cauchy_tau':         _Prop.real(csleqp.SLEQP_SETTINGS_REAL_CAUCHY_TAU),
  'cauchy_eta':         _Prop.real(csleqp.SLEQP_SETTINGS_REAL_CAUCHY_ETA),
  'linesearch_tau':     _Prop.real(csleqp.SLEQP_SETTINGS_REAL_LINESEARCH_TAU),
  'linesearch_eta':     _Prop.real(csleqp.SLEQP_SETTINGS_REAL_LINESEARCH_ETA),
  'linesearch_cutoff':  _Prop.real(csleqp.SLEQP_SETTINGS_REAL_LINESEARCH_CUTOFF),
  'feas_tol':           _Prop.real(csleqp.SLEQP_SETTINGS_REAL_FEAS_TOL),
  'slack_tol':          _Prop.real(csleqp.SLEQP_SETTINGS_REAL_SLACK_TOL),
  'stat_tol':           _Prop.real(csleqp.SLEQP_SETTINGS_REAL_STAT_TOL),
  'accepted_reduction': _Prop.real(csleqp.SLEQP_SETTINGS_REAL_ACCEPTED_REDUCTION),
  'deadpoint_bound':    _Prop.real(csleqp.SLEQP_SETTINGS_REAL_DEADPOINT_BOUND)
}

cdef class Settings:
  cdef csleqp.SleqpSettings* settings

  def __cinit__(self):
    csleqp_call(csleqp.sleqp_settings_create(&self.settings))

  def __init__(self, **values):
    for key, value in values.items():
      setattr(self, key, value)

  def props(self) -> dict:
    return {key: getattr(self, key) for key in opt_prop_map}

  def __str__(self) -> str:
    return 'Settings: {0}'.format(self.props())

  def __repr__(self) -> str:
    return 'Settings({0})'.format(repr(self.props()))

  cdef _prop(self, name):
    prop = opt_prop_map.get(name)
    if prop is None:
      raise AttributeError('Settings has no attribute \'{0}\''.format(name))
    return prop

  def __getattr__(self, name):
    prop = self._prop(name)
    prop_val = None

    if prop.prop_type == _PropType.Integer:
      prop_val = csleqp.sleqp_settings_int_value(self.settings,
                                                 prop.value)
    if prop.prop_type == _PropType.Real:
      prop_val = csleqp.sleqp_settings_real_value(self.settings,
                                                  prop.value)
    elif prop.prop_type == _PropType.Enumerated:
      prop_val = csleqp.sleqp_settings_enum_value(self.settings,
                                                  prop.value)
    else:
      assert prop.prop_type == _PropType.Bool

      prop_val = csleqp.sleqp_settings_bool_value(self.settings,
                                                  prop.value)

    return prop.to_py(prop_val)

  def __setattr__(self, name, value):
    prop = self._prop(name)
    prop_val = prop.from_py(value)

    if prop.prop_type == _PropType.Integer:
      csleqp_call(csleqp.sleqp_settings_set_int_value(self.settings,
                                                      prop.value,
                                                      prop_val))
    elif prop.prop_type == _PropType.Real:
      csleqp_call(csleqp.sleqp_settings_set_real_value(self.settings,
                                                       prop.value,
                                                       prop_val))
    elif prop.prop_type == _PropType.Enumerated:
      csleqp_call(csleqp.sleqp_settings_set_enum_value(self.settings,
                                                       prop.value,
                                                       prop_val))
    else:
      assert prop.prop_type == _PropType.Bool

      csleqp_call(csleqp.sleqp_settings_set_bool_value(self.settings,
                                                       prop.value,
                                                       prop_val))
