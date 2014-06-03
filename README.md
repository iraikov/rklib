rklib
=====

Runge-Kutta methods in Standard ML.

This library is a port of the Haskell Runge-Kutta library by Uwe
Hollerbach (http://hackage.haskell.org/package/rungekutta).  It
contains a collection of explicit Runge-Kutta methods of various
orders, some with fixed-size steps (no error estimate) and some
that provid error estimate intended for adaptive stepsize.

The library provides routines that can use arbitrary coefficients
(Butcher tableau), then instantiates them with different tableaux
based on standard Runge-Kutta methods. Adaptive step-size methods add
a row of coefficients that are used to compute the error.


* non-adaptive solvers:

  rkfe, rk3, rk4a, rk4b

* adaptive solvers:

  rkhe, rkbs, rkf45, rkck, rkdp, rkf78, rkv65

* adaptive solvers with interpolation (CERK):

  cerkdp

* auxiliary non-adaptive solvers (error estimators from the adaptive ones):

  rkhe_aux, rkbs_aux, rkf45_aux, rkck_aux, rkdp_aux, rkf78_aux, rkv65_aux

