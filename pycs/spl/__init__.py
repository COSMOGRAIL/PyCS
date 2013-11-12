"""
Simultaneous spline fit technique. This subpackage contains functions to optimize lightcurves (their shifts and ML) so that they fit to a common spline.

* ``multiopt`` contains "builing blocks" for such optimizers, and
* ``topopt`` assembles these building blocks into big wrappers.

"""

__all__ = ["multiopt", "topopt"]
import multiopt
import topopt
