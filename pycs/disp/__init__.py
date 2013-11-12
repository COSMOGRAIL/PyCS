"""
This subpackage defines the "dispersion-inspired" family of curve shifting methods.
It contains functions that optimize lightcurves (shifts, ML) so to minimize a given dispersion between them.
Dispersion functions are defined in the module :py:mod:`pycs.disp.disps`.

* multiopt defines the building blocks of these optimizers.
* topopt assembles thoes blocks into wrapper functions.

"""

__all__ = ["disps", "multiopt", "topopt"]
import disps
import multiopt
import topopt
