"""
This subpackage contains general classes and functions to manipulate and display lightcurves and microlensing representations.
"gen" stands for general : all this is as independent as possible from the actual curve shifting techniques.

The ``lc`` module is the central module of PyCS, it defines the ``lightcurve`` class. 

"""

__all__ = ["lc", "util", "polyml", "spl", "splml", "sea", "stat", "mrg"]

import lc
import util
import polyml
import spl
import splml
import sea
import stat
import mrg
