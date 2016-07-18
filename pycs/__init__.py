"""
This is PyCS, a python toolbox to estimate time delays between COSMOGRAIL light curves.

:License: GPLv3 (see README.md)
:Website: www.cosmograil.org

"""
# Keeping this docstring short, to not be redundant with the website / documentation.

__version__ = "2.0dev"

__all__ = ["gen", "disp", "spl", "sim", "tdc", "spldiff"]


import gen
import disp
import spl
import sim
import tdc
import spldiff
#import regdiff # would require pymc, but we don't want to add weird dependencies.
#import regdiff2
