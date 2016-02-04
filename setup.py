#!/usr/bin/env python

from distutils.core import setup

setup(
	name='PyCS',
	version='2.0dev',
	description='Python Curve Shifting ',
	long_description=open('README.md').read(),
	license='GPLv3',
	url='http://www.cosmograil.org',
	packages=["pycs", "pycs.gen", "pycs.sim", "pycs.disp", "pycs.spl", "pycs.regdiff", "pycs.regdiff2", "pycs.spldiff", "pycs.tdc"],
	package_data={'pycs.gen': ['epfl.png']}
)

