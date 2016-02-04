#!/usr/bin/env python

from distutils.core import setup

setup(
	name='PyCS',
	version='2.0dev',
	description='Python Curve Shifting ',
	long_description=open('README.txt').read(),
	license='GPLv3',
	url='http://www.cosmograil.org',
	packages=["pycs", "pycs.disp", "pycs.gen", "pycs.regdiff", "pycs.sim", "pycs.spl"],
	package_data={'pycs.gen': ['epfl.png']}
)

