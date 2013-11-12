#!/usr/bin/env python

from distutils.core import setup

setup(
	name='PyCS',
	version='1.0dev',
	description='Python Curve Shifting ',
	long_description=open('README.txt').read(),
	author='Malte Tewes',
	license='GPLv3',
	author_email='malte.tewes[at]epfl.ch',
	url='http://www.cosmograil.org',
	packages=["pycs", "pycs.disp", "pycs.gen", "pycs.regdiff", "pycs.sim", "pycs.spl"],
	package_data={'pycs.gen': ['epfl.png']}
)

