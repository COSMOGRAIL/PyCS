Download & Installation
=======================


Dependencies
------------

PyCS is developed using python 2.7 and might be fine with older versions, but not with python 3.

It requires ``numpy``, ``scipy`` and ``matplotlib``. Those are all you need to run the free knot splines.

Some other estimators (e.g., the "regression difference" technique) might have furter dependencies:

* `pymc <https://github.com/pymc-devs/pymc>`_
* `scikit-learn <http://scikit-learn.org>`_


Download
--------

Get the latest PyCS by cloning it from `GitHub <https://github.com/COSMOGRAIL/PyCS>`_::

	git clone https://github.com/COSMOGRAIL/PyCS.git


If for some reason you want the exact version of PyCS that was used in the original paper, you can still find it here::

	svn export https://svn.epfl.ch/svn/pycs/tags/PyCS-1.0 ./PyCS-1.0
	


Installation
------------

PyCS is not exactly "production level code", and you might want to update or tweak the sources.
So to keep this simple, we suggest to just add your cloned repository to your ``PYTHONPATH``.

For **tcsh**, add for instance something like this to your .tcshrc ::

	#setenv PYTHONPATH /path/to/PyCS
	# or, to add to the existing stuff :
	
	setenv PYTHONPATH ${PYTHONPATH}:/path/to/PyCS
	
or, if you use **bash**, to your .bash_profile or equivalent::

	export PYTHONPATH=${PYTHONPATH}:/path/to/PyCS/

... and then from any python you can simply ``import pycs``.




If you don't plan to tweak the code, you can also simply

::

	python setup.py install

or maybe

::

	python setup.py install --user

... if you don't have write access to the global site-packages directory of your machine.
