Download & Installation
=======================


Dependencies
------------

PyCS is currently developed using python 2.7.3 (but might be fine with older versions).

It uses numpy, scipy and matplotlib. If you're new to python, you could get all this from `here <http://www.enthought.com/products/epd_free.php>`_.

If and only if you want to use the "regression difference" technique, you'll also have to install `pymc <http://pypi.python.org/pypi/pymc/>`_.


Download
--------

Get PyCS via anonymous SVN checkout::
	
	svn checkout https://svn.epfl.ch/svn/pycs/trunk/PyCS ./PyCS

Without surprise, ``pycs`` is work in progress, and we are actively committing to this repository. If you really want the exact version of PyCS that was used in the original paper, you can::

	svn export https://svn.epfl.ch/svn/pycs/tags/PyCS-1.0 ./PyCS-1.0
	

For any questions, and help/advice in case of problems, I encourage you to drop me (``mtewes at astro.uni-bonn dot de``) a line. I would be happy to interact with you about curve shifting !



.. You *could* download the version of the code that we used for our `paper <http://arxiv.org/abs/1208.5598>`_ from here :  `pycs_1.0.tgz <http://www.cosmograil.org/>`_.
.. But clearly, ``pycs`` is work in progress... To get the latest version and easier help/advice in case of problems, I encourage you to drop me (``malte.tewes at epfl dot ch``) a line, and I would be glad to give you access to our `SVN repository <https://svn.epfl.ch/svn/pycs/>`_::
.. I'm busy with my thesis until end of September, but afterwards I'm happy to help !

Installation
------------

The easy answer, if you do not plan to tweak the code : cd into the PyCS directory, and::

	python setup.py install

or maybe::

	python setup.py install --user

... if you don't have write access to the global site-packages directory of your machine. More info about this type of installation : `<http://docs.python.org/install/>`_. So these commands will copy ``pycs`` somewhere in your python path. Once done, you *could* in principle delete the PyCS directory that you downloaded, as it is not required anymore. By the way, good to know, to locate any package from within python, just::

	import pycs
	print pycs.__file__	

This tells where the ``pycs`` files that you import are located. Remove them to uninstall...

Now, the longer answer: PyCS is not really "production level code", and you will probably want to tweak/read the code and do frequent commits or updates. As PyCS is pure python, the "installation" is in fact trivial : you just have to tell python were the ``pycs`` directory is.
Hence an alternative solution, especially if you work with SVN, is to add the directory containing ``pycs`` to your ``PYTHONPATH``. That way you can put PyCS in any convenient place on your system, and modify it on the fly, without having to ``python setup.py install`` after any modification.

For tcsh, add for instance something like this to your .tcshrc ::


	#setenv PYTHONPATH /path/to/PyCS
	# or, to add to the existing stuff :
	
	setenv PYTHONPATH ${PYTHONPATH}:/path/to/PyCS
	
or, if you use bash, to your .bash_profile or equivalent::

	export PYTHONPATH=${PYTHONPATH}:/path/to/PyCS/

... and then from any python you can simply ``import pycs``.


Yet another very simple solution is to explicitly start your scripts with something like :

.. code-block:: python
	
	import sys, os
	sys.path.append("/path/to/PyCS")
	
This might be an option if you intend to work with several versions.

