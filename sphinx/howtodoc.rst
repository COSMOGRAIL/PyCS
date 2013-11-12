Building this documentation
===========================

We have now switched to Sphinx : `<http://sphinx.pocoo.org/>`_.

To build the documentation by yourself (assuming you have the latest version of sphinx) :

::

	cd sphinx
	make clean
	make apidoc
	make html 

	
You can then open ``index.html`` that you will find in ``_build/html/`` with your browser...