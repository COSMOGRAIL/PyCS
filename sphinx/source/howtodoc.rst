About this documentation
========================

The documentation is build using Sphinx : `<http://sphinx.pocoo.org/>`_.

Some of it (e.g., tutorials) is written in form of plain .rst files, and other parts are created by ``sphinx-apidoc``


Try to stick to their recommendations about documenting code in your docstrings. And make use of the described way to add links to the documentation



To build the documentation by yourself (assuming you have the latest version of sphinx) :

::

	cd sphinx
	make clean
	make apidoc
	make html 

	
You can then open ``index.html`` that you will find in ``_build/html/`` with your browser...