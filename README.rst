predsim
=======

|Build-Status| |PyPI-Status| |License| |DOI-URI|

``predsim`` is a command-line tool for simulating predictive
datasets from `MrBayes <http://mrbayes.sourceforge.net>`_ output files. 
Datasets can be simulated under the GTR+G+I substitution model or any of
its nested variants available in MrBayes (JC69, HKY85 etc.). The script 
uses `Seq-Gen <http://tree.bio.ed.ac.uk/software/seqgen/>`_ for 
simulating the DNA-sequences and builds on the third-party libraries 
`DendroPy <http://dendropy.org>`_ and `pandas <http://pandas.pydata.org>`_.

The code has been tested with Python 2.7, 3.3, 3.4 and 3.5.

Source repository: `<https://github.com/jmenglund/predsim>`_

--------------------------------

.. contents:: Table of contents
   :backlinks: top
   :local:


Requirements
------------

`Seq-Gen <http://tree.bio.ed.ac.uk/software/seqgen/>`_ must be installed on
your system.


Installation
------------

For most users, the easiest way is probably to install the latest version 
hosted on `PyPI <https://pypi.python.org/>`_:

.. code-block::

    $ pip install predsim

The project is hosted at https://github.com/jmenglund/predsim and 
can also be installed using git:

.. code-block::

    $ git clone https://github.com/jmenglund/predsim.git
    $ cd predsim
    $ python setup.py install


You may consider installing ``predsim`` and its required Python packages 
within a virtual environment in order to avoid cluttering your system's 
Python path. See for example the environment management system 
`conda <http://conda.pydata.org>`_ or the package 
`virtualenv <https://virtualenv.pypa.io/en/latest/>`_.


Usage
-----

.. code-block::
    
    $ predsim --help
    usage: predsim [-h] [-V] [-l INT] [-g INT] [-c FILE] [-s INT] [-p FILE]
                       pfile tfile [outfile]
    
    A command-line utility that reads posterior output of MrBayes and simulates
    predictive datasets with Seq-Gen.
    
    positional arguments:
      pfile                 path to a MrBayes p-file
      tfile                 path to a MrBayes t-file
      outfile               path to output file (default: <stdout>)
    
    optional arguments:
      -h, --help            show this help message and exit
      -V, --version         show program's version number and exit
      -l INT, --length INT  sequence lenght (default: 1000)
      -g INT, --gamma-cats INT
                            number of gamma rate categories (default: continuous)
      -c FILE, --commands-file FILE
                            path to output file with used Seq-Gen commands
      -s INT, --skip INT    number of records (trees) to skip at the beginning of
                            the sample (default: 0)
      -p FILE, --seqgen-path FILE
                            path to a Seq-Gen executable (default: "seq-gen")

It is strongly recommended that you use the ``-c FILE`` option to check the 
commands run by Seq-Gen.

Depending on your Python version, you might need to specify the full path to 
your Seq-Gen executable with the ``-p FILE`` option.


Running tests
-------------

Testing is carried out with `pytest <http://pytest.org>`_. Here is an 
example on how to run the test suite and generating a coverage report:

.. code-block::

    $ cd predsim
    $ pip install pytest pytest-cov pytest-pep8
    $ py.test -v --cov-report term-missing --cov predsim.py --pep8


License
-------

``predsim`` is distributed under the 
`MIT license <https://opensource.org/licenses/MIT>`_.


Citing
------

If you use results produced with this package in a scientific 
publication, please just mention the package name in the text and 
cite the Zenodo DOI of this project:

|DOI-URI|

You can select a citation style from the dropdown menu in the 
*"Cite as"* section on the Zenodo page.

``predsim`` relies on other software that also should be cited. Below are 
suggested citations for Seq-Gen, DendroPy and pandas, 
respectively:

* Rambaut A., Grassly NC 1997. Seq-Gen: an application for the Monte 
  Carlo simulation of DNA sequence evolution along phylogenetic trees. 
  Comput. Appl. Biosci. 13:235–238.

* Sukumaran J, Holder MT. 2010. DendroPy: a Python library for 
  phylogenetic computing. Bioinformatics 26:1569–1571.

* McKinney W. 2010. Data structures for statistical computing in python.
  *In* Proceedings of the 9th Python in Science Conference 
  (van der Walt S, Millman J, editors), pages 51–56.


Author
------

Markus Englund, `orcid.org/0000-0003-1688-7112 <http://orcid.org/0000-0003-1688-7112>`_

.. |Build-Status| image:: https://travis-ci.org/jmenglund/predsim.svg?branch=master
   :target: https://travis-ci.org/jmenglund/predsim
.. |PyPI-Status| image:: https://img.shields.io/pypi/v/predsim.svg
   :target: https://pypi.python.org/pypi/predsim
.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://raw.githubusercontent.com/jmenglund/predsim/master/LICENSE.txt
.. |DOI-URI| image:: https://zenodo.org/badge/23107/jmenglund/predsim.svg
      :target: https://zenodo.org/badge/latestdoi/23107/jmenglund/predsim
