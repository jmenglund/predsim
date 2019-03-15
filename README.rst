predsim
=======

|Build-Status| |Coverage-Status| |PyPI-Status| |License| |DOI-URI|

**predsim** is a simple command-line tool for simulating predictive
datasets from `MrBayes <http://mrbayes.sourceforge.net>`_ output files. 
Datasets can be simulated under the GTR+G+I substitution model or any nested 
variant available in MrBayes (JC69, HKY85 etc.). The code is contained
within a single module that can be imported using Python's import mechanism.
It uses `Seq-Gen <http://tree.bio.ed.ac.uk/software/seqgen/>`_ for simulating
the DNA-sequences and builds on the third-party library
`DendroPy <http://dendropy.org>`_.

The code has been tested with Python 3.3 and 3.6.

Source repository: `<https://github.com/jmenglund/predsim>`_

--------------------------------

.. contents:: Table of contents
   :local:
   :backlinks: none


Prerequisites
-------------

* Python 3.3+
* The Python library `DendroPy <http://dendropy.org>`_ (version 4.0 or higher)
* The command-line tool `Seq-Gen <http://tree.bio.ed.ac.uk/software/seqgen/>`_

An easy way to get Python working on your computer is to install the free
`Anaconda distribution <http://anaconda.com/download)>`_.


Installation
------------

For most users, the easiest way is probably to install the latest version 
hosted on `PyPI <https://pypi.org/>`_:

.. code-block::

    $ pip install predsim

The project is hosted at `<https://github.com/jmenglund/predsim>` and 
can also be installed using git:

.. code-block::

    $ git clone https://github.com/jmenglund/predsim.git
    $ cd predsim
    $ python setup.py install


You may consider installing predsim and its required Python packages 
within a virtual environment in order to avoid cluttering your system's 
Python path.


Usage
-----

.. code-block::
    
    predsim --help
    usage: predsim [-h] [-V] [-l N] [-g N] [-s N] [-n N] [-o {nexus,phylip}]
                   [-p FILE] [--seeds-file FILE] [--commands-file FILE]
                   [--trees-file FILE]
                   pfile tfile

    A command-line utility that reads posterior output of MrBayes and simulates
    predictive datasets with Seq-Gen.

    positional arguments:
      pfile                 path to a MrBayes p-file
      tfile                 path to a MrBayes t-file

    optional arguments:
      -h, --help            show this help message and exit
      -V, --version         show program's version number and exit
      -l N, --length N      sequence lenght (default: 1000)
      -g N, --gamma-cats N  number of gamma rate categories (default: continuous)
      -s N, --skip N        number of records (trees) to skip at the beginning of
                            the sample (default: 0)
      -n N, --num-records N
                            number of records (trees) to use in the simulation
      -o {nexus,phylip}, --output-format {nexus,phylip}
                            output format (default: "nexus")
      -p FILE, --seqgen-path FILE
                            path to a Seq-Gen executable (default: "seq-gen")
      --seeds-file FILE     path to file with seed numbers to pass to Seq-Gen
      --commands-file FILE  path to output file with commands used by Seq-Gen
      --trees-file FILE     path to output file with trees used by Seq-Gen


* It is strongly recommended that you use the ``--commands-file`` and the 
  ``--trees-file`` options to check the input given to Seq-Gen.


Running the tests
-----------------

Testing is carried out with `pytest <https://docs.pytest.org/>`_:

.. code-block::

    $ pytest test_predsim.py

Test coverage can be calculated with `Coverage.py
<https://coverage.readthedocs.io/>`_ using the following commands:

.. code-block::

    $ coverage run -m pytest test_predsim.py
    $ coverage report -m predsim.py

The code follow style conventions in `PEP 8
<https://www.python.org/dev/peps/pep-0008/>`_, which can be checked
with `pycodestyle <http://pycodestyle.pycqa.org>`_:

.. code-block::

    $ pycodestyle predsim.py test_predsim.py setup.py


License
-------

predsim is distributed under the 
`MIT license <https://opensource.org/licenses/MIT>`_.


Citing
------

If you use results produced with this package in a scientific 
publication, please just mention the package name in the text and 
cite the Zenodo DOI of this project:

|DOI-URI|

You can select your preferred citation style in the "Cite as" section
on the Zenodo page.

predsim relies on other software that also should be cited. Below are 
suggested citations for Seq-Gen and DendroPy:

* Rambaut A, Grassly NC. 1997. Seq-Gen: an application for the Monte 
  Carlo simulation of DNA sequence evolution along phylogenetic trees. 
  Comput. Appl. Biosci. 13:235–238. DOI:
  `10.1093/bioinformatics/13.3.235 <https://doi.org/10.1093/bioinformatics/13.3.235>`_

* Sukumaran J, Holder MT. 2010. DendroPy: a Python library for 
  phylogenetic computing. Bioinformatics 26:1569–1571. DOI:
  `10.1093/bioinformatics/btq228 <https://doi.org/10.1093/bioinformatics/btq228>`_


Author
------

Markus Englund, `orcid.org/0000-0003-1688-7112 <http://orcid.org/0000-0003-1688-7112>`_

.. |Build-Status| image:: https://travis-ci.org/jmenglund/predsim.svg?branch=master
   :target: https://travis-ci.org/jmenglund/predsim
.. |Coverage-Status| image:: https://codecov.io/gh/jmenglund/predsim/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/jmenglund/predsim
.. |PyPI-Status| image:: https://img.shields.io/pypi/v/predsim.svg
   :target: https://pypi.python.org/pypi/predsim
.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://raw.githubusercontent.com/jmenglund/predsim/master/LICENSE.txt
.. |DOI-URI| image:: https://zenodo.org/badge/23107/jmenglund/predsim.svg
   :target: https://zenodo.org/badge/latestdoi/23107/jmenglund/predsim
