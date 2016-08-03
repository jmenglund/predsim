predsim
=======

|Build-Status| |License|

``predsim`` is a command-line tool for simulating predictive
datasets from `MrBayes <http://mrbayes.sourceforge.net>`_ output files. 
Datasets can be simulated under the GTR+G+I substitution model or any of
its nested variants available in MrBayes (JC69, HKY85 etc.).

The script uses `Seq-Gen <http://tree.bio.ed.ac.uk/software/seqgen/>`_ for 
simulating the DNA-sequences and builds on the third-party libraries 
`DendroPy <http://dendropy.org>`_ and `pandas <http://pandas.pydata.org>`_.

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

.. code-block:: bash

    $ pip install predsim

The project is hosted at https://github.com/jmenglund/predsim and 
can also be installed using git:

.. code-block:: bash

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

.. code-block:: console
    
    $ predsim -h
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


License
-------

``predsim`` is distributed under the 
`MIT license <https://opensource.org/licenses/MIT>`_.


Citing
------

If you use results produced with this package in a scientific 
publication, please just mention the package name in the text and 
cite the Zenodo DOI of this project:

[A Zenodo DOI will be inserted here]

You can select a citation style from the dropdown menu in the 
*"Cite as"* section on the Zenodo page.


Author
------

Markus Englund, `orcid.org/0000-0003-1688-7112 <http://orcid.org/0000-0003-1688-7112>`_

.. |Build-Status| image:: https://travis-ci.org/jmenglund/predsim.svg?branch=master
   :target: https://travis-ci.org/jmenglund/predsim
.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://raw.githubusercontent.com/jmenglund/predsim/master/LICENSE.txt
