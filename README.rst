predsim
=======

|License|

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

`SeqGen <http://tree.bio.ed.ac.uk/software/seqgen/>`_ must be installed on
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
    usage: predsim [-h] [-V] [-l INT] [-g INT] [-c PATH] [-s INT] [-p PATH]
                       pfile tfile [outfile]
    
    A command-line utility that reads posterior output of MrBayes and simulates
    predictive datasets with SeqGen.
    
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
      -c PATH, --commands-file PATH
                            path to output file with used SeqGen commands
      -s INT, --skip INT    number of records (trees) to skip at the beginning of
                            the sample (default: 0)
      -p PATH, --seqgen-path PATH
                            path to a SeqGen executable (default: "seq-gen")


License
-------

``predsim`` is distributed under the 
`MIT license <https://opensource.org/licenses/MIT>`_.


Author
------

Markus Englund, `orcid.org/0000-0003-1688-7112 <http://orcid.org/0000-0003-1688-7112>`_


.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://raw.githubusercontent.com/jmenglund/predsim/master/LICENSE.txt
