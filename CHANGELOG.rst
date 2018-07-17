Change Log
==========

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.


v0.3.0 - 2018-07-17
-------------------

Changed
~~~~~~~

* Output data is now written in simple Nexus format, i.e. without a separate
  taxa block.
* Some refactoring in order to make predsim easier to use as a library (i.e.
  without its command-line interface).
  
`View commits <https://github.com/jmenglund/predsim/compare/v0.2.1...v0.3.0>`_


v0.2.1 - 2018-05-19
-------------------

Changed
~~~~~~~

* Minor updates to the documentation in ``README.rst``.
  
`View commits <https://github.com/jmenglund/predsim/compare/v0.2.0...v0.2.1>`_


v0.2.0 - 2018-05-19
-------------------

Added
~~~~~

* The command-line option ``--seeds-file`` for passing seed numbers 
  to Seq-Gen. This option allows the user to exactly repeat simulations.
* The command-line option ``--commands-file`` for outputting Seq-Gen commands 
  to a file (replaces the ``-c`` option).

Changed
~~~~~~~

* The `pandas <http://pandas.pydata.org>`_ library is no longer required.
* The principles of `Semantic Versioning <http://semver.org/>`_ will be 
  followed for new releases.

Removed
~~~~~~~

* The command-line ``-c`` option for writing Seq-Gen commands to a file 
  (replaced by the ``--commands-file`` option).


`View commits <https://github.com/jmenglund/predsim/compare/v0.1.1...v0.2.0>`_


v0.1.1 - 2016-08-11
-------------------

Changed
~~~~~~~

* Updates to the documentation in ``README.rst``.

`View commits <https://github.com/jmenglund/predsim/compare/v0.1.0...v0.1.1>`_


v0.1.0 - 2016-08-05
-------------------

Initial release.
