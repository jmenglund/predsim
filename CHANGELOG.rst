Change Log
==========

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.


v0.7.0 - 2019-07-11
-------------------

* Added flags ``-f`` and ``--freqs`` to manually set base frequences.
  As a consequence, the flags ``-o`` and ``--out-format`` have to be
  used when specifying the output format.

`View commits <https://github.com/jmenglund/predsim/compare/v0.6.0...v0.7.0>`_


v0.6.0 - 2019-05-05
-------------------

Added
~~~~~

* Ability to write used trees to a file with the ``--trees-file`` flag.
* Simple test data for various substitution models.


Changed
~~~~~~~

* Command-line flags ``-o`` and ``--output-format`` replaced with ``-f`` and
  ``--format``, respectively.
* The file with seed numbers may now contain more numbers than being used.
* Content is written to output after completing each simulation (instead
  of writing everything to output at the end).
* NEXUS output is written with separate taxa blocks.
* Missing state frequences result in equal frequences instead of raising an error.


Removed
~~~~~~~

* Support for Python versions prior to 3.3.
* Two test files.

`View commits <https://github.com/jmenglund/predsim/compare/v0.5.0...v0.6.0>`_


v0.5.0 - 2019-03-10
-------------------

Added
~~~~~

* Several tests in ``test_predsim.py``.
* A directory, ``test_files/``, for files used by the test suite.


Changed
~~~~~~~

* By default, simulated data is written to standard output (stdout)
  and it is no longer possible to name an output file (but output
  can still be redirected to a file!).
* Records are processed iteratively, and output written after the
  processing of an individual record has been finished.
* Opening and closing input files when using the command-line interface.
* Documentation in ``README.rst``.

`View commits <https://github.com/jmenglund/predsim/compare/v0.4.0...v0.5.0>`_


v0.4.0 - 2019-03-09
-------------------

Added
~~~~~

* The command-line flags ``-o`` and ``--output-format`` now allow 
  writing output to either Phylip or Nexus format. 


Changed
~~~~~~~

* The path to the t-file is now correct when using the command-line interface.
* If using number of gamma rate categories, gamma shape is now also required.

`View commits <https://github.com/jmenglund/predsim/compare/v0.3.1...v0.4.0>`_


v0.3.1 - 2019-03-09
-------------------

Changed
~~~~~~~

* Small refactoring of the function for reading the p-file.
* Fixed codecov installation in Travis-CI configuration file.
  
`View commits <https://github.com/jmenglund/predsim/compare/v0.3.0...v0.3.1>`_


v0.3.0 - 2018-07-18
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
