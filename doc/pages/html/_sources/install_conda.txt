.. -*- coding: utf-8 -*-

.. _conda: http://conda.pydata.org/
.. _pip: https://pip.pypa.io/en/stable/


.. _conda-installation:

============================================
 ``conda``-based installation of MDAnalysis
============================================

Get the appropriate `miniconda installer
<http://conda.pydata.org/miniconda.html>`_ for Python 2.7; in the example we
are using the Linux 64 bit one. Common choices:

- Linux x86_64 (64 bit): https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
- Mac OS X (64 bit): https://repo.continuum.io/miniconda/Miniconda-latest-MacOSX-x86_64.sh

Run the installer to install the :program:`conda` package manager and the
necessary packages:

.. code-block:: bash

   # example for Linux x86_64 
   wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
 
   # for Mac OS X, uncomment the following line and comment the preceding one
   # wget https://repo.continuum.io/miniconda/Miniconda-latest-MacOSX-x86_64.sh

   # from here on use the same commands for Linux/Mac OS X
   chmod +x miniconda.sh
   ./miniconda.sh -b
   export PATH=${HOME}/miniconda/bin:$PATH
   conda update --yes conda
   conda create --yes -n mdaenv python=2.7 numpy=1.9.2 scipy=0.16 nose=1.3.7 ipython
   source activate mdaenv
   conda install --yes python=2.7 cython biopython matplotlib networkx netcdf4

   # install the latest release of MDAnalysis (≥ 0.11.0)
   pip install --upgrade MDAnalysis 
   pip install --no-cache-dir --upgrade MDAnalysisTests

- The installation is performed in the `virtual environment
  <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_ named **mdaenv**,
  which must be activated for use in the each shell session:

  .. code-block:: bash

     source  activate mdaenv 

  (For more technical details see `virtualenv
  <https://virtualenv.pypa.io/en/latest/>`_.)
- Add the line:

  .. code-block:: bash

     export PATH=${HOME}/miniconda/bin:$PATH

  to your shell start-up file (e.g. ``~/.bashrc``) so that the ``activate``
  script and other commands are found.
- The ``--no-cache-dir`` option for :program:`pip` may be necessary to avoid a
  ``MemoryError`` in low-memory environments such as our virtual machines; with
  lots of memory you may omit it (see `pip issue #2984
  <https://github.com/pypa/pip/issues/2984>`_).
- :program:`conda` will also install the HDF5 and netcdf libraries for you so
  you will have a *full feature* installation of MDAnalysis
