.. -*- coding: utf-8 -*-

.. _MDAnalysis: http://www.mdanalysis.org
.. _MDAnalysisTests: http://wiki.mdanalysis.org/UnitTests
.. _conda: http://conda.pydata.org/
.. _pip: https://pip.pypa.io/en/stable/
.. _macports: https://www.macports.org/

.. _chapter-installing-mdanalysis:

======================================
 MDAnalysis installation instructions
======================================

You will need a working installation of MDAnalysis_ on your **Linux** or **Mac
OS X** laptop. Windows is not supported at this time but Windows users can try
to use a `virtual machine with MDAnalysis pre-installed`_.

For this tutorial, you will need at least version |MDAnalysis_version| of
MDAnalysis (which should be automatically installed when you follow the
instructions below).

The installation recipes will install *MDAnalysis*, *MDAnalysisTests* (test
cases) [#otherdownloads]_, and all required libraries and Python packages.

.. note:: If you have problems at any stage, ask for :ref:`help`.



Installation methods
====================

We tested two approaches to install MDAnalysis and all its dependencies. The
:ref:`first approach<conda-installation>` is based on the `conda`_ package
management system and is recommended for anyone new to scientific computing
with Python (or if the :ref:`second approach<distribution-pip-installation>`
fails).


.. _conda-installation:

Conda-based
-----------

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


.. _distribution-pip-installation: 

Pip and distribution package manager
------------------------------------

We use the distribution's package manager for most of the pre-requisites and
install any remaining packages (and MDAnalysis) with pip_.

Linux
~~~~~

.. Note:: 

   All installations of pre-requisites are performed in system directories and
   require root access (via the ``sudo`` command). MDAnalysis can also be
   installed in a user directory by providing the ``--user`` option to pip_.

Ubuntu 14.04
````````````

Recent Ubuntu distributions:

.. code-block:: bash

   sudo apt-get update
   sudo apt-get install -y build-essential python-dev python-setuptools python-pip
   sudo apt-get install -y python-numpy python-scipy python-matplotlib python-biopython python-networkx ipython
   sudo apt-get install -y libhdf5-serial-dev libnetcdf-dev
   
   sudo pip install netCDF4
   sudo pip install MDAnalysis MDAnalysisTests


Debian 7.6, 7.8 Wheezy
``````````````````````

Most recent Debian distributions should all work with the following:

.. code-block:: bash

   sudo apt-get update
   sudo apt-get install -y build-essential python-dev python-setuptools python-pip
   sudo apt-get install -y python-numpy python-scipy python-matplotlib python-biopython python-networkx ipython
   sudo apt-get install -y libhdf5-serial-dev libnetcdf-dev
   
   sudo pip install netCDF4
   sudo pip install MDAnalysis MDAnalysisTests



Mac OS X (≥ 10.6.8)
~~~~~~~~~~~~~~~~~~~

Macports
````````

Using macports_:

.. code-block:: bash

   sudo port install py27-numpy  py27-cython
   sudo port install py27-scipy  py27-matplotlib py27-biopython py27-ipython+notebook
   sudo port install hdf5 netcdf+dap+netcdf4
   
   sudo pip install netCDF4
   sudo pip install MDAnalysis MDAnalysisTests


.. _help:

Help!
=====

If there are problems then please have a closer look at the `installation
notes`_ and the `installation recipes`_; in particular, `installing the netcdf
library`_ can become more involved.

If you need help with installation issues, please do not hesitate to ask on the
`user discussion group`_ or via the `issue tracker`_.

If nothing else works you can also use a complete installation inside a
`virtual machine with MDAnalysis pre-installed`_.

.. _installation notes: http://wiki.mdanalysis.org/Install
.. _installation recipes: http://wiki.mdanalysis.org/InstallRecipes
.. _installing the netcdf library: http://wiki.mdanalysis.org/netcdf
.. _user discussion group: http://groups.google.com/group/mdnalysis-discussion
.. _tutorial git repository: https://github.com/MDAnalysis/MDAnalysisTutorial
.. _`vm/README.rst`: https://github.com/MDAnalysis/MDAnalysisTutorial/tree/master/vm
.. _`virtual machine with MDAnalysis pre-installed`: 
   http://www.mdanalysis.org/MDAnalysisTutorial/installation.html#virtual-machine
.. _issue tracker: https://github.com/MDAnalysis/mdanalysis/issues


.. rubric:: Footnotes

.. [#otherdownloads] For the workshop you will also need to additionally :ref:`download
                     example trajectories<chapter-datadownload>`.
