.. -*- coding: utf-8 -*-

.. _pip: https://pip.pypa.io/en/stable/
.. _macports: https://www.macports.org/


.. _distribution-pip-installation: 

============================================================================
 ``pip`` and native distribution package manager installation of MDAnalysis
============================================================================

We use the native distribution's package manager for most of the
pre-requisites and install any remaining packages (and MDAnalysis)
with pip_.

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


