.. -*- coding: utf-8 -*-

.. _chapter-installing-mdanalysis:

=======================
 Installing MDAnalysis
=======================

Before you can start the tutorial, you need a working installation of
MDAnalysis_ (with its test suite, MDAnalysisTests_) on your Linux or Mac
OS X machine. The following contains hints and links to further
documentation to facilitate installation of the package. In an ideal
setting (development tools (gcc tool chain), python development
headers, full scientific Python stack with numpy, scipy, matplotlib
already installed, working netcdf/HDF5 library), installation will be
as simple as ::

   pip install --user MDAnalysis MDAnalysisTests

In a less ideal setting, one typically has to install additional
packages through the distribution's package manager as described in
the links under :ref:`local-installation`. Alternatively, for this
tutorial one can also set up a Linux :ref:`virtual-machine` that
installs itself with everything needed and the most recent release of
MDAnalysis.

.. Note:: For this tutorial, you will need at least version
          |MDAnalysis_version| of MDAnalysis.

.. _MDAnalysis: http://www.mdanalysis.org
.. _MDAnalysisTests: http://wiki.mdanalysis.org/UnitTests


Installation methods
====================

.. _local-installation:

Local installation
------------------

The latest release of MDAnalysis_ (and the full test suite) can be
installed from the python package index with pip_ ::

  pip install --user MDAnalysis MDAnalysisTests

(Installation of the test suite is required for this tutorial because
we will use some of the data files that are part of the tests.)

If there are problems then please have a closer look at the
`installation notes`_ and the `installation recipes`_; in particular,
`installing the netcdf library`_ can become more involved.

If you need help with installation issues, please do not hesitate to
ask on the `user discussion group`_.

For this tutorial you can also alternatively use a complete
installation inside a :ref:`virtual-machine`.

.. _pip: http://www.pip-installer.org/en/latest/index.html
.. _installation notes: http://wiki.mdanalysis.org/Install
.. _installation recipes: http://wiki.mdanalysis.org/InstallRecipes
.. _installing the netcdf library: http://wiki.mdanalysis.org/netcdf
.. _user discussion group: http://groups.google.com/group/mdnalysis-discussion
.. _tutorial git repository: https://github.com/MDAnalysis/MDAnalysisTutorial
.. _`vm/README.rst`: https://github.com/MDAnalysis/MDAnalysisTutorial/tree/master/vm

.. _virtual-machine:

Virtual machine
---------------

You will first need to clone the `tutorial git repository`_ with
:program:`git`::

  git clone https://github.com/MDAnalysis/MDAnalysisTutorial.git 
  cd MDAnalysisTutorial

The directory ``vm`` contains configuration files for `vagrant`_
virtual machines (VM) (using `VirtualBox`_). The file `vm/README.rst`_
describes setup in more detail but provided that `vagrant`_ and
`VirtualBox`_ are installed, the following should provide you with a
working VM::

  cd vm/Ubuntu/14.04
  vagrant up
  vagrant ssh

When you are done, just exit the ssh session (``exit`` or Control-D)
and halt the VM::

  vagrant halt

You can access your real home directory (:envvar:`$HOME`) in the virtual
machine at the path ``/myhome``. Anything that you do in this
directory will be reflected in your real home directory, including
deletion of files!

.. _Vagrant: https://www.vagrantup.com/
.. _VirtualBox: https://www.virtualbox.org/



Testing the installation
========================

.. _test cases: http://wiki.mdanalysis.org/UnitTests

MDAnalysis comes with over 2500 `test cases`_ that check its
functionality. These test cases can be run with the command ::

  python -c 'from MDAnalysis.tests import test; test(label="full", verbose=3, extra_argv=["--exe"])'

This can take a few minutes. Ideally, you should only get passing
tests ("ok" or just a single dot "." when using :code:`verbose=1`) or
"KnownFailures".

.. Note:: 
   The test suite consumes a considerable amount of memory (> 4GB) and
   thus it might fail or become very slow on machines with
   insufficient memory such as a :ref:`virtual-machine`. (This is a known
   problem with the test suite and will be addressed in the future.)


