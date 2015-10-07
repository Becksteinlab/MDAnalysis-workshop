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

.. toctree::
   :maxdepth: 1

   install_conda
   install_native_pip


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
