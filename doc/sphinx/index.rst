.. -*- mode: reST; coding: utf-8 -*-
.. MDAnalysis Tutorial documentation master file, created by
   sphinx-quickstart on Thu Oct 30 00:40:26 2014.

MDAnalysis Tutorial
===================

:MDAnalysis version: ≥ |MDAnalysis_version|
:Tutorial release: |release|
:Last updated: |today|


MDAnalysis_ is an open source Python library that helps you to quickly
write your own analysis algorithm for studying trajectories produced
by the most popular simulation packages [Michaud-Agrawal2011]_. 

This tutorial serves as an **introduction to MDAnalysis**. It starts
out with *installing the library* and then introduces the key components
of the library. It will show you 

* how to load a structure or a MD trajectory; 
* how to select parts of your system;
* how to work with atoms, residues and molecules through the object-oriented
  interface of MDAnalysis;
* how to analyze MD trajectories;
* how to write out modified trajectories;
* how to use algorithms in the  `MDAnalysis.analysis`_ module
  (intermediate level of difficulty)

The tutorial contains many links to the `online documentation`_ ,
which you can use to learn more about the functions, classes, an
methods that are discussed. The online help together with the
interactive Python documentation (``help(...)`` or ``...?`` in
:program:`ipython`) should help you while you are using the library.

If you have **questions or suggestions** please post them in the
`MDAnalysis User Discussion Group`_.

.. _MDAnalysis: http://www.mdanalysis.org
.. _online documentation: http://pythonhosted.org/MDAnalysis/
.. _MDAnalysis.analysis: https://pythonhosted.org/MDAnalysis/documentation_pages/analysis_modules.html
.. _MDAnalysis User Discussion Group: http://groups.google.com/group/mdnalysis-discussion


Contents
--------

.. toctree::
   :maxdepth: 1

   howto
   installation
   preparations
   basics
   atomgroups   
   trajectories
   writing
   analysismodule	


References
----------

.. [Michaud-Agrawal2011] N. Michaud-Agrawal, E. J. Denning,
   T. B. Woolf, and O. Beckstein.  MDAnalysis: A Toolkit for the
   Analysis of Molecular Dynamics Simulations.
   *J. Comput. Chem.* **32** (2011), 2319–2327, doi:`10.1002/jcc.21787`_
   PMCID:`PMC3144279`_

.. [Beckstein2009] O Beckstein. EJ Denning, JR Perilla, and TB
   Woolf. Zipping and Unzipping of Adenylate Kinase: Atomistic
   Insights into the Ensemble of Open/Closed Transitions. *J Mol Biol*
   **394** (2009), 160–176. doi:`10.1016/j.jmb.2009.09.009`_

.. [Seyler2014] L. Seyler and O. Beckstein. Sampling of large
   conformational transitions: Adenylate kinase as a testing
   ground. *Molec. Simul.* **40** (2014), 855–877. doi:
   `10.1080/08927022.2014.919497`_.

.. _10.1002/jcc.21787: http://doi.org/10.1002/jcc.21787
.. _PMC3144279: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3144279

.. _10.1016/j.jmb.2009.09.009: http://doi.org/10.1016/j.jmb.2009.09.009
.. _10.1080/08927022.2014.919497: http://doi.org/10.1080/08927022.2014.919497


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

