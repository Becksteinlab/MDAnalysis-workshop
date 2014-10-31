.. -*- encoding: utf-8 -*-

=========================
 Working with AtomGroups
=========================

A :class:`~MDAnalysis.core.AtomGroup.AtomGroup` has a large number of
methods attributes defined that provide information about the atoms
such as names, indices, or the coordinates in the
:attr:`~MDAnalysis.core.AtomGroup.AtomGroup.positions` attribute::

  >>> CA = u.selectAtoms("protein and name CA")
  >>> r = CA.positions
  >>> r.shape
  (214, 3)

The resulting output is a :class:`numpy.ndarray`. The main purpose of
MDAnalysis is to get trajectory data into numpy arrays!


Important methods and attributes of AtomGroup
=============================================

The coordinates :attr:`~MDAnalysis.core.AtomGroup.AtomGroup.positions`
attribute is probably the most important information that you can get
from an :class:`~MDAnalysis.core.AtomGroup.AtomGroup`.

Other quantities that can be easily calculated for a
:class:`~MDAnalysis.core.AtomGroup.AtomGroup` are

* the center of mass
  :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.centerOfMass` and the
  center of geoemtry (or centroid)
  :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.centerOfGeometry`
  (equivalent to
  :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.centroid`);

* the total mass
  :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.totalMass`;

* the total charge
  :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.totalCharge` (if partial
  charges are defined in the topology);

* the radius of gyration

  .. math:: 

        R_\mathrm{gyr} = \sqrt{\frac{1}{M}\sum_{i=1}^{N} m_i(\mathbf{r}_i - \mathbf{R})^2}

  with :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.radiusOfGyration`;

* the principal axes :math:`\mathbf{p}_1, \mathbf{p}_2, \mathbf{p}_1`
  from :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.principalAxes` via
  a diagonalization of the tensor of inertia
  :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.momentOfInertia`,

  .. math::

       \Lambda = U^T I U, \quad \text{with}\  U=(\mathbf{p}_1,
       \mathbf{p}_2, \mathbf{p}_3)

  where :math:`U` is a rotation matrix whose columns are the
  eigenvectors that form the principal axes, :math:`\Lambda` is the
  diagonal matrix of eigenvalues (sorted from largest to smallest)
  known as the principal moments of inertia, and :math:`I =
  \sum_{i=1}^{N} m_i [(\mathbf{r}_i\cdot\mathbf{r}_i)\sum_{\alpha=1}^3
  \mathbf{e}_\alpha \otimes \mathbf{e}_\alpha - \mathbf{r}_i
  \otimes \mathbf{r}_i]` is the tensor of inertia.


.. _AdK-domains:

Exercises 3
===========

.. image:: /figs/angle_defs.*
   :width: 40%
   :align: right

AdK consists of three domains:

* *CORE* residues 1-29, 60-121, 160-214 (gray)
* *NMP* residues 30-59 (blue)
* *LID* residues 122-159 (yellow)

1. Calculate the center of mass and the center of geometry for each of
   the three domains. 

   * What are the distances between the centers of mass? 

     (Hint: you can use :func:`numpy.linalg.norm` or use a function
     like :func:`veclength` that you defined previously)

   * Does it matter to use center of mass vs center of geometry?

.. image:: /figs/6b_angle_def_open.*
   :width: 40%
   :align: right

AdK undergoes a conformational transition during which CORE and LID
move relative to each other. The movement can be characterized by two
angles, :math:`\theta_\text{NMP}` and :math:`\theta_\text{LID}`, which
are defined between the *centers of geometry* of the *backbone and*
:math:`\text{C}_\beta` atoms between groups of residues
[Beckstein2009]_:

definition of :math:`\theta_\text{NMP}`
   A: 115-125, B: 90-100, C: 35-55

definition of :math:`\theta_\text{LID}`
  A: 179-185, B: 115-125, C: 125-153

The angle between vectors :math:`\vec{BA}` and :math:`\vec{BC}` is 

.. math::

   \theta = \arccos\left(\frac{\vec{BA}\cdot\vec{BC}}{|\vec{BA}||\vec{BC}|}\right)

2. Write a function :func:`theta_NMP` that takes a :class:`Universe`
   as an argument and computes :math:`\theta_\text{NMP}`:

   .. function:: theta_NMP(u)
 
      Calculate the NMP-CORE angle for E. coli AdK in degrees from
      :class:`~MDAnalysis.core.AtomGroup.Universe` *u*      

   Use the following **incomplete** code as a starting point::

     import numpy as np
     from np.linalg import norm

     def theta_NMP(u):
        """Calculate the NMP-CORE angle for E. coli AdK in degrees"""
	A = u.selectAtoms("resid 115:125 and (backbone or name CB)").centerOfGeometry()
	B = 
	C = 
	BA = A - B
	BC = 
        theta = np.arccos( 
        return np.rad2deg(theta)

   Write the function in a file :file:`adk.py` and use
   :program:`ipython` :code:`%run adk.py` to load the function while
   working on it.

   Test it on the AdK simulation (actually, the first frame)::
     
     >>> theta_NMP(u)
     44.124821
      
3. Add the corresponding function :func:`theta_LID` to :file:`adk.py`.

   Test it::
  
     >>> theta_LID(u)
     107.00881


.. _processing-atomgroups:

Processing AtomGroups
=====================

You can directly write a :class:`~MDAnalysis.core.AtomGroup.AtomGroup`
to a file with the :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.write`
method::

   CORE = u.selectAtoms("resid 1:29 or resid 60:121 or resid 160:214")
   CORE.write("AdK_CORE.pdb")

(The extension determines the file type.)

You can do fairly complicated things on the fly, such as writing the
hydration shell around a protein to a file ::

   u.selectAtoms("byres (name OW and around 4.0 protein)").write("hydration_shell.pdb")

for further analysis or visualization.

You can also write Gromacs_ index files (in case you don't like
:program:`make_ndx`...) with the
:meth:`~MDAnalysis.core.AtomGroup.AtomGroup.write_selection` method::

  CORE.write_selection("CORE.ndx", name="CORE")
  

.. SeeAlso::  The lists of supported

   * `coordinate file formats`_
   * `selection formats`_


.. _coordinate file formats: 
   http://pythonhosted.org/MDAnalysis/documentation_pages/coordinates/init.html#id1
.. _selection formats:
   http://pythonhosted.org/MDAnalysis/documentation_pages/selections_modules.html#selection-exporters
.. _Gromacs: http://www.gromacs.org
