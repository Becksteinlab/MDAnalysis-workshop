.. -*- coding: utf-8 -*-

=========================
 Working with AtomGroups
=========================

A :class:`~MDAnalysis.core.AtomGroup.AtomGroup` has a large number of
methods attributes defined that provide information about the atoms
such as names, indices, or the coordinates in the
:attr:`~MDAnalysis.core.AtomGroup.AtomGroup.positions` attribute::

  >>> CA = u.select_atoms("protein and name CA")
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
  :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.center_of_mass` and the
  center of geoemtry (or centroid)
  :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.center_of_geometry`
  (equivalent to
  :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.centroid`);

* the total mass
  :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.total_mass`;

* the total charge
  :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.total_charge` (if partial
  charges are defined in the topology);

* the radius of gyration

  .. math:: 

        R_\mathrm{gyr} = \sqrt{\frac{1}{M}\sum_{i=1}^{N} m_i(\mathbf{r}_i - \mathbf{R})^2}

  with :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.radius_of_gyration`;

* the principal axes :math:`\mathbf{p}_1, \mathbf{p}_2, \mathbf{p}_1`
  from :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.principal_axes` via
  a diagonalization of the tensor of inertia
  :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.moment_of_inertia`,

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

        >>> domains = {
        >>>   'CORE': u.select_atoms("protein and (resid 1-29 or resid 60-121 or resid 160-214)"),
        >>>   'NMP': u.select_atoms("protein and resid 30-59"),
	>>>   'LID': u.select_atoms("protein and resid 122-159")
        >>>   }
        >>> cg = dict((name, dom.centroid()) for name,dom in domains.items())
        >>> cm = dict((name, dom.center_of_mass()) for name,dom in domains.items())
        >>> print(cg)
        {'LID': array([-15.16074944,   2.11599636,  -4.37305355], dtype=float32), 
         'CORE': array([ 4.43884087,  2.05389476,  1.63895261], dtype=float32), 
         'NMP': array([ -2.99990702, -13.62531662,  -2.93235731], dtype=float32)}
        >>> print(cm)
        {'LID': array([-15.11337499,   2.12292226,  -4.40910485]), 
         'CORE': array([ 4.564116  ,  2.08700105,  1.54992649]), 
         'NMP': array([ -3.20330174, -13.60247613,  -3.06221538])}

   * What are the distances between the centers of mass? 

     (Hint: you can use :func:`numpy.linalg.norm` or calculate it
     manually.) ::

        >>> from numpy.linalg import norm
        >>> print(norm(cm['CORE'] - cm['NMP']))
        18.1042626244
        >>> print(norm(cm['CORE'] - cm['LID']))
        20.5600339602
        >>> print(norm(cm['NMP'] - cm['LID']))
        19.7725089609
	
   * Does it matter to use center of mass vs center of geometry? 

        >>> print(norm(cg['CORE'] - cg['NMP']))
        17.9463
        >>> print(norm(cg['CORE'] - cg['LID']))
        20.501
        >>> print(norm(cg['NMP'] - cg['LID']))
        19.9437


.. image:: /figs/6b_angle_def_open.*
   :width: 40%
   :align: right

AdK undergoes a conformational transition during which the NMP and LID domain
move relative to the CORE domain. The movement can be characterized by two
angles, :math:`\theta_\text{NMP}` and :math:`\theta_\text{LID}`, which are
defined between the *centers of geometry* of the *backbone and*
:math:`\text{C}_\beta` atoms between groups of residues [Beckstein2009]_:

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

   Use the following *incomplete* code as a starting point::

     import numpy as np
     from numpy.linalg import norm

     def theta_NMP(u):
        """Calculate the NMP-CORE angle for E. coli AdK in degrees"""
	A = u.select_atoms("resid 115:125 and (backbone or name CB)").center_of_geometry()
	B = 
	C = 
	BA = A - B
	BC = 
        theta = np.arccos( 
        return np.rad2deg(theta)

   Write the function to file :file:`adk.py` and inside :program:`ipython` run
   the file with :code:`%run adk.py` to load the function while working on it.

   Test it on the AdK simulation (actually, the first frame)::
     
     >>> theta_NMP(u)
     44.124821
      
3. Add the corresponding function :func:`theta_LID` to :file:`adk.py`.

   Test it::
  
     >>> theta_LID(u)
     107.00881

(See below for a solution.)

.. rubric:: Calculation of the domain angles of AdK

.. literalinclude:: /code/adk.py
   :linenos:
   :lines: 1-23
   
   
   


.. _processing-atomgroups:

Processing AtomGroups
=====================

You can directly write a :class:`~MDAnalysis.core.AtomGroup.AtomGroup`
to a file with the :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.write`
method::

   CORE = u.select_atoms("resid 1:29 or resid 60:121 or resid 160:214")
   CORE.write("AdK_CORE.pdb")

(The extension determines the file type.)

You can do fairly complicated things on the fly, such as writing the
hydration shell around a protein to a file ::

   u.select_atoms("byres (name OW and around 4.0 protein)").write("hydration_shell.pdb")

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
