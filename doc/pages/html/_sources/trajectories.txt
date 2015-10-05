.. -*- coding: utf-8 -*-

=====================
 Trajectory analysis
=====================

The :class:`~MDAnalysis.core.AtomGroup.Universe` binds together the
static **topology** (which atoms, how are they connected, what
un-changing properties do the atoms possess (such as partial charge),
...) and the changing coordinate information, which is stored in the
**trajectory**.

The length of a trajectory (number of frames) is ::

  len(u.trajectory)

The standard way to assess each time step (or frame) in a trajectory
is to *iterate* over the :attr:`Universe.trajectory` attribute (which
is an instance of :class:`~MDAnalysis.coordinates.base.Reader` class)::

  for ts in u.trajectory:
      print("Frame: %5d, Time: %8.3f ps" % (ts.frame, u.trajectory.time))
      print("Rgyr: %g A" % (u.atoms.radius_of_gyration(), ))

The :attr:`~MDAnalysis.coordinates.base.Reader.time` attribute
contains the *current time step*. The
:class:`~MDAnalysis.coordinates.base.Reader` **only contains
information about one time step**: imagine a cursor or pointer moving
along the trajectory file. Where the cursor points, there's you
current coordinates, frame number, and time.

Normally you will collect the data in a list or array, e.g. ::

  Rgyr = []
  protein = u.select_atoms("protein")
  for ts in u.trajectory:
     Rgyr.append((u.trajectory.time, protein.radius_of_gyration()))
  Rgyr = np.array(Rgyr)

.. Note:: 

   It is important to note that the coordinates and related properties
   calculated from the coordinates such as the radius of gyration
   *change* while selections_ on static properties (such as
   :code:`protein` in the example) do not change when moving through a
   trajectory: You can define the selection *once* and then
   recalculate the property of interest for each frame of the
   trajectory.

   However, if selections contain distance-dependent queries (such as
   ``around`` or ``point``, see `selection keywords`_ for more
   details) then one might have to recalculate the selection for each
   time step and one would put it inside the loop over frames.

.. _selections: 
   http://pythonhosted.org/MDAnalysis/documentation_pages/selections.html

.. _selection keywords:
   http://pythonhosted.org/MDAnalysis/documentation_pages/selections.html#selection-keywords

The data can be plotted to give the graph below::

  # quick plot
  import matplotlib.pyplot as plt
  ax = plt.subplot(111)
  ax.plot(Rgyr[:,0], Rgyr[:,1], 'r--', lw=2, label=r"$R_G$")
  ax.set_xlabel("time (ps)")
  ax.set_ylabel(r"radius of gyration $R_G$ ($\AA$)")
  ax.figure.savefig("Rgyr.pdf")
  plt.draw()

The :math:`R_G(t)` increases over time, indicating an opening up of
the AdK enzyme.

.. image:: /figs/Rgyr.*
   :width: 40%
   :align: center



Exercises 4
===========

1. Take the functions to calculate :math:`\theta_\text{NMP}` and
   :math:`\theta_\text{LID}` ::

      import numpy as np
      from numpy.linalg import norm

      def theta_NMP(u):
	  """Calculate the NMP-CORE angle for E. coli AdK in degrees"""
	  C = u.select_atoms("resid 115:125 and (backbone or name CB)").center_of_geometry()
	  B = u.select_atoms("resid 90:100 and (backbone or name CB)").center_of_geometry()
	  A = u.select_atoms("resid 35:55 and (backbone or name CB)").center_of_geometry()
	  BA = A - B
	  BC = C - B
	  theta = np.arccos(np.dot(BA, BC)/(norm(BA)*norm(BC)))
	  return np.rad2deg(theta)

      def theta_LID(u):
	  """Calculate the LID-CORE angle for E. coli AdK in degrees"""
	  C = u.select_atoms("resid 179:185 and (backbone or name CB)").center_of_geometry()
	  B = u.select_atoms("resid 115:125 and (backbone or name CB)").center_of_geometry()
	  A = u.select_atoms("resid 125:153 and (backbone or name CB)").center_of_geometry()
	  BA = A - B
	  BC = C - B
	  theta = np.arccos(np.dot(BA, BC)/(norm(BA)*norm(BC)))
	  return np.rad2deg(theta)

   and calculate the time series :math:`\theta_\text{NMP}(t)` and
   :math:`\theta_\text{LID}(t)`.

   Plot them together in one plot.

2. Plot :math:`\theta_\text{NMP}(t)` against
   :math:`\theta_\text{LID}(t)`. 

   What does the plot show?

   Why could such a plot be useful?

.. image:: /figs/NMP_LID_angle_projection.*
   :width: 50%
   :align: center

The code to generate the figure contains :func:`theta_LID` and
:func:`theta_NMP`. 

.. literalinclude:: /code/adk.py
   :emphasize-lines: 30-32
   :linenos:

Note that one would normally write the code more efficiently and
generate the atom groups once and then pass them to a simple function
to calculate the angle ::

   def theta(A, B, C):
       """Calculate the angle between BA and BC for AtomGroups A, B, C"""
       B_center = B.centroid()
       BA = A.centroid() - B_center
       BC = C.centroid() - B_center
       theta = np.arccos(np.dot(BA, BC)/(norm(BA)*norm(BC)))
       return np.rad2deg(theta)
  


Bells and whistles
==================

.. rubric:: Quick data acquisition

Especially useful for interactive analysis in :program:`ipython
--pylab` using list comprehensions (implicit for loops)::

  protein = u.select_atoms("protein")
  data = np.array([(u.trajectory.time, protein.radius_of_gyration()) for ts in u.trajectory])
  time, RG = data.T
  plot(time, RG)  

.. rubric:: More on the trajectory iterator

One can directly **jump to a frame** by using "indexing syntax":

  >>> u.trajectory[50]
  < Timestep 51 with unit cell dimensions array([  0.,   0.,   0.,  90.,  90.,  90.], dtype=float32) >
  >>> ts.frame
  50

You can also **slice trajectories**, e.g. if you want to start at the 10th
frame and go to 10th before the end, and only use every 5th frame::

  for ts in u.trajectory[9:-10:5]:
      print(ts.frame)
      ...

.. Note:: 

   Trajectory indexing and slicing uses 0-based indices (as in
   standard Python) and MDAnalysis also numbers frames starting
   with 0. Thus the "tenth frame" in a trjectory has ``ts.frame ==
   9``.

   
.. Note::

   Not all trajectory readers support direct access and arbitrary
   slices, although many commonly used ones such as DCD, XTC/TRR, and
   Amber NETCDF do.

.. SeeAlso:: One can iterate through multiple trajectories in parallel
             with the help of :func:`itertools.izip` from the
             :mod:`itertools` module, which also provide other
             interesting ways to work with trajectory iterators.
