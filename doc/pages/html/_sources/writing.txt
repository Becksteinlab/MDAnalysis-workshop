.. -*- coding: utf-8 -*-

=====================
 Writing coordinates
=====================

MDAnalysis also supports writing of data in a range of file formats
(see the `Table of supported coordinate formats`_ for
details). MDAnalysis supports both *single frame* writers (such as a
simple PDB or GRO file) and *trajectory* writers (e.g. XTC, DCD, but
also multi-frame PDB files).

.. _`Table of supported coordinate formats`:
   https://pythonhosted.org/MDAnalysis/documentation_pages/coordinates/init.html#id1

.. _writing-single-frames:

Single frames
=============

The most straightforward way to write to a file that can only hold a
single frame is to use the
:meth:`~MDAnalysis.core.AtomGroup.AtomGroup.write` method of any
:class:`~MDAnalysis.core.AtomGroup.AtomGroup` as already also shown under
:ref:`processing-atomgroups`. For instance, to only write out the
protein without solvent to a file in GRO format::

   u = MDAnalysis.Universe(PDB)
   protein = u.select_atoms("protein")
   protein.write("protein.gro")

MDAnalysis uses the file suffix to determine the output file format
(unless the *format* keyword is specified) and will raise an exception
if it is not suitable for single frame writing.


.. _writing-trajectories:

Trajectories
============

The typical use pattern is to

#. Get a trajectory writer with :func:`MDAnalysis.Writer` (which is
   the same as :func:`MDAnalysis.coordinates.core.writer`), typically
   specifying in advance how many atoms a frame will contain.
#. Use the :meth:`~MDAnalysis.coordinates.base.Writer.write` method to
   write a new time step to the trajectory.
#. Close the trajectory with
   :meth:`~MDAnalysis.coordinates.base.Writer.close` (although it is
   recommended to simply use the writer with the :keyword:`with`
   statement and have the context manager close the file automatically).

Example: Protein-only trajectory
--------------------------------

In practice, the second step is typically repeated in a loop as in the
example below::

  import MDAnalysis
  from MDAnalysis.tests.datafiles import PDB, XTC
  
  u = MDAnalysis.Universe(PDB, XTC)
  protein = u.select_atoms("protein")
  with MDAnalysis.Writer("protein.xtc", protein.n_atoms) as W:
      for ts in u.trajectory:
          W.write(protein)

The loop steps through the input trajectory frame by frame. The
coordinates of the selection (the
:class:`~MDAnalysis.core.AtomGroup.AtomGroup` ``protein``) change
accordingly and are then written as a new frame into the output
trajectory.

The output trajectory only contains the coordinates of the
protein. For this trajectory to be useful, a protein-only topology
file also has to be stored, as in the example under
:ref:`writing-single-frames`.


Example: Saving dynamic per-atom properties in B-factor
-------------------------------------------------------

It is often very useful to project per-atom properties on the
structure. A common approach is to save scalar values in the B-factor
field of a PDB file and then color atoms by B-factor (also known as
temperature factor or just "beta").

The following example computes the shift of each atom in AdK relative
to a reference structure (line 23). We take as reference the closed
conformation (after a structural superposition on the CORE domain with
:func:`~MDAnalysis.analysis.align.alignto`). The shifts are written
into the B-factor with the
:meth:`~MDAnalysis.core.AtomGroup.AtomGroup.set_bfactor` method of
:class:`~MDAnalysis.core.AtomGroup.AtomGroup`. Each frame is written
out as part of a multi-frame PDB file:

.. literalinclude:: /code/bfacmovie.py
   :linenos:
   :emphasize-lines: 22,23,25,26

To visualize in VMD_, use the :ref:`pdbbfactor Tcl script <pdbbfactor-tcl-script>` below on
the VMD Tcl commandline:

.. code-block:: tcl

  source pdbbfactor.tcl
  pdbbfactor adk_distance_bfac.pdb

Rendered snapshots from the beginning, middle, and end of the
trajectroy are shown below. Note that the tip of the LID domain moves
by almost 25 Å, which provides some justification for calling the AdK
closed/open transition a "large conformational change" [Seyler2014]_.

+----------------------------------------+----------------------------------------+----------------------------------------+
| .. image:: /figs/AdK_distance_0001.*   | .. image:: /figs/AdK_distance_0040.*   | .. image:: /figs/AdK_distance_0097.*   |
|    :scale: 30%                         |    :scale: 30%                         |    :scale: 30%                         |
+========================================+========================================+========================================+
|   AdK closed conformation.             | AdK intermediate conformation, atoms   | AdK open conformation, atoms colored   |
|                                        | colored by displacement from the closed| by displacement from the closed        |
|                                        | conformation. Color scale ranges from  | conformation. Color scale ranges from  |
|                                        | 0 Å (blue) to 25 Å (red).              | 0 Å (blue) to 25 Å (red).              |
+----------------------------------------+----------------------------------------+----------------------------------------+

.. _pdbbfactor-tcl-script:

.. rubric:: pdbbfactor Tcl script
 
`pdbbfactor`_ was originally written by Justin Gullingsrud (2004) and
slightly modified for this tutorial:

.. literalinclude:: /code/pdbbfactor.tcl
   :language: tcl
   :linenos:

.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _`pdbbfactor`: http://www.ks.uiuc.edu/Research/vmd/script_library/scripts/pdbbfactor/

.. Note:: We could have also directly loaded the scalar data into the
          ``User`` field in VMD_; this is demonstrated in
          `vmduser.py`_.

.. _`vmduser.py`: 
   https://github.com/MDAnalysis/MDAnalysisTutorial/blob/master/doc/sphinx/code/vmduser.py

