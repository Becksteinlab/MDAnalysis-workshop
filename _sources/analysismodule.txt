.. -*- coding: utf-8 -*-

=======================================
 Using the MDAnalysis.analysis modules
=======================================

MDAnalysis comes with a number of existing analysis code in the
`MDAnalysis.analysis`_ module and `example scripts`_ (see also the
`Examples`_ on the MDAnalysis wiki).


RMSD
====

As an example we will use the :func:`MDAnalysis.analysis.rms.rmsd`
function from the :mod:`MDAnalysis.analysis.rms` module. It computes
the coordinate root mean square distance between two sets of
coordinates. For example for the AdK trajectory the backbone RMSD
between first and last frame is ::

    >>> import MDAnalysis.analysis.rms
    >>> u = MDAnalysis.Universe(PSF, DCD)
    >>> bb = u.select_atoms('backbone')
    >>> A = bb.positions  # coordinates of first frame
    >>> u.trajectory[-1]      # forward to last frame
    >>> B = bb.positions  # coordinates of last frame
    >>> MDAnalysis.analysis.rms.rmsd(A,B)
    6.8342494129169804


Superposition of structure
==========================

In order to superimpose two structures in a way that minimizes the
RMSD we have functions in the :mod:`MDAnalysis.analysis.align` module.

The example uses files provided as part of the MDAnalysis test suite
(in the variables :data:`~MDAnalysis.tests.datafiles.PSF`,
:data:`~MDAnalysis.tests.datafiles.DCD`, and
:data:`~MDAnalysis.tests.datafiles.PDB_small`). For all further
examples execute first ::

   >>> import MDAnalysis
   >>> from MDAnalysis.analysis import align, rms
   >>> from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small


In the simplest case, we can simply calculate the C-alpha RMSD between
two structures, using :func:`~MDAnalysis.analysis.rms.rmsd`::

   >>> ref = MDAnalysis.Universe(PDB_small)
   >>> mobile = MDAnalysis.Universe(PSF,DCD)
   >>> rms.rmsd(mobile.atoms.CA.positions, ref.atoms.CA.positions)
   18.858259026820352

Note that in this example translations have not been removed. In order
to look at the pure rotation one needs to superimpose the centres of
mass (or geometry) first:

   >>> ref0 =  ref.atoms.CA.positions - ref.atoms.CA.center_of_mass()
   >>> mobile0 =  mobile.atoms.CA.positions - mobile.atoms.CA.center_of_mass()
   >>> rms.rmsd(mobile0, ref0)
    6.8093965864717951

The rotation matrix that superimposes *mobile* on *ref* while
minimizing the CA-RMSD is obtained with the
:func:`~MDAnalysis.analysis.align.rotation_matrix` function ::

   >>> R, rmsd = align.rotation_matrix(mobile0, ref0)
   >>> print rmsd
   6.8093965864717951
   >>> print R
   [[ 0.14514539 -0.27259113  0.95111876]
    [ 0.88652593  0.46267112 -0.00268642]
    [-0.43932289  0.84358136  0.30881368]]

Putting all this together one can superimpose all of *mobile* onto *ref*::

   >>> mobile.atoms.translate(-mobile.atoms.CA.center_of_mass())
   >>> mobile.atoms.rotate(R)
   >>> mobile.atoms.translate(ref.atoms.CA.center_of_mass())
   >>> mobile.atoms.write("mobile_on_ref.pdb")



Exercise 5
==========

Use the above in order to investigate how rigid the :ref:`CORE, NMP,
and LID domains <AdK-domains>` are during the transition: Compute time
series of the CA RMSD of each domain relative to its own starting
structure, when superimposed on the starting structure.

*  You will need to make a copy of the starting *reference*
   coordinates that are needed for the shifts, e.g. ::

     NMP = u.select_atoms("resid 30:59")
     u.trajectory[0]   # make sure to be on initial frame
     ref_com = NMP.select_atoms("name CA").center_of_mass()
     ref0 = NMP.positions - ref_com

   which is then used instead of ``ref.atoms.CA.center_of_mass()``
   (which would *change* for each time step).

* I suggest writing a function that does the superposition for a given
  time step, reference, and mobile :class:`AtomGroup` to make the code
  more manageable (or use :func:`MDAnalysis.analysis.align.alignto`)

.. rubric:: Possible solution

.. image:: /figs/AdK_domain_rigidity.*
   :width: 50%
   :align: center

The code contains a function :func:`superpose` and :func:`rmsd`. The
latter is marginally faster because we only need the calculated RMSD
and not the full rotation matrix. (We are calling the lower-level
function :func:`MDAnalysis.core.qcprot.CalcRMSDRotationalMatrix`
directly, which has somewhat non-intuitive calling
conventions). :func:`superpose` also does the superposition of the
mobile group to the references (but
:func:`~MDAnalysis.analysis.align.alignto` is actually a more flexible
tool for doing this). Otherwise it is mostly book-keeping, which is
solved by organizing everything in dictionaries with keys "CORE",
"NMP", "LID".

.. literalinclude:: /code/domrigid.py
   :linenos:

.. _MDAnalysis.analysis: http://pythonhosted.org/MDAnalysis/documentation_pages/analysis_modules.html
.. _Examples: 
   http://wiki.mdanalysis.org/Examples
.. _example scripts:
   https://github.com/MDAnalysis/mdanalysis/tree/develop/package/examples


