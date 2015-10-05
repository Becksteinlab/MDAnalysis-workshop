.. -*- coding: utf-8 -*-

.. _chapter-basics:

========
 Basics
========

We discuss the fundamental objects in MDAnalysis, the
:ref:`universe-and-atomgroup`, and the facilities for
:ref:`selections` of atoms. These selections themselves are again an
:class:`~MDAnalysis.core.AtomGroup.AtomGroup`.


.. _universe-and-atomgroup:

Universe and AtomGroup
======================

MDAnalysis is **object oriented**. Molecular systems consist of
:class:`~MDAnalysis.core.AtomGroup.Atom` objects (instances of the
class :class:`MDAnalysis.core.AtomGroup.Atom`), which are grouped in
:class:`~MDAnalysis.core.AtomGroup.AtomGroup` instances. You build the
:class:`~MDAnalysis.core.AtomGroup.AtomGroup` of your system by
loading a **topology** (list of atoms and possibly their connectivity)
together with a **trajectory** (coordinate information) into the
central data structure, the
:class:`~MDAnalysis.core.AtomGroup.Universe` object::

  >>> u = MDAnalysis.Universe(PSF, DCD)
  >>> print(u)
  <Universe with 3341 atoms>

The atoms are stored in the attribute
:attr:`~MDAnalysis.core.AtomGroup.Universe.atoms` of the
:class:`MDAnalysis.core.AtomGroup.Universe`::

  >>> print(u.atoms)
  <AtomGroup with 3341 atoms>
  >>> list(u.atoms[:5])
  [< Atom 1: name 'N' of type '56' of resname 'MET', resid 1 and segid '4AKE'>,
   < Atom 2: name 'HT1' of type '2' of resname 'MET', resid 1 and segid '4AKE'>,
   < Atom 3: name 'HT2' of type '2' of resname 'MET', resid 1 and segid '4AKE'>,
   < Atom 4: name 'HT3' of type '2' of resname 'MET', resid 1 and segid '4AKE'>,
   < Atom 5: name 'CA' of type '22' of resname 'MET', resid 1 and segid '4AKE'>]

Any :class:`~MDAnalysis.core.AtomGroup.AtomGroup` knows the residues
that the atoms belong to via the attribute
:attr:`~MDAnalysis.core.AtomGroup.AtomGroup.residues`, which produces a
:class:`~MDAnalysis.core.AtomGroup.ResidueGroup`. A
:class:`~MDAnalysis.core.AtomGroup.ResidueGroup` acts like a list of
:class:`~MDAnalysis.core.AtomGroup.Residue` objects::

  >>> u.atoms[100:130].residues
  <ResidueGroup [<Residue 'LEU', 6>, <Residue 'GLY', 7>, <Residue 'ALA', 8>]>

Larger organizational units are
:class:`~MDAnalysis.core.AtomGroup.Segment` instances, for example one
protein or all the solvent molecules or simply the whole
system. :class:`~MDAnalysis.core.AtomGroup.Atom`,
:class:`~MDAnalysis.core.AtomGroup.AtomGroup`,
:class:`~MDAnalysis.core.AtomGroup.Residue`, and
:class:`~MDAnalysis.core.AtomGroup.ResidueGroup` have an
attribute :attr:`~MDAnalysis.core.AtomGroup.AtomGroup.segments` that
will list the segment IDs ("segids") as a
:class:`~MDAnalysis.core.AtomGroup.SegmentGroup`::

  >>> u.atoms.segments
  <SegmentGroup [<Segment '4AKE'>]>

The converse is also true: each "higher" level in the hierarchy also
know about the :class:`~MDAnalysis.core.AtomGroup.Residue` and
:class:`~MDAnalysis.core.AtomGroup.Atom` instances it contains. For
example, to list the atoms of the
:class:`~MDAnalysis.core.AtomGroup.ResidueGroup` we had before::

  >>> r = u.atoms[100:130].residues
  >>> r.atoms
  <AtomGroup with 36 atoms>


Exercise 1
----------

1. What residue ("resname") does the last atom belong to in the above
   example? ::

    >>> r = u.atoms[100:130].residues
    >>> r.atoms[-1]
    < Atom 136: name 'O' of type '70' of resname 'ALA', resid 8 and segid '4AKE'>   
 

2. Why does the expression ::

     len(u.atoms[100:130]) == len(u.atoms[100:130].residues.atoms)
   
   return ``False``?

   Because the complete residues contain more atoms than the arbitrary
   slice of atoms.

3. How many residues are in the
   :class:`~MDAnalysis.core.AtomGroup.AtomGroup.Universe` ``u``? ::

     >>> len(u.atoms.residues)
     >>> u.atoms.n_residues
     214

   How do you get a list of the residue names (such as ``["Ala",
   "Gly", "Gly", "Asp", ...]``) and residue numbers ("resid") for
   atoms 1000 to 1300? And as a list of tuples ``(resname, resid)``
   (Hint: :func:`zip`)?::

     >>> resnames = u.atoms[999:1300].residues.resnames
     >>> resids = u.atoms[999:1300].residues.resids
     >>> zip(resnames, resids)

   How do you obtain the resid and the resname for the 100th residue?
   (Hint: investigate the :class:`~MDAnalysis.core.AtomGroup.Residue`
   object interactively with :kbd:`TAB` completion) ::

     >>> r100 = u.atoms.residues[99]
     >>> print(r100.id, r100.name)
     100 GLY


4. How many segments are there?  ::

     >>> len(u.segments)
     >>> len(u.atoms.segments)
     >>> u.atoms.n_segments
     1

   What is the segment identifier of the first
   :class:`~MDAnalysis.core.AtomGroup.Segment`? ::

     >>> s1 = u.segments[0]
     >>> s1.id
     '4AKE'
   

.. SeeAlso:: 

   Methods of :class:`~MDAnalysis.core.AtomGroup.AtomGroup`,
   :class:`~MDAnalysis.core.AtomGroup.ResidueGroup`, and
   :class:`~MDAnalysis.core.AtomGroup.SegmentGroup`
           
   * :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.numberOfResidues` and 
     :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.numberOfAtoms`
   * :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.resids`
   * :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.resnames`


.. _selections:

Selections
==========

.. TODO: named selections

MDAnalysis comes with a fairly complete `atom selection`_
facility. Primarily, one uses the method
:meth:`~MDAnalysis.core.AtomGroup.Universe.select_atoms` of a
:class:`~MDAnalysis.core.AtomGroup.Universe`::

  >>> CA = u.select_atoms("protein and name CA")
  >>> CA
  >>> <AtomGroup with 214 atoms>

but really any :class:`~MDAnalysis.core.AtomGroup.AtomGroup` has a
:meth:`~MDAnalysis.core.AtomGroup.AtomGroup.select_atoms` method::

  >>> acidic = CA.select_atoms("resname ASP or resname GLU")
  >>> acidic
  >>> <AtomGroup with 35 atoms>
  >>> acidic.residues
  <ResidueGroup [<Residue 'GLU', 22>, <Residue 'ASP', 33>, <Residue 'GLU', 44>, <Residue 'ASP', 51>, <Residue 'ASP', 54>, <Residue 'ASP', 61>, <Residue 'GLU', 62>, <Residue 'GLU', 70>, <Residue 'GLU', 75>, <Residue 'ASP', 76>, <Residue 'ASP', 84>, <Residue 'ASP', 94>, <Residue 'GLU', 98>, <Residue 'ASP', 104>, <Residue 'GLU', 108>, <Residue 'ASP', 110>, <Residue 'ASP', 113>, <Residue 'GLU', 114>, <Residue 'ASP', 118>, <Residue 'GLU', 143>, <Residue 'ASP', 146>, <Residue 'ASP', 147>, <Residue 'GLU', 151>, <Residue 'GLU', 152>, <Residue 'ASP', 158>, <Residue 'ASP', 159>, <Residue 'GLU', 161>, <Residue 'GLU', 162>, <Residue 'GLU', 170>, <Residue 'GLU', 185>, <Residue 'GLU', 187>, <Residue 'ASP', 197>, <Residue 'GLU', 204>, <Residue 'ASP', 208>, <Residue 'GLU', 210>]>
  
.. SeeAlso:: All the `selection keywords`_ are described in the documentation.

Selections can be combined with boolean expression and it is also
possible to select by geometric criteria, e.g. with the :samp:`around
{distance} {selection}` keyword::

  u.select_atoms("((resname ASP or resname GLU) and not (backbone or name CB or name CG)) \
                   and around 4.0 ((resname LYS or resname ARG) \
                                    and not (backbone or name CB or name CG))").residues

This selection will find atoms potentially involved in salt bridges
between acidic and basic residues.


Exercises 2
-----------

1. Select the range of resids 100 to 200 ("100-200") with a
   selection. Compare the result to what you get by slicing the
   :attr:`u.atoms.residues` appropriately.

   Which approach would you prefer to use in a analysis script?

   Solution::

      >>> u.select_atoms("resid 100-200")
      <AtomGroup with 1609 atoms>

   Compare to the slicing solution (doing an element-wise comparison,
   i.e. residue by residue in each :func:`list`)::

      >>> list(u.select_atoms("resid 100-200").residues) == list(u.atoms.residues[99:200])

   If one wants to get specific residues in scripts one typically uses
   selections instead of slicing because the index in the slice might
   not correspond to the actual residue ids (minus 1): If a number of
   residues (e.g. 150-160) are missing from the structure then the
   selection will simply give you residues 100-149 and 151-200 but the
   slice 99:200 would give you residues 100-149 and *161-209*.

2. Select all residues that do not contain a :math:`\mathrm{C}_\beta`
   ("CB") atom. How many are there? What residue names did you find? 

   Solution::

      >>> sel = u.select_atoms("(byres name CA) and not (byres name CB)").residues
      >>> len(sel)
      20

   These are all Glycines, as can be seen by comparing the residue
   groups element-wise::

      >>> glycines = u.select_atoms("resname GLY")
      >>> list(sel) == list(glycines.residues)
      True


.. _atom selection: 
   http://pythonhosted.org/MDAnalysis/documentation_pages/selections.html

.. _selection keywords:
   http://pythonhosted.org/MDAnalysis/documentation_pages/selections.html#selection-keywords

