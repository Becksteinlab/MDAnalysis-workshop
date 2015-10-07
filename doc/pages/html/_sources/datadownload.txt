.. -*- mode: reST; coding: utf-8 -*-

.. _chapter-datadownload:

===============
 Data Download
===============


Download a set of trajectories for the tutorial (equilibrium trajectories and
DIMS trajectories) from dropbox `CECAM_Workshop/MDAnalysis`_. You can do all of
it from the commandline with `curl`_::

  curl -o mdatrj.zip -L 'https://www.dropbox.com/sh/am6y00kac8myihe/AABDiQI28fWnRZueQTT7W2s1a?dl=1'
  unzip mdatrj.zip && rm mdatrj.zip

You should now have two directories named ``equilibrium/`` and ``dims/``
containing PSF and DCD files.

Note that this is about 318 MiB so **download well in advance of the workshop**
where (and when) you have a good internet connection. You can also use
Dropbox_, add the shared `CECAM_Workshop/MDAnalysis`_ folder to your own
Dropbox account and then work from there during the tutorial.


.. _`CECAM_Workshop/MDAnalysis`:
   https://www.dropbox.com/sh/ln0klc9j7mhvxkg/AAB0gMcPPsrDhdVrM2PWmopXa?dl=0
.. _Dropbox: https://www.dropbox.com
.. _curl: http://curl.haxx.se/
