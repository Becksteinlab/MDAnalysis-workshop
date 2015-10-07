.. -*- coding: utf-8 -*-

Testing the installation
========================

MDAnalysis comes with over 2500 `test cases`_ that check its
functionality. These test cases can be run with the command ::

  python -c 'from MDAnalysis.tests import test; test(label="full",  \
             extra_argv=["--exe", "--processes=4", "--process-timeout=120"])'

Depending on how many parallel processes (4 in the example) you can use, this
can take from a minute to 20 minutes. Ideally, you should only get
passing tests ("ok" or just a single dot ".") or "KnownFailures" (K).

.. Note:: 

   The test suite consumes a considerable amount of memory (> 4GB) and thus it
   might fail or become very slow on machines with insufficient memory such as
   a virtual machine. (This is a known problem with the test suite and will be
   addressed in the future.)

.. _test cases: http://wiki.mdanalysis.org/UnitTests
