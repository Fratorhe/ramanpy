Post-processing Classes
========================

We can either process a single file, or multiple files.

Multiple files
---------------

To read multiple files, we use the class :meth:`ramanpy.ResultsDataFrames`. This class, internally
uses :meth:`ramanpy.ReadResultParamsFit` to read each result file.
The idea is to have here everything needed to compute for example the intensity ratio.

.. autoclass:: ramanpy.ResultsDataFrames
    :members:

Single file
-------------

To read the results of a single file, we use the class :meth:`ramanpy.ReadResultParamsFit`. This should work for both
Raman and XRD.
However, most of the time we will process several files at once so this is not really important.

.. autoclass:: ramanpy.ReadResultParamsFit
    :members: