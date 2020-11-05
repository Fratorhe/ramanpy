Runners
=====================

Instead of going step by step, we can also define runners.
These functions wrap the complexity of a typical case, so that user only needs a couple of lines of code.
The list of all the runners is provided below.
First, we describe here the one runner for Raman carbon material.

In this case, only the file to analyze and the peaks file are required. If the file peaks has any problem, the default
one will be read instead.
Then, the code will read the peaks file, and the extra data provided in other_data.
This will instantiate an object of the class :meth:`ramanpy.RamanFit`. Following, some methods will be run on the data
(see source code for information). Then, the fit will be performed, the figure and report also saved.

.. code-block:: python

    import ramanpy as rp
    rp.runners.raman_fit_carbon(file_to_analyze, file_peaks)

.. automodule:: ramanpy.runners
    :members: