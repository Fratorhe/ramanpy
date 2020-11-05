================
Fitting Classes
================
The main class is the ramanpy.GenericFit, which implements all the required methods to fit spectra.

.. autoclass:: ramanpy.GenericFit
    :members:

However, other classes inherit from it: RamanFit and XRDFit. The main changes are the reading of the experimental
files and the bounds used in the lorentzians.

.. autoclass:: ramanpy.RamanFit
    :members:

.. autoclass:: ramanpy.XRDFit
    :members:

