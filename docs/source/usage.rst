Step by step
=====================
We will use an example of fitting Raman spectra, but XRD should be analogous.

Start by importing ramanpy as rp

.. code-block:: python

    import ramanpy as rp

Get the peaks file and the extra data if it applies.
The peaks file should looks like this:

.. literalinclude:: ../../ramanpy/peak_files/raman_linear_carbon.ini

As you can see, in this case, the extra data is also included there **[other data]**.

To read these files, we do the following:

.. code-block:: python

    file_peaks = 'raman_linear_carbon.ini'

    peaks = RamanFit.read_peaks_configfile(file_peaks)
    other_data = RamanFit.read_otherdata_configfile(file_peaks)


Then, instantiate an object of the class :meth:`ramanpy.Ramanfit`.

.. code-block:: python

    raman_carbon = RamanFit(file_to_analyze=file_to_analyze, peaks=peaks, other_data=other_data)

We can apply different data treatment such as:

.. code-block:: python

    raman_carbon.apply_smoothing()
    raman_carbon.apply_normalize()

Then, we can set the tolerances for each variable (peak center, amplitude and sigma), and build the model with the
different peaks.
Info here: `LMFIT Lorentzian model <https://lmfit.github.io/lmfit-py/builtin_models.html#lorentzianmodel>`_

.. code-block:: python

    raman_carbon.set_tolerances_fit()
    raman_carbon.build_fitting_model_peaks()

Once the model is built, we can fit it:

.. code-block:: python

        raman_carbon.run_fit_model()

And finally, plot and save the results:

.. code-block:: python

    raman_carbon.plot_results()
    raman_carbon.save_results()