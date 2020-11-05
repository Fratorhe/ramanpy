Post-processing Example
========================

We do this directly for multiple files.

We first prepare the inputs. We need the files that have been processes, the sample names, and the peak names.

.. code-block:: python

    file_names = ['TC4_1', 'TC4_2', 'TC4_3']
    sample_names = ['1', '2', '3']
    peaks_names = ['D','G', "G'", 'D+G']

Then, the output files. We call them like this:

.. code-block:: python

    # output files
    intensity_ratios_file = 'intensity_ratios.cfg'
    results_table = 'table_results.csv'
    equivalent_La_file = 'equivalent_La.cfg'

We instantiate the results class, which takes care of reading the results files:

.. code-block:: python

    ## create the results object
    results = ResultsDataFrames(file_names, peaks_names=peaks_names, sample_names=sample_names)

Finally, we perform the calculations we need:

.. code-block:: python

    # compute the intensity ratio D/G for each sample and spit it in a file
    intensity_ratios = results.compute_intensity_ratio_each_sample(file_to_save=intensity_ratios_file)
    # get access to the dataframe of results
    data_pandas = results.data_pandas
    results.add_xypositions() # adds the position where the sample was taken.
    # write results to a file
    results.to_csv(results_table)

    # compute the equivalent La for each sample
    results.compute_equivalent_La(Lambda=532, file_to_save=equivalent_La_file)