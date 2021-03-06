from ramanpy import RamanFit, XRDFit


def raman_fit_carbon(file_to_analyze, file_peaks):
    """
    Runner for carbon typical Raman case.

    Parameters
    ------------
    file_to_analyze: str
        name of the file to analyze
    file_peaks: str
        name of the peaks file, if not provided, use the default one.
    """
    # get already the default peaks file in case...
    default_peaks_file = 'raman_linear_carbon.ini'

    peaks = RamanFit.read_peaks_configfile(file_peaks, default_peaks_file=default_peaks_file)
    other_data = RamanFit.read_otherdata_configfile(file_peaks)
    raman_carbon = RamanFit(file_to_analyze=file_to_analyze, peaks=peaks, other_data=other_data)

    raman_carbon.apply_smoothing()
    raman_carbon.apply_normalize()
    raman_carbon.set_tolerances_fit()
    raman_carbon.build_fitting_model_peaks()
    raman_carbon.run_fit_model()
    raman_carbon.plot_results()
    raman_carbon.save_results()


def xrd_fit_carbon(file_to_analyze, file_peaks):
    """
    Runner for carbon typical XRD case.

    Parameters
    ------------
    file_to_analyze: str
        name of the file to analyze
    file_peaks: str
        name of the peaks file, if not provided, use the default one.
    """
    # get already the default peaks file in case...
    default_peaks_file = 'XRD_linear_carbon.ini'

    peaks = XRDFit.read_peaks_configfile(file_peaks, default_peaks_file=default_peaks_file)
    other_data = XRDFit.read_otherdata_configfile(file_peaks,default_config_file=default_peaks_file)
    xrd_carbon = XRDFit(file_to_analyze=file_to_analyze, peaks=peaks, other_data=other_data)

    xrd_carbon.apply_smoothing()
    xrd_carbon.apply_normalize()
    xrd_carbon.set_tolerances_fit()
    xrd_carbon.build_fitting_model_peaks()
    xrd_carbon.run_fit_model()
    xrd_carbon.plot_results()
    xrd_carbon.save_results()


