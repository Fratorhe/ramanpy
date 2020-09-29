import os
from itertools import takewhile

import pandas as pd
from .generic_fit_class import GenericFit
from .tools import cleanup_header
from scipy.signal import detrend


class RamanFit(GenericFit):
    """
    A class to fit raman spectra

    ...

    Attributes
    ----------
    data : df
        experimental data
    metadata : list
        metadata from the experiments
    peaks : list
        list of peaks to be fit
    other_data : configObj
        as it says...other data that could be parsed, so far, almost empty
    folder_out: str
        folder to save the reports from the fit
    var_x: str
        variable name for x, usually wavenumber, for plot only
    var_y: str
        variable name for y, usually intensity, for plot only
    x: list
        data x, in this case, wavenumber
    y: list
        data y, in this case, intensity

    Methods
    -------
    read_data_raman(file_to_analyze):
        reads the data and dumps in pandas dataframe
    fit_lorentzians():
        fits the lorentzians (from add_peak) to the data. Returns plot, print and report file.
    """

    def __init__(self, file_to_analyze, peaks, other_data=None, folder_out=None):
        experimental_data, metadata = self.read_data_raman(file_to_analyze)
        super().__init__(experimental_data=experimental_data, peaks=peaks, other_data=other_data, folder_out=folder_out)

        self.var_x = 'Wavenumber, cm$^{-1}$'  # for plots
        self.var_y = 'Intensity, -'
        self.x = self.experimental_data['wavenumber'].values
        self.y = self.experimental_data['intensity'].values
        self.filename = file_to_analyze.split(".")[0]  # remove the extension

    @staticmethod
    def read_data_raman(file_to_analyze):
        """
        read data and put in a dataframe, and metadata
        :param file_to_analyze: filename
        :return: pandas df and metadata
        """
        data, metadata = RamanFit.reader_single_point(file_to_analyze)
        return data, metadata

    @staticmethod
    def read_header(filename):
        '''
        Function to read the header of a file whose comments start with #
        :param filename: str name of file
        :return: header data cleaned up
        '''
        with open(filename, 'r', errors='ignore') as myfile:
            headiter = takewhile(lambda s: s.startswith('#'), myfile)
            header = cleanup_header(headiter)
        header = dict(element.split('=', 1) for element in header)
        return header


    @staticmethod
    def reader_single_point(filename, normalize=False, remove_offset=False):
        """
        Reader for raman spectrometry files of a single location.
        :param filename: str with filename
        :param normalize: bool. Normalize to maximum (max of the counts will be 1)
        :param remove_offset: bool. Remove linear offset using detrend from scipy. Mostly for visualization
        :return:
            :data: pandas dataframe with two columns: wavenumber and intensity
            :header: metadata from the measurement
        """
        data = pd.read_csv(filename, comment='#', sep='\t', index_col=False, names=['wavenumber', 'intensity'])
        header = RamanFit.read_header(filename=filename)
        if normalize:
            data.intensity = data.intensity.apply(lambda x: x / data.intensity.max())

        if remove_offset:
            data.intensity = detrend(data.intensity, type='linear')

        return data, header

    def set_tolerances_fit(self):
        min_max_amplitude = self.try_get_other_data(self.other_data, 'min_max_amplitude', default_value=(0, 200))
        min_max_sigma = self.try_get_other_data(self.other_data, 'min_max_sigma', default_value=(0, 200))
        tolerance_center = self.try_get_other_data(self.other_data, 'peak_center_tolerance', default_value=(10,))[0]
        amplitude = self.try_get_other_data(self.other_data, 'amplitude', default_value=(10,))[0]
        sigma = self.try_get_other_data(self.other_data, 'sigma', default_value=(10,))[0]

        self.dict_tolerances_fit = {
            'min_max_amplitude': min_max_amplitude,
            'min_max_sigma': min_max_sigma,
            'tolerance_center': tolerance_center,
            'amplitude': amplitude,
            'sigma': sigma
        }


class XRDFit(GenericFit):
    """
    A class to fit XRD spectra

    ...

    Attributes
    ----------
    data : df
        experimental data
    metadata : list
        metadata from the experiments
    peaks : list
        list of peaks to be fit
    other_data : configObj
        as it says...other data that could be parsed, so far, almost empty
    folder_out: str
        folder to save the reports from the fit
    var_x: str
        variable name for x, usually angle, for plot only
    var_y: str
        variable name for y, usually intensity, for plot only
    x: list
        data x, in this case, angle
    y: list
        data y, in this case, intensity

    Methods
    -------
    read_data_XRD(file_to_analyze):
        reads the data and dumps in pandas dataframe
    fit_lorentzians():
        fits the lorentzians (from add_peak) to the data. Returns plot, print and report file.
    """

    def __init__(self, file_to_analyze, peaks, other_data=None, folder_out=None):
        experimental_data, metadata = self.read_data_xrd(file_to_analyze)
        super().__init__(experimental_data=experimental_data, peaks=peaks, other_data=other_data, folder_out=folder_out)

        self.var_x = '$2-\\theta$, deg'
        self.var_y = 'Intensity, -'
        self.x = self.experimental_data['angle'].values
        self.y = self.experimental_data['intensity'].values
        self.filename = file_to_analyze.split(".")[0]  # remove the extension

    @staticmethod
    def read_data_xrd(filename, normalize=False):
        """
        Reader for XRD files.
        :param normalize: bool
        :param filename: str with filename
        :return:
            :data: pandas dataframe with two columns: angle and intensity
            :header: metadata from the measurement
        """
        data = pd.read_csv(filename, skiprows=1, comment='#', delim_whitespace=True,
                           index_col=False, names=['angle', 'intensity'])
        header = os.path.splitext(filename)[0]

        if normalize:
            data.intensity = data.intensity.apply(lambda x: x / data.intensity.max())

        return data, header

    def set_tolerances_fit(self):
        min_max_amplitude = self.try_get_other_data(self.other_data, 'min_max_amplitude', default_value=(0, 10))
        min_max_sigma = self.try_get_other_data(self.other_data, 'min_max_sigma', default_value=(0, 10))
        tolerance_center = self.try_get_other_data(self.other_data, 'peak_center_tolerance', default_value=(5,))[0]
        amplitude = self.try_get_other_data(self.other_data, 'peak_center_tolerance', default_value=(10,))[0]
        sigma = self.try_get_other_data(self.other_data, 'peak_center_tolerance', default_value=(10,))[0]

        self.dict_tolerances_fit = {
            'min_max_amplitude': min_max_amplitude,
            'min_max_sigma': min_max_sigma,
            'tolerance_center': tolerance_center,
            'amplitude': amplitude,
            'sigma': sigma
        }