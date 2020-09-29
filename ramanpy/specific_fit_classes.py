import os
from itertools import takewhile

import pandas as pd
from ramanpy.generic_fit_class import GenericFit
from ramanpy.tools import cleanup_header
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

        self.var_x = 'Wavenumber, cm^{-1}'  # for plots
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
    def read_xyz(filename):
        '''
        Function to read the xyz positions from the header.
        :param filename: str name of file
        :return: list
        '''
        header = RamanFit.read_header(filename=filename)
        axis = ['X', 'Y', 'Z']
        positions = []
        for element in axis:
            position = float(
                header.get(f'{element}(µm)', 'nan'))  # get the position in micrometers, otherwise put a nan.
            positions.append(position)
        return positions

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
