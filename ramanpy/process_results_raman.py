import string

import numpy as np
import pandas as pd
from configobj import ConfigObj

from . import RamanFit


class ReadResultParamsFit:
    """
    A class to read a result file from fitting of raman/xrd spectra
    This class is used in ResultsDataFrame to read the parameters, but can also be used independently
    ...

    Attributes
    ----------
    dict_results : dict
        dictionary with the results from the fit, read from the _params.txt file
    number_of_lorentzians: int
        number of lorentzian peaks
    peaks_names: list
        names for the peaks (typically D, G, G')
    params_of_interest: list
        parameters to be studied fwhm, center, height
    lorentzians: dict
        dict with peaks_names and values from dict_results
    lorentzians_stderr: dict
        same but the stderr, not really used.

    Methods
    -------
    better_structure_params():
        builds the dict of lorentzians and lorentzians_std with the names.
    method1():
        fits the lorentzians (from add_peak) to the data. Returns plot, print and report file.
    """

    def __init__(self, params_file, peaks_names=None):
        """

        :param params_file: filename with the results of the fit
        :param peaks_names: names for the peaks (D, G, etc)
        """
        # print(params_file)
        self.dict_results = ConfigObj(params_file, file_error=True)

        # find number of lorentzians
        keys = self.dict_results.keys()
        if peaks_names is None: # if the name of the peaks is not provided, generate some
            self.number_of_lorentzians = max(set([get_num(key) for key in keys]))
            self.peaks_names = list(string.ascii_lowercase)[0:self.number_of_lorentzians]
        else:
            self.number_of_lorentzians = len(peaks_names)
            self.peaks_names = peaks_names

        self.params_of_interest = ['fwhm', 'center', 'height']

        self.lorentzians, self.lorentzians_stderr = self.better_structure_params()

    def better_structure_params(self):
        """
        This is used to separate the average and the stderr from the reading of the lorentzians (_params.txt file)
        :return:
        """
        lorentzians = dict((key, {}) for key in self.peaks_names)
        lorentzians_stderr = dict((key, {}) for key in self.peaks_names)

        # build word
        for parameter in self.params_of_interest:
            for name, num in zip(self.peaks_names, range(1, self.number_of_lorentzians + 1)):
                key_dict = 'lz' + str(num) + parameter
                value_in_dict_results = self.dict_results[key_dict]

                lorentzians[name][parameter] = float(value_in_dict_results)

        return lorentzians, lorentzians_stderr


class ResultsDataFrame:
    """
        A class to read multiple results from fitting of raman/xrd spectra
        Performs also typical calculations to obtain the equivalent La, intensity ratios, etc.

        ...

        Attributes
        ----------
        file_names_params : list
            names of the _params.txt files
        file_names_experiment: list
            name of the base files (experimental raw data)
        sample_names: list
            same as file_names_experiment but without extension
        data_dict: dict of obj ReadResultsRamanFit
            data from all the params files read
        data_pandas: pandas
            same as data_dict, but in a pandas dataframe
        peak_names: list
            list of names for the peaks

        Methods
        -------
        to_csv():
            print the data_pandas to a csv file
        compute_statistics():
            applies a statistic (avg, std) to a given column of the dataframe defined in _apply_statistics.
            prints data to a file
        _apply_statistics:
            applies a statistic (avg, std)
        compute_intensity_ratio_each_sample:
            computes the D/G ratio for each sample
        compute_equivalent_La:
            computes the La as in Cancado et al
        add_xypositions:
            adds the xy position of the point to the data_pandas attribute
        _read_xyz:
            reads the coordinates in the experiment file to pass it to add_xypositions.
        """
    def __init__(self, file_names, peaks_names=None, sample_names=None):
        """

        :param file_names: file names without "_params.txt" this will be added inside
        :param peaks_names: for the lorentzians, if not given, it will be letters
        :param sample_names: for the dataframe, if not given, it will equal to file_names
        """

        self.file_names_params = [file_to_process + '_params.txt' for file_to_process in file_names]
        self.file_names_experiment = [file_to_process + '.txt' for file_to_process in file_names]
        if sample_names is None:
            self.sample_names = file_names
        else:
            self.sample_names = sample_names

        self.peak_names = peaks_names

        self.data_dict = {}

        for file_name_param, sample_name in zip(self.file_names_params, self.sample_names):
            print(sample_name)
            self.data_dict[sample_name] = ReadResultParamsFit(file_name_param, peaks_names=self.peak_names).lorentzians

        self.data_pandas = pd.concat({k: pd.DataFrame(v).T for k, v in self.data_dict.items()}, axis=0)

        # ensure we get the right peak_names back from the ReadResultParamsFit.
        # if they are provided it is a repetition, otherwise it will set them here, because they were set to None
        self.peak_names = self.data_pandas.index.levels[1].to_list()
        # Get the column names ie. variables to be analyzed
        self.cols_dataframe = self.data_pandas.columns.to_list()

        self.unstacked = False  # dirty trick to unstack only once

    def to_csv(self, filename):
        """
        passes the current data_pandas to a dataframe. the use of unstacked is due to the inclusion of x and y positions
         (see add_xypositions)
        :param filename: file to dump the data. typically table_results.csv
        """
        if self.unstacked:
            self.data_pandas.to_csv(filename)
        else:
            self.data_pandas.unstack().to_csv(filename)

        # TODO: re-impliment save in json instead of csv, it's just a headache.

    def compute_statistics(self, filename=None):
        """
        computes the statistics defined in another function to each
        column and each peak for the different samples
        it will print them to a file and return them here as a dictionary
        :param filename: file to be saved
        :return: dict of data
        """

        # create the idx variable for better indexing in multiindex pandas
        idx = pd.IndexSlice  # this helps for the multiindexing

        self.dict_stats = {}
        # compute for each peak
        for peak in self.peak_names:
            # compute for each variable of interest (center, height, fwhm)
            dict_stats_data = {}
            for column in self.cols_dataframe:
                # grab the data
                data = self.data_pandas.loc[idx[:, peak], column].values

                # compute statistics and put them in a dict
                dict_stats_data[column] = self._apply_statistics(data)

            # combine all the stats into one
            self.dict_stats[peak] = dict_stats_data

        # To dump in a ConfigObj
        if filename is not None:
            dump_file = ConfigObj(indent_type='\t')
            dump_file.filename = filename

            for peak in self.peak_names:
                # compute for each variable of interest (center, height, fwhm)
                dump_file[peak] = {}
                for column in self.cols_dataframe:
                    dump_file[peak][column] = self.dict_stats[peak][column]
            dump_file.write()

        return self.dict_stats

    @staticmethod
    def _apply_statistics(data):
        """
        Apply any stadistics to the data.
        :param data: dataframe with intensities for example.
        :return: dictionary with the results.
        """
        average = np.average(data)
        std = np.std(data)

        dict_stats_data = {'average': average, 'std': std}
        return dict_stats_data

    def compute_intensity_ratio_each_sample(self, file_to_save=None):
        """
        computes the intensity ratio as D/G
        :param file_to_save: filename to save the dataframe (table_results.csv)
        :return: self.intensity ratios. A set with the ratio values.
        """

        # create the idx variable for better indexing in multiindex pandas
        idx = pd.IndexSlice  # this helps for the multiindexing

        self.intensity_ratios = {}
        for sample in self.sample_names:
            # grab the data
            D_band = self.data_pandas.loc[idx[sample, 'D'], 'height']
            G_band = self.data_pandas.loc[idx[sample, 'G'], 'height']

            self.intensity_ratios[sample] = D_band / G_band

        if file_to_save is not None:
            # TODO: move this somewhere else
            dump_file = ConfigObj(indent_type='\t')
            dump_file.filename = file_to_save

            for sample in self.sample_names:
                # compute for each variable of interest (center, height, fwhm)
                dump_file[sample] = self.intensity_ratios[sample]
            dump_file.write()

        return self.intensity_ratios

    def compute_equivalent_La(self, Lambda=532, file_to_save=None):
        """
        Computes the equivalent La from the formula of Cancado 2006
        La = (2.4*10**(-10))*Lambda**4*(I_D/I_G)**-1
        the laser wavelength is typically 532 nm
        The result is also in nm.
        *WARNING*: this method considers that this one compute_intensity_ratio_each_sample has been run before.
        :return: La: the set of computed equivalent La
        :param Lambda: wavelength of the laser
        :param file_to_save: filename of to save to
        """
        La = {}
        for sample in self.sample_names:
            La[sample] = (2.4 * 10 ** (-10)) * Lambda ** 4 * (self.intensity_ratios[sample]) ** -1

        if file_to_save is not None:
            dump_file = ConfigObj(indent_type='\t')
            dump_file.filename = file_to_save

            for sample in self.sample_names:
                # compute for each variable of interest (center, height, fwhm)
                dump_file[sample] = La[sample]
            dump_file.write()
        return La

    def add_xypositions(self):
        """
        adds the x, y position of the raman measurement to the dataframe in order to be output in table_results
        """
        x_positions = []
        y_positions = []
        self.data_pandas = self.data_pandas.unstack()
        self.unstacked = True
        for filename in self.file_names_experiment:
            positions = self._read_xyz(filename=filename)
            x_positions.append(positions[0])
            y_positions.append(positions[1])

        self.data_pandas['x(um)'] = x_positions
        self.data_pandas['y(um)'] = y_positions

    @staticmethod
    def _read_xyz(filename):
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


def get_num(string_with_number):
    """
    get numbers in a string
    :param string_with_number:
    :return:
    """
    try:
        return int(''.join(ele for ele in string_with_number if ele.isdigit()))  # get digits
    except ValueError:  # if no digits in the string, just assign -1
        return -1
