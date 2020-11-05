import os
from abc import ABC, abstractmethod
from pathlib import Path

from configobj import ConfigObj
from lmfit.models import LorentzianModel, QuadraticModel, LinearModel, ConstantModel, PolynomialModel
from matplotlib import pyplot as plt
from scipy.signal import savgol_filter

try:
    from plot_python_vki import apply_style

    apply_style()
except ImportError:
    pass


class GenericFit(ABC):
    """
    Generic Fit class
    This abstract class is used as base for Raman and XRD specialized classes.


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

    """


    def __init__(self, experimental_data=None, peaks=None, other_data=None, folder_out=None):
        """

        Parameters
        ----------
        experimental_data: df with experimental data
        peaks: list of peaks to be retrieved
        other_data: if needed
        folder_out: str folder where report will be saved.
        """

        if peaks is None:
            self.peaks = []
        else:
            self.peaks = peaks

        if other_data is not None:
            self.other_data = other_data
        else:  # here we can add some default values in a dictionary
            self.other_data = dict()
            self.other_data['_normalize_data'] = True
            self.other_data['bkg'] = 'quadratic'

        if folder_out is not None:
            self.folder_out = Path(folder_out)
        else:  # output folder for the fitting report, which is not used normally
            self.folder_out = Path('out_report')

        # create the out folder. if it exists, just pass
        os.makedirs(self.folder_out, exist_ok=True)

        self.experimental_data = experimental_data

        # extract the experimental data into two variables. Gets extended in inheritance
        self.var_x = None  # name of the varible x
        self.var_y = None
        self.x = None  # values of x
        self.y = None
        self.model = None
        self.params = None
        self.filename = None
        self.dict_tolerances_fit = None

    def apply_normalize(self):
        """
        performs the normalization.
        """
        self.y = self._normalize_data(self.y)

    def apply_smoothing(self):
        """
        performs smoothing using _sav_gol filter.
        :param inplace: bool. if False, returns a new column for the df called smoothed.
        """

        win_size = self._try_get_other_data(self.other_data, 'window_size', default_value=(15,))[0]
        poly_order = self._try_get_other_data(self.other_data, 'poly_order', default_value=(3,))[0]

        self.y = self._sav_gol(self.y, win_size=win_size, poly_order=poly_order)

    @abstractmethod
    def set_tolerances_fit(self):
        pass

    def build_fitting_model_peaks(self):
        """
        Builds the fitting model with parameters.
        It uses a quadraticModel to remove background noise, even though it is not the most important.
        :return:

        """

        poly_type = self.other_data.get('poly_type')
        model, params = self.create_bkg_model()
        for i, cen in enumerate(self.peaks):
            peak, pars = self._add_peak('lz%d' % (i + 1), cen, amplitude=self.dict_tolerances_fit['amplitude'],
                                        sigma=self.dict_tolerances_fit['sigma'],
                                        tolerance_center=self.dict_tolerances_fit['tolerance_center'],
                                        min_max_amplitude=self.dict_tolerances_fit['min_max_amplitude'],
                                        min_max_sigma=self.dict_tolerances_fit['min_max_sigma'])
            model = model + peak
            params.update(pars)

        self.model = model
        self.params = params

    def run_fit_model(self):
        """
        Perform the fit
        """
        result, components = self._fit_lorentzians(self.x, self.y, self.model, self.params)
        self.result = result
        self.components = components

    def save_results(self):
        """
        Saves 2 types of files:
            report file : with a lot of data
            params file : with the actual paramters and their std.

        """
        # save fit report to a file:
        with open(f'{self.folder_out / self.filename}_report', 'w') as fh:
            fh.write(self.result.fit_report())

        with open(f'{self.filename}_params.txt', 'w') as fh:
            for key in self.result.params:
                fh.write(key + " = " + str(self.result.params[key].value) + '\n')
                fh.write(key + '_stderr = ' + str(self.result.params[key].stderr) + '\n')

    def plot_results(self):
        """
        Plots the results of the fit.
        """
        plt.plot(self.x, self.y, label='data')
        plt.plot(self.x, self.result.best_fit, label='best fit')
        for name, component in self.components.items():
            if isinstance(component, float):
                plt.axhline(component, linestyle='--', label=name)
            else:
                plt.plot(self.x, component, linestyle='--', label=name)
        plt.xlabel(self.var_x)
        plt.ylabel(self.var_y)
        plt.legend(loc='upper right')
        plt.savefig(self.filename + '.png')
        plt.close()

    def create_bkg_model(self):
        """
        Creates a bkg model for removing the background from the signals.
        Gets the data from the other_data part of the input file. Otherwise it will assign quadratic.

        :return: model lmfit  for the bkg function.
        :return: params lmfit parameters to be adjusted.
        """
        bkg_model = self._choose_bkg_model(self.other_data['poly_type'])
        model = bkg_model[0](**bkg_model[1])
        params = model.make_params(bkg_model[2])

        return model, params

    #########
    # the static methods in the following are basically the ones doing the tasks.
    # The underscore is to treat them as private
    #########
    @staticmethod
    def _add_peak(prefix, center, amplitude, sigma, tolerance_center,
                  min_max_amplitude, min_max_sigma):
        """
        adds a peak using a LorentzianModel from lmfit. Peaks can be summed as a linear combination


        :param prefix: str
                name of the peak
        :param center: float
                center location
        :param amplitude: float
                amplitude of the peak
        :param sigma: float
                controls shape of the peak
        :param min_max_amplitude: tuple
                for the amplitude of the peak
        :param tolerance_center: float
                plus minus this quantity for the peak center location
        :param min_max_sigma: tuple
                for the sigma of the peak
        :return: peak lmfit model with the peak and its properties.
        :return: pars lmfit parameters to be adjusted.
        """
        peak = LorentzianModel(prefix=prefix)  # created a lorentzian function
        pars = peak.make_params()

        pars[prefix + 'center'].set(center, min=center - tolerance_center, max=center + tolerance_center)
        pars[prefix + 'amplitude'].set(amplitude, min=min_max_amplitude[0], max=min_max_amplitude[1])
        pars[prefix + 'sigma'].set(sigma, min=min_max_sigma[0], max=min_max_sigma[1])
        return peak, pars

    @staticmethod
    def _fit_lorentzians(x, y, model, params):
        """
        Fits the lorentzians to the experimental data.
        It uses a quadraticModel to remove background noise, even though it is not the most important.
        :param x: 1D array like
                with the x values, namely 2theta or raman displacement
        :param y: 1D array like
                with intensity counts
        :param model: lmfit model
                to be fit
        :param params: lmfit params
                to be adjusted
        :return:
        """

        init = model.eval(params, x=x)
        result = model.fit(y, params, x=x)
        components = result.eval_components()

        return result, components

    @staticmethod
    def _choose_bkg_model(poly_type):
        """
        Selects a bkg model for the fit. If not available, it will use the default quadratic.

        :param poly_type: str
                Type of bkg: linear, quadratic, constant, cubic.
        :return: lmfit model
                the bkg model to be added in the fitting.
        """
        poly_type = poly_type.lower()  # to avoid typos

        poly_type_dict = {
            'quadratic': (QuadraticModel, {'prefix': 'bkg'}, {'a': 0, 'b': 0, 'c': 0}),
            'linear': (LinearModel, {'prefix': 'bkg'}, {'intercept': 0, 'slope': 0}),
            'constant': (ConstantModel, {'prefix': 'bkg'}, {'c': 0}),
            'cubic': (PolynomialModel, {'prefix': 'bkg', 'degree': 3}, {'c': 0}),
        }

        try:
            bkg_model = poly_type_dict[poly_type]
        except KeyError:
            print('Background model not available, using quadratic')
            bkg_model = poly_type_dict['quadratic']

        return bkg_model

    @staticmethod
    def _try_get_other_data(other_data, string_to_find, default_value):
        """
        This method tries to get the default data for a given property. If it does not find it, the value returned
        will be the default one.

        :param other_data: dict
                dictionary with extra data passed
        :param string_to_find: str
                parameter to find
        :param default_value: tuple or float or else
                default value if the string is not found
        :return: list_numbers
                either a list of numbers, or float, or else, corresponding to the values specififed for the quantity.

        """
        try:
            data_requested = other_data[string_to_find]
        except KeyError:
            list_numbers = default_value
            print(f'{string_to_find} range not found, set to default: {default_value}')
            return list_numbers

        if isinstance(data_requested, str):
            data_requested = (data_requested,)

        list_numbers = tuple(map(float, data_requested))

        return list_numbers

    @staticmethod
    def _sav_gol(intensity_data, win_size=11, poly_order=4):
        """
        applies the savgol_filter for a 1D data. set as static method for convenience.
        
        :param intensity_data: 
            1D array with the original data
        :return: 1D array
            with data smoothed
        """
        data_smoothed = savgol_filter(intensity_data, window_length=int(win_size), polyorder=int(poly_order), axis=0)
        return data_smoothed

    @staticmethod
    def _normalize_data(intensity_data):
        """
        Here we normalize as z = z - min(x)/(max(x)-min(x)).
        :param intensity_data
            1D array with the original data
        :return: intensity_data_scaled:
            scaled intensity data
        """

        min_intensity = min(intensity_data)
        max_intensity = max(intensity_data)

        intensity_data_scaled = (intensity_data - min_intensity) / (max_intensity - min_intensity)
        return intensity_data_scaled

    @staticmethod
    def read_otherdata_configfile(config_file):
        """
        Read if there is any extra data in the configfile

        :param config_file: str
                name of the file with the extra data
        :return: dict
                with the other data.
        """

        config = ConfigObj(config_file)
        other_data = config.get('other data', None)

        return other_data

    @staticmethod
    def read_peaks_configfile(config_file, default_peaks_file, default_folder=None):
        """
        Alternate constructor from configobj file

        :param file_to_analyze: filename
        :param config_file: configobj file
        :return:
        """
        if default_folder is None:
            default_folder = Path(os.path.dirname(__file__))

        try:
            config = ConfigObj(config_file)
            # get the peaks, transform them to floats, and put them in a list, then sort the list
            peaks = list(map(float, config['peaks']))
            peaks.sort()
        except:
            print(
                'Data peaks not found or corrupted, using the one, which is in the ramanpy folder')
            config = ConfigObj(str(default_folder / default_peaks_file))
            # get the peaks, transform them to floats, and put them in a list, then sort the list
            peaks = list(map(float, config['peaks']))
            peaks.sort()

        return peaks
