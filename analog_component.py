'''
Analog component class definition.
'''

import numpy as np

class AnalogComponent(object):

    def __init__(self, name = "", comp_type = "", reference_file = "", skiprows = -1, gain_row = -1, units = "-1"):
        self.name = name
        self.comp_type = comp_type

        self.gain_data = [-1., -1]
        self.noise_data = [-1., -1]

        self.filled_data = False
        self.filled_noise_data = False

        self.OIP3_data = np.zeros((2, 0)) + -999.
        self.OIP3_filled = False
        self.comp_data = np.zeros((2, 0)) + -999.
        self.output_compression_filled = False


        if reference_file != "":
            self.fill_data_array(reference_file, skiprows, gain_row, units)
            if comp_type == "filter":
                self.fill_filter_noise_figure()
            elif comp_type == "attenuator":
                self.fill_attn_noise_figure()
            

    def fill_amplifier_array(self, reference_file, skiprows, comp_row, OIP3_row, units = "-1"):

        data = np.loadtxt(reference_file, skiprows=skiprows, unpack=1)

        if units == "Hz":
            self.comp_data = np.array([data[0], data[comp_row]])
            self.OIP3_data = np.array([data[0], data[OIP3_row]])
        elif units == "MHz":
            self.comp_data = np.array([data[0] * 1.e6, data[comp_row]])
            self.OIP3_data = np.array([data[0] * 1.e6, data[OIP3_row]])

        self.OIP3_filled = True
        self.output_compression_filled = True

    def fill_data_array(self, reference_file, skiprows, gain_row, units):
        '''
        Fill data array with information from a reference file.
        '''
                
        data = np.loadtxt(reference_file, skiprows=skiprows, unpack=1)
        if units == "Hz":
            self.gain_data = np.array([data[0], data[gain_row]])
        elif units == "MHz":
            self.gain_data = np.array([data[0] * 1.e6, data[gain_row]])

        self.filled_data = True

    def set_data_array(self, f, g):
        '''
        Fill data array manually.
        params: 
            f: frequency array
            g: gain array
        '''
        self.filled_data = True 
        self.gain_data = [np.array(f), np.array(g)]

    def fill_attn_noise_figure(self):
        '''
        Automattically fill noise figure data array.
        For an attenuator, NF = attenuation = -1 * gain
        '''
        self.noise_data[0] = np.array(self.gain_data[0])
        self.noise_data[1] = np.array(-self.gain_data[1])
        self.filled_noise_data = True

    def fill_filter_noise_figure(self):
        '''
        Automattically fill noise figure data array.
        For a purely capacitative filter, NF = 0.
        '''
        self.noise_data[0] = np.array(self.gain_data[0])
        self.noise_data[1] = np.array(np.zeros(len(self.gain_data[0])))
        self.filled_noise_data = True        

    def fill_noise_figure(self, reference_file, skiprows, noise_figure_row, units):
        '''
        Fill data array with information from a reference file.
        '''
                
        data = np.loadtxt(reference_file, skiprows=skiprows, unpack=1)
        if units == "Hz":
            self.noise_data = np.array([data[0], data[noise_figure_row]])
        elif units == "MHz":
            self.noise_data = np.array([data[0] * 1.e6, data[noise_figure_row]])

        self.filled_noise_data = True


    def set_noise_figure(self, f, noise_figure):
        '''
        Manually set noise figure data array
        '''
        self.noise_data = [np.array(f), np.array(noise_figure)]
        self.filled_noise_data = True

        
    def get_gain(self, f):
        '''
        Return gain of analog component at frequency f
        '''
        return np.interp(f, self.gain_data[0], self.gain_data[1])

    def get_OIP3(self, f):
        '''
        Return gain of analog component at frequency f
        '''
        if self.comp_type == "amplifier":
            return np.interp(f, self.OIP3_data[0], self.OIP3_data[1])        
        else:
            return 99.

    def get_compression_point(self, f):
        '''
        Return compression point of analog component at frequency f
        '''
        if self.comp_type == "amplifier":
            return np.interp(f, self.comp_data[0], self.comp_data[1])    
        else:
            return 999

    def get_NF(self, f):  
        '''
        Return the noise figure of a component at frequency f
        '''

        return np.interp(f, self.noise_data[0], self.noise_data[1])

    def get_noise_temperature(self, f):
        
        NF = np.interp(f, self.noise_data[0], self.noise_data[1])
        T = ((10.**(NF / 10.)) - 1. ) * (273.15 + 16.85)

        return T




class Cable(AnalogComponent):

    def __init__(self, **kwargs):
        
        self.k1 = kwargs.pop('k1')
        self.k2 = kwargs.pop('k2')        
        self.length = kwargs.pop('length') #meters
        super(Cable, self).__init__(**kwargs)
        
        test_freq = np.linspace(0.1e6, 3200.e6, 32000)
        self.set_data_array(test_freq, -self.cable_atten(test_freq) * self.length / 100.) 
        self.fill_attn_noise_figure()


    def cable_atten(self, freq):  
        '''
        Return attenuation in dB / 100 m
        '''
        unit_conv = 3.2808399 #ft per m

        atten = self.k1 * np.sqrt(freq / 1.e6) + self.k2 * freq / 1.e6 #dB / 100 ft
        return atten * unit_conv # (1/ft) * (ft/m) = 1 / m

class Amplifier(AnalogComponent):
    '''
    Amplifier subclass of AnalogComponent
    '''

    def __init__(self, **kwargs):

        self.k = 1


class Attenuator(AnalogComponent):
    '''
    Attenuator subclass of AnalogComponent
    '''

    def __init__(self, **kwargs):

        self.k = 1        