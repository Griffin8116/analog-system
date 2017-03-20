'''
Analog component class definition.
'''

import numpy as np

class AnalogComponent(object):
    '''
    Analogue component base class.
    '''

    def __init__(self, **kwargs):

        self.name = kwargs.pop('name')
        self.comp_type = kwargs.pop('comp_type')

        self.gain_data = [-1., -1]
        self.noise_data = [-1., -1]

        self.filled_data = False
        self.filled_noise_data = False
        
        if kwargs.get('gain_reference_file', "") != "":
            self.fill_gain_array(**kwargs)
            
    def fill_gain_array(self, **kwargs):#reference_file, skiprows, gain_row, units):
        '''
        Fill data array with information from a reference file.
        Assumes that in the file the first column is frequency 
        and the gain row is passed as a parameter.
        '''

        reference_file = kwargs['gain_reference_file']
        skiprows = kwargs['gain_file_skiprows']
        gain_row = kwargs['gain_row']
        units = kwargs['gain_file_units']
                
        data = np.loadtxt(reference_file, skiprows=skiprows, unpack=1)

        if units == "MHz":
            data[0] = data[0] * 1.e6

        self.gain_data = np.array([data[0], data[gain_row]])

        self.filled_data = True

    def set_noise_array(self, freq, noise_figure):
        '''
        Fill data array manually.
        params: 
            freq: frequency array
            noise_figure: noise figure array
        '''
        self.filled_noise_data = True 
        self.noise_data = np.array([freq, noise_figure])

    def set_data_array(self, freq, gain):
        '''
        Fill data array manually.
        params: 
            freq: frequency array
            gain: gain array
        '''
        self.filled_data = True 
        self.gain_data = np.array([freq, gain])


    def fill_noise_figure(self, **kwargs):#reference_file, skiprows, noise_figure_row, units):
        '''
        Fill data array with information from a reference file.
        '''

        reference_file = kwargs['noise_reference_file']
        skiprows = kwargs['noise_file_skiprows']
        noise_figure_row = kwargs['noise_row']
        units = kwargs['noise_file_units']
                
        data = np.loadtxt(reference_file, skiprows=skiprows, unpack=1)

        if units == "MHz":
            data[0] = data[0] * 1.e6

        self.noise_data = np.array([data[0], data[noise_figure_row]])
        self.filled_noise_data = True


    def set_noise_figure(self, freq, noise_figure):
        '''
        Manually set noise figure data array
        '''
        self.noise_data = [np.array(freq), np.array(noise_figure)]
        self.filled_noise_data = True

        
    def get_gain(self, freq):
        '''
        Return gain of analog component at frequency f
        '''
        return np.interp(freq, self.gain_data[0], self.gain_data[1])

    def get_NF(self, freq):  
        '''
        Return the noise figure of a component at a specified frequency.
        '''

        return np.interp(freq, self.noise_data[0], self.noise_data[1])

    def get_noise_temperature(self, freq):

        '''
        Calculate the noise temperature of a component at a given frequency.
        '''
        
        noise_figure = self.get_NF(freq) 
        temperature = ((10.**(noise_figure / 10.)) - 1.) * (273.15 + 16.85)

        return temperature




class Amplifier(AnalogComponent):
    '''
    Amplifier subclass of AnalogComponent
    '''

    def __init__(self, **kwargs):

        kwargs['comp_type'] = "amplifier"
        super(Amplifier, self).__init__(**kwargs)      

        self.OIP3_data = np.zeros((2, 0)) + -999.
        self.OIP3_filled = False

        self.comp_data = np.zeros((2, 0)) + -999.
        self.output_compression_filled = False  

        self.fill_amplifier_array(**kwargs)

    def fill_amplifier_array(self, **kwargs):
        '''
        Fill the noise figure and OIP3 data for an amplifier.
        '''
        reference_file = kwargs['amplifier_reference_file']
        skiprows = kwargs['amplifier_file_skiprows'] 
        units = kwargs['amplifier_units']
        
        NF_row = kwargs['NF_row']
        comp_row = kwargs['comp_row'] 
        OIP3_row = kwargs['OIP3_row']

        data = np.loadtxt(reference_file, skiprows=skiprows, unpack=1)

        if units == "MHz":
            data[0] = data[0] * 1.e6
            
        self.noise_data = np.array([data[0], data[NF_row]])
        self.comp_data = np.array([data[0], data[comp_row]])
        self.OIP3_data = np.array([data[0], data[OIP3_row]])

        self.OIP3_filled = True
        self.output_compression_filled = True
        self.filled_noise_data = True


    def get_OIP3(self, freq):
        '''
        Return gain of analog component at frequency freq
        '''

        if self.comp_type == "amplifier":
            return np.interp(freq, self.OIP3_data[0], self.OIP3_data[1])        
        else:
            return None

    def get_compression_point(self, freq):
        '''
        Return compression point of analog component at frequency f
        '''

        if self.comp_type == "amplifier":
            return np.interp(freq, self.comp_data[0], self.comp_data[1])    
        else:
            return None


class Attenuator(AnalogComponent):
    '''
    Attenuator subclass of AnalogComponent
    '''

    def __init__(self, **kwargs):

        kwargs['comp_type'] = kwargs.get('comp_type', "attenuator")
        super(Attenuator, self).__init__(**kwargs)   
        test_freq = np.linspace(0.1e6, 3200.e6, 32000)
        self.set_data_array(test_freq, np.zeros(len(test_freq)) - kwargs.get('attenuation', 0.))
        self.fill_attn_noise_figure()

    def fill_attn_noise_figure(self):
        '''
        Automattically fill noise figure data array.
        For an attenuator, NF = attenuation = -1 * gain
        '''
        self.noise_data[0] = np.array(self.gain_data[0])
        self.noise_data[1] = np.array(-self.gain_data[1])
        self.filled_noise_data = True

class Cable(Attenuator):
    '''
    Cable subclass of Attenuator
    '''

    def __init__(self, **kwargs):
        
        self.cable_k1 = kwargs.pop('k1')
        self.cable_k2 = kwargs.pop('k2')        
        self.length = kwargs.pop('length') #meters
        kwargs['comp_type'] = "cable"
        kwargs['gain_reference_file'] = ""

        super(Cable, self).__init__(**kwargs)
        
        test_freq = np.linspace(0.1e6, 3200.e6, 32000)
        self.set_data_array(test_freq, -self.cable_atten(test_freq) * self.length / 100.) 
        self.fill_attn_noise_figure()

    def cable_atten(self, freq):  
        '''
        Return attenuation in dB / 100 m
        '''
        unit_conv = 3.2808399 #ft per m

        atten = self.cable_k1 * np.sqrt(freq / 1.e6) + self.cable_k2 * freq / 1.e6 #dB / 100 ft
        return atten * unit_conv # (1/ft) * (ft/m) = 1 / m



class Filter(AnalogComponent):
    '''
    Filter subclass of AnalogComponent. 
    Is similar to an attenduator but has no noise figure.
    '''
    
    def __init__(self, **kwargs):
        kwargs['comp_type'] = "filter"        
        super(Filter, self).__init__(**kwargs)        
        self.fill_filter_noise_figure()

    def fill_filter_noise_figure(self):
        '''
        Automattically fill noise figure data array.
        For a purely capacitative filter, NF = 0.
        '''

        self.noise_data[0] = np.array(self.gain_data[0])
        self.noise_data[1] = np.array(np.zeros(len(self.gain_data[0])))
        self.filled_noise_data = True      
