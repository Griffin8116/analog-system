import os
import numpy as np

import analog_component
from analog_component import AnalogComponent
from analog_component import Cable
from analog_component import Amplifier
from analog_component import Attenuator
from analog_component import Filter


component_path = os.path.split(os.path.abspath(analog_component.__file__))[0]

def lin_gain(gain_db):
    '''
    Convert a gain in dB to a linear gain.
    '''
    return pow(10., np.array(gain_db) / 10.)

def calc_gain_before_component(chain, component_index, freq, linear=False):
    '''
    Calculate the total gain of components in a chain before a specified component.
    '''

    gain = 0
    for component in chain[0:component_index]:
        gain += component.get_gain(freq)
    
    if linear:
        gain = pow(10., gain / 10.)
        
    return gain

def calc_total_gain(freq, chain, linear=False):
    '''
    Calculate total gain of a component chain at frequency freq.
    '''
  
    gain = 0
    for comp in chain:
        gain = gain + comp.get_gain(freq)
        
    if linear:
        return lin_gain(gain)
    else:
        return gain

def refer_component_temp_to_input(chain, component_index, freq):
    '''
    Refer the noise temperature of a single component in a chain to the input of the chain.
    '''

    component_temperature = chain[component_index].get_noise_temperature(freq)
    
    gain = calc_gain_before_component(chain, component_index, freq, True)
    
    return component_temperature / gain        


def calc_noise_power_before_component(chain, component_index, freq):
    '''
    Calculate noise power at the input to a component. 
    '''
    BOLTZMANN_CONSTANT = 1.380648522000 * pow(10.,-23.) # J / K

    pre_gain = calc_gain_before_component(chain, component_index, freq, True)
    
    if component_index > 0:
        T_before_index_at_input = np.zeros(len(freq))
        for i in xrange(component_index):
            T_before_index_at_input += refer_component_temp_to_input(chain, i, freq)
    else:
        T_before_index_at_input = np.array([0.])
    
    noise_power_at_component = T_before_index_at_input * BOLTZMANN_CONSTANT * pre_gain

    return noise_power_at_component



def calc_Tchain(temperatures=[3], gains=[3, 5], verbose=False):
    '''
    Calculate the reciever noise temperature T_R 
    (the referenced noise temperature due to amplifiers and such)
    at the input of the first component in the chain.
    '''
    
    if len(temperatures) != len(gains):
        print "ERROR! Temperature and gain arrays have mismatched lengths. Cannot compute."
        return -999
    
    temperature_total = 0  
    
    gains = lin_gain(gains)
    
    for (index, temperature) in enumerate(temperatures):
        #print "Temperature index: %i" % index
        if index == 0:
            temperature_total += temperature
        else:
            referenced_temperature = temperature
            total_gain = 1 
            for gain in gains[0:index]:
                total_gain = total_gain * gain
            
            if verbose: 
                print "Total gain for stage: %.2f (ratio)" % total_gain
                print "Total gain for stage: %.2f (dB)" % (10 * np.log10(total_gain))
            
            referenced_temperature = referenced_temperature / total_gain
            
            temperature_total += referenced_temperature
            
    return temperature_total


def calc_chain_temperature(chain, freq):
    (gains, temperatures) = np.array(zip(*[[i.get_gain(freq), i.get_noise_temperature(freq)] for i in chain]))
   
    return np.array([calc_Tchain(temperatures[:,i], gains[:,i], False) for i in range(len(freq))])


def calc_T_at_component(temperatures=[3], gains=[3, 5], component = 0, verbose = False):
    '''
    Calculate the noise temperature of the system referenced to a point at index 'component'. 
    '''
    
    linear_gains = lin_gain(np.array(gains))
    
    T0 = calc_Tchain(temperatures, gains, verbose)
    
    if component > len(temperatures):
        print "ERROR: Component out of bounds."
        return -998
    
    if component == 0:
        return T0
    
    Tdex = T0
    for i in np.arange(0, component):
        Tdex = Tdex * linear_gains[i]
                 
    return Tdex


def calc_T_at_component_at_freq(chain, component_index, freq):
    '''
    Calculate the noise temperature of the system referenced to a point at index 'component'. 
    '''
    
    component_temp = np.empty(len(freq))
    
    for i, frequency in enumerate(freq):
        #temperatures = np.array([])
        #gains = np.array([])
        temperatures = np.empty(len(chain))
        gains = np.empty(len(chain))

        for j, component in enumerate(chain):
            temperatures[j] = component.get_noise_temperature(frequency)
            gains[j] = component.get_gain(frequency)

        component_temp[i] = calc_T_at_component(temperatures, gains, component_index)

    return component_temp

        
def print_chain(chain):
    '''
    Print the names of the components in an AnalogComponent chain.
    '''

    print "Chain: "
    components = str(chain[0].name)
    
    for compoennt in chain[1:]:
        components += str(" --> " + compoennt.name)
    
    print components

def attenuator(attenuation):
    '''
    Build an attenuator with attenuation of specified value.
    '''

    if attenuation - int(attenuation) == 0:
        name = "atten-%i" % (int(attenuation))
    else:
        name = "atten-%.1f" % (attenuation)
    return Attenuator(name=name, attenuation=attenuation)



def build_LNA():

    spec_dict = {}
    spec_dict['name'] = "LNA"
    spec_dict['gain_reference_file'] = component_path + "/data/lna/ZX60-P103LN+_5V_Plus25DegC_Unit1.S2P"
    spec_dict['gain_file_units'] = "Hz"
    spec_dict['gain_file_skiprows'] = 13
    spec_dict['gain_row'] = 3

    spec_dict['amplifier_reference_file'] = component_path + "/data/lna/LNA_noise.txt"
    spec_dict['amplifier_file_skiprows'] = 2    
    spec_dict['amplifier_units'] = "MHz"
    spec_dict['NF_row'] = -1
    spec_dict['comp_row'] = -2
    spec_dict['OIP3_row'] = -3

    lna = Amplifier(**spec_dict)

    return lna


def build_AMP():
    
    spec_dict = {}
    spec_dict['name'] = "AMP"
    spec_dict['gain_reference_file'] = component_path + "/data/amp/amp.S2P"

    spec_dict['gain_file_units'] = "MHz"
    spec_dict['gain_file_skiprows'] = 5
    spec_dict['gain_row'] = 3

    spec_dict['amplifier_reference_file'] = component_path + "/data/amp/ampNOISE.txt"
    spec_dict['amplifier_file_skiprows'] = 2    
    spec_dict['amplifier_units'] = "MHz"
    spec_dict['NF_row'] = 6
    spec_dict['comp_row'] = -1
    spec_dict['OIP3_row'] = -3


    lna = Amplifier(**spec_dict)


    return lna    
'''
def build_LMR400(length):
    cable = AnalogComponent("LMR-400", "attenuator")
    #LMR-400; length
    freq = 1.e6 * np.array([30, 50, 150, 220, 450, 900, 1500, 1800, 2000, 2500, 5800])
    gain = -1*    np.array([2.2, 2.9, 5.0, 6.1, 8.9, 12.8, 16.8, 18.6, 19.6, 22.2, 35.5])
    cable.set_data_array(freq, (length / 100.)*gain)
    cable.fill_attn_noise_figure()
    return cable
'''

def build_chime_filter():
    '''
    Build a CHIME filter object.
    '''
    
    spec_dict = {}
    spec_dict['name'] = "CHIME_filter"
    spec_dict['gain_reference_file'] = component_path + "/data/chime_filter/chimefilter.S2P"

    spec_dict['gain_file_units'] = "MHz"
    spec_dict['gain_file_skiprows'] = 11
    spec_dict['gain_row'] = 3

    chime_filter = Filter(**spec_dict)
    return chime_filter

def build_low340_filter():

    spec_dict = {}
    spec_dict['name'] = "low340_filter"
    spec_dict['gain_reference_file'] = component_path + "/data/low340_filter/low_filter.S2P"

    spec_dict['gain_file_units'] = "MHz"
    spec_dict['gain_file_skiprows'] = 19
    spec_dict['gain_row'] = 3

    return Filter(**spec_dict)


def build_high875_filter():

    spec_dict = {}
    spec_dict['name'] = "high875_filter"
    spec_dict['gain_reference_file'] = component_path + "/data/high875_filter/high875_filter.txt"

    spec_dict['gain_file_units'] = "MHz"
    spec_dict['gain_file_skiprows'] = 1
    spec_dict['gain_row'] = 2

    return Filter(**spec_dict)
    

def build_FM_filter():
    fm_filter = 
    
    return Filter(name="FM_filter",
                       gain_reference_file=component_path + "/data/fm_filter/FM_filter_25deg.S2P", 
                       gain_file_skiprows=18, 
                       gain_row=3,
                       gain_file_units="MHz")    


def build_chime_ADC():
    #boltzmann_constant = 1.38064852*10**-23
    #ADC_quantization_noise_temperature = 
    # 0.5V, 100 Ohms, 400 MHz 
    #((((0.5 / pow(2., 8)) / np.sqrt(12))**2)/100) / 400.e6/ boltzmann_constant 

    spec_dict = {}
    spec_dict['name'] = "ADC"
    spec_dict['comp_type'] = "ADC"
    spec_dict['gain_reference_file'] =\
     "/home/sean/work/cosmology/suit/analog_files/data/ADC/ADC.txt"

    spec_dict['gain_file_units'] = "Hz"
    spec_dict['gain_file_skiprows'] = 1
    spec_dict['gain_row'] = 1

    spec_dict['noise_reference_file'] = spec_dict['gain_reference_file']
    spec_dict['noise_file_units'] = "Hz"
    spec_dict['noise_file_skiprows'] = 1
    spec_dict['noise_row'] = -1
    
    ADC = Attenuator(**spec_dict)
    ADC.fill_gain_array(**spec_dict)
    ADC.fill_noise_figure(**spec_dict)

    return ADC


def build_cable(name, k1, k2, length):

    '''
    Build a cable of a specified length 
    '''

    #kwargs = {}
    #kwargs['name'] = name
    #kwargs['k1'] = k1
    #kwargs['k2'] = k2
    #kwargs['length'] = length
    #cable = Cable(**kwargs)

    return Cable(name=name, k1=k1, k2=k2, length=length)



LNA = build_LNA()
AMP = build_AMP()
LMR400 = build_cable("LMR-400", 0.122290, 0.000260, 30.)
CHIMEFILTER = build_chime_filter()
low340_filter = build_low340_filter()
high875_filter = build_high875_filter()
ADC = build_chime_ADC()


TUNABLE0 = Filter(name="tunable0", 
                  gain_reference_file=\
                  "/home/sean/work/cosmology/suit/analog_files/data/tuned_655/FILTER0.s2p",
                  gain_file_skiprows=15,
                  gain_row=3,
                  gain_file_units="Hz"
                 )

TUNABLE1 = Filter(name="tunable1", 
                  gain_reference_file=\
                  "/home/sean/work/cosmology/suit/analog_files/data/tuned_655/FILTER1.s2p",
                  gain_file_skiprows=15,
                  gain_row=3,
                  gain_file_units="Hz"
                 )

TUNABLE2 = Filter(name="tunable2", 
                  gain_reference_file=\
                  "/home/sean/work/cosmology/suit/analog_files/data/tuned_655/FILTER2.s2p",
                  gain_file_skiprows=15,
                  gain_row=3,
                  gain_file_units="Hz"
                 )

TUNABLE3 = Filter(name="tunable3", 
                  gain_reference_file=\
                  "/home/sean/work/cosmology/suit/analog_files/data/tuned_655/FILTER3.s2p",
                  gain_file_skiprows=15,
                  gain_row=3,
                  gain_file_units="Hz"
                 )


TUNABLE0_ENC = Filter(name="tunable0_enc", 
                      gain_reference_file=\
                      "/home/sean/work/cosmology/suit/analog_files/data/tuned_655/TE0.s2p",
                      gain_file_skiprows=15,
                      gain_row=3,
                      gain_file_units="Hz"
                     )

TUNABLE1_ENC = Filter(name="tunable1_enc", 
                      gain_reference_file=\
                      "/home/sean/work/cosmology/suit/analog_files/data/tuned_655/TE1.s2p",
                      gain_file_skiprows=15,
                      gain_row=3,
                      gain_file_units="Hz"
                     )

TUNABLE2_ENC = Filter(name="tunable2_enc", 
                      gain_reference_file=\
                      "/home/sean/work/cosmology/suit/analog_files/data/tuned_655/TE2.s2p",
                      gain_file_skiprows=15,
                      gain_row=3,
                      gain_file_units="Hz"
                     )

TUNABLE3_ENC = Filter(name="tunable3_enc", 
                      gain_reference_file=\
                      "/home/sean/work/cosmology/suit/analog_files/data/tuned_655/TE3.s2p",
                      gain_file_skiprows=15,
                      gain_row=3,
                      gain_file_units="Hz"
                     )


FM_FILTER = Filter(name="FM_filter",
                   gain_reference_file=component_path + "/data/fm_filter/FM_filter_25deg.S2P", 
                   gain_file_skiprows=18, 
                   gain_row=3,
                   gain_file_units="MHz")
   

#TEST_FREQ = np.linspace(0.1e6, 3200.e6, 3200)
TEST_AMP = Amplifier(name="test_amp",
                     gain_reference_file=component_path + "/data/test_amp/ampNOISE.txt",
                     gain_file_skiprows=1,
                     gain_row=1,
                     gain_file_units="Hz",
                     amplifier_reference_file=component_path + "/data/test_amp/ampNOISE.txt",
                     amplifier_file_skiprows=1,
                     amplifier_units="Hz",
                     NF_row=-1,
                     comp_row=-1,
                     OIP3_row=-1
                    )

BP_filters = []

for i in [1, 2, 3, 5, 6, 8, 9]:

  BP = Filter(name="mnt_tuned_{0:d}".format(i), 
              gain_reference_file=\
              "/home/sean/work/cosmology/suit/analog_files/data/655_bandpass/F{0:d}.s2p".format(i),
              gain_file_skiprows=15,
              gain_row=3,
              gain_file_units="Hz"
             )

  BP_filters.append(BP)

'''
BP_1 = Filter(name="mnt_tuned_1", 
              gain_reference_file=\
              "/home/sean/work/cosmology/suit/analog_files/data/655_bandpass/P1.s2p",
              gain_file_skiprows=15,
              gain_row=3,
              gain_file_units="Hz"
             )

BP_2 = Filter(name="mnt_tuned_2", 
              gain_reference_file=\
              "/home/sean/work/cosmology/suit/analog_files/data/655_bandpass/P2.s2p",
              gain_file_skiprows=15,
              gain_row=3,
              gain_file_units="Hz"
             )

BP_3 = Filter(name="mnt_tuned_3", 
              gain_reference_file=\
              "/home/sean/work/cosmology/suit/analog_files/data/655_bandpass/P3.s2p",
              gain_file_skiprows=15,
              gain_row=3,
              gain_file_units="Hz"
             )

BP_4 = Filter(name="mnt_tuned_4", 
              gain_reference_file=\
              "/home/sean/work/cosmology/suit/analog_files/data/655_bandpass/P4.s2p",
              gain_file_skiprows=15,
              gain_row=3,
              gain_file_units="Hz"
             )

'''
'''

    spec_dict = {}
    spec_dict['name'] = "AMP"
    spec_dict['gain_reference_file'] = component_path + "/data/amp/amp.S2P"

    spec_dict['gain_file_units'] = "MHz"
    spec_dict['gain_file_skiprows'] = 5
    spec_dict['gain_row'] = 3

    spec_dict['amplifier_reference_file'] = component_path + "/data/amp/ampNOISE.txt"
    spec_dict['amplifier_file_skiprows'] = 2    
    spec_dict['amplifier_units'] = "MHz"
    spec_dict['NF_row'] = 6
    spec_dict['comp_row'] = -1
    spec_dict['OIP3_row'] = -3


    lna = Amplifier(**spec_dict)


TEST_FREQ = np.linspace(0.1e6, 3200.e6, 32000)
TEST_AMP = AnalogComponent("test_amp (10dB)", "amplifier")
TEST_AMP.set_data_array(TEST_FREQ, np.zeros(len(TEST_FREQ)) + 10.)
TEST_AMP.set_noise_figure(TEST_FREQ, np.zeros(len(TEST_FREQ)) + 10 * np.log10(1 + 290. / 290.))

TEST_FILTER = AnalogComponent("test_filter (-3dB)", "filter")
TEST_FILTER.set_data_array(TEST_FREQ, np.zeros(len(TEST_FREQ)) - 3.)
TEST_FILTER.fill_filter_noise_figure()

TEST_ATTENUATOR = AnalogComponent("test_atten (-1dB)", "attenuator")
TEST_ATTENUATOR.set_data_array(TEST_FREQ, np.zeros(len(TEST_FREQ)) -1.)
TEST_ATTENUATOR.fill_attn_noise_figure()

TEST_NOTCH = AnalogComponent("test_notch (600MHz)", "filter")
NOTCH_GAIN = np.zeros(len(TEST_FREQ)) - 10.
NOTCH_GAIN[(TEST_FREQ >= 590.e6) & (TEST_FREQ <= 610.e6)] = -3.
TEST_NOTCH.set_data_array(TEST_FREQ, NOTCH_GAIN)
TEST_NOTCH.fill_filter_noise_figure()
'''