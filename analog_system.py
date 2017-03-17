import os
import numpy as np

import analog_component
from analog_component import AnalogComponent
from analog_component import Cable


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
    kB = 1.380648522000 * pow(10.,-23.) # J / K

    pre_gain = calc_gain_before_component(chain, component_index, freq, True)
    
    if component_index > 0:
        T_before_index_at_input = np.zeros(len(freq))
        for i in xrange(component_index):
            T_before_index_at_input += refer_component_temp_to_input(chain, i, freq)
    else:
        T_before_index_at_input = np.array([0.])
    
    noise_power_at_component = T_before_index_at_input * kB * pre_gain

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
            #temperatures = np.append(temperatures, component.get_noise_temperature(frequency))
            #gains = np.append(gains, component.get_gain(frequency))
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
    Build an attenuator with attenuation of specified value
    '''

    if attenuation - int(attenuation) == 0:
        name = "atten-%i" % (int(attenuation))
    else:
        name = "atten-%.1f" % (attenuation)
    atten = AnalogComponent(name, "attenuator")
    
    freq = np.logspace(-2, 3, 1000.)
    attenuation = np.zeros(len(freq)) - attenuation
           
    atten.set_data_array(freq, attenuation)
    atten.fill_attn_noise_figure()
    
    return atten    


def build_LNA():
    #lna = AC("LNA","amplifier","/home/sean/work/cosmology/suit/analog_files/data/lna/LNA.S2P", skiprows=13, gain_row=3, units="Hz")
    lna = AnalogComponent("LNA","amplifier",component_path + "/data/lna/LNA_noise.txt", skiprows=2, gain_row=1, units="MHz")
    lna.fill_noise_figure(component_path + "/data/lna/LNA_noise.txt", skiprows=2, noise_figure_row=-1, units="MHz")
    lna.fill_amplifier_array(component_path + "/data/lna/LNA_noise.txt", skiprows=2, comp_row=-2, OIP3_row=-3, units="MHz")
    return lna

def build_AMP():
    amp = AnalogComponent("AMP","amplifier",component_path + "/data/amp/amp.S2P", skiprows= 5, gain_row=3, units="MHz")
    amp.fill_noise_figure(component_path + "/data/amp/ampNOISE.txt", skiprows=1, noise_figure_row=6, units="MHz")
    amp.fill_amplifier_array(component_path + "/data/amp/ampNOISE.txt", skiprows=2, comp_row=-1, OIP3_row=-3, units="MHz")
    return amp

def build_LMR400(length):
    cable = AnalogComponent("LMR-400", "attenuator")
    #LMR-400; length
    freq = 1.e6 * np.array([30, 50, 150, 220, 450, 900, 1500, 1800, 2000, 2500, 5800])
    gain = -1*    np.array([2.2, 2.9, 5.0, 6.1, 8.9, 12.8, 16.8, 18.6, 19.6, 22.2, 35.5])
    cable.set_data_array(freq, (length / 100.)*gain)
    cable.fill_attn_noise_figure()
    return cable


def build_chime_filter():
    '''
    Build a CHIME filter object.
    '''
    chime_filter = AnalogComponent("CHIMEFILTER", "filter",
                                   component_path + "/data/chime_filter/chimefilter.S2P", 
                                   skiprows=11, 
                                   gain_row=3,
                                   units="MHz")
    chime_filter.fill_filter_noise_figure()
    return chime_filter

def build_low340_filter():
    low340_filter = AnalogComponent("low340_filter", "filter",
                                    component_path + "/data/low340_filter/low_filter.S2P", 
                                    skiprows=19, 
                                    gain_row=3,
                                    units="MHz")
    low340_filter.fill_filter_noise_figure()
    return low340_filter

def build_high875_filter():
    high875_filter = AnalogComponent("high875_filter", "filter",
                                     component_path + "/data/high875_filter/high875_filter.txt", 
                                     skiprows=1, 
                                     gain_row=2,
                                     units="MHz")
    high875_filter.gain_data[1] = -1*high875_filter.gain_data[1] #correct sign of insertion loss
    high875_filter.fill_filter_noise_figure()
    return high875_filter

def build_FM_filter():
    fm_filter = AnalogComponent("FM_filter", "filter",
                                component_path + "/data/fm_filter/FM_filter_25deg.S2P", 
                                skiprows=18, 
                                gain_row=3,
                                units="MHz")
    
    fm_filter.fill_filter_noise_figure()
    return fm_filter    

#RBP650_filter = AnalogComponent("RBP650_filter", "filter", 
#                                component_path + "/data/RBP-650/RBP-650+_Plus25degC.S2P",
#                                skiprows=11,
#                                gain_row=5,
#                                units="MHz")
#RBP650_filter.fill_filter_noise_figure()

#EDU1063_filter = AnalogComponent("EDU1063_filter", "filter", 
#                                 component_path + "/data/EDU1063/EDU1063.txt",
#                                 skiprows=1,
#                                 gain_row=1,
#                                 units="MHz")
#EDU1063_filter.fill_filter_noise_figure()

def build_chime_ADC():
    boltzmann_constant = 1.38064852*10**-23
    ADC_quantization_noise_temperature = ((((0.5 / pow(2., 8)) / np.sqrt(12))**2)/100) / 400.e6/ boltzmann_constant # 0.5V, 100 Ohms, 400 MHz
    ADC = AnalogComponent("ADC", "attenuator",
                          "/home/sean/work/cosmology/suit/analog_files/data/ADC/ADC.txt",
                          skiprows=1,
                          gain_row=1,
                          units="Hz")
    ADC.fill_noise_figure("/home/sean/work/cosmology/suit/analog_files/data/ADC/ADC.txt",
                          skiprows=1,
                          noise_figure_row=-1,
                          units="Hz")    

    return ADC

def build_cable(name, k1, k2, length):
    kwargs = {}
    kwargs['name'] = name
    kwargs['k1'] = k1
    kwargs['k2'] = k2
    kwargs['length'] = length
    cable = Cable(**kwargs)

    return cable



LNA = build_LNA()
AMP = build_AMP()
LMR400 = build_LMR400(30.)
CHIMEFILTER = build_chime_filter()
low340_filter = build_low340_filter()
high875_filter = build_high875_filter()
ADC = build_chime_ADC()


TUNABLE0 = AnalogComponent("tunable0", "filter", 
                "/home/sean/work/cosmology/suit/analog_files/data/tuned_655/FILTER0.s2p",
                 skiprows=15,
                 gain_row=3,
                 units="Hz")
TUNABLE0.fill_filter_noise_figure()

TUNABLE1 = AnalogComponent("tunable1", "filter", 
                "/home/sean/work/cosmology/suit/analog_files/data/tuned_655/FILTER1.s2p",
                 skiprows=15,
                 gain_row=3,
                 units="Hz")
TUNABLE1.fill_filter_noise_figure()

TUNABLE2 = AnalogComponent("tunable2", "filter", 
                "/home/sean/work/cosmology/suit/analog_files/data/tuned_655/FILTER2.s2p",
                 skiprows=15,
                 gain_row=3,
                 units="Hz")
TUNABLE2.fill_filter_noise_figure()

TUNABLE3 = AnalogComponent("tunable3", "filter", 
                "/home/sean/work/cosmology/suit/analog_files/data/tuned_655/FILTER3.s2p",
                 skiprows=15,
                 gain_row=3,
                 units="Hz")
TUNABLE3.fill_filter_noise_figure()


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
