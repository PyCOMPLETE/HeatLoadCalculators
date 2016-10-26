from scipy.constants import e as qe
from scipy.constants import m_p, epsilon_0, c
import numpy as np

lhc_bending_radius = 2803.95 # design report 
lhc_circumference = 26658.883 # design report
lhc_arc_length = 171.4*2 + 23*2*53.45 #Q7toQ7


class HeatLoadCalculatorSynchrotronRadiationLHCArc(object):
    def __init__(self):
        pass
        
        
    def calculate_P_Wm(self, ppb, sigma_t, energy_eV=0., n_bunches=None):
        if n_bunches is None:
            if not hasattr(ppb, '__iter__'):
                raise ValueError('N_bunches has to be specified somehow!')
            total_intensity = np.sum(ppb, axis=1)
                
        else:
            if hasattr(ppb, '__iter__'):
                raise ValueError('N_bunches should not be specified if a ppb array is provided!')
            total_intensity = n_bunches*ppb
            
        beam_gamma = energy_eV * qe / (m_p*c**2)
        energy_loss_turn_beam = (qe**2 * beam_gamma**4) / (3. * epsilon_0*lhc_bending_radius) * total_intensity
        
        power_loss_beam = energy_loss_turn_beam*c/lhc_circumference
        
        # we assume it is spread over the arcs (0 in the straight sections)
        power_arc_m = power_loss_beam/(8*lhc_arc_length)
        
        return power_arc_m
            
        
