import os
import numpy as np
from scipy.constants import e as qe, c, m_p, epsilon_0
from scipy.interpolate import interp1d
from scipy.special import gamma as gammafunc




lhc_measurement_tools_path = os.path.dirname(os.path.realpath(__file__))
default_rho_vs_T_file = lhc_measurement_tools_path + '/copper_rho_vs_T.txt'

lhc_circumference = 26658.883
lhc_arc_chamb_radius = 18.4e-3 
lhc_arc_weld_thickness = 2e-3
rho_stainless_steel = 6e-7 # Elias

Z0_vac = 376.73031

data_Cu = np.loadtxt(default_rho_vs_T_file)

class HeatLoadCalculatorImpedance(object):
    def __init__(self, temperature_K, circumference_m, chamb_radius_m,
             weld_Thickness_m=None, weld_Rho_Ohm_m=None, rho_vs_T = data_Cu):
        # load resistivity
        if type(rho_vs_T) is str:
            data = np.loadtxt(rho_vs_T)
        else:
            data = rho_vs_T  
        self.copper_rho_Ohm_m = interp1d(data[:,0], data[:,1]*1e-8)
        
        self.temperature_K = temperature_K
        self.circumference_m = circumference_m
        self.chamb_radius_m = chamb_radius_m
        self.weld_Thickness_m = weld_Thickness_m
        self.weld_Rho_Ohm_m = weld_Rho_Ohm_m
        
        
    def calculate_P_Wm(self, ppb, sigma_t, b_field=0., n_bunches=None):
    
            
        """
        If n_bunches is not specified, ppb has to be either a numpy array with dimensions 
        time_steps x n_bunches or n_bunches.
        If n_bunches is specified, ppb can be either a number or an arry with dimensions of timesteps
        sigma_t is either a float or a np array of dimensions time_steps or timp_steps * n_bunches
        Do not specify n_bunches explicitly and implicitly at the same time!
        """
        
        if n_bunches is None:
            if type(ppb) is float or type(ppb) is int:
                raise ValueError('N_bunches has to be specified somehow!')
            else:
                n_bunches = 1.

        rho_B0_T = self.copper_rho_Ohm_m(self.temperature_K) # 0.014 *1e-8 for 20 K
        
        if not hasattr(b_field, '__iter__') and b_field == 0.:
            rho_B_T = rho_B0_T
        else:
            rho_B0_273 = self.copper_rho_Ohm_m(273.)
            #elias
            rho_B_T = rho_B_T = rho_B0_T*(1 + 10**(1.055*np.log10(b_field*(rho_B0_273/rho_B0_T))- 2.69))
           
            #philipp
            #try:
            #   rho_B0_4 = self.copper_rho_Ohm_m(4.)
            #except ValueError:
            #    rho_B0_4 = 3.*self.copper_rho_Ohm_m(4.2) - 2.*self.copper_rho_Ohm_m(4.3)            
            #rho_B_T = rho_B0_T * (1. + 10.**(-2.69) * (b_field*rho_B0_273/rho_B0_4)**1.055)
            
            #benoit/nicolas
            #rho_B_T = rho_B0_T *(1.0048 + 0.0038*b_field*self.copper_rho_Ohm_m(300.)/self.copper_rho_Ohm_m(20.))
            #For 20K!!!
            #rho_B_T = 2.4e-10 *(1.0048 + 0.0038*b_field*70.)

        #print rho_B_T
        per_bunch_factor = ppb**2 * sigma_t**(-1.5)

        if hasattr(ppb, '__iter__') and len(ppb.shape) == 2:
            per_bunch_factor = np.sum(per_bunch_factor, axis=1)

        P_no_weld = 1./self.circumference_m*gammafunc(0.75)*n_bunches/self.chamb_radius_m * (qe/2./np.pi)**2* np.sqrt(c*rho_B_T*Z0_vac/2.) * per_bunch_factor

        if self.weld_Thickness_m is not None:
            weld_Factor = 1. + np.sqrt(self.weld_Rho_Ohm_m/rho_B_T)*self.weld_Thickness_m/(2.*np.pi*self.chamb_radius_m)
            #weld_Factor = 1.4 # benoit assumption
        else:
            weld_Factor = 1.

        return P_no_weld*weld_Factor
        

                
def B_lhc_arc_dipole(energy_eV):
    return 8.33*energy_eV/7e12
def grad_lhc_arc_quad(energy_eV):
    return 200.*energy_eV/7e12
                

class HeatLoadCalculatorImpedanceLHCArc(object):
    
    def __init__(self, lhc_arc_bs_temperature = 20.):
        self.hl_beam_screen = HeatLoadCalculatorImpedance(temperature_K = lhc_arc_bs_temperature, 
            circumference_m=lhc_circumference, chamb_radius_m=lhc_arc_chamb_radius,
            weld_Thickness_m=lhc_arc_weld_thickness, weld_Rho_Ohm_m=rho_stainless_steel)
    
    def calculate_P_Wm_dipole(self, ppb, sigma_t, energy_eV=0., n_bunches=None):
        return self.hl_beam_screen.calculate_P_Wm(ppb, sigma_t, b_field=B_lhc_arc_dipole(energy_eV), n_bunches=n_bunches)
        
    def calculate_P_Wm_quad(self, ppb, sigma_t, energy_eV=0., n_bunches=None):
        return self.hl_beam_screen.calculate_P_Wm(ppb, sigma_t, 
            b_field=grad_lhc_arc_quad(energy_eV)*self.hl_beam_screen.chamb_radius_m, n_bunches=n_bunches)
    
    def calculate_P_Wm_drift(self, ppb, sigma_t, energy_eV=0., n_bunches=None):
        return self.hl_beam_screen.calculate_P_Wm(ppb, sigma_t, 
            b_field=0., n_bunches=n_bunches)
            
    def calculate_P_Wm(self, ppb, sigma_t, energy_eV=0., n_bunches=None):
        #return self.calculate_P_Wm_dipole(ppb, sigma_t, energy_eV, n_bunches)#benoit assumption
         return (self.calculate_P_Wm_dipole(ppb, sigma_t, energy_eV, n_bunches)*14.3*3+\
                 self.calculate_P_Wm_quad(ppb, sigma_t, energy_eV, n_bunches)*3.1+\
                 self.calculate_P_Wm_drift(ppb, sigma_t, energy_eV, n_bunches)*7.)/53.
                        
                


  
