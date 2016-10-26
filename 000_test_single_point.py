import sys, os
BIN = os.path.expanduser("../")
sys.path.append(BIN)



import impedance_heatload as ihl




#test
ppb_test = 1.15e11
sigma_t_test = 1e-9/4.
energy_eV_test=7000e9
n_bunches_test = 2808

hl_calculator  = ihl.HeatLoadCalculatorImpedanceLHCArc()
hl_imped_dip = hl_calculator.calculate_P_Wm_dipole(ppb=ppb_test, sigma_t=sigma_t_test,  energy_eV=energy_eV_test, n_bunches=n_bunches_test)
print 'Impedance heating in the dipoles: %.1f mW\m'%(hl_imped_dip*1000)
hl_imped_quad = hl_calculator.calculate_P_Wm_quad(ppb=ppb_test, sigma_t=sigma_t_test,  energy_eV=energy_eV_test, n_bunches=n_bunches_test)  
print 'Impedance heating in the quadrupoles: %.1f mW\m'%(hl_imped_quad*1000)
hl_imped_drift = hl_calculator.calculate_P_Wm_drift(ppb=ppb_test, sigma_t=sigma_t_test,  energy_eV=energy_eV_test, n_bunches=n_bunches_test)
print 'Impedance heating in the drifts: %.1f mW\m'%(hl_imped_drift*1000)
hl_imped_hcell = hl_calculator.calculate_P_Wm(ppb=ppb_test, sigma_t=sigma_t_test,  energy_eV=energy_eV_test, n_bunches=n_bunches_test)  
print 'Impedance heating average half-cell: %.1f mW\m'%(hl_imped_hcell*1000)







from LHCMeasurementTools.TimberManager import parse_timber_file
import FillCalculator as fc

filln = 5219

fill_dict = {}
fill_dict.update(parse_timber_file('../fill_basic_data_csvs/basic_data_fill_%d.csv' % filln, verbose=False))
fill_dict.update(parse_timber_file('../fill_bunchbybunch_data_csvs/bunchbybunch_data_fill_%d.csv' % filln, verbose=False))
    
hl_imped_fill = fc.HeatLoad_calculated_fill(fill_dict, hl_calculator)

import pylab as pl
pl.close('all')
pl.plot(hl_imped_fill.t_stamps, hl_imped_fill.heat_load_calculated_per_beam_Wm['beam_1']*53., 'b')
pl.plot(hl_imped_fill.t_stamps, hl_imped_fill.heat_load_calculated_per_beam_Wm['beam_2']*53, 'r')

pl.show()
