import sys, os
BIN = os.path.expanduser("../")
sys.path.append(BIN)



import impedance_heatload as ihl
import synchrotron_radiation_heatload as srhl


ppb_test = 1.15e11
sigma_t_test = 1e-9/4.
energy_eV_test=7000e9
n_bunches_test = 2808

hli_calculator  = ihl.HeatLoadCalculatorImpedanceLHCArc()
hlsr_calculator  = srhl.HeatLoadCalculatorSynchrotronRadiationLHCArc()


hl_imped_dip = hli_calculator.calculate_P_Wm_dipole(ppb=ppb_test, sigma_t=sigma_t_test,  energy_eV=energy_eV_test, n_bunches=n_bunches_test)
print('Impedance load in the dipoles: %.1f mW/m/beam'%(hl_imped_dip*1000))
hl_imped_quad = hli_calculator.calculate_P_Wm_quad(ppb=ppb_test, sigma_t=sigma_t_test,  energy_eV=energy_eV_test, n_bunches=n_bunches_test)  
print('Impedance load in the quadrupoles: %.1f mW/m/beam'%(hl_imped_quad*1000))
hl_imped_drift = hli_calculator.calculate_P_Wm_drift(ppb=ppb_test, sigma_t=sigma_t_test,  energy_eV=energy_eV_test, n_bunches=n_bunches_test)
print('Impedance load in the drifts: %.1f mW/m/beam'%(hl_imped_drift*1000))
hl_imped_hcell = hli_calculator.calculate_P_Wm(ppb=ppb_test, sigma_t=sigma_t_test,  energy_eV=energy_eV_test, n_bunches=n_bunches_test)  
print('Impedance load average half-cell: %.1f mW/m/beam'%(hl_imped_hcell*1000))


hl_imped_sr = hlsr_calculator.calculate_P_Wm(ppb=ppb_test, sigma_t=sigma_t_test,  energy_eV=energy_eV_test, n_bunches=n_bunches_test)
print('Synchrotron radiation load average half-cell: %.1f mW/m/beam'%(hl_imped_sr*1000))
