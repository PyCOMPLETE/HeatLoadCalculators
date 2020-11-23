import sys, os
BIN = os.path.expanduser("../")
sys.path.append(BIN)



# Load fills
import LHCMeasurementTools.TimberManager as tm

filln = 6740

fill_dict = {}
fill_dict.update(tm.CalsVariables_from_h5(
    ('../../LHC_followup_download_scripts/fill_basic_data_h5s'
    '/basic_data_fill_%d.h5' % filln)))
fill_dict.update(tm.CalsVariables_from_h5(
    ('../../LHC_followup_download_scripts/fill_bunchbybunch_data_h5s/'
    'bunchbybunch_data_fill_%d.h5' % filln)))



# Build heat load calculators
import impedance_heatload as ihl
import synchrotron_radiation_heatload as srhl

hli_calculator  = ihl.HeatLoadCalculatorImpedanceLHCArc()
hlsr_calculator  = srhl.HeatLoadCalculatorSynchrotronRadiationLHCArc()

# Use fill calculator
import FillCalculator as fc
    
hl_imped_fill = fc.HeatLoad_calculated_fill(fill_dict, hli_calculator)
hl_sr_fill = fc.HeatLoad_calculated_fill(fill_dict, hlsr_calculator)


# Plot
import pylab as pl
pl.close('all')
pl.plot((hl_imped_fill.t_stamps-hl_imped_fill.t_stamps[0])/3600., hl_imped_fill.heat_load_calculated_per_beam_Wm['beam_1']*53., 'b')
pl.plot((hl_imped_fill.t_stamps-hl_imped_fill.t_stamps[0])/3600., hl_imped_fill.heat_load_calculated_per_beam_Wm['beam_2']*53, 'r')
pl.plot((hl_imped_fill.t_stamps-hl_imped_fill.t_stamps[0])/3600., hl_sr_fill.heat_load_calculated_per_beam_Wm['beam_1']*53., '--b')
pl.plot((hl_imped_fill.t_stamps-hl_imped_fill.t_stamps[0])/3600., hl_sr_fill.heat_load_calculated_per_beam_Wm['beam_2']*53, '--r')

pl.ylabel('W/hcell/beam')
pl.xlabel('Time [h]')

pl.show()
