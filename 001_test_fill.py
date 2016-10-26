import sys, os
BIN = os.path.expanduser("../")
sys.path.append(BIN)



# Load fills
from LHCMeasurementTools.TimberManager import parse_timber_file

filln = 5219

fill_dict = {}
fill_dict.update(parse_timber_file('../fill_basic_data_csvs/basic_data_fill_%d.csv' % filln, verbose=False))
fill_dict.update(parse_timber_file('../fill_bunchbybunch_data_csvs/bunchbybunch_data_fill_%d.csv' % filln, verbose=False))



# Build heat load calculators
import impedance_heatload as ihl
import synchrotron_radiation_heatload as srhl 

hli_calculator  = ihl.HeatLoadCalculatorImpedanceLHCArc()
hlsr_calculator  = srhl.HeatLoadCalculatorSynchrotronRadiationLHCArc()

# Use fill calculators
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
