import numpy as np



class HeatLoad_calculated_fill(object):
    def __init__(self, fill_dict, heat_load_calculator, Dt_calc = 60.):
 
        """
        Returns the half cell heat load for a fill dict, which has to consist of basic and bunchbybunch data
        """
        from LHCMeasurementTools.LHC_BCT import BCT
        from LHCMeasurementTools.LHC_FBCT import FBCT
        from LHCMeasurementTools.LHC_BQM import filled_buckets, blength
        from LHCMeasurementTools.LHC_Heatloads import magnet_length
        from LHCMeasurementTools.LHC_Energy import energy
        
        
        self.heat_load_calculator = heat_load_calculator
 
        ene = energy(fill_dict, beam=1)
 
        # Different for both beams
        self.heat_load_calculated_per_beam_Wm = {}
        for beam_ctr in (1,2):

            bct = BCT(fill_dict, beam=beam_ctr)
            fbct = FBCT(fill_dict, beam_ctr)
            bunch_length = blength(fill_dict, beam = beam_ctr)

            
            if beam_ctr == 1:
                self.t_stamps = np.arange(bct.t_stamps[0], bct.t_stamps[-1], Dt_calc)
                
            ppb = []
            energy_eV = []
            sigma_t = []
            
            for tt in self.t_stamps:
                fbct_trace = fbct.nearest_older_sample(tt)

                summed = np.sum(fbct_trace)
                if summed != 0.:
                    fbct_trace *= bct.nearest_older_sample(tt)/summed
                bl_trace = bunch_length.nearest_older_sample(tt)
                
                mask_no_bunch = bl_trace<1e-15
                bl_trace[mask_no_bunch] = 1.
                fbct_trace[mask_no_bunch] = 0.
    
                ppb.append(fbct_trace)
                sigma_t.append(bl_trace/4)
                energy_eV.append(ene.nearest_older_sample(tt)*1e9)
                
                
                
            ppb = np.array(ppb)
            sigma_t = np.array(sigma_t)
            energy_eV = np.array(energy_eV)
            
            self.heat_load_calculated_per_beam_Wm['beam_%d'%beam_ctr] = self.heat_load_calculator.calculate_P_Wm(ppb, sigma_t, energy_eV=energy_eV)
        
        self.heat_load_calculated_total = self.heat_load_calculated_per_beam_Wm['beam_1'] + self.heat_load_calculated_per_beam_Wm['beam_2']
               
