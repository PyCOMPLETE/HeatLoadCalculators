import mydownselscan as mds
import myfilemanager as mfm 
from scipy import interpolate
import numpy as np 


class HeatLoadCalculatorElectronCloud(object):
    def __init__(self, fname_sim_db, sey, varname_hl='hl', rescaling_factor=1. ):
        self.varname_hl = varname_hl
        self.rescaling_factor = rescaling_factor
        
        
        ob_all = mfm.myloadmat_to_obj(fname_sim_db, squeeze=False)
        ob_all = mds.downsel_scan(ob_all, {}, varname_list = [varname_hl])        
        self.database = mds.downsel_scan(ob_all, {'sey':sey}, varname_list = [varname_hl])
        self.matrix = getattr(self.database, self.varname_hl+'_sqz')
        
        
        # Build heat load functions (and change units)
        self.interpolating_function = interpolate.RectBivariateSpline(self.database.int*1e11, 
            self.database.bl*1e-9/4,self.matrix*rescaling_factor, kx=2, ky=2).ev
        
    
    def patch_missing(self):
        
        mat_to_patch = self.matrix
        
        for ii in range(len(mat_to_patch)):
            current_row = mat_to_patch[ii,:]
            mask_missing = (current_row==0)
            current_row[mask_missing] = np.interp(self.database.bl[mask_missing],
                            self.database.bl[~mask_missing],current_row[~mask_missing])
            mat_to_patch[ii,:] = current_row
            
        self.interpolating_function = interpolate.RectBivariateSpline(self.database.int*1e11, 
            self.database.bl*1e-9/4,self.matrix*self.rescaling_factor, kx=2, ky=2).ev
    
    def calculate_P_Wm(self, ppb, sigma_t, energy_eV=0., n_bunches=None):
        
        return self.interpolating_function(ppb, sigma_t)*n_bunches
        
        
