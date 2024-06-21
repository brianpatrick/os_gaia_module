# PY SCRIPT TO PERFORM COMMON FUNCTIONS TO PROCESS GAIA DATA
# 
# 
#
#
# ZACK REEVES
# CREATED: 2024
#
# VERSIONS:
#  1.1  JUNE 2024 CREATE PY FILE


import pandas as pd
import numpy as np
import sys
import collections

import astropy.units as u
import astropy.coordinates
from astropy.table import Table, join, vstack

sys.path.insert(0, '..')
from common import file_functions, calculations





#setting dcalc and distance
#requires r_med_geo and r_med_photogeo
def set_bj_distance(data:Table):
    #setting dcalc based on r_med_geo (if>500pc and photogeo exists, we choose photogeo and set dcalc to 1, else geo and dcalc to 2)
    data['dcalc'] = [1 if((not(np.ma.is_masked(data['r_med_photogeo'][i])))and(data['r_med_geo'][i]>500)) else 2 for i in range(len(data))]
    
    #setting metadata for dcalc
    data['dcalc'] = data.Column(data['dcalc'],
                                meta=collections.OrderedDict([('ucd', 'meta.dcalc')]),
                                description='Distance Indicator: 1 indicates a Bailer-Jones photogeometric distance; 2 indicates a Bailer-Jones geometric distance')
    
    #Choosing distance based on dcalc
    data['bj_distance'] = [data['r_med_photogeo'][i] if data['dcalc'][i]==1 else data['r_med_geo'][i] for i in range(len(data))]
    data['bj_distance'].unit=u.pc
    
    #Choosing and calculating distance error based on the distance we chose
    data['e_bj_dist'] = [((data['r_hi_photogeo'][i]-data['r_lo_photogeo'][i])/2)*u.pc if((not(np.ma.is_masked(data['r_med_photogeo'][i])))and(data['r_med_geo'][i]>500)) else ((data['r_hi_geo'][i]-data['r_lo_geo'][i])/2)*u.pc for i in range(len(data))]



#calculating absolute magnitudes and setting a column for apparent magnitudes
#currently uses Gaia green band magnitude
#requires phot_g_mean_mag and dist_pc
def get_magnitudes(data:Table):
    data['appmag'] = data.MaskedColumn(data=data['phot_g_mean_mag'],
                                       unit=u.mag,
                                       meta=collections.OrderedDict([('ucd', 'phot.mag;em.opt.G')]),
                                       format='{:.6f}',
                                       description='Apparent magnitude in Gaia G-band')
    data['absmag'] = data.MaskedColumn(data=[data['appmag'][i]+5-5*np.log10(data['dist_pc'][i]) for i in range(len(data))],
                                 unit=u.mag,
                                 meta=collections.OrderedDict([('ucd', 'phot.magAbs;em.opt.G')]),
                                 format='{:.6f}',
                                 description='Absolute magnitude in Gaia G-band')



#calculate luminosity based on absolute magnitude
def get_luminosity(data:Table):
    data['lum'] = [10**(1.89 - 0.4*data['absmag'][i]) for i in range(len(data))]
    small_luminosities = np.where((data['lum']>0.0) & (data['lum']<0.001))[0]
    data['lum'][small_luminosities] = [0.001]*len(small_luminosities)
    
    data['lum'] = data.MaskedColumn(data=data['lum'],
                                    unit=u.solLum,
                                    meta=collections.OrderedDict([('ucd', 'phys.luminosity')]),
                                    format='{:.6f}',
                                    description='Stellar Luminosity')



#setting color
def get_bp_g_color(data:Table):
    data['color'] = data.MaskedColumn(data=data['bp_g'],
                                      unit=u.solLum,
                                      meta=collections.OrderedDict([('ucd', 'phys.color')]),
                                      format='{:.2f}',
                                      description='Gaia BP-G color')