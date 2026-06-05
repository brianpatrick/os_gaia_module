# Process the Gaia Catalog of Nearby Stars:
# https://cdsarc.cds.unistra.fr/viz-bin/cat/J/MNRAS/508/3877#/browse

# Zack Reeves
# Created: 2024
# Edited by Cade Mohrhardt: 2026

# Versions:
#  1.1  FEB 2024 CREATE JUPYTER NOTEBOOK
#  1.02  MAY 2026 MODERNIZE CODE TO MATCH THE DIGITAL UNIVERSE

#  Python 3.12.12

import pandas as pd
import numpy as np
import sys
import os
import collections
import astropy.units as u
import astropy.coordinates
from astropy.table import Table
from astroquery.gaia import Gaia
from matplotlib import pyplot as plt, colors

sys.path.insert(0, '..')
from common import file_functions, calculations, asset_creation


GENERATE_SPECK = True
GENERATE_ASSET_FILE = True
READ_LOCAL_CATALOG = True


# Define the metadata for the data set. 
def generate_metadata():
    metadata = {}

    metadata['project'] = 'Digital Universe Atlas Gaia Subsets'
    metadata['sub_project'] = 'Gaia Catalog of Nearby Stars'

    metadata['catalog'] = 'The Gaia Catalogue of Nearby Stars (Gaia Collaboration, 2021)'
    metadata['catalog_author'] = 'Gaia Collaboration'
    metadata['catalog_year'] = '2021'
    metadata['catalog_doi'] = 'doi:10.1051/0004-6361/202039498'
    metadata['catalog_bibcode'] = '2021A&A...649A...6G'
    metadata['prepared_by'] = 'Brian Abbott, Zack Reeves, Cade Mohrhardt'
    metadata['version'] = '1.02'

    metadata['dir'] = metadata['sub_project'].replace(' ', '_').lower()
    metadata['raw_data_dir'] = ''

    metadata['data_group_title'] = 'GCNS'
    metadata['data_group_desc'] = 'Nearby stars in the Milky Way mapped by Gaia'
    metadata['data_group_desc_long'] = 'Have you ever wondered what is out there in space? Now, thanks to Gaia EDR3, the solar neighbourhood has been mapped with great precision out to 100 pc (326 light years)'
    metadata['fileroot'] = 'gcns'

    metadata['OS_identifier'] = metadata['sub_project']
    metadata['OS_gui_name'] = metadata['data_group_title']
    metadata['OS_gui_path'] = '/Milky Way/Stars'
    metadata['OS_gui_description'] = 'Thanks to Gaia EDR3, the solar neighbourhood has been mapped with great precision out to 100 pc (326 light years)'

    return metadata


metadata = generate_metadata()

file_functions.generate_license_file(metadata)
#file_functions.generate_asset_file(metadata)



if READ_LOCAL_CATALOG == False:
    #reading in the data
    #log in to Gaia Server - Can change to different credentials file for a different user
    Gaia.login(credentials_file='../common/gaia_credentials.txt')
        
    #get username from credentials file for query
    with open('../common/gaia_credentials.txt', 'r') as file:
        username = file.readline()

    #Query Gaia DR3 source for parallaxes
    #add distances and error
    job = Gaia.launch_job_async("select a.source_id, a.ra, a.dec, a.dist_50, a.dist_16, a.dist_84, a.pmra, a.pmdec, a.adoptedrv as radial_velocity, "
                                "bj.r_med_geo, bj.r_hi_geo, bj.r_lo_geo, bj.r_med_photogeo, bj.r_hi_photogeo, bj.r_lo_photogeo, "
                                "c.phot_g_mean_mag, c.bp_g, c.teff_gspphot "
                                "from external.gaiaedr3_gcns_main_1 a inner join external.gaiaedr3_distance bj on a.source_id = bj.source_id "
                                "inner join gaiadr3.gaia_source c on a.source_id = c.source_id",
                                dump_to_file=False)

    #Put the resulting table into a Table
    data = job.get_results()
    
    #Download the query results to a csv
    download = data.to_pandas()
    download.to_csv('raw_data/gcnsquery.csv', index=False)

    #Deleting job from Gaia ESA server so we don't clog the memory
    Gaia.remove_jobs(job.jobid)

    Gaia.logout()

if READ_LOCAL_CATALOG == True:
    data = pd.read_csv('raw_data/gcnsquery.csv')
    data = Table.from_pandas(data)



#setting units and metadata for important columns
data['ra'] = data.MaskedColumn(data=data['ra'], 
                                      unit=u.deg,
                                      meta = collections.OrderedDict([('ucd', 'pos.eq.ra')]),
                                      format='{:.6f}', 
                                      description='Right Ascension')

data['dec'] = data.MaskedColumn(data=data['dec'], 
                                      unit=u.deg,
                                      meta = collections.OrderedDict([('ucd', 'pos.eq.dec')]),
                                      format='{:.6f}', 
                                      description='Declination')

data['pmra'] = data.MaskedColumn(data=data['pmra'], 
                                       unit=u.mas/u.yr,
                                       meta = collections.OrderedDict([('ucd', 'pos.eq.ra')]),
                                       format='{:.6f}', 
                                       description='Proper Motion of RA')

data['pmdec'] = data.MaskedColumn(data=data['pmdec'], 
                                       unit=u.mas/u.yr,
                                       meta = collections.OrderedDict([('ucd', 'pos.eq.dec')]),
                                       format='{:.6f}', 
                                       description='Proper Motion of DEC')

data['radial_velocity'] = data.MaskedColumn(data=data['radial_velocity'], 
                                       unit=u.km/u.s,
                                       meta = collections.OrderedDict([('ucd', 'pos.eq.dec')]),
                                       format='{:.6f}', 
                                       description='Radial Velocity')


len(data[data['teff_gspphot']>0])

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

data.remove_rows(np.where(data['bj_distance']>500)[0])

#calculating distance in light years and parsecs
calculations.get_distance(data, dist='bj_distance', use='distance')

#calculating cartesian coordinates
calculations.get_cartesian(data, ra='ra', dec='dec', pmra='pmra', pmde='pmdec', radial_velocity='radial_velocity', frame='icrs')


#calculating absolute magnitudes
#calculate absolute V mag based on apparent magnitude and distance
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
data['lum'] = [10**(1.89 - 0.4*data['absmag'][i]) for i in range(len(data))]
small_luminosities = np.where((data['lum']>0.0) & (data['lum']<0.001))[0]
data['lum'][small_luminosities] = [0.001]*len(small_luminosities)

data['lum'] = data.MaskedColumn(data=data['lum'],
                             unit=u.solLum,
                             meta=collections.OrderedDict([('ucd', 'phys.luminosity')]),
                             format='{:.6f}',
                             description='Stellar Luminosity')

#setting color and visualizing
data['color'] = data.MaskedColumn(data=data['bp_g'],
                             unit=u.solLum,
                             meta=collections.OrderedDict([('ucd', 'phys.color')]),
                             format='{:.2f}',
                             description='Gaia BP-G color')
plt.hist(data['color'], bins=250)


# #2D Visualization
# fig, ax = plt.subplots(1, 2)

# #XY Plane
# ax[0].scatter(data['x'], data['y'])
# ax[0].set_title('XY Plane')

# #XZ Plane
# ax[1].scatter(data['x'], data['z'])
# ax[1].set_title('XZ Plane')

# #set good spacing
# fig.tight_layout()
# fig.set_size_inches(10, 4, forward=True)
# plt.show

# #2D Density Visualization
# fig, ax = plt.subplots(1, 2)

# #XY Plane
# ax[0].hist2d(data['x'], data['y'], 
#            bins = 200,  
#            norm = colors.LogNorm(),  
#            cmap = "RdYlGn_r",) 
# ax[0].set_title('XY Plane')

# #XZ Plane
# ax[1].hist2d(data['x'], data['z'], 
#            bins = 200,  
#            norm = colors.LogNorm(),  
#            cmap = "RdYlGn_r",) 
# ax[1].set_title('XZ Plane')

# #set good spacing
# fig.tight_layout()
# fig.set_size_inches(10, 4, forward=True)
# #plt.show


#construct a speck comment column
data['speck_label'] = data.Column(data=['#  '+str(name) for name in data['source_id']], 
                                  meta=collections.OrderedDict([('ucd', 'meta.id')]),
                                  description='Gaia EDR3 Source ID')

#construct a label column
data['label'] = ['GaiaEDR3_'+ str(source) for source in data['source_id']]  #leaving for now in case we want to add other labels

#setting texture number column
data['texnum'] = data.Column(data=[1]*len(data), 
                                  meta=collections.OrderedDict([('ucd', 'meta.texnum')]),
                                  description='Texture Number')


#Getting the column metadata
columns = file_functions.get_metadata(data, columns=['x', 'y', 'z', 'color', 'lum', 'absmag', 'appmag', 'texnum', 'dist_ly', 'dcalc', 'u', 'v', 'w', 'speed', 'speck_label'])



if GENERATE_SPECK:
    # Print the speck file using the to_speck function in file_functions
    file_functions.to_speck(metadata, Table.to_pandas(data), columns)


# Print the csv file using the to_csv function in file_functions
file_functions.to_csv(metadata, Table.to_pandas(data), columns)

# Print the label file using the to_label function in file_functions
file_functions.to_label(metadata, Table.to_pandas(data))


df = Table.to_pandas(data)
file_functions.generate_plot_pdf(df[columns['name']], metadata)



def asset_main():
    """Generate the asset file for stars"""

    metadata = generate_metadata()
    datainfo = {
        "renderable": "RenderableStars",
        "filename": metadata['fileroot'],
		"asset_dir": "",
        "local_modules": True,
        "data": {
            "File": metadata['fileroot']+".speck",
            "Name": metadata['data_group_title']+" Speck Files",
            "Identifier": "gaia_"+metadata['fileroot']+"_speck",
            "Version": 6
            },
        "Texture": {
            "Glare": "halo.png",
            "Core": "glare.png",
            "Name": "Stars Textures",
            "Identifier": "stars_textures",
            "Version": 1
            },
        "ColorMap":{
            "ColorMap": "colorbv.cmap",
            "OtherDataColorMap": "viridis.cmap",
            "Name": "Stars Color Table",
            "Identifier": "stars_colormap",
            "Version": 3
            },
        "Identifier": metadata['fileroot'],
        "Bv_column": "color",
        "Luminance_column": "lum",
        "AbsoluteMagnitude_column": "absmag",
        "ApparentMagnitude_column": "appmag",
        "Vx_column": "u",
        "Vy_column": "v",
        "Vz_column": "w",
        "Speed_column": "speed",
        "GUI": {
            "Name": metadata['OS_gui_name'],
            "Path": metadata['OS_gui_path'],
            "Description": metadata['OS_gui_description']
            },
        "meta_name": metadata['sub_project'],
		"author": metadata['prepared_by']
    }
    asset_creation.write_asset(datainfo)


if __name__ == "__main__" and GENERATE_ASSET_FILE:
    asset_main()
    print("Asset file for "+metadata['data_group_title']+" generated successfully.")