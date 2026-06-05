# PROCESS THE GAIA DR3 RADIAL VELOCITY CATALOG:
# https://gea.esac.esa.int/archive/

# Zack Reeves
# Created: 2024
# Edited by Cade Mohrhardt: 2026

# Versions:
#  1.1  MAR 2024 CREATE JUPYTER NOTEBOOK
#  1.02  MAY 2026 MODERNIZE CODE TO MATCH THE DIGITAL UNIVERSE

#  Python 3.12.12


# **STEPS TO RUN THIS CODE**

# There are two queries you can run to generate this catalog
# The first is a query to grab every star in Gaia with a radial velocity (and a gmag) ~33 million stars
# The second is a query to grab every Gaia star with a radial velocity AND which pass some data quality tests (parallax error) ~29 million stars

# First choose a query, then run it.  
# If you have a decent internet connection and can run the code for hours, the code will execute the query and download the table
# If the code times out, the query results should still be available on your gaia archive account

# Download the data either from the code or from the gaia archive, then run the rest of the processing code.  This can take several hours
# Recommended to download as a vot.gz - this is the smallest file size the gaia server offers
# If any errors occur, consider slicing the data down to the first 1000 rows (data[:1000]) to debug
# Can also add "select TOP 1000" to the query to grab 1000 stars for testing

#This code pulls a catalogue from Gaia DR3 of each star that has a reliable radial velocity
#Based on https://www.aanda.org/articles/aa/full_html/2023/06/aa44220-22/aa44220-22.html#S14,
#We correct radial velocities for stars of high magnitude grvs_mag>11 by Katz et al.
#We also correct stars of high effective temperature 14500>rv_template_teff>8500 and 6>grvs_mag>12 by Blomme et al. as recommended by Katz^

import pandas as pd
import numpy as np
import sys
import os
import collections
import astropy.units as u
import astropy.coordinates
from astropy.table import Table
from astropy.io import ascii
from astropy.io import fits
from astroquery.gaia import Gaia
from matplotlib import pyplot as plt, colors

sys.path.insert(0, '..')
from common import file_functions, calculations, gaia_functions, asset_creation


GENERATE_SPECK = True
GENERATE_ASSET_FILE = True
READ_LOCAL_CATALOG = True


# Define the metadata for the data set. 
def generate_metadata():
    metadata = {}

    metadata['project'] = 'Digital Universe Atlas Gaia Subsets'
    metadata['sub_project'] = 'Gaia DR3 Radial Velocities'

    metadata['catalog'] = 'Gaia Data Release 3: Properties and validation of the radial velocities (Katz et al., 2023)'
    metadata['catalog_author'] = 'Katz et al.'
    metadata['catalog_year'] = '2023'
    metadata['prepared_by'] = 'Zack Reeves (AMNH), Cade Mohrhardt (AMNH)'
    metadata['version'] = '1.1'

    metadata['dir'] = metadata['sub_project'].replace(' ', '_').lower()
    metadata['raw_data_dir'] = ''

    metadata['data_group_title'] = 'Radial Velocity Stars'
    metadata['data_group_desc'] = 'Gaia DR3 Radial Velocity'
    metadata['data_group_desc_long'] = 'Gaia DR3 Stars with Radial Velocity'
    metadata['fileroot'] = 'gdr3rv'

    metadata['OS_identifier'] = metadata['sub_project']
    metadata['OS_gui_name'] = metadata['data_group_title']
    metadata['OS_gui_path'] = '/Milky Way/Stars'
    metadata['OS_gui_description'] = metadata['data_group_desc_long']

    return metadata


metadata = generate_metadata()

file_functions.generate_license_file(metadata)
#file_functions.generate_asset_file(metadata)



#query the catalogue from Gaia
#https://gea.esac.esa.int/archive/
#The query pulls the source id, positional data, velocity data as well as teff and magnitude for correction purposes
#corrective data included in the query was informed by https://www.aanda.org/articles/aa/full_html/2023/06/aa44220-22/aa44220-22.html#S14
# if READ_LOCAL_CATALOG == False:
#     #QUERY #1 - 33 million stars

#     #log in to Gaia Server - Can change to different credentials file for a different user
#     #query runs in a little over an hour
#     #file is 3.2 gigabytes, 33,653,049 objects
#     Gaia.login(credentials_file='../common/gaia_credentials.txt')

#     #Query Gaia DR3 source for parallaxes
#     job = Gaia.launch_job_async("select a.source_id, a.ra, a.dec, a.pmra, a.pmdec, a.parallax, a.parallax_error, a.phot_g_mean_mag, a.bp_g, a.radial_velocity, a.radial_velocity_error, a.grvs_mag, a.rv_template_teff, "
#                                 "bj.r_med_geo, bj.r_hi_geo, bj.r_lo_geo, bj.r_med_photogeo, bj.r_hi_photogeo, bj.r_lo_photogeo "
#                                 "from gaiadr3.gaia_source a left join external.gaiaedr3_distance bj on a.source_id = bj.source_id "
#                                 "where a.radial_velocity is not null and a.phot_g_mean_mag > 0 and parallax > 0",
#                                 dump_to_file=False)

#     #Put the resulting table into a Table
#     data = job.get_results()

#     # Uncomment to download the query results to a csv
#     #download = data.to_pandas()
#     #download.to_csv('raw_data/dr3rvquery.csv', index=False)
    
#     Gaia.remove_jobs(job.jobid)

#     Gaia.logout()

if READ_LOCAL_CATALOG == False:
    #QUERY #2 - 29 million stars

    #log in to Gaia Server - Can change to different credentials file for a different user
    #query runs in a little over an hour
    #file is 2.8 gigabytes, 29,946,388 objects
    Gaia.login(credentials_file='../common/gaia_credentials.txt')

    #Query Gaia DR3 source for parallaxes
    job = Gaia.launch_job_async("select a.source_id, a.ra, a.dec, a.pmra, a.pmdec, a.parallax, a.parallax_error, a.phot_g_mean_mag, a.bp_g, a.radial_velocity, a.radial_velocity_error, a.grvs_mag, a.rv_template_teff, "
                                "bj.r_med_geo, bj.r_hi_geo, bj.r_lo_geo, bj.r_med_photogeo, bj.r_hi_photogeo, bj.r_lo_photogeo "
                                "from gaiadr3.gaia_source a left join external.gaiaedr3_distance bj on a.source_id = bj.source_id "
                                "where a.radial_velocity is not null and a.phot_g_mean_mag > 0 and parallax > 0 and a.parallax / a.parallax_error > 5",
                                dump_to_file=False)

    #Put the resulting table into a Table
    data = job.get_results()
        
    #Download the query results to a csv
    download = data.to_pandas()
    download.to_csv('raw_data/dr3rvquery.csv', index=False)
    
    Gaia.remove_jobs(job.jobid)

    Gaia.logout()

if READ_LOCAL_CATALOG == True:
    data = pd.read_csv('raw_data/dr3rvquery.csv')
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


gaia_functions.set_bj_distance(data)

#calculating distance in light years and parsecs
calculations.get_distance(data, dist='bj_distance', use='distance')

# data quality check
data.remove_rows(np.where(data['dist_pc']<=0)[0])

gaia_functions.get_magnitudes(data)

gaia_functions.get_luminosity(data)

gaia_functions.get_bp_g_color(data) #may want to change to bp_rp

# data check on the G mag
x = data['phot_g_mean_mag']
q25, q75 = np.percentile(x, [25, 75])
bin_width = 2 * (q75 - q25) * len(x) ** (-1/3)
bins = round((x.max() - x.min()) / bin_width)
print("Freedman–Diaconis number of bins:", bins)
plt.hist(x, bins=bins);

#applying corrections from papers
data['radial_velocity_correction'] = [0.0]*len(data)

#Katz correction
katz_indexes = np.where(data['grvs_mag']>11)[0]
data['radial_velocity_correction'][katz_indexes] = [(0.02755*data['grvs_mag'][i]**2 - 0.55863*data['grvs_mag'][i] + 2.81129) for i in katz_indexes]

#Blomme correction
blomme_indexes = np.where((data['grvs_mag']>11)&(data['rv_template_teff']>8500)&(data['rv_template_teff']<14500))[0]
data['radial_velocity_correction'][blomme_indexes] = [7.98 - 1.135*data['grvs_mag'][i] for i in blomme_indexes]

data['corrected_radial_velocity'] = np.subtract(data['radial_velocity'], data['radial_velocity_correction'])
data['corrected_radial_velocity'].unit=u.km/u.s

#calculating cartesian coordinates
calculations.get_cartesian(data, ra='ra', dec='dec', pmra='pmra', pmde='pmdec', radial_velocity='corrected_radial_velocity', frame='icrs')


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
                                  description='Gaia DR3 Source ID')

#construct a label column
data['label'] = ['GaiaDR3_'+ str(source) for source in data['source_id']]  #leaving for now in case we want to add other labels

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