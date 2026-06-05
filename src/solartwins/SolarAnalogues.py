# Process the Gaia catalogue of Solar Analogues:
# Gaia has astrophysical parameters which can be used to find stars similar to the Sun
# Based on https://ui.adsabs.harvard.edu/abs/2023A%26A...674A..39G/abstract , we use Teff, logg, and M/H to find solar analogues

# Zack Reeves
# Created: 2024
# Edited by Cade Mohrhardt: 2026
#
# Versions:
#  1.1  JUN 2024 CREATE JUPYTER NOTEBOOK
#  1.02  MAY 2026 MODERNIZE CODE TO MATCH THE DIGITAL UNIVERSE

#  Python 3.12.12 OCT 2025

import pandas as pd
import numpy as np
import sys
import os
import collections
import astropy.units as u
import astropy.coordinates
from astropy.table import Table, join, vstack
from astropy.io import ascii
from astroquery.gaia import Gaia
from matplotlib import pyplot as plt, colors

sys.path.insert(0, '..')
from common import file_functions, calculations, gaia_functions, asset_creation


GENERATE_SPECK = True
GENERATE_ASSET_FILE = True
READ_LOCAL_CATALOG = True


# Define the metadata for the data set.
#https://ui.adsabs.harvard.edu/abs/2023A%26A...674A..39G/abstract
def generate_metadata():
    metadata = {}

    metadata['project'] = 'Digital Universe Atlas Gaia Subsets'
    metadata['sub_project'] = 'Solar Twins'

    metadata['catalog'] = 'Gaia Data Release 3. A golden sample of astrophysical parameters (Gaia Collaboration, 2023)'
    metadata['catalog_author'] = 'Gaia Collaboration'
    metadata['catalog_year'] = '2023'
    metadata['catalog_doi'] = 'doi:10.1051/0004-6361/202243800'
    metadata['catalog_bibcode'] = '2023A&A...674A..39G'

    metadata['prepared_by'] = 'Brian Abbott, Zack Reeves, Cade Mohrhardt'
    metadata['version'] = '1.1'

    metadata['dir'] = metadata['sub_project'].replace(' ', '_').lower()
    metadata['raw_data_dir'] = ''

    metadata['data_group_title'] = 'Solar Twins'
    metadata['data_group_desc'] = 'Stars similar to the Sun in the Milky Way mapped by Gaia'
    metadata['data_group_desc_long'] = 'The Sun is the reference point in much of stellar astronomy and astrophysics. Solar twins are stars that resemble the Sun in terms of a restricted set of parameters. In contrast to the Sun, they can be observed in the night sky and with the very same instruments used to study stars in the Milky Way.'
    metadata['fileroot'] = 'solar_twins'

    metadata['OS_identifier'] = metadata['sub_project']
    metadata['OS_gui_name'] = metadata['data_group_title']
    metadata['OS_gui_path'] = '/Milky Way/Stars'
    metadata['OS_gui_description'] = metadata['data_group_desc_long']

    return metadata


metadata = generate_metadata()

file_functions.generate_license_file(metadata)
#file_functions.generate_asset_file(metadata)



#query Gaia

#based on the paper, we want to find stars with teff, logg, and M/H within reasonable error of the Sun's
#we first select candiates based on their data reliability.  We constrain our selection to stars with apparent magnitude brighter than G<16.  We also use parallax/parallax error > 20 to ensure reliable distance metrics
#we then choose candidates based on the following GSP Spec criteria:
# - Teff must be within 100K of 5772K
# - logg must be within 0.25 of 4.44
# - [M/H] must be within 0.2 of 0.0
# - good gspspec flag: 0 for characters 1-13 except for 8 which can be 0 or 1
#we then futher thresh our sample by checking the FLAME parameters:
# - stellar mass measured by FLAME must be within 0.05 M. of 1 M.
# - stellar radius measured by FLAME must be within 0.2 R. of 1 R.

#we don't actually have to craft this query ourselves, this list is published in the Gaia archive
if READ_LOCAL_CATALOG == False:
    #log in to Gaia Server - Can change to different credentials file for a different user
    #query runs in less than a minute
    #file is small, 5683 objects
    Gaia.login(credentials_file='../common/gaia_credentials.txt')

    #Query Gaia DR3 **we can add more params later, maybe get metallicity

    job = Gaia.launch_job_async("select a.source_id, "
                                "b.ra, b.dec, b.pmra, b.pmdec, b.parallax, b.parallax_error, b.phot_g_mean_mag, b.bp_g, b.radial_velocity, b.radial_velocity_error, b.grvs_mag, b.rv_template_teff, "
                                "bj.r_med_geo, bj.r_hi_geo, bj.r_lo_geo, bj.r_med_photogeo, bj.r_hi_photogeo, bj.r_lo_photogeo "
                                "from gaiadr3.gold_sample_solar_analogues a inner join gaiadr3.gaia_source b on a.source_id = b.source_id "
                                "left join external.gaiaedr3_distance bj on a.source_id = bj.source_id",
                                dump_to_file=False)

    #Put the resulting table into a Table
    data = job.get_results()

    #Download the query results to a csv
    download = data.to_pandas()
    download.to_csv('raw_data/solaranaloguesquery.csv', index=False)

    Gaia.remove_jobs(job.jobid)

    Gaia.logout()

if READ_LOCAL_CATALOG == True:
    data = pd.read_csv('raw_data/solaranaloguesquery.csv')
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

calculations.get_distance(data, dist='bj_distance', use='distance')

calculations.get_cartesian(data, ra='ra', dec='dec', pmra='pmra', pmde='pmdec', radial_velocity='radial_velocity', frame='icrs')

gaia_functions.get_magnitudes(data)

gaia_functions.get_luminosity(data)

gaia_functions.get_bp_g_color(data)


#2D Visualization
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
#plt.show


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