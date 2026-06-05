# PROCESS THE GAIA CATALOG OF RGB Stars:
# https://ui.adsabs.harvard.edu/abs/2023ApJS..267....8A/abstract
# https://zenodo.org/records/7945154

# Zack Reeves
# Created: 2024
# Edited by Cade Mohrhardt: 2026

# Versions:
#  1.1  JUN 2024 CREATE JUPYTER NOTEBOOK
#  1.02  MAY 2026 MODERNIZE CODE TO MATCH THE DIGITAL UNIVERSE

#  Python 3.12.12

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


# Define the metadata for the data set.
#https://ui.adsabs.harvard.edu/abs/2023A%26A...674A..39G/abstract
def generate_metadata():
    metadata = {}

    metadata['project'] = 'Digital Universe Atlas Gaia Subsets'
    metadata['sub_project'] = 'Red Giant Branch Stars'

    metadata['catalog'] = 'Robust Data-driven Metallicities for 175 Million Stars from Gaia XP Spectra (Andrae, 2023)'
    metadata['catalog_author'] = 'Andrae+'
    metadata['catalog_year'] = '2023'
    metadata['catalog_doi'] = 'doi:10.3847/1538-4365/acd53e'
    metadata['catalog_bibcode'] = '2023ApJS..267....8A'

    metadata['prepared_by'] = 'Brian Abbott, Zack Reeves, Cade Mohrhardt'
    metadata['version'] = '1.1'

    metadata['dir'] = metadata['sub_project'].replace(' ', '_').lower()
    metadata['raw_data_dir'] = ''

    metadata['data_group_title'] = 'Giants'
    metadata['data_group_desc'] = 'Red Giant Branch Stars' 
    metadata['data_group_desc_long'] = 'The Sun is the reference point in much of stellar astronomy and astrophysics. Solar analogues are stars that resemble the Sun in terms of a restricted set of parameters. In contrast to the Sun, they can be observed in the night sky and with the very same instruments used to study stars in the Milky Way.'
    metadata['fileroot'] = 'giants'

    metadata['OS_identifier'] = metadata['sub_project']
    metadata['OS_gui_name'] = metadata['data_group_title']
    metadata['OS_gui_path'] = '/Milky Way/Stars'
    metadata['OS_gui_description'] = metadata['data_group_desc_long']

    return metadata


metadata = generate_metadata()

file_functions.generate_license_file(metadata)
#file_functions.generate_asset_file(metadata)



#download the data from https://zenodo.org/records/7945154 
#~12 million stars
data = Table.read('raw_data/table_2_catwise.fits.gz')



#calculating distance in light years and parsecs
#this dataset only uses gaia parallaxes to calculate distance to avoid the cpmutational expense of uploading >3 million stars to grab BJ distances
data['parallax'].unit=u.mas
calculations.get_distance(data, parallax='parallax')

#setting metadata for dcalc
data['dcalc'] = data.Column([3]*len(data),
                            meta=collections.OrderedDict([('ucd', 'meta.dcalc')]),
                            description='Distance Indicator: 1 indicates a Bailer-Jones photogeometric distance; 2 indicates a Bailer-Jones geometric distance; 3 indicates a Gaia parallax-based distance')


#setting necessary units and calculating galactic cartesian XYZ
data['ra'].unit=u.deg
data['dec'].unit=u.deg
data['pmra'].unit=u.mas/u.yr
data['pmdec'].unit=u.mas/u.yr
data['radial_velocity'].unit=u.km/u.s

calculations.get_cartesian(data, ra='ra', dec='dec', pmra='pmra', pmde='pmdec', radial_velocity='radial_velocity', frame='icrs')


#setting necessary units
data['phot_g_mean_mag'].unit=u.mag
data['phot_bp_mean_mag'].unit=u.mag
data['phot_rp_mean_mag'].unit=u.mag

#calculating absolute and apparent magnitudes, luminosity, and color
gaia_functions.get_magnitudes(data)
gaia_functions.get_luminosity(data)
data['bp_rp'] = [data['phot_bp_mean_mag'][i]-data['phot_rp_mean_mag'][i] for i in range(len(data))]
gaia_functions.get_bp_g_color(data, color='bp_rp')


# plt.hist(data['bp_rp'], bins=250);

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