# PROCESS THE STELLAR CLUSTER CATAOLOG:
# https://arxiv.org/abs/2308.04546
# https://zenodo.org/records/10042028
#
#
# ZACK REEVES
# CREATED: 2023
# CADE MOHRHARDT
# UPDATED: 2025
#
# VERSIONS:
#  1.1  OCT 2023 CREATE JUPYTER NOTEBOOK
#  Python 3.12.12 OCT 2025

import pandas as pd
import numpy as np
import sys
import collections

from astropy.io import ascii
import astropy.units as u
import astropy.coordinates
from astropy.coordinates import Angle
from astropy.table import Table
import os

from astroquery.vizier import Vizier

sys.path.insert(0, '..')
from common import file_functions, calculations, asset_creation

import matplotlib.pyplot as plt
import pandas as pd
import sys
import collections
import astropy.units as u
from astropy.table import Table
from astropy.io import fits
import astropy.cosmology.units as cu

PRINT_CENSUS = True
GENERATE_ASSET_FILE = True

# Define the metadata for the data set. 
def generate_metadata():
    metadata = {}

    metadata['project'] = 'Digital Universe Atlas'
    metadata['sub_project'] = 'Stellar Clusters'

    metadata['catalog'] = 'The Unified Cluster Catalogue: towards a comprehensive and homogeneous data base of stellar clusters (Perren+, 2023)'
    metadata['catalog_author'] = 'Perren+'
    metadata['prepared_by'] = 'Zack Reeves (AMNH), Cade Mohrhardt (AMNH)'
    metadata['catalog_year'] = '2023'
    metadata['version'] = '1.1'

    metadata['dir'] = metadata['sub_project'].replace(' ', '_').lower()
    metadata['raw_data_dir'] = ''

    metadata['data_group_title'] = 'Stellar Clusters'
    metadata['data_group_desc'] = 'Stellar Cluster catalog'
    metadata['data_group_desc_long'] = 'This is a large catalog containing "Milky Way open clusters" using Gaia DR3 data.'
    metadata['fileroot'] = 'sc'
    return metadata

metadata = generate_metadata()

file_functions.generate_license_file(metadata)
#file_functions.generate_asset_file(metadata)
#file_functions.generate_asset_file(metadata, RenderableType="RenderableGaiaStars")


#Reading data into Astropy Table

#The catalog is downloaded from https://zenodo.org/records/10042028
#The catalog that should be downloaded to view the 14,000 clusters is UCC_cat.csv.gz, 
#the other catalog is the 1,300,000 stellar members

data = Table.from_pandas(pd.read_csv('UCC_cat.csv.gz'))

#view the data
#print(data.columns)

#setting units and metadata for important columns (ID, RA, DEC, parallax, N_50, r_50)

data['ID'] = data.MaskedColumn(data=data['ID'], 
                               meta = collections.OrderedDict([('ucd', 'meta.id')]),
                               description='Median parallax estimated using the selected members')

data['RA_ICRS_m'] = data.MaskedColumn(data=data['RA_ICRS_m'], 
                                      unit=u.deg,
                                      meta = collections.OrderedDict([('ucd', 'pos.eq.ra')]),
                                      format='{:.6f}', 
                                      description='Median RA estimated using the selected members')

data['DE_ICRS_m'] = data.MaskedColumn(data=data['DE_ICRS_m'], 
                                      unit=u.deg,
                                      meta = collections.OrderedDict([('ucd', 'pos.eq.de')]),
                                      format='{:.6f}', 
                                      description='Median DE estimated using the selected members')

data['plx_m'] = data.MaskedColumn(data=data['plx_m'], 
                                  unit=u.mas,
                                  meta = collections.OrderedDict([('ucd', 'pos.parallax.trig')]),
                                  format='{:.6f}', 
                                  description='Median parallax estimated using the selected members')

data['N_50'] = data.MaskedColumn(data=data['N_50'], 
                                 dtype=int,
                                 meta = collections.OrderedDict([('ucd', 'meta.number')]),
                                 description='Number of estimated members with P>0.5')

data['r_50'] = data.MaskedColumn(data=data['r_50'], 
                                 unit=u.arcmin,
                                 meta = collections.OrderedDict([('ucd', 'phys.ang.size')]),
                                 format='{:.6f}',
                                 description='Radius that contains half the members (in arcmin)')

#Thresh data based on parallax
#some rows have plx_m <= 0 so we start there

data.remove_rows(np.where(data['plx_m']<=0.0)[0])

#calculating distance in light years and parsecs
calculations.get_distance(data, parallax='plx_m', use='parallax')

#calculating cartesian coordinates
calculations.get_cartesian(data, ra='RA_ICRS_m', dec='DE_ICRS_m')

#playing around with threshing on distance
data.remove_rows(np.where(data['dist_pc']>20000)[0])

#2D Visualization
fig, ax = plt.subplots(1, 2)

#XY Plane
ax[0].scatter(data['x'], data['y'])
ax[0].set_title('XY Plane')

#XZ Plane
ax[1].scatter(data['x'], data['z'])
ax[1].set_title('XZ Plane')

#set good spacing
fig.tight_layout()
fig.set_size_inches(10, 4, forward=True)
plt.show

#construct a speck comment column
data['speck_label'] = data.Column(data=['#__'+name for name in data['ID']], 
                                  meta=collections.OrderedDict([('ucd', 'meta.id')]),
                                 description='Object ID')

#construct a label column
data['label'] = data['ID']  #leaving for now in case we want to add other labels

#construct a metadata table
columns = file_functions.get_metadata(data, columns=['x', 'y', 'z', 'dist_ly', 'N_50', 'r_50', 'speck_label'])
columns

# Print the speck file using the to_speck function in file_functions
file_functions.to_speck(metadata, Table.to_pandas(data), columns)

# Print the label file using the to_label function in file_functions
file_functions.to_label(metadata, Table.to_pandas(data))

file_functions.to_csv(metadata, Table.to_pandas(data), columns)

def asset_main():
    """Generate the asset file for Gaia Stellar Clusters."""

    datainfo = {
        # Core script requirements
        "renderable": "RenderablePolygonCloud",
        "filename": metadata['fileroot'],
        "Identifier": "GaiaStellarClusters",
        "local_modules": True,
        "asset_dir": "",
        "Enabled": "true",
        "name": "Gaia Stellar Clusters",
        "data": {
            "File": "sc.speck",
            "Name": "Stellar Clusters Speck Files",
            "Identifier": "gaia_stellarclusters_speck",
            "Version": 3
            },
        "Texture": {
            "File": "stars_textures"
        },
        
        "local colormaps": """asset.resource({
    Name = "Stars Color Table",
    Type = "HttpSynchronization",
    Identifier = "stars_colormap",
    Version = 3
})""",

        "local textures": """asset.resource({
    Name = "Stars Textures",
    Type = "HttpSynchronization",
    Identifier = "stars_textures",
    Version = 1
})""",

        # Renderable Parameters
        "opacity": 0.9,
        "PolygonSides": "12",
        "Unit": "pc",
        "FixedColor": "0.0, 0.50, 0.50",
        
        # Label Settings
        "labels": {
            "file": "sc.label",
            "color": "{ 0.0, 0.36, 0.14 }",
            "size": 15.5,
            "min_max": "{ 4, 30 }"
        },
        
        # Size Settings
        "ScaleExponent": "15.7",
        "MaxSize": "23.0",
        "enable_max_size": "true",

        # GUI & Meta
        "GUI": {
        "Path": "/Milky Way/Star Clusters",
        "Name": "Gaia Stellar Clusters",
        "Description": "This is a large catalog containing 'Milky Way Stellar Clusters' using Gaia DR3 data.",
        "url": "https://www.amnh.org/research/hayden-planetarium/digital-universe",
        "license": "AMNH Gaia"
        },
        "meta_name": metadata['sub_project'],
		"author": metadata['prepared_by']
    }

    # Generate the file using your internal utility
    asset_creation.write_asset(datainfo)

if __name__ == "__main__":
    asset_main()
    print("Stellar clusters asset generated successfully.")